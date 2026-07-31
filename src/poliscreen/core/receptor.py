"""Preparación de receptores: descarga del PDB, inspección y limpieza tipo DockPrep.

La protonación usa PDBFixer porque conserva la numeración de autor, la que aparece en la literatura
y la que reporta PLIP. PDBFixer elimina todo lo no proteico, así que los cofactores a conservar se
extraen, se protonan y se vuelven a unir aparte.
"""
from __future__ import annotations

import subprocess
import tempfile
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

WATERS = {"HOH", "WAT", "H2O", "DOD"}
RCSB_URL = "https://files.rcsb.org/download/{}.pdb"


class ReceptorError(RuntimeError):
    pass


@dataclass(frozen=True)
class Het:
    """Un grupo heteroatomico: ligando, cofactor, ion o agua."""
    resname: str
    chain: str
    resseq: str
    n_atoms: int

    @property
    def key(self) -> str:
        return f"{self.resname}|{self.chain}|{self.resseq}"

    @property
    def label(self) -> str:
        return f"{self.resname} {self.chain}:{self.resseq} ({self.n_atoms} atomos)"

    @property
    def is_water(self) -> bool:
        return self.resname in WATERS


@dataclass(frozen=True)
class Modificado:
    """Un residuo modificado de una cadena: fosfotirosina, selenometionina, fosfoserina...

    Aunque el PDB lo declare como heteroátomo, forma parte de la cadena y su modificación suele ser
    funcional. Se distingue para decidir si se conserva o se sustituye por el aminoácido de origen.
    """
    resname: str
    chain: str
    resseq: str
    estandar: str

    @property
    def key(self) -> str:
        return f"{self.resname}|{self.chain}|{self.resseq}"

    @property
    def label(self) -> str:
        return f"{self.resname} {self.chain}:{self.resseq} (deriva de {self.estandar})"


def modified_residues(pdb) -> list:
    """Residuos modificados de las cadenas, detectados con la tabla de PDBFixer.

    Se usa su misma detección, no una lista propia, para que lo ofrecido al usuario coincida con lo
    que ocurriría sin intervenir. El PDB los declara también en MODRES, pero ese registro falta a veces.
    """
    try:
        from pdbfixer import PDBFixer
    except Exception:
        return []
    try:
        f = PDBFixer(filename=str(pdb))
        f.findNonstandardResidues()
        return [Modificado(r.name, r.chain.id, str(r.id).strip(), nuevo)
                for r, nuevo in f.nonstandardResidues]
    except Exception:
        return []


@dataclass
class Structure:
    path: Path
    chains: list = field(default_factory=list)
    het: list = field(default_factory=list)     # sin aguas
    n_waters: int = 0
    n_atoms: int = 0

    def find(self, key: str) -> Optional[Het]:
        return next((h for h in self.het if h.key == key), None)

    def summary(self) -> str:
        lines = [f"{self.path.name}: {self.n_atoms} atomos, cadenas {','.join(self.chains) or '-'}, "
                 f"{self.n_waters} aguas"]
        lines += [f"  {h.label}" for h in self.het]
        return "\n".join(lines)


def fetch_pdb(pdb_id: str, out_dir) -> Path:
    """Descarga un PDB por su identificador. Reutiliza el archivo si ya existe."""
    pdb_id = pdb_id.strip().upper()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / f"{pdb_id}.pdb"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    try:
        urllib.request.urlretrieve(RCSB_URL.format(pdb_id), dest)
    except Exception as e:
        raise ReceptorError(f"No pude descargar {pdb_id}: {e}. "
                            "Causa probable: identificador invalido o sin conexion.")
    return dest


def inspect(pdb) -> Structure:
    """Lee cadenas, heterogrupos y aguas sin modificar nada."""
    pdb = Path(pdb)
    chains, counts, n_atoms, n_waters = [], {}, 0, 0
    for line in pdb.read_text(errors="ignore").splitlines():
        if line.startswith("ATOM"):
            n_atoms += 1
            ch = (line[21].strip() or "_")
            if ch not in chains:
                chains.append(ch)
        elif line.startswith("HETATM"):
            n_atoms += 1
            rn = line[17:20].strip()
            if rn in WATERS:
                n_waters += 1
                continue
            k = (rn, (line[21].strip() or "_"), line[22:26].strip())
            counts[k] = counts.get(k, 0) + 1
    het = [Het(rn, ch, rs, n) for (rn, ch, rs), n in sorted(counts.items())]
    return Structure(path=pdb, chains=chains, het=het, n_waters=n_waters, n_atoms=n_atoms)


def _het_lines(pdb, het: Het) -> list:
    out = []
    for l in Path(pdb).read_text(errors="ignore").splitlines():
        if (l.startswith("HETATM") and l[17:20].strip() == het.resname
                and (l[21].strip() or "_") == het.chain and l[22:26].strip() == het.resseq):
            out.append(l)
    return out


def extract_ligand(pdb, het: Het, out_path, ph: float = 7.4, smiles: Optional[str] = None) -> Path:
    """Extrae un heterogrupo como molécula independiente, útil como control de referencia.

    Con `smiles` se corrigen los ordenes de enlace desde una plantilla: el PDB no los guarda
    y sin ellos algunos ligandos quedan con valencias imposibles.
    """
    lines = _het_lines(pdb, het)
    if not lines:
        raise ReceptorError(f"No encontre el grupo {het.label} en {Path(pdb).name}.")
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        raw = Path(td) / "raw.pdb"
        raw.write_text("\n".join(lines) + "\nEND\n")
        if smiles:
            try:
                from rdkit import Chem, RDLogger
                from rdkit.Chem import AllChem
                RDLogger.DisableLog("rdApp.*")
                tmpl = Chem.MolFromSmiles(smiles)
                m = Chem.MolFromPDBFile(str(raw), removeHs=True, sanitize=False)
                if tmpl is not None and m is not None:
                    m.UpdatePropertyCache(strict=False)
                    m = AllChem.AssignBondOrdersFromTemplate(tmpl, m)
                    m = Chem.AddHs(m, addCoords=True)
                    w = Chem.SDWriter(str(out_path))
                    m.SetProp("_Name", f"{het.resname}_{het.chain}{het.resseq}")
                    w.write(m)
                    w.close()
                    if out_path.exists() and out_path.stat().st_size > 0:
                        return out_path
            except Exception:
                pass  # si la plantilla no aplica, se cae a obabel
        # -r: se queda con el fragmento conectado mayor. Extraer por geometría (sin plantilla SMILES)
        # puede fragmentar o duplicar el ligando; el control debe ser una sola molécula limpia.
        subprocess.run(["obabel", str(raw), "-O", str(out_path), "-p", str(ph), "-r"], capture_output=True)
    if not out_path.exists() or out_path.stat().st_size == 0:
        raise ReceptorError(f"No pude escribir el ligando {het.label}.")
    return out_path


def extract_chain(pdb, chain: str, out_path, solo_polipeptido: bool = True,
                  on_aviso=None) -> Path:
    """Extrae una cadena como control de referencia. Devuelve la ruta escrita.

    Un péptido cocristalizado es una cadena del modelo, no un heteroátomo, así que no figura entre
    los cofactores extraíbles. Sin esta vía, un cribado de péptidos tendría que compararse con una
    molécula pequeña —otro motor, otro quimiotipo— en vez de con el propio péptido del cristal.

    solo_polipeptido descarta los heteroátomos de la cadena. Una cadena puede llevar grupos
    conjugados (un análogo de nucleótido unido por un conector, en un inhibidor bisustrato) que el
    motor peptídico no reproduce; conservarlos haría a referencia y pose moléculas distintas y el
    RMSD de validación no existiría. Se escribe en PDB, no SDF, para conservar nombres de residuo y
    átomo, que son los que permiten reconocer la secuencia y derivar el control al motor peptídico.
    """
    pdb, out_path = Path(pdb), Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    de_la_cadena = [l for l in pdb.read_text(errors="ignore").splitlines()
                    if l.startswith(("ATOM", "HETATM")) and l[21] == chain
                    and l[17:20].strip() not in WATERS]
    if not de_la_cadena:
        raise ReceptorError(f"La cadena {chain} no tiene atomos en {pdb.name}.")

    lineas = [l for l in de_la_cadena if l.startswith("ATOM")] if solo_polipeptido else de_la_cadena
    if not lineas:
        lineas = de_la_cadena
    fuera = {l[17:20].strip() for l in de_la_cadena if l not in lineas}
    if fuera and on_aviso:
        on_aviso(f"La cadena {chain} lleva grupos no peptídicos ({', '.join(sorted(fuera))}) que se "
                 f"han excluido del control: el acoplamiento de péptidos reproduce solo la parte "
                 f"peptídica, y conservarlos haría que la referencia y la pose fuesen moléculas "
                 f"distintas.")
    out_path.write_text("\n".join(lineas) + "\nEND\n")
    return out_path


def prepare(pdb, out_path, keep_chains: Optional[Sequence[str]] = None,
            keep_het: Sequence[str] = (), ph: float = 7.4, add_hydrogens: bool = True,
            keep_modified: Sequence[str] = (), on_aviso=None,
            add_missing_residues: bool = False) -> Path:
    """Deja el receptor listo para acoplar: sin aguas, con hidrogenos y con lo que se pida conservar.

    keep_chains           cadenas a conservar (None = todas)
    keep_het              claves de heterogrupos a conservar, p. ej. un cofactor
    add_missing_residues  reconstruir lazos ausentes; por defecto no, para no inventar
                          geometría cerca del sitio de unión
    """
    pdb, out_path = Path(pdb), Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        from openmm.app import PDBFile
        from pdbfixer import PDBFixer
    except Exception as e:
        raise ReceptorError(f"Falta PDBFixer/OpenMM: {e}. Instalalos con conda para preparar receptores.")

    fixer = PDBFixer(filename=str(pdb))
    if keep_chains:
        keep = {c.strip() for c in keep_chains}
        drop = [ch.id for ch in fixer.topology.chains() if ch.id not in keep]
        if drop:
            fixer.removeChains(chainIds=drop)
    fixer.findMissingResidues()
    if not add_missing_residues:
        fixer.missingResidues = {}
    # Residuos modificados (fosfotirosina, selenometionina...). Por defecto se sustituyen por su
    # aminoácido estándar, que es lo que piden los campos de fuerza; pero eso ELIMINA la
    # modificación, y a menudo la modificación es la función (un bucle de activación fosforilado deja
    # de estarlo). Lo decide quien prepara, y para eso los ve: se listan aparte de los cofactores.
    fixer.findNonstandardResidues()
    todos_mod = list(getattr(fixer, "nonstandardResidues", []))

    def _clave(res):
        return f"{res.name}|{res.chain.id}|{str(res.id).strip()}"

    conservar = set(keep_modified or ())
    fixer.nonstandardResidues = [(r, n) for r, n in todos_mod if _clave(r) not in conservar]
    # Todos los modificados se excluyen del paso de cofactores: ya están en la cadena, y reañadirlos
    # pondría dos átomos en la misma posición, lo que rompe el cálculo de ángulos y duplica contactos.
    reemplazados = set()
    for res, _n in todos_mod:
        try:
            reemplazados.add((res.chain.id, str(res.id).strip()))
        except Exception:
            continue
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    if add_hydrogens:
        fixer.addMissingHydrogens(ph)
    with open(out_path, "w") as fh:
        # keepIds imprescindible: sin él OpenMM renumera desde 1 y TYR157 pasa a VAL157, un fallo
        # silencioso que arruina el análisis farmacofórico.
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)

    # Los modificados conservados hay que reponerlos: PDBFixer no los sustituye, pero su limpieza de
    # heteroátomos los borra por no ser residuos estándar y deja un hueco. Se reponen con sus líneas
    # originales (traen la modificación completa) y no se duplican porque ya no están en el archivo.
    repuestos = []
    for res, _n in todos_mod:
        clave = _clave(res)
        if clave not in conservar:
            continue
        lineas = [l for l in Path(pdb).read_text(errors="ignore").splitlines()
                  if l.startswith(("ATOM", "HETATM")) and l[21] == res.chain.id
                  and l[22:27].strip() == str(res.id).strip()
                  and l[17:20].strip() == res.name]
        if not lineas:
            continue
        cuerpo = "".join(l for l in out_path.read_text().splitlines(keepends=True)
                         if not l.strip().startswith(("END", "ENDMDL")))
        out_path.write_text(cuerpo + "".join(l + "\n" for l in lineas) + "END\n")
        repuestos.append(f"{res.name} {res.chain.id}:{str(res.id).strip()}")

    omitidos = []
    for key in keep_het:
        st = inspect(pdb)
        het = st.find(key)
        if het is None:
            raise ReceptorError(f"No existe el heterogrupo {key} en {pdb.name}.")
        if (het.chain, str(het.resseq).strip()) in reemplazados:
            # Residuo modificado de la cadena, no cofactor: ya está; reañadirlo duplicaría sus átomos.
            omitidos.append(het.label)
            continue
        # Líneas originales tal cual: protonar el fragmento con obabel lo renombraría y PLIP dejaría
        # de distinguir el cofactor del ligando acoplado. Sus hidrógenos los pone PLIP internamente.
        body = "".join(l for l in out_path.read_text().splitlines(keepends=True)
                       if not l.strip().startswith(("END", "ENDMDL")))
        extra = "".join(l + "\n" for l in _het_lines(pdb, het))
        out_path.write_text(body + extra + "END\n")

    n_dup = _quitar_superpuestos(out_path)
    if on_aviso:
        for lab in omitidos:
            on_aviso(f"{lab} es un residuo modificado de la cadena, no un cofactor: se trata en su "
                     f"propia sección y no se añade por separado.")
        sustituidos = [f"{r.name} {r.chain.id}:{str(r.id).strip()}"
                       for r, _n in todos_mod if _clave(r) not in conservar]
        if sustituidos:
            on_aviso("Sustituidos por su aminoácido estándar (se pierde la modificación): "
                     + ", ".join(sustituidos))
        if repuestos:
            on_aviso("Conservados con su modificación: " + ", ".join(repuestos))
        if n_dup:
            on_aviso(f"Se eliminaron {n_dup} átomo(s) superpuestos del receptor preparado.")
    return out_path


def _quitar_superpuestos(pdb: Path, tol: float = 1e-3) -> int:
    """Elimina átomos que ocupan la misma posición. Devuelve cuántos quitó.

    Dos átomos en el mismo punto hacen fallar el cálculo de ángulos que asigna tipos e hibridación
    (ADFRsuite aborta con división por cero) y duplican cada contacto. Siempre es patológico.
    """
    vistos, salida, quitados = set(), [], 0
    for l in pdb.read_text(errors="ignore").splitlines(keepends=True):
        if l.startswith(("ATOM", "HETATM")):
            try:
                k = (round(float(l[30:38]) / tol), round(float(l[38:46]) / tol),
                     round(float(l[46:54]) / tol))
            except ValueError:
                salida.append(l)
                continue
            if k in vistos:
                quitados += 1
                continue
            vistos.add(k)
        salida.append(l)
    if quitados:
        pdb.write_text("".join(salida))
    return quitados
