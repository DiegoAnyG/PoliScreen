"""Complejos receptor-pose e interacciones con PLIP.

Un solo motor de interacciones. PLIP protona internamente y conserva la numeración de autor
del PDB, que es la que aparece en la literatura y en los visores. El ligando acoplado se
escribe como HETATM/LIG para que PLIP lo identifique y no confunda cofactores con el ligando.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Optional, Sequence

import pandas as pd

from .screening import compound_from_pose_name, receptor_from_name

LIGAND_RESNAME = "LIG"
LIGAND_CHAIN = "Z"
LIGAND_RESSEQ = 900
WATERS = {"HOH", "WAT", "H2O", "DOD"}

# etiqueta del XML de PLIP -> tipo de interacción que usamos como columna
PLIP_TAGS = {
    "hbond": ["hydrogen_bond"],
    "hydrophobic": ["hydrophobic_interaction"],
    "saltbridge": ["salt_bridge"],
    "water": ["water_bridge"],
    "halogen": ["halogen_bond"],
    "pistack": ["pi_stack"],
    "pication": ["pi_cation_interaction"],
    "metal": ["metal_complex"],
}


class InteractionsError(RuntimeError):
    pass


# ---------------------------------------------------------------- fusión
def receptor_text(pdb, cache_dir=None) -> str:
    """Cuerpo del receptor sin las líneas de cierre, convirtiendo a PDB si hace falta."""
    pdb = Path(pdb)
    if pdb.suffix.lower() != ".pdb":
        cache = Path(cache_dir or pdb.parent)
        cache.mkdir(parents=True, exist_ok=True)
        conv = cache / f"{pdb.stem}_rec.pdb"
        if not conv.exists():
            subprocess.run(["obabel", str(pdb), "-O", str(conv)], capture_output=True)
        pdb = conv
    return "".join(l for l in pdb.read_text(errors="ignore").splitlines(keepends=True)
                   if not l.strip().startswith(("END", "ENDMDL")))


def ligand_text(pose_pdb) -> str:
    """Reescribe la pose como HETATM/LIG en su propia cadena. Conserva los CONECT."""
    out = []
    for l in Path(pose_pdb).read_text(errors="ignore").splitlines():
        if l.startswith(("ATOM", "HETATM")):
            try:
                serial = l[6:11].strip()
                name = l[12:16]
                x, y, z = float(l[30:38]), float(l[38:46]), float(l[46:54])
                el = (l[76:78].strip() or name.strip()[0])
                out.append(f"{'HETATM':<6}{serial:>5} {name:<4} {LIGAND_RESNAME:>3} {LIGAND_CHAIN}{LIGAND_RESSEQ:>4}"
                           f"    {x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00          {el:>2}\n")
            except Exception:
                continue
        elif l.startswith("CONECT"):
            out.append(l + "\n")
    return "".join(out)


def fuse(receptor, pose_pdb, out_pdb, cache_dir=None) -> Path:
    out_pdb = Path(out_pdb)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    body = receptor_text(receptor, cache_dir)
    if not body.endswith("\n"):
        body += "\n"
    out_pdb.write_text(body + ligand_text(pose_pdb) + "END\n")
    return out_pdb


def fuse_batch(receptors: Sequence, poses_dir, out_dir, cache_dir=None, on_progress=None,
               stem_to_file=None) -> list:
    """Fusiona cada pose con su receptor. Salta las que ya existen.

    El receptor se deduce del nombre de la pose. `stem_to_file` (id_sitio -> archivo) permite que
    varios sitios de docking apunten al mismo PDB (docking hibrido); si no se pasa, cada receptor es
    su propio sitio.
    """
    by_stem = dict(stem_to_file) if stem_to_file else {Path(r).stem: r for r in receptors}
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    models = sorted(Path(poses_dir).glob("*-model*.pdb"))
    made = []
    for i, pose in enumerate(models, 1):
        out = out_dir / f"Complejo_{pose.stem}.pdb"
        if not out.exists():
            rec = by_stem.get(receptor_from_name(pose.name))
            if rec is None:
                continue
            fuse(rec, pose, out, cache_dir)
            made.append(out)
        if on_progress:
            on_progress(i, len(models), out.name)
    return made


# ---------------------------------------------------------------- PLIP
def sanitize_pdb(src, dst) -> bool:
    """Reescribe el complejo con columnas limpias; PLIP es estricto con el formato."""
    clean = []
    for line in Path(src).read_text(errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                name = line[12:16]
                res = line[17:20].strip()
                chain = line[21] if len(line) > 21 else "A"
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                serial, resseq = line[6:11].strip(), line[22:26].strip()
                el = (line[76:78].strip() or name.strip()[0])
                clean.append(f"{line[0:6]:<6}{serial:>5} {name:<4} {res:>3} {chain}{resseq:>4}"
                             f"    {x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00 20.00          {el:>2}\n")
            except Exception:
                continue
        elif line.startswith("TER"):
            clean.append("TER\n")
    if not clean:
        return False
    Path(dst).parent.mkdir(parents=True, exist_ok=True)
    Path(dst).write_text("".join(clean) + "END\n")
    return True


def parse_plip_xml(xml_path, ligand_resname: str = LIGAND_RESNAME) -> dict:
    """Devuelve {Residuo+numero_tipo: conteo} solo del ligando acoplado; ignora cofactores."""
    feats = {}
    try:
        root = ET.parse(xml_path).getroot()
    except Exception:
        return feats
    for site in root.iter("bindingsite"):
        hetid = (site.findtext(".//hetid") or "").strip().upper()
        if hetid and hetid != ligand_resname.upper():
            continue
        for itype, tags in PLIP_TAGS.items():
            for tag in tags:
                for it in site.iter(tag):
                    rn = (it.findtext("resnr") or "").strip()
                    rt = (it.findtext("restype") or "").strip()
                    if rn and rt:
                        k = f"{rt.capitalize()}{rn}_{itype}"
                        feats[k] = feats.get(k, 0) + 1
    return feats


def plip_available() -> bool:
    return shutil.which("plip") is not None


def run_plip(complex_pdb, xml_dir, san_dir) -> Optional[Path]:
    """Corre PLIP sobre un complejo y deja su XML. Reutiliza el que ya exista."""
    complex_pdb = Path(complex_pdb)
    xml_dir, san_dir = Path(xml_dir), Path(san_dir)
    xml_dir.mkdir(parents=True, exist_ok=True)
    san_dir.mkdir(parents=True, exist_ok=True)
    xml = xml_dir / f"{complex_pdb.stem}.xml"
    if xml.exists():
        return xml
    san = san_dir / f"fx_{complex_pdb.stem}.pdb"
    if not sanitize_pdb(complex_pdb, san):
        return None
    tmp = xml_dir / f"_tmp_{complex_pdb.stem}"
    tmp.mkdir(exist_ok=True)
    subprocess.run(["plip", "-f", str(san), "-x", "--nopdb", "-o", str(tmp)], capture_output=True, text=True)
    found = list(tmp.glob("*.xml"))
    if found:
        found[0].rename(xml)
    shutil.rmtree(tmp, ignore_errors=True)
    return xml if xml.exists() else None


def crystal_fingerprint(receptor, control_file, work_dir, cache_dir=None) -> dict:
    """Huella PLIP del ligando Cristalográfico en su pose real (la del archivo del control, no una
    pose dockeada). Fusiona receptor + control y perfila con PLIP. Es la línea base objetiva: las
    interacciones de los compuestos se comparan contra esta, no contra el docking del control."""
    work = Path(work_dir); work.mkdir(parents=True, exist_ok=True)
    cf = Path(control_file)
    if cf.suffix.lower() != ".pdb":                  # PLIP necesita el ligando en PDB para fusionarlo
        pdbc = work / f"{cf.stem}_xtal.pdb"
        if not pdbc.exists():
            # -r: fragmento mayor. Si el control salio fragmentado/duplicado, la huella de referencia
            # no debe incluir esos restos espurios.
            subprocess.run(["obabel", str(cf), "-O", str(pdbc), "-r"], capture_output=True)
        cf = pdbc
    if not cf.exists() or cf.stat().st_size == 0:
        return {}
    comp = work / f"Complejo_xtal_{Path(receptor).stem}.pdb"
    fuse(receptor, cf, comp, cache_dir)
    xml = run_plip(comp, work / "xml", work / "san")
    return parse_plip_xml(xml) if xml else {}


def hetero_fingerprint(receptor_pdb, box, work_dir, min_atoms: int = 6) -> tuple:
    """Huella PLIP de un cofactor/ligando co-cristalizado que quede DENTRO de la caja de un sitio.
    Sirve de referencia objetiva para un bolsillo secundario del docking hibrido. Devuelve (feats, etiqueta)
    o ({}, None) si no hay ninguno de tamaño suficiente (los iones de 1 átomo no cuentan).

    Se corre PLIP sobre el receptor tal cual (que ya contiene el cofactor) y se lee su sitio de unión,
    en vez de re-anexarlo, para no contar el cofactor dos veces."""
    receptor_pdb = Path(receptor_pdb)
    hx, hy, hz = box.sx / 2.0, box.sy / 2.0, box.sz / 2.0
    grupos = {}
    for l in receptor_pdb.read_text(errors="ignore").splitlines():
        if not l.startswith("HETATM"):
            continue
        rn = l[17:20].strip()
        if rn in WATERS:
            continue
        try:
            x, y, z = float(l[30:38]), float(l[38:46]), float(l[46:54])
        except ValueError:
            continue
        if abs(x - box.cx) <= hx and abs(y - box.cy) <= hy and abs(z - box.cz) <= hz:
            grupos.setdefault((rn, (l[21].strip() or "_"), l[22:26].strip()), []).append(l)
    if not grupos:
        return {}, None
    (rn, ch, seq), atoms = max(grupos.items(), key=lambda kv: len(kv[1]))
    if len(atoms) < min_atoms:
        return {}, None
    work = Path(work_dir); work.mkdir(parents=True, exist_ok=True)
    xml = run_plip(receptor_pdb, work / "xml", work / "san")
    if not xml:
        return {}, None
    return parse_plip_xml(xml, ligand_resname=rn), f"{rn} {ch}:{seq}"


def crystal_fingerprints(receptors: Sequence, controls: Sequence, control_assign: dict,
                         work_dir, cache_dir=None) -> dict:
    """{stem_receptor: feats} tomando, para cada receptor, el control que tiene asignado.
    Con un solo control se usa para el único receptor. Sirve de referencia objetiva del ranking."""
    from .screening import normalize_key
    by_ck = {normalize_key(Path(c).stem): c for c in controls}
    out = {}
    for r in receptors:
        stem = Path(r).stem
        ck = next((k for k, rc in control_assign.items() if rc == stem), None)
        cfile = by_ck.get(ck) if ck else None
        if cfile is None and len(controls) == 1:
            cfile = controls[0]
        if cfile:
            fp = crystal_fingerprint(r, cfile, work_dir, cache_dir)
            if fp:
                out[stem] = fp
    return out


def plip_batch(complexes: Sequence, work_dir, workers: int = 0, force: bool = False,
               on_progress=None) -> pd.DataFrame:
    """Perfila todos los complejos y devuelve la matriz pose x interacción.

    Reutiliza interacciones.csv salvo que force=True. Escribe el CSV en work_dir.
    """
    work = Path(work_dir)
    csv = work / "interacciones.csv"
    xml_dir, san_dir = work / "xml_plip", work / "san"
    esperados = {Path(c).stem for c in complexes}
    # La cache solo vale si cubre TODOS los complejos actuales. Devolver un CSV que no incluye
    # los compuestos recien acoplados deja su huella vacia y falsea el ranking entero.
    if csv.exists() and not force:
        previo = pd.read_csv(csv)
        if esperados and esperados.issubset(set(previo.get("name", []))):
            return previo
    if force:
        for d in (xml_dir, san_dir):
            shutil.rmtree(d, ignore_errors=True)
    if not plip_available():
        raise InteractionsError(
            "PLIP no esta instalado. Instalalo sin compilar OpenBabel: "
            "conda install -c conda-forge -c bioconda plip"
        )
    complexes = [Path(c) for c in complexes]
    if not complexes:
        raise InteractionsError("No hay complejos que analizar. Corre antes la fusion.")
    n = workers if workers > 0 else max(1, (os.cpu_count() or 2) - 2)
    with ThreadPoolExecutor(max_workers=n) as ex:
        for i, _ in enumerate(ex.map(lambda c: run_plip(c, xml_dir, san_dir), complexes), 1):
            if on_progress:
                on_progress(i, len(complexes))
    # Solo los complejos de ESTA corrida: así el CSV nunca arrastra compuestos de corridas viejas.
    rows = []
    for c in complexes:
        x = xml_dir / f"{c.stem}.xml"
        if not x.exists():
            continue
        row = parse_plip_xml(x)
        row["name"] = x.stem
        rows.append(row)
    if not rows:
        raise InteractionsError("PLIP no produjo resultados. Revisa que los complejos tengan el ligando como LIG.")
    df = pd.DataFrame(rows).fillna(0)
    cols = ["name"] + sorted(c for c in df.columns if c != "name")
    df = df[cols]
    for c in cols[1:]:
        df[c] = df[c].astype(int)
    work.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv, index=False)
    return df
