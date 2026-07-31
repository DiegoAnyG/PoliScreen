"""Generación combinatoria de péptidos y sus propiedades.

Segunda vía de diseño, independiente de la síntesis por reacción: aquí los ligandos no se
construyen uniendo un nucleo a un reactivo, sino enumerando secuencias de aminoacidos bajo
reglas. Pensada para péptidos antimicrobianos, donde la actividad depende sobre todo de la
carga neta positiva y de la anfipaticidad, no de un farmacoforo puntual.

La generación es determinista: misma semilla y mismas reglas producen la misma biblioteca.
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

# Un control cocristalizado más largo que esto no es un péptido de diseño sino una cadena proteica,
# y tratarlo como ligando no tendría sentido.
MAX_LARGO = 60

MIN_LONGITUD = 1
MAX_LONGITUD = 20

# Aminoacido: nombre, clases a las que pertenece, carga a pH 7.4 e hidropatía (Kyte-Doolittle).
# Las clases son las que el usuario combina para restringir el alfabeto.
AMINOACIDOS = {
    "A": ("Alanine",       ("hidrofobico", "alifatico", "no_esencial"),            0.0,  1.8),
    "R": ("Arginine",      ("hidrofilico", "cargado_pos", "no_esencial"),          1.0, -4.5),
    "N": ("Asparagine",    ("hidrofilico", "polar", "no_esencial"),                0.0, -3.5),
    "D": ("Aspartate",     ("hidrofilico", "cargado_neg", "no_esencial"),         -1.0, -3.5),
    "C": ("Cysteine",      ("hidrofobico", "azufrado", "no_esencial"),             0.0,  2.5),
    "E": ("Glutamate",     ("hidrofilico", "cargado_neg", "no_esencial"),         -1.0, -3.5),
    "Q": ("Glutamine",     ("hidrofilico", "polar", "no_esencial"),                0.0, -3.5),
    "G": ("Glycine",       ("especial", "no_esencial"),                            0.0, -0.4),
    "H": ("Histidine",     ("hidrofilico", "aromatico", "cargado_pos", "esencial"), 0.1, -3.2),
    "I": ("Isoleucine",    ("hidrofobico", "alifatico", "esencial"),               0.0,  4.5),
    "L": ("Leucine",       ("hidrofobico", "alifatico", "esencial"),               0.0,  3.8),
    "K": ("Lysine",        ("hidrofilico", "cargado_pos", "esencial"),             1.0, -3.9),
    "M": ("Methionine",    ("hidrofobico", "azufrado", "esencial"),                0.0,  1.9),
    "F": ("Phenylalanine", ("hidrofobico", "aromatico", "esencial"),               0.0,  2.8),
    "P": ("Proline",       ("especial", "no_esencial"),                            0.0, -1.6),
    "S": ("Serine",        ("hidrofilico", "polar", "no_esencial"),                0.0, -0.8),
    "T": ("Threonine",     ("hidrofilico", "polar", "esencial"),                   0.0, -0.7),
    "W": ("Tryptophan",    ("hidrofobico", "aromatico", "esencial"),               0.0, -0.9),
    "Y": ("Tyrosine",      ("hidrofilico", "aromatico", "no_esencial"),            0.0, -1.3),
    "V": ("Valine",        ("hidrofobico", "alifatico", "esencial"),               0.0,  4.2),
}

CLASES = {
    "esencial":     "Essential (not synthesized by the body)",
    "no_esencial":  "Non-essential",
    "hidrofobico":  "Hydrophobic",
    "hidrofilico":  "Hydrophilic",
    "aromatico":    "Aromatic (F, W, Y, H)",
    "alifatico":    "Aliphatic (A, I, L, V)",
    "polar":        "Uncharged polar",
    "cargado_pos":  "Positive charge (K, R, H)",
    "cargado_neg":  "Negative charge (D, E)",
    "azufrado":     "Sulfur-containing (C, M)",
    "especial":     "Special (G, P)",
}

# Hidrofobicidad consenso de Eisenberg: se usa para el momento hidrofóbico, que mide la
# anfipaticidad de una hélice y es el descriptor que mejor separa péptidos antimicrobianos
# activos de inactivos.
_EISENBERG = {
    "A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29, "E": -0.74, "Q": -0.85,
    "G": 0.48, "H": -0.40, "I": 1.38, "L": 1.06, "K": -1.50, "M": 0.64, "F": 1.19,
    "P": 0.12, "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08,
}
# Índice de Boman: energía de unión a proteínas. Alto = tiende a unirse a muchas proteínas
# (potencial promiscuidad); bajo = más selectivo.
_BOMAN = {
    "A": -0.5, "R": 2.5, "N": 0.2, "D": 3.0, "C": -1.0, "E": 3.0, "Q": 0.2, "G": 0.0,
    "H": -0.5, "I": -1.8, "L": -1.8, "K": 3.0, "M": -1.3, "F": -2.5, "P": 0.0,
    "S": 0.3, "T": -0.4, "W": -3.4, "Y": -2.3, "V": -1.5,
}

ANGULO_HELICE = 100.0    # grados por residuo en una hélice alfa


def alfabeto(incluir: Sequence[str] = (), excluir_clases: Sequence[str] = (),
             excluir_residuos: Sequence[str] = ()) -> list:
    """Alfabeto de aminoacidos permitido. `incluir` son clases (unión); vacio = los 20."""
    if incluir:
        sel = {a for a, (_n, cls, _c, _h) in AMINOACIDOS.items() if set(cls) & set(incluir)}
    else:
        sel = set(AMINOACIDOS)
    if excluir_clases:
        sel -= {a for a, (_n, cls, _c, _h) in AMINOACIDOS.items() if set(cls) & set(excluir_clases)}
    sel -= set(x.upper() for x in excluir_residuos)
    return sorted(sel)


@dataclass
class Reglas:
    """Restricciones de la biblioteca. Todas son opcionales y se combinan."""
    longitud: int = 7
    alfabeto: Sequence[str] = field(default_factory=lambda: sorted(AMINOACIDOS))
    sin_repetir: bool = False              # cada aminoacido aparece como mucho una vez
    max_consecutivos: int = 0              # 0 = sin límite; 2 prohíbe AAA pero permite AA
    max_por_residuo: int = 0               # 0 = sin límite; veces que puede aparecer un mismo residuo
    prefijo: str = ""                      # la secuencia empieza así
    sufijo: str = ""                       # y termina así
    carga_min: Optional[float] = None      # carga neta a pH 7.4
    carga_max: Optional[float] = None
    gravy_min: Optional[float] = None      # hidropatía media (Kyte-Doolittle)
    gravy_max: Optional[float] = None

    def validar(self) -> list:
        """Avisos legibles; lista vacia = reglas coherentes."""
        avisos = []
        if not (MIN_LONGITUD <= self.longitud <= MAX_LONGITUD):
            avisos.append(f"Length must be between {MIN_LONGITUD} and {MAX_LONGITUD}.")
        if not self.alfabeto:
            avisos.append("The alphabet is empty: no selected class leaves any amino acids.")
        fijos = len(self.prefijo) + len(self.sufijo)
        if fijos > self.longitud:
            avisos.append(f"Prefix and suffix add up to {fijos} residues and the length is {self.longitud}.")
        if self.sin_repetir and self.longitud > len(self.alfabeto):
            avisos.append(f"Without repeats, {self.longitud} residues do not fit in an alphabet of "
                          f"{len(self.alfabeto)}.")
        for patron, nombre in ((self.prefijo, "prefix"), (self.sufijo, "suffix")):
            fuera = set(patron.upper()) - set(AMINOACIDOS)
            if fuera:
                avisos.append(f"The {nombre} contains symbols that are not amino acids: {', '.join(sorted(fuera))}.")
        # Conflictos entre reglas: se detectan aquí para explicar por que no saldría ninguna
        # secuencia, en vez de dejar que la generación devuelva una lista vacia sin motivo.
        extremos = (self.prefijo + self.sufijo).upper()
        if self.sin_repetir and len(set(extremos)) != len(extremos):
            avisos.append("«No repeats» conflicts with the prefix or suffix, which already repeat a residue.")
        if self.max_consecutivos and self.max_consecutivos > 0:
            racha, previo = 1, ""
            for a in extremos:
                racha = racha + 1 if a == previo else 1
                previo = a
                if racha > self.max_consecutivos:
                    avisos.append(f"The prefix or suffix repeats {racha} residues in a row and the maximum "
                                  f"consecutive is {self.max_consecutivos}.")
                    break
        fuera_alfabeto = set(extremos) - set(self.alfabeto)
        if fuera_alfabeto:
            avisos.append(f"The prefix or suffix uses residues outside the chosen alphabet: "
                          f"{', '.join(sorted(fuera_alfabeto))}.")
        return avisos

    def espacio(self) -> float:
        """Cota superior del número de secuencias posibles. Sirve para avisar cuando se pide
        más biblioteca de la que el espacio combinatorio permite."""
        libres = max(0, self.longitud - len(self.prefijo) - len(self.sufijo))
        n = len(self.alfabeto)
        if libres == 0:
            return 1.0
        if self.sin_repetir:
            total, disp = 1.0, n
            for _ in range(libres):
                total *= max(disp, 0); disp -= 1
            return total
        return float(n) ** libres


def _cumple(seq: str, r: Reglas) -> bool:
    if r.sin_repetir and len(set(seq)) != len(seq):
        return False
    if r.max_por_residuo:
        for a in set(seq):
            if seq.count(a) > r.max_por_residuo:
                return False
    if r.max_consecutivos:
        racha, previo = 1, ""
        for a in seq:
            racha = racha + 1 if a == previo else 1
            if racha > r.max_consecutivos:
                return False
            previo = a
    if r.carga_min is not None or r.carga_max is not None:
        q = carga_neta(seq)
        if r.carga_min is not None and q < r.carga_min:
            return False
        if r.carga_max is not None and q > r.carga_max:
            return False
    if r.gravy_min is not None or r.gravy_max is not None:
        g = gravy(seq)
        if r.gravy_min is not None and g < r.gravy_min:
            return False
        if r.gravy_max is not None and g > r.gravy_max:
            return False
    return True


def generate(reglas: Reglas, n: int, seed: int = 42, max_intentos: int = 200) -> tuple:
    """Genera hasta `n` secuencias únicas que cumplen las reglas. Devuelve (secuencias, aviso).

    Muestreo aleatorio con semilla fija en vez de enumeracion exhaustiva: el espacio crece como
    20^longitud y enumerarlo es inviable salvo en péptidos muy cortos. Se para cuando junta n
    secuencias o cuando deja de encontrar nuevas, de modo que pedir más de las que existen no
    cuelga la aplicación.
    """
    rng = random.Random(seed)
    pre, suf = reglas.prefijo.upper(), reglas.sufijo.upper()
    libres = reglas.longitud - len(pre) - len(suf)
    if libres < 0 or not reglas.alfabeto:
        return [], "Las reglas no permiten construir ninguna secuencia."

    vistas, salida, fallos = set(), [], 0
    limite = max_intentos * max(n, 1)
    intentos = 0
    while len(salida) < n and intentos < limite and fallos < max_intentos * 20:
        intentos += 1
        if reglas.sin_repetir:
            usados = set(pre) | set(suf)
            libres_pool = [a for a in reglas.alfabeto if a not in usados]
            if len(libres_pool) < libres:
                break
            medio = "".join(rng.sample(libres_pool, libres))
        else:
            medio = "".join(rng.choice(reglas.alfabeto) for _ in range(libres))
        seq = pre + medio + suf
        if seq in vistas:
            fallos += 1
            continue
        if not _cumple(seq, reglas):
            fallos += 1
            continue
        vistas.add(seq); salida.append(seq); fallos = 0

    aviso = ""
    if len(salida) < n:
        esp = reglas.espacio()
        aviso = (f"{len(salida)} of {n} sequences generated: with these rules the available "
                 f"space is about {esp:.0f} and the filters discard the rest.")
    return salida, aviso


# ---------------------------------------------------------------- propiedades
def carga_neta(seq: str, ph: float = 7.4, c_amida: bool = False,
               n_acetil: bool = False, ciclico: bool = False) -> float:
    """Carga neta aproximada a pH fisiologico. En péptidos antimicrobianos es el descriptor
    más asociado a la actividad: la membrana bacteriana es anionica.

    La química de los extremos entra aquí porque cambia la carga, que es justamente el descriptor
    que se usa para decidir: amidar el carboxilo suprime una carga negativa (+1 neto), acetilar el
    amino suprime una positiva (-1 neto) y el cierre cabeza-cola consume ambos extremos, de modo
    que solo quedan las cargas de las cadenas laterales.
    """
    q = sum(AMINOACIDOS[a][2] for a in seq.upper() if a in AMINOACIDOS)
    if ciclico:
        return round(q, 2)
    if not n_acetil:
        q += 1.0                          # amino terminal protonado
    if not c_amida:
        q -= 1.0                          # carboxilo terminal desprotonado
    return round(q, 2)


def gravy(seq: str) -> float:
    """Hidropatía media (Kyte-Doolittle). Positivo = hidrofóbico global."""
    vals = [AMINOACIDOS[a][3] for a in seq.upper() if a in AMINOACIDOS]
    return round(sum(vals) / len(vals), 3) if vals else 0.0


def momento_hidrofobico(seq: str, angulo: float = ANGULO_HELICE) -> float:
    """Momento hidrofóbico normalizado de Eisenberg. Mide la ANFIPATICIDAD: si al plegarse en
    hélice los residuos hidrofóbicos quedan en una cara y los polares en la otra. Es lo que
    permite a un péptido insertarse en la membrana."""
    s = [a for a in seq.upper() if a in _EISENBERG]
    if not s:
        return 0.0
    rad = math.radians(angulo)
    sx = sum(_EISENBERG[a] * math.cos(i * rad) for i, a in enumerate(s))
    sy = sum(_EISENBERG[a] * math.sin(i * rad) for i, a in enumerate(s))
    return round(math.hypot(sx, sy) / len(s), 3)


def indice_boman(seq: str) -> float:
    """Potencial de unión a otras proteínas (kcal/mol). Por encima de ~2.5 se considera que el
    péptido es promiscuo y puede tener más efectos inespecificos."""
    vals = [_BOMAN[a] for a in seq.upper() if a in _BOMAN]
    return round(sum(vals) / len(vals), 2) if vals else 0.0


def fraccion_hidrofobica(seq: str) -> float:
    s = [a for a in seq.upper() if a in AMINOACIDOS]
    if not s:
        return 0.0
    h = sum(1 for a in s if "hidrofobico" in AMINOACIDOS[a][1])
    return round(h / len(s), 3)


_TRES_A_UNA = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
               "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
               "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
               "TYR": "Y", "VAL": "V"}


def secuencia_de_estructura(ruta) -> Optional[tuple]:
    """(secuencia, ciclico) si el archivo contiene un polipeptido; None si no lo es.

    Permite reconocer un péptido sin depender de ningun metadato: hace falta para los controles
    cocristalizados, que se extraen del cristal y no pasan por la tabla de ligandos. Se exige
    esqueleto completo (N, CA y C) en cada residuo para no confundir con un polipeptido cualquier
    molécula que contenga fragmentos de aminoacido.
    """
    p = Path(ruta)
    if p.suffix.lower() not in (".pdb", ".pdbqt", ".ent"):
        return None
    res, esq = {}, {}
    try:
        for l in p.read_text(errors="ignore").splitlines():
            if not l.startswith(("ATOM", "HETATM")):
                continue
            nom3, clave = l[17:20].strip().upper(), (l[21], l[22:27].strip())
            if nom3 not in _TRES_A_UNA:
                continue
            res[clave] = nom3
            esq.setdefault(clave, set()).add(l[12:16].strip())
    except Exception:
        return None
    completos = [k for k, v in res.items() if {"N", "CA", "C"} <= esq.get(k, set())]
    if not 1 <= len(completos) <= MAX_LARGO:
        return None
    orden = sorted(completos, key=lambda k: (k[0], _num(k[1])))
    seq = "".join(_TRES_A_UNA[res[k]] for k in orden)

    # Ciclado cabeza-cola: el nitrogeno del primer residuo esta a distancia de enlace del carbono
    # del último. Se mide sobre las coordenadas porque un PDB no declara los enlaces.
    ciclo = False
    xyz = {}
    for l in p.read_text(errors="ignore").splitlines():
        if l.startswith(("ATOM", "HETATM")):
            k = (l[21], l[22:27].strip())
            if k in (orden[0], orden[-1]):
                xyz[(k, l[12:16].strip())] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    a, b = xyz.get((orden[0], "N")), xyz.get((orden[-1], "C"))
    if a and b:
        ciclo = sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5 < 1.8
    return seq, ciclo


def _num(s: str) -> int:
    d = "".join(c for c in s if c.isdigit() or c == "-")
    return int(d) if d and d != "-" else 0


def etiqueta(seq: str, n_acetil: bool = False, c_amida: bool = False,
             ciclico: bool = False) -> str:
    """Nombre que hace visible la química de los extremos: 'Ac-KWKLF-NH2', 'ciclo-KWKLF'.

    Sin él, dos moléculas distintas comparten nombre en tablas y archivos de poses, y un cribado con
    extremos protegidos sería indistinguible de uno sin proteger. Sin caracteres que compliquen el
    nombre de archivo.
    """
    seq = seq.upper()
    if ciclico:
        return f"cyclo-{seq}"
    return ("Ac-" if n_acetil else "") + seq + ("-NH2" if c_amida else "")


def propiedades(seq: str, c_amida: bool = False, n_acetil: bool = False,
                ciclico: bool = False) -> dict:
    """Descriptores de la secuencia, sin construir la molécula 3D."""
    seq = seq.upper()
    return {
        "secuencia": seq,
        "nombre": etiqueta(seq, n_acetil=n_acetil, c_amida=c_amida, ciclico=ciclico),
        "longitud": len(seq),
        "carga_neta": carga_neta(seq, c_amida=c_amida, n_acetil=n_acetil, ciclico=ciclico),
        "gravy": gravy(seq),
        "momento_hidrofobico": momento_hidrofobico(seq),
        "fraccion_hidrofobica": fraccion_hidrofobica(seq),
        "indice_boman": indice_boman(seq),
    }


# ---------------------------------------------------------------- estructura
def to_smiles(seq: str, n_acetil: bool = False, c_amida: bool = False,
              ciclico: bool = False) -> Optional[str]:
    """Secuencia -> SMILES, con los aminoácidos en configuración L (la natural).

    n_acetil / c_amida protegen los extremos, habitual en péptidos antimicrobianos: amidar el
    carboxilo quita una carga negativa y sube la carga neta positiva, lo que favorece la interacción
    con la membrana bacteriana. ciclico cierra cabeza-cola, rigidizando y resistiendo proteasas.
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    seq = "".join(a for a in seq.upper() if a in AMINOACIDOS)
    if not seq:
        return None
    m = Chem.MolFromSequence(seq)
    if m is None:
        return None
    n_idx, c_idx, o_idx = _extremos(m)

    if ciclico:
        # Cierre cabeza-cola: se enlaza el N terminal con el C del carboxilo y se elimina su OH
        # (condensación, se pierde una agua). El enlace va ANTES de borrar el oxígeno porque
        # RemoveAtom reindexa; la valencia intermedia no es válida, de ahí sanitizar solo al final.
        if None in (n_idx, c_idx, o_idx):
            return None
        rw = Chem.RWMol(m)
        rw.AddBond(n_idx, c_idx, Chem.BondType.SINGLE)
        rw.RemoveAtom(o_idx)
        try:
            cerrado = rw.GetMol()
            Chem.SanitizeMol(cerrado)
            return Chem.MolToSmiles(cerrado)
        except Exception:
            return None

    # Los extremos se modifican sobre la misma molécula, sin pasar por SMILES: cada ida y vuelta
    # pierde la información de residuo y reordena los átomos, invalidando los índices localizados. Se
    # amida antes de acetilar porque ReplaceAtom no reindexa y AddAtom solo añade al final.
    rw = Chem.RWMol(m)
    if c_amida:
        if o_idx is None:
            return None
        rw.ReplaceAtom(o_idx, Chem.Atom(7))      # el OH del carboxilo pasa a NH2
    if n_acetil:
        if n_idx is None:
            return None
        c = rw.AddAtom(Chem.Atom(6)); o = rw.AddAtom(Chem.Atom(8)); me = rw.AddAtom(Chem.Atom(6))
        rw.AddBond(n_idx, c, Chem.BondType.SINGLE)
        rw.AddBond(c, o, Chem.BondType.DOUBLE)
        rw.AddBond(c, me, Chem.BondType.SINGLE)
    try:
        out = rw.GetMol()
        Chem.SanitizeMol(out)
        return Chem.MolToSmiles(out)
    except Exception:
        return None


def _extremos(mol) -> tuple:
    """(N del extremo N, C del carboxilo terminal, O de su OH), por información de residuo.

    Identificarlos por su forma química falla en silencio: la lisina aporta una amina primaria libre
    y el aspartato/glutamato un ácido carboxílico, así que quedarse con el primero que encaja cierra
    cualquier péptido con D o E por la cadena lateral —un lactama con otro anillo, carga y forma, pero
    con el nombre de ciclo cabeza-cola. MolFromSequence sí conserva residuo y nombre de átomo, que
    son inequívocos: el esqueleto se llama N y C en todos los aminoácidos.
    """
    info = [(a.GetIdx(), a.GetPDBResidueInfo()) for a in mol.GetAtoms()]
    nums = [i.GetResidueNumber() for _, i in info if i is not None]
    if not nums:
        return None, None, None
    primero, ultimo = min(nums), max(nums)
    n_idx = c_idx = o_idx = None
    for idx, i in info:
        if i is None:
            continue
        nombre = i.GetName().strip()
        if i.GetResidueNumber() == primero and nombre == "N":
            n_idx = idx
        elif i.GetResidueNumber() == ultimo and nombre == "C":
            c_idx = idx
    if c_idx is not None:
        o_idx = next((v.GetIdx() for v in mol.GetAtomWithIdx(c_idx).GetNeighbors()
                      if v.GetSymbol() == "O" and v.GetTotalNumHs() == 1), None)
    return n_idx, c_idx, o_idx


def viabilidad_docking(longitud: int, n_peptidos: int = 1, hay_adcp: bool = False) -> tuple:
    """(nivel, mensaje) sobre el coste y la fiabilidad de acoplar péptidos de esa longitud.

    Los umbrales NO son teoricos. Con AutoDock Vina, medido sobre saFtsZ (caja de 23 A,
    exhaustividad 8, un hilo), enlaces rotables y tiempo por péptido:
        3 residuos -> 15 rotables ->  ~98 s
        5 residuos -> 23 rotables -> más de 2 min
        8 residuos -> 39 rotables -> más de 2 min
    El coste crece con los grados de libertad, y la fiabilidad cae por el mismo motivo: Vina trata
    el ligando como un arbol de torsiones independientes, sin termino de energía conformacional, así
    que con muchas torsiones el muestreo deja de cubrir el espacio y la pose pierde sentido.

    ADCP no comparte esa limitacion, porque muestrea la conformación con rotámeros en vez de
    enumerar torsiones: sobre 8HTB, un octapeptido con 250.000 pasos x 10 replicas tarda 35 s con
    seis hilos. A cambio necesita al menos cinco residuos. Los dos programas se reparten por tanto
    el intervalo, y cual este disponible cambia el consejo.
    """
    if hay_adcp and longitud >= 5:
        minutos = max(1, round(n_peptidos * 0.6))
        coste = (f" Estimate for {n_peptidos} peptides: on the order of {minutos} minutes."
                 if n_peptidos > 1 else "")
        return "bueno", ("Length within ADCP's range, which generates the conformation during "
                         "docking instead of starting from a rigid structure. Raise steps and "
                         "replicas if the energy still improves." + coste)

    minutos = n_peptidos * (1.6 if longitud <= 3 else 3.0 if longitud <= 6 else 6.0)
    coste = (f" Estimate for {n_peptidos} peptides: on the order of {minutos:.0f} minutes on one "
             f"thread; reduce the time by raising «dockings in parallel»." if n_peptidos > 1 else "")
    if longitud <= 3:
        return "bueno", ("Flexibility comparable to a small molecule. Even so, each docking "
                         "costs about a minute and a half." + coste)
    if longitud <= 6:
        return "medio", ("High flexibility for Vina: docking is slow and the exact pose is "
                         "unreliable, though the relative order stays informative. Tighten the box "
                         "to the site, consider cyclizing the peptide to rigidify it, or install ADCP, "
                         "which covers this length without that problem." + coste)
    return "malo", ("Above 6 residues, rigid docking with Vina stops being practical: "
                    "the number of torsions exceeds what the algorithm samples reasonably. "
                    "This is exactly the length ADCP exists for; install it with "
                    "scripts/get_adcp.sh instead of forcing Vina here." + coste)
