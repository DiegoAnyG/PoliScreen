"""Combinatorial peptide generation and their properties.

Second design route, independent of reaction-based synthesis: here the ligands are not built by
joining a core to a reagent, but by enumerating amino-acid sequences under rules. Aimed at
antimicrobial peptides, where activity depends mostly on positive net charge and amphipathicity,
not on a point pharmacophore.

Generation is deterministic: the same seed and rules produce the same library.
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

# Longer than this it is a protein chain, not a design peptide.
MAX_CHAIN_LENGTH = 60

MIN_LENGTH = 1
MAX_LENGTH = 20

# Amino acid: name, classes, charge at pH 7.4 and Kyte-Doolittle hydropathy.
AMINO_ACIDS = {
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

CLASSES = {
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

_EISENBERG = {
    "A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29, "E": -0.74, "Q": -0.85,
    "G": 0.48, "H": -0.40, "I": 1.38, "L": 1.06, "K": -1.50, "M": 0.64, "F": 1.19,
    "P": 0.12, "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08,
}
# Boman index: protein-binding energy; high means promiscuous.
_BOMAN = {
    "A": -0.5, "R": 2.5, "N": 0.2, "D": 3.0, "C": -1.0, "E": 3.0, "Q": 0.2, "G": 0.0,
    "H": -0.5, "I": -1.8, "L": -1.8, "K": 3.0, "M": -1.3, "F": -2.5, "P": 0.0,
    "S": 0.3, "T": -0.4, "W": -3.4, "Y": -2.3, "V": -1.5,
}

HELIX_ANGLE = 100.0


def alphabet(include: Sequence[str] = (), exclude_classes: Sequence[str] = (),
             exclude_residues: Sequence[str] = ()) -> list:
    """Allowed amino-acid alphabet. `include` are classes (union); empty = all 20."""
    if include:
        sel = {a for a, (_n, cls, _c, _h) in AMINO_ACIDS.items() if set(cls) & set(include)}
    else:
        sel = set(AMINO_ACIDS)
    if exclude_classes:
        sel -= {a for a, (_n, cls, _c, _h) in AMINO_ACIDS.items() if set(cls) & set(exclude_classes)}
    sel -= set(x.upper() for x in exclude_residues)
    return sorted(sel)


@dataclass
class Rules:
    """Library constraints. All are optional and combine."""
    length_: int = 7
    alphabet: Sequence[str] = field(default_factory=lambda: sorted(AMINO_ACIDS))
    no_repeats: bool = False
    max_consecutive: int = 0
    max_per_residue: int = 0
    prefix_: str = ""
    suffix_: str = ""
    charge_min: Optional[float] = None
    charge_max: Optional[float] = None
    gravy_min: Optional[float] = None
    gravy_max: Optional[float] = None

    def validate(self) -> list:
        """Readable warnings; empty list = consistent rules."""
        notices = []
        if not (MIN_LENGTH <= self.length_ <= MAX_LENGTH):
            notices.append(f"Length must be between {MIN_LENGTH} and {MAX_LENGTH}.")
        if not self.alphabet:
            notices.append("The alphabet is empty: no selected class leaves any amino acids.")
        fijos = len(self.prefix_) + len(self.suffix_)
        if fijos > self.length_:
            notices.append(f"Prefix and suffix add up to {fijos} residues and the length is {self.length_}.")
        if self.no_repeats and self.length_ > len(self.alphabet):
            notices.append(f"Without repeats, {self.length_} residues do not fit in an alphabet of "
                          f"{len(self.alphabet)}.")
        for patron, name_ in ((self.prefix_, "prefix"), (self.suffix_, "suffix")):
            outside = set(patron.upper()) - set(AMINO_ACIDS)
            if outside:
                notices.append(f"The {name_} contains symbols that are not amino acids: {', '.join(sorted(outside))}.")
        termini = (self.prefix_ + self.suffix_).upper()
        if self.no_repeats and len(set(termini)) != len(termini):
            notices.append("«No repeats» conflicts with the prefix or suffix, which already repeat a residue.")
        if self.max_consecutive and self.max_consecutive > 0:
            streak, previous = 1, ""
            for a in termini:
                streak = streak + 1 if a == previous else 1
                previous = a
                if streak > self.max_consecutive:
                    notices.append(f"The prefix or suffix repeats {streak} residues in a row and the maximum "
                                  f"consecutive is {self.max_consecutive}.")
                    break
        outside_alphabet = set(termini) - set(self.alphabet)
        if outside_alphabet:
            notices.append(f"The prefix or suffix uses residues outside the chosen alphabet: "
                          f"{', '.join(sorted(outside_alphabet))}.")
        return notices

    def space(self) -> float:
        """Upper bound on the number of possible sequences. Used to warn when more library is
        requested than the combinatorial space allows."""
        free_slots = max(0, self.length_ - len(self.prefix_) - len(self.suffix_))
        n = len(self.alphabet)
        if free_slots == 0:
            return 1.0
        if self.no_repeats:
            total, disp = 1.0, n
            for _ in range(free_slots):
                total *= max(disp, 0); disp -= 1
            return total
        return float(n) ** free_slots


def _cumple(seq: str, r: Rules) -> bool:
    if r.no_repeats and len(set(seq)) != len(seq):
        return False
    if r.max_per_residue:
        for a in set(seq):
            if seq.count(a) > r.max_per_residue:
                return False
    if r.max_consecutive:
        streak, previous = 1, ""
        for a in seq:
            streak = streak + 1 if a == previous else 1
            if streak > r.max_consecutive:
                return False
            previous = a
    if r.charge_min is not None or r.charge_max is not None:
        q = carga_neta(seq)
        if r.charge_min is not None and q < r.charge_min:
            return False
        if r.charge_max is not None and q > r.charge_max:
            return False
    if r.gravy_min is not None or r.gravy_max is not None:
        g = gravy(seq)
        if r.gravy_min is not None and g < r.gravy_min:
            return False
        if r.gravy_max is not None and g > r.gravy_max:
            return False
    return True


def generate(rules: Rules, n: int, seed: int = 42, max_attempts: int = 200) -> tuple:
    """Generates up to `n` unique sequences that satisfy the rules. Returns (sequences, notice).

    Random sampling with a fixed seed instead of exhaustive enumeration: the space grows as
    20^length and enumerating it is infeasible except for very short peptides. It stops when it has
    n sequences or when it stops finding new ones, so asking for more than exist does not hang the app.
    """
    rng = random.Random(seed)
    pre, suf = rules.prefix_.upper(), rules.suffix_.upper()
    free_slots = rules.length_ - len(pre) - len(suf)
    if free_slots < 0 or not rules.alphabet:
        return [], "The rules do not allow building any sequence."

    seen_, output, failures = set(), [], 0
    limite = max_attempts * max(n, 1)
    attempts = 0
    while len(output) < n and attempts < limite and failures < max_attempts * 20:
        attempts += 1
        if rules.no_repeats:
            used_set = set(pre) | set(suf)
            free_pool = [a for a in rules.alphabet if a not in used_set]
            if len(free_pool) < free_slots:
                break
            mid = "".join(rng.sample(free_pool, free_slots))
        else:
            mid = "".join(rng.choice(rules.alphabet) for _ in range(free_slots))
        seq = pre + mid + suf
        if seq in seen_:
            failures += 1
            continue
        if not _cumple(seq, rules):
            failures += 1
            continue
        seen_.add(seq); output.append(seq); failures = 0

    notice = ""
    if len(output) < n:
        esp = rules.space()
        notice = (f"{len(output)} of {n} sequences generated: with these rules the available "
                 f"space is about {esp:.0f} and the filters discard the rest.")
    return output, notice


def carga_neta(seq: str, ph: float = 7.4, c_amida: bool = False,
               n_acetil: bool = False, cyclic: bool = False) -> float:
    """Approximate net charge at physiological pH. In antimicrobial peptides it is the descriptor
    most associated with activity: the bacterial membrane is anionic.

    Terminus chemistry enters here because it changes the charge, which is precisely the descriptor
    used to decide: amidating the carboxyl removes a negative charge (+1 net), acetylating the amino
    removes a positive one (-1 net) and head-to-tail closure consumes both termini, so only the
    side-chain charges remain.
    """
    q = sum(AMINO_ACIDS[a][2] for a in seq.upper() if a in AMINO_ACIDS)
    if cyclic:
        return round(q, 2)
    if not n_acetil:
        q += 1.0
    if not c_amida:
        q -= 1.0
    return round(q, 2)


def gravy(seq: str) -> float:
    """Mean hydropathy (Kyte-Doolittle). Positive = overall hydrophobic."""
    vals = [AMINO_ACIDS[a][3] for a in seq.upper() if a in AMINO_ACIDS]
    return round(sum(vals) / len(vals), 3) if vals else 0.0


def momento_hidrofobico(seq: str, angulo: float = HELIX_ANGLE) -> float:
    """Eisenberg normalized hydrophobic moment. Measures AMPHIPATHICITY: whether, folded into a
    helix, the hydrophobic residues sit on one face and the polar ones on the other. This is what
    lets a peptide insert into the membrane."""
    s = [a for a in seq.upper() if a in _EISENBERG]
    if not s:
        return 0.0
    rad = math.radians(angulo)
    sx = sum(_EISENBERG[a] * math.cos(i * rad) for i, a in enumerate(s))
    sy = sum(_EISENBERG[a] * math.sin(i * rad) for i, a in enumerate(s))
    return round(math.hypot(sx, sy) / len(s), 3)


def indice_boman(seq: str) -> float:
    """Potential to bind other proteins (kcal/mol). Above ~2.5 the peptide is considered promiscuous
    and may have more nonspecific effects."""
    vals = [_BOMAN[a] for a in seq.upper() if a in _BOMAN]
    return round(sum(vals) / len(vals), 2) if vals else 0.0


def fraccion_hidrofobica(seq: str) -> float:
    s = [a for a in seq.upper() if a in AMINO_ACIDS]
    if not s:
        return 0.0
    h = sum(1 for a in s if "hidrofobico" in AMINO_ACIDS[a][1])
    return round(h / len(s), 3)


_THREE_TO_ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
               "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
               "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
               "TYR": "Y", "VAL": "V"}


def sequence_from_structure(path_) -> Optional[tuple]:
    """(sequence, cyclic) if the file contains a polypeptide; None if it is not one.

    Lets a peptide be recognized without depending on any metadata: needed for the co-crystallized
    controls, which are extracted from the crystal and do not pass through the ligand table. A
    complete backbone (N, CA and C) is required at each residue so as not to mistake for a
    polypeptide any molecule containing amino-acid fragments.
    """
    p = Path(path_)
    if p.suffix.lower() not in (".pdb", ".pdbqt", ".ent"):
        return None
    res, esq = {}, {}
    try:
        for l in p.read_text(errors="ignore").splitlines():
            if not l.startswith(("ATOM", "HETATM")):
                continue
            nom3, key_ = l[17:20].strip().upper(), (l[21], l[22:27].strip())
            if nom3 not in _THREE_TO_ONE:
                continue
            res[key_] = nom3
            esq.setdefault(key_, set()).add(l[12:16].strip())
    except Exception:
        return None
    completos = [k for k, v in res.items() if {"N", "CA", "C"} <= esq.get(k, set())]
    if not 1 <= len(completos) <= MAX_CHAIN_LENGTH:
        return None
    order_ = sorted(completos, key=lambda k: (k[0], _num(k[1])))
    seq = "".join(_THREE_TO_ONE[res[k]] for k in order_)

    # Head-to-tail cyclization: amide between the first nitrogen and the last carbon.
    ciclo = False
    xyz = {}
    for l in p.read_text(errors="ignore").splitlines():
        if l.startswith(("ATOM", "HETATM")):
            k = (l[21], l[22:27].strip())
            if k in (order_[0], order_[-1]):
                xyz[(k, l[12:16].strip())] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    a, b = xyz.get((order_[0], "N")), xyz.get((order_[-1], "C"))
    if a and b:
        ciclo = sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5 < 1.8
    return seq, ciclo


def _num(s: str) -> int:
    d = "".join(c for c in s if c.isdigit() or c == "-")
    return int(d) if d and d != "-" else 0


def label_for(seq: str, n_acetil: bool = False, c_amida: bool = False,
             cyclic: bool = False) -> str:
    """Name that makes the terminus chemistry visible: 'Ac-KWKLF-NH2', 'cyclo-KWKLF'.

    Without it, two different molecules share a name in tables and pose files, and a screening with
    protected termini would be indistinguishable from one without protection. No characters that
    complicate the file name.
    """
    seq = seq.upper()
    if cyclic:
        return f"cyclo-{seq}"
    return ("Ac-" if n_acetil else "") + seq + ("-NH2" if c_amida else "")


def properties(seq: str, c_amida: bool = False, n_acetil: bool = False,
                cyclic: bool = False) -> dict:
    """Sequence descriptors, without building the 3D molecule."""
    seq = seq.upper()
    return {
        "sequence": seq,
        "name": label_for(seq, n_acetil=n_acetil, c_amida=c_amida, cyclic=cyclic),
        "length": len(seq),
        "net_charge": carga_neta(seq, c_amida=c_amida, n_acetil=n_acetil, cyclic=cyclic),
        "gravy": gravy(seq),
        "hydrophobic_moment": momento_hidrofobico(seq),
        "hydrophobic_fraction": fraccion_hidrofobica(seq),
        "boman_index": indice_boman(seq),
    }


def to_smiles(seq: str, n_acetil: bool = False, c_amida: bool = False,
              cyclic: bool = False) -> Optional[str]:
    """Sequence -> SMILES, with the amino acids in L configuration (the natural one).

    n_acetil / c_amida protect the termini, common in antimicrobial peptides: amidating the carboxyl
    removes a negative charge and raises the positive net charge, which favors interaction with the
    bacterial membrane. ciclico closes head-to-tail, rigidifying and resisting proteases.
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    seq = "".join(a for a in seq.upper() if a in AMINO_ACIDS)
    if not seq:
        return None
    m = Chem.MolFromSequence(seq)
    if m is None:
        return None
    n_idx, c_idx, o_idx = _extremos(m)

    if cyclic:
        # Condensation loses one water; the bond is made before deleting the oxygen so the indices stay valid.
        if None in (n_idx, c_idx, o_idx):
            return None
        rw = Chem.RWMol(m)
        rw.AddBond(n_idx, c_idx, Chem.BondType.SINGLE)
        rw.RemoveAtom(o_idx)
        try:
            closed = rw.GetMol()
            Chem.SanitizeMol(closed)
            return Chem.MolToSmiles(closed)
        except Exception:
            return None

    # The termini are modified on the molecule: a SMILES round trip loses the residue information.
    rw = Chem.RWMol(m)
    if c_amida:
        if o_idx is None:
            return None
        rw.ReplaceAtom(o_idx, Chem.Atom(7))
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
    """(N of the N-terminus, C of the terminal carboxyl, O of its OH), by residue information.

    Identifying them by their chemical form fails silently: lysine contributes a free primary amine
    and aspartate/glutamate a carboxylic acid, so keeping the first match closes any peptide with D or
    E through the side chain —a lactam with another ring, charge and shape, but named a head-to-tail
    cycle. MolFromSequence does keep residue and atom name, which are unambiguous: the backbone is
    called N and C in every amino acid.
    """
    info = [(a.GetIdx(), a.GetPDBResidueInfo()) for a in mol.GetAtoms()]
    nums = [i.GetResidueNumber() for _, i in info if i is not None]
    if not nums:
        return None, None, None
    first_, last_ = min(nums), max(nums)
    n_idx = c_idx = o_idx = None
    for idx, i in info:
        if i is None:
            continue
        name_ = i.GetName().strip()
        if i.GetResidueNumber() == first_ and name_ == "N":
            n_idx = idx
        elif i.GetResidueNumber() == last_ and name_ == "C":
            c_idx = idx
    if c_idx is not None:
        o_idx = next((v.GetIdx() for v in mol.GetAtomWithIdx(c_idx).GetNeighbors()
                      if v.GetSymbol() == "O" and v.GetTotalNumHs() == 1), None)
    return n_idx, c_idx, o_idx


def docking_feasibility(length_: int, n_peptides: int = 1, has_adcp: bool = False) -> tuple:
    """(level, message) about the cost and reliability of docking peptides of that length.

    The thresholds are NOT theoretical. With AutoDock Vina, measured on saFtsZ (23 A box,
    exhaustiveness 8, one thread), rotatable bonds and time per peptide:
        3 residues -> 15 rotatable ->  ~98 s
        5 residues -> 23 rotatable -> more than 2 min
        8 residues -> 39 rotatable -> more than 2 min
    Cost grows with the degrees of freedom, and reliability falls for the same reason: Vina treats
    the ligand as an independent torsion tree, with no conformational-energy term, so with many
    torsions the sampling stops covering the space and the pose loses meaning.

    ADCP does not share that limitation, because it samples the conformation with rotamers instead of
    enumerating torsions: on 8HTB, an octapeptide with 250,000 steps x 10 replicas takes 35 s on six
    threads. In exchange it needs at least five residues. The two programs thus split the range, and
    which one is available changes the advice.
    """
    if has_adcp and length_ >= 5:
        minutes = max(1, round(n_peptides * 0.6))
        coste = (f" Estimate for {n_peptides} peptides: on the order of {minutes} minutes."
                 if n_peptides > 1 else "")
        return "good", ("Length within ADCP's range, which generates the conformation during "
                         "docking instead of starting from a rigid structure. Raise steps and "
                         "replicas if the energy still improves." + coste)

    minutes = n_peptides * (1.6 if length_ <= 3 else 3.0 if length_ <= 6 else 6.0)
    coste = (f" Estimate for {n_peptides} peptides: on the order of {minutes:.0f} minutes on one "
             f"thread; reduce the time by raising «dockings in parallel»." if n_peptides > 1 else "")
    if length_ <= 3:
        return "good", ("Flexibility comparable to a small molecule. Even so, each docking "
                         "costs about a minute and a half." + coste)
    if length_ <= 6:
        return "mid", ("High flexibility for Vina: docking is slow and the exact pose is "
                         "unreliable, though the relative order stays informative. Tighten the box "
                         "to the site, consider cyclizing the peptide to rigidify it, or install ADCP, "
                         "which covers this length without that problem." + coste)
    return "bad", ("Above 6 residues, rigid docking with Vina stops being practical: "
                    "the number of torsions exceeds what the algorithm samples reasonably. "
                    "This is exactly the length ADCP exists for; install it with "
                    "scripts/get_adcp.sh instead of forcing Vina here." + coste)
