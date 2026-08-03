"""Receptor preparation: PDB download, inspection and DockPrep-style cleanup.

Protonation uses PDBFixer because it keeps the author numbering, the one that appears in the
literature and that PLIP reports. PDBFixer strips everything non-protein, so the cofactors to keep are
extracted, protonated and reattached separately.
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
    """A hetero group: ligand, cofactor, ion or water."""
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
class ModifiedResidue:
    """A modified residue of a chain: phosphotyrosine, selenomethionine, phosphoserine...

    Even though the PDB declares it as a heteroatom, it is part of the chain and its modification is
    usually functional. It is distinguished to decide whether to keep it or replace it with the parent
    amino acid.
    """
    resname: str
    chain: str
    resseq: str
    standard: str

    @property
    def key(self) -> str:
        return f"{self.resname}|{self.chain}|{self.resseq}"

    @property
    def label(self) -> str:
        return f"{self.resname} {self.chain}:{self.resseq} (derives from {self.standard})"


def modified_residues(pdb) -> list:
    """Modified residues of the chains, detected with PDBFixer's table.

    Its own detection is used, not a separate list, so what is offered to the user matches what would
    happen without intervening. The PDB also declares them in MODRES, but that record is sometimes missing.
    """
    try:
        from pdbfixer import PDBFixer
    except Exception:
        return []
    try:
        f = PDBFixer(filename=str(pdb))
        f.findNonstandardResidues()
        return [ModifiedResidue(r.name, r.chain.id, str(r.id).strip(), new_)
                for r, new_ in f.nonstandardResidues]
    except Exception:
        return []


@dataclass
class Structure:
    path: Path
    chains: list = field(default_factory=list)
    het: list = field(default_factory=list)
    n_waters: int = 0
    n_atoms: int = 0

    def find(self, key: str) -> Optional[Het]:
        return next((h for h in self.het if h.key == key), None)

    def summary(self) -> str:
        lines = [f"{self.path.name}: {self.n_atoms} atoms, chains {','.join(self.chains) or '-'}, "
                 f"{self.n_waters} waters"]
        lines += [f"  {h.label}" for h in self.het]
        return "\n".join(lines)


def fetch_pdb(pdb_id: str, out_dir) -> Path:
    """Downloads a PDB by its identifier. Reuses the file if it already exists."""
    pdb_id = pdb_id.strip().upper()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / f"{pdb_id}.pdb"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    try:
        urllib.request.urlretrieve(RCSB_URL.format(pdb_id), dest)
    except Exception as e:
        raise ReceptorError(f"Could not download {pdb_id}: {e}. "
                            "Probable cause: invalid identifier or no connection.")
    return dest


def inspect(pdb) -> Structure:
    """Reads chains, hetero groups and waters without modifying anything."""
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
    """Extracts a hetero group as a standalone molecule, useful as a reference control.

    With `smiles` the bond orders are corrected from a template: the PDB does not store them and
    without them some ligands end up with impossible valences.
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
                pass
        # Keeps the largest fragment: extracting by geometry can split or duplicate the ligand.
        subprocess.run(["obabel", str(raw), "-O", str(out_path), "-p", str(ph), "-r"], capture_output=True)
    if not out_path.exists() or out_path.stat().st_size == 0:
        raise ReceptorError(f"Could not write the ligand {het.label}.")
    return out_path


def extract_chain(pdb, chain: str, out_path, polypeptide_only: bool = True,
                  on_notice=None) -> Path:
    """Extracts a chain as a reference control. Returns the written path.

    A co-crystallized peptide is a chain of the model, not a heteroatom, so it does not appear among
    the extractable cofactors. Without this route, a peptide screening would have to be compared with a
    small molecule —another engine, another chemotype— instead of with the crystal's own peptide.

    solo_polipeptido discards the chain's heteroatoms. A chain may carry conjugated groups (a
    nucleotide analogue linked by a connector, in a bisubstrate inhibitor) that the peptide engine does
    not reproduce; keeping them would make reference and pose different molecules and the validation
    RMSD would not exist. It is written in PDB, not SDF, to keep residue and atom names, which are what
    allow recognizing the sequence and routing the control to the peptide engine.
    """
    pdb, out_path = Path(pdb), Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    chain_lines = [l for l in pdb.read_text(errors="ignore").splitlines()
                    if l.startswith(("ATOM", "HETATM")) and l[21] == chain
                    and l[17:20].strip() not in WATERS]
    if not chain_lines:
        raise ReceptorError(f"Chain {chain} has no atoms in {pdb.name}.")

    lines_ = [l for l in chain_lines if l.startswith("ATOM")] if polypeptide_only else chain_lines
    if not lines_:
        lines_ = chain_lines
    outside = {l[17:20].strip() for l in chain_lines if l not in lines_}
    if outside and on_notice:
        on_notice(f"Chain {chain} carries non-peptide groups ({', '.join(sorted(outside))}) that have "
                 f"been excluded from the control: peptide docking reproduces only the peptide "
                 f"part, and keeping them would make the reference and the pose different "
                 f"molecules.")
    out_path.write_text("\n".join(lines_) + "\nEND\n")
    return out_path


def prepare(pdb, out_path, keep_chains: Optional[Sequence[str]] = None,
            keep_het: Sequence[str] = (), ph: float = 7.4, add_hydrogens: bool = True,
            keep_modified: Sequence[str] = (), on_notice=None,
            add_missing_residues: bool = False) -> Path:
    """Leaves the receptor ready to dock: no waters, with hydrogens and with whatever is kept.

    keep_chains           chains to keep (None = all)
    keep_het              hetero-group keys to keep, e.g. a cofactor
    add_missing_residues  rebuild missing loops; off by default, so as not to invent
                          geometry near the binding site
    """
    pdb, out_path = Path(pdb), Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        from openmm.app import PDBFile
        from pdbfixer import PDBFixer
    except Exception as e:
        raise ReceptorError(f"PDBFixer/OpenMM missing: {e}. Install them with conda to prepare receptors.")

    fixer = PDBFixer(filename=str(pdb))
    if keep_chains:
        keep = {c.strip() for c in keep_chains}
        drop = [ch.id for ch in fixer.topology.chains() if ch.id not in keep]
        if drop:
            fixer.removeChains(chainIds=drop)
    fixer.findMissingResidues()
    if not add_missing_residues:
        fixer.missingResidues = {}
    # PDBFixer replaces a modified residue with its standard amino acid, losing the modification.
    fixer.findNonstandardResidues()
    all_modified = list(getattr(fixer, "nonstandardResidues", []))

    def _clave(res):
        return f"{res.name}|{res.chain.id}|{str(res.id).strip()}"

    keep_keys = set(keep_modified or ())
    fixer.nonstandardResidues = [(r, n) for r, n in all_modified if _clave(r) not in keep_keys]
    # A modified residue is already in the chain: adding it as a cofactor would duplicate its atoms.
    replaced = set()
    for res, _n in all_modified:
        try:
            replaced.add((res.chain.id, str(res.id).strip()))
        except Exception:
            continue
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    if add_hydrogens:
        fixer.addMissingHydrogens(ph)
    with open(out_path, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)

    kept_modified = []
    for res, _n in all_modified:
        key_ = _clave(res)
        if key_ not in keep_keys:
            continue
        lines_ = [l for l in Path(pdb).read_text(errors="ignore").splitlines()
                  if l.startswith(("ATOM", "HETATM")) and l[21] == res.chain.id
                  and l[22:27].strip() == str(res.id).strip()
                  and l[17:20].strip() == res.name]
        if not lines_:
            continue
        body_ = "".join(l for l in out_path.read_text().splitlines(keepends=True)
                         if not l.strip().startswith(("END", "ENDMDL")))
        out_path.write_text(body_ + "".join(l + "\n" for l in lines_) + "END\n")
        kept_modified.append(f"{res.name} {res.chain.id}:{str(res.id).strip()}")

    skipped = []
    for key in keep_het:
        st = inspect(pdb)
        het = st.find(key)
        if het is None:
            raise ReceptorError(f"No existe el heterogrupo {key} en {pdb.name}.")
        if (het.chain, str(het.resseq).strip()) in replaced:
            skipped.append(het.label)
            continue
        # The cofactor keeps its original lines: protonating it renames the residue and PLIP would confuse it with the ligand.
        body = "".join(l for l in out_path.read_text().splitlines(keepends=True)
                       if not l.strip().startswith(("END", "ENDMDL")))
        extra = "".join(l + "\n" for l in _het_lines(pdb, het))
        out_path.write_text(body + extra + "END\n")

    n_dup = _remove_overlapping(out_path)
    if on_notice:
        for lab in skipped:
            on_notice(f"{lab} is a modified residue of the chain, not a cofactor: it is handled in its "
                     f"own section and not added separately.")
        replaced_res = [f"{r.name} {r.chain.id}:{str(r.id).strip()}"
                       for r, _n in all_modified if _clave(r) not in keep_keys]
        if replaced_res:
            on_notice("Replaced by their standard amino acid (the modification is lost): "
                     + ", ".join(replaced_res))
        if kept_modified:
            on_notice("Kept with their modification: " + ", ".join(kept_modified))
        if n_dup:
            on_notice(f"Removed {n_dup} overlapping atom(s) from the prepared receptor.")
    return out_path


def _remove_overlapping(pdb: Path, tol: float = 1e-3) -> int:
    """Removes atoms occupying the same position. Returns how many were removed.

    Two atoms at the same point make the angle calculation that assigns types and hybridization fail
    (ADFRsuite aborts with a division by zero) and double every contact. It is always pathological.
    """
    seen_items, output, removed = set(), [], 0
    for l in pdb.read_text(errors="ignore").splitlines(keepends=True):
        if l.startswith(("ATOM", "HETATM")):
            try:
                k = (round(float(l[30:38]) / tol), round(float(l[38:46]) / tol),
                     round(float(l[46:54]) / tol))
            except ValueError:
                output.append(l)
                continue
            if k in seen_items:
                removed += 1
                continue
            seen_items.add(k)
        output.append(l)
    if removed:
        pdb.write_text("".join(output))
    return removed
