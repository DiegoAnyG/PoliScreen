"""Core tests: chemistry, geometry and routing. They need no external binaries.

Each test corresponds to a failure found in real use; the comment says which, because a test without
that context ends up deleted when it gets in the way.
"""
from pathlib import Path

import pytest

from poliscreen.core import docking as dk
from poliscreen.core import peptides as pp
from poliscreen.core import session as ss


# --------------------------------------------------------------- hand-typed paths
@pytest.mark.parametrize("entry,expected", [
    (r"\\wsl.localhost\Ubuntu-24.04\home\u\proy", "/home/u/proy"),
    (r"\\wsl$\Ubuntu\home\u\proy", "/home/u/proy"),
    ("/home/u/proy", "/home/u/proy"),
    ('"/home/u/proy"', "/home/u/proy"),
])
def test_windows_paths_are_translated(entry, expected):
    """A Windows path created ONE folder with backslashes in the name."""
    assert str(ss.normalize_path(entry)[0]) == expected


def test_an_empty_path_gives_a_valid_destination():
    p, _ = ss.normalize_path("")
    assert p.is_absolute()


def test_the_default_root_follows_the_container_volume(monkeypatch):
    """In Docker the home directory is not persistent: results written there would be lost
    with the container, so the image points this at its mounted volume."""
    monkeypatch.setenv("POLISCREEN_PROJECTS", "/data")
    assert ss.default_root() == Path("/data")
    # Which folder inside it belongs to test_project_root.py; here only that it is inside.
    assert ss.normalize_path("")[0].parent == Path("/data")
    monkeypatch.delenv("POLISCREEN_PROJECTS")
    assert ss.default_root().parent == Path.home()


# --------------------------------------------------------------- peptide chemistry
def test_cyclization_is_head_to_tail_with_reactive_side_chains():
    """Looking for the COOH by its shape closed through aspartate or glutamate."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, rdMolDescriptors as rmd
    RDLogger.DisableLog("rdApp.*")
    for seq in ("SDCEFGQ", "ELQGRAK", "KWKLF", "ICIWDDS"):
        lin = Chem.MolFromSmiles(pp.to_smiles(seq))
        cic = Chem.MolFromSmiles(pp.to_smiles(seq, cyclic=True))
        macro = max((len(r) for r in cic.GetRingInfo().AtomRings() if len(r) > 8), default=0)
        assert macro == 3 * len(seq), f"{seq}: ring of {macro}, expected {3 * len(seq)}"
        loss = Descriptors.MolWt(lin) - Descriptors.MolWt(cic)
        assert abs(loss - 18.02) < 0.1, f"{seq}: lost {loss:.2f}, expected one water"
        assert rmd.CalcNumRings(cic) == rmd.CalcNumRings(lin) + 1


def test_only_the_termini_are_modified():
    """Acetylating touched the lysine and amidating the aspartate."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    AMINE = Chem.MolFromSmarts("[NX3;H2;!$(N-C=O)]")
    ACID = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")

    def count_of(smi, patron):
        return len(Chem.MolFromSmiles(smi).GetSubstructMatches(patron))

    assert count_of(pp.to_smiles("KKKGG"), AMINE) - count_of(pp.to_smiles("KKKGG", n_acetil=True), AMINE) == 1
    assert count_of(pp.to_smiles("DDEGG"), ACID) - count_of(pp.to_smiles("DDEGG", c_amida=True), ACID) == 1


def test_charge_reflects_the_termini_chemistry():
    """The table ignored acetylation and cyclization when computing the charge."""
    base = pp.properties("KWKLF")["net_charge"]
    assert pp.properties("KWKLF", n_acetil=True)["net_charge"] == base - 1
    assert pp.properties("KWKLF", c_amida=True)["net_charge"] == base + 1
    assert pp.properties("KWKLF", cyclic=True)["net_charge"] == base


def test_names_tell_the_variants_apart():
    """Two different molecules shared a name in the tables and in the pose files."""
    names_ = {pp.label_for("KWKLF", **kw) for kw in
               ({}, {"n_acetil": True}, {"c_amida": True},
                {"n_acetil": True, "c_amida": True}, {"cyclic": True})}
    assert len(names_) == 5


def test_smiles_fails_instead_of_returning_the_linear_peptide():
    """Returning the open peptide with a cyclic name hid the error."""
    assert pp.to_smiles("") is None


# --------------------------------------------------------------- resource allocation
def test_parallelism_drops_with_flexible_ligands(tmp_path):
    """Splitting only by cores exhausted the memory and the system killed the process."""
    box_ = dk.Box(0, 0, 0, 24, 24, 24)
    assert dk.memory_cost_gb(box_, 25) > dk.memory_cost_gb(box_, 5)
    big_box = dk.Box(0, 0, 0, 60, 60, 60)
    assert dk.memory_cost_gb(big_box, 5) > dk.memory_cost_gb(box_, 5)
    assert dk.safe_parallelism([box_], []) >= 1


def test_parallelism_leaves_the_machine_room(monkeypatch):
    """Vina takes whatever it is given. Filling every core throttles a laptop, and a throttled
    run is slower than fewer jobs at full clock."""
    monkeypatch.setattr(dk.os, "cpu_count", lambda: 16)
    monkeypatch.setattr(dk, "available_memory_gb", lambda: 64.0)
    assert dk.safe_parallelism([dk.Box(0, 0, 0, 20, 20, 20)], []) == 8


def test_torsdof_is_read_from_the_pdbqt(tmp_path):
    f = tmp_path / "l.pdbqt"
    f.write_text("ATOM      1  C   LIG A   1       0.0   0.0   0.0\nTORSDOF 22\n")
    assert dk.torsdof(f) == 22
    assert dk.torsdof(tmp_path / "no_existe.pdbqt") == 0


# --------------------------------------------------------------- recognize peptides
def test_a_small_molecule_is_not_taken_for_a_peptide(tmp_path):
    """Recognition decides the engine: a false positive would send a ligand to the wrong engine."""
    f = tmp_path / "lig.pdb"
    f.write_text("HETATM    1  C1  LIG A   1       0.0   0.0   0.0  1.00  0.00           C\nEND\n")
    assert pp.sequence_from_structure(f) is None
    assert pp.sequence_from_structure(tmp_path / "x.sdf") is None


def test_a_peptide_with_a_full_backbone_is_recognized(tmp_path):
    f = tmp_path / "pep.pdb"
    rows_ = []
    n = 1
    for i, (res, ) in enumerate([("ALA",), ("GLY",), ("LYS",), ("TRP",), ("PHE",)], start=1):
        for at in ("N", "CA", "C"):
            rows_.append("ATOM  %5d  %-3s %3s B%4d    %8.3f%8.3f%8.3f  1.00  0.00"
                         % (n, at, res, i, i * 3.0, 0.0, 0.0))
            n += 1
    f.write_text("\n".join(rows_) + "\nEND\n")
    seq, ciclo = pp.sequence_from_structure(f)
    assert seq == "AGKWF"
    assert ciclo is False


# --------------------------------------------------------------- export
def test_the_package_is_built_in_memory(tmp_path):
    """The export left loose copies in the results folder."""
    (tmp_path / "ranking.csv").write_text("compound,best_dock\na,-8.0\n")
    datos, included_items = ss.package_bytes(tmp_path, ["results_csv"], methods_text="# m")
    assert datos[:2] == b"PK"
    assert any("ranking.csv" in i for i in included_items)
    assert not list(tmp_path.glob("export_*"))


def test_the_catalog_marks_what_is_recommended(tmp_path):
    (tmp_path / "ranking.csv").write_text("a,b\n1,2\n")
    cat = ss.catalog(tmp_path)
    assert cat["results_csv"]["has"] and cat["results_csv"]["reason"]
    assert not cat["poses_zip"]["has"]


# --------------------------------------------------------------- modified residues
def test_preparation_does_not_duplicate_atoms(tmp_path):
    """A modified residue requested as a cofactor was added on top of the one already in the chain."""
    from poliscreen.core import receptor as rc
    f = tmp_path / "r.pdb"
    f.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       1.500   0.000   0.000  1.00  0.00           C\n"
        "END\n")
    assert rc._remove_overlapping(f) == 1
    assert rc._remove_overlapping(f) == 0        # idempotent


def test_conect_only_when_the_ring_closes(tmp_path):
    """The viewer drew a cyclic peptide open; it is closed only if the geometry justifies it."""
    from poliscreen.core import adcp

    def pose(dist):
        return (f"ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
                f"ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00  0.00           C\n"
                f"ATOM      3  C   ALA A   2       {dist:.3f}   0.000   0.000  1.00  0.00           C\n"
                f"END\n")
    closed_pose = tmp_path / "cerrada.pdb"
    closed_pose.write_text(pose(1.5))                 # N(res1)-C(res2) at 1.5 A
    d = adcp._close_ring(closed_pose)
    assert d is not None and "CONECT" in closed_pose.read_text()

    abierta = tmp_path / "abierta.pdb"
    abierta.write_text(pose(4.0))                 # at 4 A: no bond
    adcp._close_ring(abierta)
    assert "CONECT" not in abierta.read_text()


def test_it_explains_why_there_is_no_rmsd(tmp_path):
    """'not computable' said nothing; the usual cause is comparing different molecules."""
    from poliscreen.core import validation as vl
    a, b = tmp_path / "a.pdb", tmp_path / "b.pdb"
    a.write_text("ATOM      1  C1  LIG A   1       0.0   0.0   0.0  1.00  0.00           C\n"
                 "ATOM      2  P1  LIG A   1       1.0   0.0   0.0  1.00  0.00           P\n")
    b.write_text("ATOM      1  C1  LIG A   1       0.0   0.0   0.0  1.00  0.00           C\n")
    assert "are not the same molecule" in vl._reason_no_rmsd(a, b)
    assert vl._reason_no_rmsd(a, a) == "not computable"


def test_a_long_peptide_reaches_a_3d_structure(tmp_path):
    """Two 13-mers from a library of ten disappeared with no message: the embedding failed."""
    from poliscreen.core import ligands as lg
    seqs = ["DHITYAVHVQIRW", "WMHSPRFKIVVKW"]
    smis = [pp.to_smiles(s, c_amida=True) for s in seqs]
    assert all(smis)
    made = lg.materialize(smis, tmp_path, names=seqs)
    assert len(made) == len(seqs), "a ligand was lost while generating the 3D structure"


def test_the_gnina_wrapper_is_preferred(tmp_path, monkeypatch):
    """The standalone gnina binary is not self-contained: choosing it forced exporting a variable."""
    tools = tmp_path / "poliscreen_tools"
    tools.mkdir()
    for name_ in ("gnina", "gnina-run"):
        f = tools / name_
        f.write_text("#!/bin/sh\n")
        f.chmod(0o755)
    monkeypatch.delenv("POLISCREEN_GNINA", raising=False)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    assert dk.gnina_exe().endswith("gnina-run")


def test_modified_residues_are_detected_without_our_own_lists(tmp_path):
    """Detection uses PDBFixer's table; with no valid structure it returns an empty list, does not fail."""
    from poliscreen.core import receptor as rc
    f = tmp_path / "vacio.pdb"
    f.write_text("END\n")
    assert rc.modified_residues(f) == []
