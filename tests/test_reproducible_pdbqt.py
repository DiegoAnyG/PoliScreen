"""The receptor that reaches Vina must be the receptor that came out of preparation.

Seeding the hydrogen placement in receptor.prepare() made the prepared PDB reproducible, and two
machines then held byte-identical files -- and still ranked the compounds differently, because the
conversion to pdbqt threw that away: obabel's -p strips every hydrogen and re-adds them at random
positions. Those are the atoms Vina types HD and PLIP reads its hydrogen bonds off.
"""
import shutil

import pytest

from poliscreen.core import docking as dk

pytestmark = pytest.mark.skipif(shutil.which("obabel") is None, reason="needs OpenBabel")

# Three residues with the polar hydrogens obabel liked to move: a backbone amide, a hydroxyl and a
# charged amine. Written the way OpenMM writes them, hydrogens included.
PROTONATED = """\
ATOM      1  N   SER A   1      -1.234   0.456   2.789  1.00  0.00           N
ATOM      2  H   SER A   1      -1.567   1.234   3.210  1.00  0.00           H
ATOM      3  CA  SER A   1       0.123   0.789   2.345  1.00  0.00           C
ATOM      4  HA  SER A   1       0.456   1.567   3.012  1.00  0.00           H
ATOM      5  CB  SER A   1       1.234   0.123   1.567  1.00  0.00           C
ATOM      6  HB2 SER A   1       1.789   0.890   1.012  1.00  0.00           H
ATOM      7  HB3 SER A   1       0.890  -0.678   0.901  1.00  0.00           H
ATOM      8  OG  SER A   1       2.345  -0.234   2.456  1.00  0.00           O
ATOM      9  HG  SER A   1       3.012   0.345   2.789  1.00  0.00           H
ATOM     10  C   SER A   1       0.567   1.234   0.987  1.00  0.00           C
ATOM     11  O   SER A   1       0.123   2.345   0.678  1.00  0.00           O
ATOM     12  N   LYS A   2       1.456   0.678   0.234  1.00  0.00           N
ATOM     13  H   LYS A   2       1.789  -0.234   0.456  1.00  0.00           H
ATOM     14  CA  LYS A   2       1.901   1.123  -1.012  1.00  0.00           C
ATOM     15  HA  LYS A   2       1.234   1.890  -1.456  1.00  0.00           H
ATOM     16  CB  LYS A   2       3.234   1.789  -0.890  1.00  0.00           C
ATOM     17  HB2 LYS A   2       3.567   2.012  -1.901  1.00  0.00           H
ATOM     18  HB3 LYS A   2       3.123   2.678  -0.345  1.00  0.00           H
ATOM     19  C   LYS A   2       2.012   0.012  -2.012  1.00  0.00           C
ATOM     20  O   LYS A   2       2.789  -0.890  -1.789  1.00  0.00           O
ATOM     21  N   ALA A   3       1.234  -0.012  -3.123  1.00  0.00           N
ATOM     22  H   ALA A   3       0.678   0.789  -3.345  1.00  0.00           H
ATOM     23  CA  ALA A   3       1.123  -1.012  -4.178  1.00  0.00           C
ATOM     24  HA  ALA A   3       1.890  -1.789  -4.012  1.00  0.00           H
ATOM     25  CB  ALA A   3      -0.234  -1.678  -4.123  1.00  0.00           C
ATOM     26  HB1 ALA A   3      -0.345  -2.345  -4.987  1.00  0.00           H
ATOM     27  HB2 ALA A   3      -1.012  -0.901  -4.178  1.00  0.00           H
ATOM     28  HB3 ALA A   3      -0.345  -2.234  -3.178  1.00  0.00           H
ATOM     29  C   ALA A   3       1.345  -0.345  -5.567  1.00  0.00           C
ATOM     30  O   ALA A   3       2.234   0.456  -5.678  1.00  0.00           O
END
"""

STRIPPED = "\n".join(ln for ln in PROTONATED.splitlines()
                     if ln[76:78].strip() != "H") + "\n"


def _convert(src, dst):
    assert dk.to_pdbqt(src, dst, receptor=True), "obabel produced nothing"
    return dst.read_text()


def test_the_same_receptor_converts_to_the_same_pdbqt(tmp_path):
    src = tmp_path / "ready.pdb"
    src.write_text(PROTONATED)
    first = _convert(src, tmp_path / "a.pdbqt")
    for i in range(3):
        assert _convert(src, tmp_path / f"b{i}.pdbqt") == first


def test_the_prepared_hydrogens_are_the_ones_that_reach_vina(tmp_path):
    """Not merely reproducible: reproducing PDBFixer's placement, not obabel's guess at it."""
    src = tmp_path / "ready.pdb"
    src.write_text(PROTONATED)
    out = _convert(src, tmp_path / "a.pdbqt")
    assert " 3.012   0.345   2.789" in out, "the serine hydroxyl hydrogen was moved"


def test_a_receptor_without_hydrogens_still_gets_them(tmp_path):
    """Skipping the protonation must not leave a bare structure bare."""
    src = tmp_path / "bare.pdb"
    src.write_text(STRIPPED)
    assert not dk.has_hydrogens(src)
    out = _convert(src, tmp_path / "a.pdbqt")
    assert any(ln.rstrip().endswith(("HD", " H")) for ln in out.splitlines())
