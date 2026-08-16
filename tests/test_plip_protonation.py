"""Who places the hydrogens decides which hydrogen bonds exist.

PLIP is built on Open Babel and, left alone, protonates the complex itself before profiling it.
That put the hydrogens that decide every hydrogen bond in the hands of the one library whose
builds differ between the machines being compared -- which is why two of them still disagreed on
the anchor residues after the ligand chemistry had been fixed. Measured on a real complex,
--nohydro lost nothing and found one more bond, so the flag is used whenever the hydrogens are
already there, and only then.
"""
from pathlib import Path

from poliscreen.core import interactions as ix

WITH_H = """HETATM    1  C1  LIG A 401      -1.000   0.000   0.000  1.00  0.00           C
HETATM    2  H1  LIG A 401      -1.500   0.500   0.000  1.00  0.00           H
END
"""
WITHOUT_H = """HETATM    1  C1  LIG A 401      -1.000   0.000   0.000  1.00  0.00           C
HETATM    2  O1  LIG A 401      -1.500   0.500   0.000  1.00  0.00           O
END
"""
NO_ELEMENT_COLUMN = """HETATM    1  C1  LIG A 401      -1.000   0.000   0.000  1.00  0.00
HETATM    2  HG  LIG A 401      -1.500   0.500   0.000  1.00  0.00
END
"""
PROTEIN_ONLY = """ATOM      1  N   HIS A  10     -15.185  -6.358  34.779  1.00  0.00           N
ATOM      2  H   HIS A  10     -15.449  -5.830  33.734  1.00  0.00           H
END
"""


def _pdb(tmp_path, text, name="c.pdb"):
    p = tmp_path / name
    p.write_text(text)
    return p


def test_a_ligand_that_has_its_hydrogens_is_recognised(tmp_path):
    assert ix.has_ligand_hydrogens(_pdb(tmp_path, WITH_H))


def test_a_ligand_stripped_of_hydrogens_still_gets_them_from_plip(tmp_path):
    """Not adding them there would silently drop every hydrogen bond the ligand makes."""
    assert not ix.has_ligand_hydrogens(_pdb(tmp_path, WITHOUT_H))


def test_the_element_column_may_be_missing(tmp_path):
    """Older writers leave columns 77-78 blank; the atom name is what is left to read."""
    assert ix.has_ligand_hydrogens(_pdb(tmp_path, NO_ELEMENT_COLUMN))


def test_the_receptor_hydrogens_do_not_count(tmp_path):
    """A protonated receptor says nothing about the ligand, and the ligand is what is profiled."""
    assert not ix.has_ligand_hydrogens(_pdb(tmp_path, PROTEIN_ONLY))
