"""Deterministic protonation: the pieces the portability fix would be built from.

PLIP protonates the complex with openbabel, which is compiled per platform, and that is where an
18.6 % disagreement between two machines came from. These two functions move the decision: the
receptor's network is rebuilt by PDB2PQR (noarch, the same file everywhere) and the ligand's
hydrogens come from the SMILES it was built from. Whether that trade is worth making is being
measured by scripts/compare_protonation.py; what is tested here is that the pieces behave.
"""
import os
from pathlib import Path

import pytest

from poliscreen.core import ligands as lg
from poliscreen.core import receptor as rc

Chem = pytest.importorskip("rdkit.Chem")

# A pose as Vina returns it: heavy atoms only, no bond orders, nothing to rebuild hydrogens from.
ETHANOL_POSE = """HETATM    1  C1  LIG A 900       0.000   0.000   0.000  1.00  0.00           C
HETATM    2  C2  LIG A 900       1.520   0.000   0.000  1.00  0.00           C
HETATM    3  O1  LIG A 900       2.000   1.330   0.000  1.00  0.00           O
END
"""


def _pose(tmp_path):
    p = tmp_path / "pose.pdb"
    p.write_text(ETHANOL_POSE)
    return p


def test_a_pose_gets_its_hydrogens_from_its_recorded_chemistry(tmp_path):
    out = lg.with_hydrogens(_pose(tmp_path), "CCO", tmp_path / "lig.pdb")
    assert out is not None
    m = Chem.MolFromPDBFile(str(out), removeHs=False)
    assert sum(1 for a in m.GetAtoms() if a.GetSymbol() == "H") == 6, "ethanol carries six"


def test_a_template_that_is_not_the_molecule_says_so(tmp_path):
    """Silently returning the pose unprotonated would lose every hydrogen bond it makes."""
    assert lg.with_hydrogens(_pose(tmp_path), "c1ccccc1", tmp_path / "lig.pdb") is None


def test_no_template_is_not_a_crash(tmp_path):
    assert lg.with_hydrogens(_pose(tmp_path), "", tmp_path / "lig.pdb") is None


def test_the_executable_can_be_pointed_at(monkeypatch, tmp_path):
    """Same rule as every other external tool: no hardcoded path, an environment variable wins."""
    fake = tmp_path / "pdb2pqr30"
    fake.write_text("")
    monkeypatch.setenv("POLISCREEN_PDB2PQR", str(fake))
    assert rc.pdb2pqr_exe() == str(fake)


def test_a_missing_pdb2pqr_is_not_fatal(monkeypatch, tmp_path):
    """It improves the preparation; it is not a requirement for it."""
    monkeypatch.delenv("POLISCREEN_PDB2PQR", raising=False)
    monkeypatch.setattr(rc.shutil, "which", lambda _n: None)
    assert rc.pdb2pqr_exe() is None
    assert rc.optimize_hydrogens(tmp_path / "nothing.pdb", tmp_path / "out.pdb") is None
