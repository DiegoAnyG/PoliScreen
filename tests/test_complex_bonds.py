"""A fused complex has to say the ligand's bonds without lying about the receptor's.

A pose numbers its atoms from 1, and so does the receptor. Fusing them wrote CONECT records that
named receptor atoms: every complex the project has ever written carries bonds that are not there,
which any viewer opening the file believes. The pipeline never saw it because it stripped CONECT
before PLIP, and stripping them is also what hid the corruption from the interaction analysis.
"""
from pathlib import Path

from poliscreen.core import interactions as ix

RECEPTOR = """ATOM      1  N   HIS A  10     -15.185  -6.358  34.779  1.00  0.00           N
ATOM      2  CA  HIS A  10     -14.609  -7.736  34.643  1.00  0.00           C
ATOM      3  C   HIS A  10     -13.100  -7.700  34.400  1.00  0.00           C
END
"""

POSE = """HETATM    1  C1  UNL A 900       0.000   0.000   0.000  1.00  0.00           C
HETATM    2  O1  UNL A 900       1.430   0.000   0.000  1.00  0.00           O
CONECT    1    2
CONECT    2    1
END
"""


def _fused(tmp_path):
    rec, pose = tmp_path / "rec.pdb", tmp_path / "pose.pdb"
    rec.write_text(RECEPTOR)
    pose.write_text(POSE)
    return ix.fuse(rec, pose, tmp_path / "complex.pdb")


def test_the_ligand_is_numbered_above_the_receptor(tmp_path):
    text = Path(_fused(tmp_path)).read_text()
    serials = [int(l[6:11]) for l in text.splitlines() if l[17:20].strip() == "LIG"]
    receptor = [int(l[6:11]) for l in text.splitlines() if l.startswith("ATOM")]
    assert min(serials) > max(receptor), "the ligand reuses receptor serial numbers"


def test_the_bonds_follow_the_renumbering(tmp_path):
    """A CONECT left pointing at the old serial is a bond drawn to a receptor atom."""
    text = Path(_fused(tmp_path)).read_text()
    conect = [l for l in text.splitlines() if l.startswith("CONECT")]
    assert conect, "the ligand's bonds were dropped"
    referenced = {int(l[i:i + 5]) for l in conect for i in range(6, len(l.rstrip()), 5)}
    receptor = {int(l[6:11]) for l in text.splitlines() if l.startswith("ATOM")}
    assert not (referenced & receptor), f"bonds point at receptor atoms: {referenced & receptor}"


def test_sanitizing_keeps_them(tmp_path):
    """The column rewriting exists to protect the residue numbering, not to discard bonds."""
    fused = _fused(tmp_path)
    out = tmp_path / "san.pdb"
    assert ix.sanitize_pdb(fused, out)
    assert "CONECT" in out.read_text()


def test_the_residue_numbering_still_survives_sanitizing(tmp_path):
    """What this function was written for: a misread column turns Asn5 into Asn6."""
    fused = _fused(tmp_path)
    out = tmp_path / "san.pdb"
    ix.sanitize_pdb(fused, out)
    numbers = {l[22:26].strip() for l in out.read_text().splitlines() if l.startswith("ATOM")}
    assert numbers == {"10"}
