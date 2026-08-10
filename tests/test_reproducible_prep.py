"""The prepared receptor has to come out the same twice.

OpenMM places new hydrogens at random positions and relaxes them with a short minimisation --
"The hydrogens were added at random positions", says its own source -- taking its numbers from
Python's global RNG. Unseeded, the same PDB produced a different receptor on every run: heavy
atoms identical to the last decimal, hydrogens up to 2.3 A apart. Nothing downstream can be
reproducible on top of that: the pdbqt changes, Vina lands elsewhere, PLIP reads its H-bonds off
those very hydrogens, and the ranking is normalised against a control that has also moved.

Seeding fixes the starting positions; pinning the Reference platform fixes the minimisation,
which on a multi-threaded platform sums in a different order every time and left ~0.1 A behind.
Both halves are needed, so this asks for the whole file back, byte for byte.
"""
from pathlib import Path

import pytest

pytest.importorskip("pdbfixer")
pytest.importorskip("openmm")

from poliscreen.core import receptor as rc  # noqa: E402

# Three residues with a hydroxyl and a charged amine: the rotatable hydrogens are the ones that
# used to land somewhere else on every run. Small enough that preparing it twice stays fast.
FRAGMENT = """\
ATOM      1  N   SER A   1      -0.518   1.363   0.000  1.00  0.00           N
ATOM      2  CA  SER A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   SER A   1       1.520   0.000   0.000  1.00  0.00           C
ATOM      4  O   SER A   1       2.147   1.055   0.000  1.00  0.00           O
ATOM      5  CB  SER A   1      -0.507  -0.751  -1.226  1.00  0.00           C
ATOM      6  OG  SER A   1      -1.923  -0.792  -1.258  1.00  0.00           O
ATOM      7  N   LYS A   2       2.117  -1.185   0.000  1.00  0.00           N
ATOM      8  CA  LYS A   2       3.570  -1.322   0.000  1.00  0.00           C
ATOM      9  C   LYS A   2       4.240  -0.011   0.386  1.00  0.00           C
ATOM     10  O   LYS A   2       3.606   0.919   0.886  1.00  0.00           O
ATOM     11  CB  LYS A   2       4.049  -1.788  -1.377  1.00  0.00           C
ATOM     12  CG  LYS A   2       5.566  -1.949  -1.446  1.00  0.00           C
ATOM     13  CD  LYS A   2       6.019  -2.437  -2.819  1.00  0.00           C
ATOM     14  CE  LYS A   2       7.535  -2.594  -2.887  1.00  0.00           C
ATOM     15  NZ  LYS A   2       7.984  -3.070  -4.220  1.00  0.00           N
ATOM     16  N   ALA A   3       5.549   0.055   0.155  1.00  0.00           N
ATOM     17  CA  ALA A   3       6.317   1.257   0.470  1.00  0.00           C
ATOM     18  C   ALA A   3       7.798   1.019   0.199  1.00  0.00           C
ATOM     19  O   ALA A   3       8.180  -0.056  -0.271  1.00  0.00           O
ATOM     20  CB  ALA A   3       5.821   2.440  -0.354  1.00  0.00           C
TER
END
"""


def _prepared(tmp_path: Path, tag: str, seed: int = 42) -> str:
    src = tmp_path / "fragment.pdb"
    src.write_text(FRAGMENT)
    out = tmp_path / f"ready_{tag}.pdb"
    rc.prepare(src, out, ph=7.4, seed=seed)
    # The REMARK carries the date OpenMM wrote it, which is not part of the structure.
    return "\n".join(l for l in out.read_text().splitlines() if not l.startswith("REMARK"))


def test_the_same_receptor_twice_is_the_same_file(tmp_path):
    assert _prepared(tmp_path, "a") == _prepared(tmp_path, "b"), \
        "the receptor changes between runs: nothing downstream of it can be reproducible"


def test_the_seed_is_what_decides_it(tmp_path):
    """If a different seed changed nothing, the seeding above would be decorative."""
    assert _prepared(tmp_path, "c", seed=42) != _prepared(tmp_path, "d", seed=7)
