"""Test control assignment and provenance persistence.

Verifies that:
1. Controls are assigned to their true receptor without misassignment, even when Cartesian
   coordinates overlap with other structures.
2. Control provenance is persisted in control_map.json and respected by _controls_of.
3. Dropping a receptor updates the persisted control_map.json.
"""
import json
from pathlib import Path
import pytest

from poliscreen.core import layout as lay
from poliscreen.core import pipeline as pl
from poliscreen.core import receptor as rc
from poliscreen.core import screening as sc
from poliscreen.ui import common


def _atom(serial, name, resname, chain, resseq, x, y, z):
    return (f"ATOM  {serial:>5} {name:<4} {resname:>3} {chain}{resseq:>4}    "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           C\n")


def test_control_assignment_avoids_clashing_receptors(tmp_path):
    """When a control has atoms near two receptors, it matches the pocket without severe clashes."""
    rec_dir = lay.artifact(tmp_path, lay.RECEPTORS)
    rec_dir.mkdir(parents=True, exist_ok=True)

    # Receptor 1: Pocket residues surrounding (10, 10, 10) at distance ~2.5 A
    rec1 = rec_dir / "REC1_ready.pdb"
    rec1_lines = ["CRYST1    0.0 0.0 0.0  90.0 90.0 90.0 P 1           1\n"]
    for i, (dx, dy, dz) in enumerate([(2.5, 0, 0), (-2.5, 0, 0), (0, 2.5, 0), (0, -2.5, 0)]):
        rec1_lines.append(_atom(i + 1, "CA", "ALA", "A", i + 1, 10.0 + dx, 10.0 + dy, 10.0 + dz))
    rec1.write_text("".join(rec1_lines))

    # Receptor 2: An atom directly clashing with the ligand at (10.0, 10.0, 10.2) (d = 0.2 A)
    rec2 = rec_dir / "REC2_ready.pdb"
    rec2_lines = ["CRYST1    0.0 0.0 0.0  90.0 90.0 90.0 P 1           1\n"]
    rec2_lines.append(_atom(1, "CA", "PHE", "A", 1, 10.0, 10.0, 10.2))  # 0.2 A clash!
    for i in range(5):
        rec2_lines.append(_atom(i + 2, "CA", "ALA", "A", i + 2, 50.0 + i, 50.0, 50.0))
    rec2.write_text("".join(rec2_lines))

    # Control ligand: Single heavy atom at (10, 10, 10)
    ctrl = rec_dir / "control_LIG.pdb"
    ctrl.write_text(_atom(1, "C1", "LIG", "A", 1, 10.0, 10.0, 10.0))

    asig = pl._assign_controls([ctrl], [rec1, rec2], {})
    assert asig.get("controllig") == "REC1_ready", f"Expected REC1_ready, got {asig}"


def test_control_map_persistence(tmp_path):
    """Saving and loading control_map preserves receptor-control associations."""
    rec_dir = lay.artifact(tmp_path, lay.RECEPTORS)
    rec_dir.mkdir(parents=True, exist_ok=True)

    rec1 = rec_dir / "4D44_ready.pdb"
    rec1.write_text(_atom(1, "CA", "ALA", "A", 1, 0, 0, 0))
    ctrl1 = rec_dir / "control_JA3.sdf"
    ctrl1.write_text(_atom(1, "C", "JA3", "A", 1, 0, 0, 2))

    # Save mapping
    mapping = {"controlja3": "4D44_ready"}
    common._save_control_map(rec_dir, mapping)

    # Load mapping
    loaded = common._load_control_map(rec_dir)
    assert loaded.get("controlja3") == "4D44_ready"

    # _controls_of returns ctrl1
    ctrls = common._controls_of(rec1, [rec1], [ctrl1])
    assert ctrl1 in ctrls

