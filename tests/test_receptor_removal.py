"""Removing a receptor from the Receptors panel.

Without this, the only way to drop a receptor was to point the project at another folder. The
control has to go with it: extracted from that structure it shares its coordinate system, and left
behind the pipeline would assign it by geometry to whichever receptor remained.
"""
from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest  # noqa: E402

from poliscreen.core import layout as lay  # noqa: E402

APP = str(Path(__file__).resolve().parent.parent / "src" / "poliscreen" / "ui" / "streamlit_app.py")


def _atom(serial, name, resname, chain, resseq, x, y, z):
    return (f"ATOM  {serial:>5} {name:<4} {resname:>3} {chain}{resseq:>4}    "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           C")


def _project(tmp_path):
    """A project with two receptors far apart, each with its own control beside it."""
    rec_dir = lay.artifact(tmp_path, lay.RECEPTORS)
    rec_dir.mkdir(parents=True, exist_ok=True)
    made = {}
    for i, (stem, x0) in enumerate((("AAAA", 0.0), ("BBBB", 500.0))):
        rec = rec_dir / f"{stem}_ready.pdb"
        rec.write_text("\n".join(
            _atom(n + 1, "CA", "ALA", "A", n + 1, x0 + n, 0.0, 0.0) for n in range(8)) + "\nEND\n")
        ctrl = rec_dir / f"control_LG{i}.pdb"
        ctrl.write_text("\n".join(
            _atom(n + 1, "C1", "LIG", "A", 1, x0 + n, 1.0, 0.0) for n in range(3)) + "\nEND\n")
        made[stem] = (rec, ctrl)
    return made


def _run(tmp_path, made):
    at = AppTest.from_file(APP, default_timeout=600)
    at.run()
    at.session_state["proj_dir"] = str(tmp_path)
    at.session_state["stage"] = "Receptors"
    at.session_state["receptors"] = [str(r) for r, _ in made.values()]
    at.session_state["controls"] = [str(c) for _, c in made.values()]
    at.run()
    assert not at.exception, "; ".join(str(getattr(e, "value", e))[:300] for e in at.exception)
    return at


def test_removing_one_receptor_takes_its_control_and_leaves_the_other(tmp_path):
    made = _project(tmp_path)
    at = _run(tmp_path, made)

    at.button(key="drop_receptor_aaaa").click().run()
    assert not at.exception, "; ".join(str(getattr(e, "value", e))[:300] for e in at.exception)

    gone_rec, gone_ctrl = made["AAAA"]
    kept_rec, kept_ctrl = made["BBBB"]
    assert not gone_rec.exists(), "the receptor file survived"
    assert not gone_ctrl.exists(), "its control was left orphaned on disk"
    assert kept_rec.exists() and kept_ctrl.exists(), "the other receptor was taken down with it"
    assert at.session_state["receptors"] == [str(kept_rec)]
    assert at.session_state["controls"] == [str(kept_ctrl)]


def test_remove_all_empties_the_folder(tmp_path):
    made = _project(tmp_path)
    at = _run(tmp_path, made)

    at.button(key="wipe_receptors").click().run()
    assert not at.exception, "; ".join(str(getattr(e, "value", e))[:300] for e in at.exception)

    for rec, ctrl in made.values():
        assert not rec.exists() and not ctrl.exists()
    assert at.session_state["receptors"] == []
    assert at.session_state["controls"] == []


def test_the_choices_made_for_a_removed_receptor_are_forgotten(tmp_path):
    made = _project(tmp_path)
    at = AppTest.from_file(APP, default_timeout=600)
    at.run()
    at.session_state["proj_dir"] = str(tmp_path)
    at.session_state["stage"] = "Receptors"
    at.session_state["receptors"] = [str(r) for r, _ in made.values()]
    at.session_state["controls"] = [str(c) for _, c in made.values()]
    at.session_state["rec_chains_aaaa"] = ["A"]
    at.session_state["_signature_prep_aaaa"] = ("whatever",)
    at.run()

    at.button(key="drop_receptor_aaaa").click().run()
    assert "rec_chains_aaaa" not in at.session_state, "the chain choice outlived its receptor"
    assert "_signature_prep_aaaa" not in at.session_state, "it would still count as prepared"
