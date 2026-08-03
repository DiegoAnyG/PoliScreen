"""Interface tests with AppTest: walks every stage and every mode looking for exceptions.

A failure in a stage only shows up when it is drawn, and the app only draws the active one. Without
this, a change in Results can break Ligands and go unnoticed until someone opens it.
"""
from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest  # noqa: E402

APP = str(Path(__file__).resolve().parent.parent / "src" / "poliscreen" / "ui" / "streamlit_app.py")
STAGES = ["Receptors", "Ligands", "Run", "Results"]
MODES = ["Build by reaction", "Generate peptides", "Upload ready ligands"]


def _app(tmp_path, **state):
    at = AppTest.from_file(APP, default_timeout=600)
    at.run()
    at.session_state["proj_dir"] = str(tmp_path)
    for k, v in state.items():
        at.session_state[k] = v
    at.run()
    return at


def _no_exception(at, what):
    if at.exception:
        detail = "; ".join(str(getattr(e, "value", e))[:200] for e in at.exception)
        pytest.fail(f"{what}: {detail}")


@pytest.mark.parametrize("stage", STAGES)
def test_every_stage_is_drawn(tmp_path, stage):
    _no_exception(_app(tmp_path, stage=stage), f"stage {stage}")


@pytest.mark.parametrize("stage", STAGES)
def test_every_stage_is_drawn_in_spanish(tmp_path, stage):
    """The Spanish interface must render every stage without exceptions, and actually translate."""
    at = _app(tmp_path, _lang_pick="Español", stage=stage)
    _no_exception(at, f"stage {stage} in Spanish")
    assert any(b.label == "Ejecutar" for b in at.button), "the stage bar did not translate to Spanish"


@pytest.mark.parametrize("mode", MODES)
def test_every_ligand_mode_is_drawn(tmp_path, mode):
    _no_exception(_app(tmp_path, stage="Ligands", ligand_mode=mode), f"mode {mode}")


def test_an_empty_folder_breaks_no_stage(tmp_path):
    """A freshly created project has no tables or poses: every view must tolerate it."""
    for stage in STAGES:
        _no_exception(_app(tmp_path, stage=stage), f"empty project in {stage}")


def test_both_viewer_views(tmp_path):
    for view in ("Summary", "3D complex"):
        _no_exception(_app(tmp_path, stage="Results", vis_res_view=view), f"view {view}")


def test_state_survives_changing_stage(tmp_path):
    """Streamlit discards widgets not drawn; without reassigning them the parameters were lost."""
    at = _app(tmp_path, stage="Ligands", ligand_mode="Generate peptides", pep_len=7)
    at.session_state["stage"] = "Results"
    at.run()
    at.session_state["stage"] = "Ligands"
    at.run()
    assert at.session_state["ligand_mode"] == "Generate peptides"
    assert at.session_state["pep_len"] == 7
    _no_exception(at, "round trip between stages")


def test_a_windows_path_creates_no_folder_with_backslashes(tmp_path):
    at = AppTest.from_file(APP, default_timeout=600)
    at.run()
    at.session_state["proj_dir"] = str(tmp_path).replace("/", "\\")
    at.run()
    _no_exception(at, "path with backslashes")
    assert not [p for p in Path.cwd().iterdir() if "\\" in p.name]


def test_the_downloads_dialog_opens(tmp_path):
    (tmp_path / "ranking.csv").write_text("compound,best_dock\na,-8.0\n")
    at = _app(tmp_path, stage="Results", _open_downloads=True)
    _no_exception(at, "downloads dialog")
    assert [c for c in at.checkbox if c.key and c.key.startswith("dl_")]
