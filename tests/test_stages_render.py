"""Every stage has to draw without raising.

The rest of the interface tests read the source: they catch a shadowed helper or a widget that
declares a default it will not get. They cannot catch a name that is never defined on the path
taken, because nothing imports that path -- a missing `import re` inside a function passed the
whole suite and broke the Run stage on sight.

Streamlit's own AppTest executes the script, which is the only thing that does.
"""
import pytest

pytest.importorskip("streamlit")

from streamlit.testing.v1 import AppTest  # noqa: E402

APP = str(__import__("pathlib").Path(__file__).resolve().parent.parent
          / "src" / "poliscreen" / "ui" / "streamlit_app.py")


@pytest.fixture(autouse=True)
def isolated_projects_dir(monkeypatch, tmp_path):
    monkeypatch.setenv("POLISCREEN_PROJECTS", str(tmp_path))
    mpl_dir = tmp_path / "mpl"
    mpl_dir.mkdir(parents=True, exist_ok=True)
    monkeypatch.setenv("MPLCONFIGDIR", str(mpl_dir))


@pytest.mark.parametrize("stage", ["Receptors", "Ligands", "Run", "Results"])
def test_the_stage_draws(stage):
    app = AppTest.from_file(APP, default_timeout=180)
    app.session_state["stage"] = stage
    app.run()
    assert not app.exception, f"{stage}: {app.exception[0].message if app.exception else ''}"


def test_both_stages_that_carry_tunnels_show_the_tab():
    """Run and Results each hold a Transport tunnels tab; a rename that loses one is silent."""
    for stage in ("Run", "Results"):
        app = AppTest.from_file(APP, default_timeout=180)
        app.session_state["stage"] = stage
        app.run()
        assert "Transport tunnels" in [s.value for s in app.get("subheader")], stage


def test_results_stage_with_screening_results(tmp_path, monkeypatch):
    """Results stage with real ranking table must render downloads, charts and tables without NameError."""
    import json
    import pandas as pd

    proj = tmp_path / "proj1"
    proj.mkdir()
    meta = {"project": "proj1", "receptors": ["rec1.pdb"]}
    (proj / "run.json").write_text(json.dumps(meta), encoding="utf-8")
    rk = pd.DataFrame([{
        "compound": "cmp1",
        "receptor": "rec1",
        "is_control": 1,
        "best_dock": -8.0,
        "best_inter": 1.0,
        "effectiveness_pct": 100.0,
        "confidence": 0.8,
        "pKi": 6.0,
        "LE": 0.3,
        "pareto_rank": 1,
        "pose": 1,
    }])
    rk.to_csv(proj / "ranking.csv", index=False)

    app = AppTest.from_file(APP, default_timeout=180)
    app.session_state["stage"] = "Results"
    app.session_state["proj_dir"] = str(proj)
    app.run()
    assert not app.exception, f"Results with screening data failed: {app.exception[0].message if app.exception else ''}"

