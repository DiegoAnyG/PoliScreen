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
