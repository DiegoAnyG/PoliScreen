"""Ticking a tunnel draws it, on that click and not the next.

This started as st.data_editor with a checkbox column. The editor keeps the user's edits as a diff
and replays it over the frame it is handed -- and that frame was rebuilt from the state the diff
had just produced. The two disagreed for one run, which is felt as having to press twice. Removing
its key fixed the box's own appearance and not the viewer, because the lag was in the mechanism
rather than in the key.

st.checkbox has no diff to replay: its value *is* session state, so what it returns is what was
just done. It is also the only one of the two that AppTest can operate, which is why the table is
built from checkboxes and why this file can exist at all.

The app is driven here rather than a helper called, because every version of this bug lived in the
gap between what the widget returned and what the next panel read.
"""
from pathlib import Path

import pytest

pytest.importorskip("streamlit")

from streamlit.testing.v1 import AppTest  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
APP = str(ROOT / "src" / "poliscreen" / "ui" / "streamlit_app.py")

# One sphere is enough: this is about the selection, not about the geometry.
CLUSTER = ("MODEL        0\n"
           "ATOM      1  H   FIL T   1      -1.000   0.000   0.000        2.00\n")


@pytest.fixture()
def app(tmp_path):
    """The Run stage with three tunnels already found."""
    clusters = tmp_path / "caver" / "out" / "data" / "clusters"
    clusters.mkdir(parents=True)
    for n in (1, 2, 3):
        (clusters / f"tun_cl_00{n}.pdb").write_text(CLUSTER)

    at = AppTest.from_file(APP, default_timeout=180)
    at.session_state["stage"] = "Run"
    at.session_state["tun_drawn"] = str(tmp_path / "caver" / "out")
    at._tunnel_root = str(tmp_path / "caver" / "out")
    return at


def drawn(at):
    """What the viewer would draw, which is the thing the user is actually judging.

    The root is passed explicitly: the viewer reads it from the running session, and these tests
    run outside that session. AppTest's session_state also has no .get(), so it is read by
    subscript.
    """
    from poliscreen.ui.streamlit_app import _tunnel_groups

    try:
        picked = at.session_state["tun_shown"]
    except (KeyError, AttributeError):
        picked = []
    return sorted(g["number"] for g in _tunnel_groups(set(picked or []), root=at._tunnel_root))


def boxes(at):
    return [c for c in at.checkbox if c.key and c.key.startswith("tun_draw_")]


def test_the_first_tunnel_is_drawn_before_anything_is_touched(app):
    app.run()
    assert not app.exception
    assert app.session_state["tun_shown"] == [1]
    assert drawn(app) == [1]


def test_one_click_draws_it(app):
    """The bug: it took two. The first press was overwritten before the viewer read it."""
    app.run()
    boxes(app)[2].check()
    app.run()
    assert app.session_state["tun_shown"] == [1, 3]
    assert drawn(app) == [1, 3]


def test_one_click_undraws_it(app):
    app.run()
    boxes(app)[0].uncheck()
    app.run()
    assert app.session_state["tun_shown"] == []
    assert drawn(app) == []


def test_unticking_everything_leaves_the_viewer_empty_and_it_stays_empty(app):
    """Tunnel 1 used to come back on the redraw, because an empty list is falsy."""
    app.run()
    boxes(app)[0].uncheck()
    app.run()
    app.run()
    assert app.session_state["tun_shown"] == []
    assert drawn(app) == []


def test_choosing_the_third_from_empty_gives_the_third(app):
    """Reported exactly this way: with none ticked, pressing 3 ticked 1."""
    app.run()
    boxes(app)[0].uncheck()
    app.run()
    boxes(app)[2].check()
    app.run()
    assert app.session_state["tun_shown"] == [3]
    assert drawn(app) == [3]


def test_each_row_keeps_its_own_answer(app):
    app.run()
    boxes(app)[1].check()
    app.run()
    boxes(app)[0].uncheck()
    app.run()
    assert app.session_state["tun_shown"] == [2]
    assert drawn(app) == [2]
