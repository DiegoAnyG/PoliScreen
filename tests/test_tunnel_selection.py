"""Unticking every tunnel has to mean none, not "back to the first one".

Two bugs met here, and the first hid the second.

`S.get("tun_shown") or numbers[:1]` treats an empty list as absent, because an empty list is
falsy. So unticking the last row put tunnel 1 back on the next redraw -- and since the frame is
rebuilt from that same state, the click on tunnel 3 was overwritten before it was ever read. The
user saw tunnel 1 tick itself while the row they pressed stayed empty.

Underneath, st.data_editor was given a key. With one it stores the edits as a diff and replays
them over the frame it is handed, which is the frame we rebuild from state: the two fight and a
change only lands on the second press. The editor is now uncontrolled by key and controlled by the
frame, which is one source of truth instead of two.
"""
from poliscreen.ui.streamlit_app import _selected_tunnels

AVAILABLE = [1, 2, 3, 4, 5, 6]


def test_nothing_chosen_yet_shows_the_first():
    """None is 'the panel has just opened', and a blank viewer would look broken."""
    assert _selected_tunnels(AVAILABLE, None) == [1]


def test_unticking_everything_is_obeyed():
    """The bug: an empty list is falsy, so tunnel 1 came back and stole the next click."""
    assert _selected_tunnels(AVAILABLE, []) == []


def test_a_choice_is_kept_exactly():
    assert _selected_tunnels(AVAILABLE, [3]) == [3]
    assert _selected_tunnels(AVAILABLE, [2, 5]) == [2, 5]


def test_choosing_the_third_when_none_are_chosen_gives_the_third():
    """Stated the way it was reported: pick 3 from empty and 3 is what is drawn."""
    after_clearing = _selected_tunnels(AVAILABLE, [])
    assert _selected_tunnels(AVAILABLE, after_clearing + [3]) == [3]


def test_tunnels_from_a_previous_receptor_are_dropped():
    """A different structure has different clusters; a remembered 9 must not survive into it."""
    assert _selected_tunnels([1, 2], [2, 9]) == [2]
    assert _selected_tunnels([1, 2], [9]) == []


def test_the_editor_is_not_given_a_key():
    """A key makes st.data_editor keep the edits as a diff and replay them over the frame we
    rebuild from state. That is the two-clicks-to-untick bug, and no test that does not render
    can see it -- so the call is checked instead.

    Read with the parser, not by slicing the text: the first attempt split the source on the first
    ")" and only ever looked at `pd.DataFrame(rows)`, so putting the key back went unnoticed.
    """
    import ast
    import inspect
    import textwrap

    from poliscreen.ui import streamlit_app

    tree = ast.parse(textwrap.dedent(inspect.getsource(streamlit_app._tunnel_table)))
    calls = [n for n in ast.walk(tree)
             if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
             and n.func.attr == "data_editor"]
    assert calls, "the tunnel table no longer uses st.data_editor; check this still applies"
    for call in calls:
        assert "key" not in {k.arg for k in call.keywords}, "the tunnel editor took a key again"
