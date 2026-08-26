"""Interface rules that fail silently.

Streamlit does not raise for either of these: it logs a warning nobody reads and
carries on with the wrong value, or it strips markup and shows nothing.
"""
import ast
from pathlib import Path

import pytest

pytest.importorskip("streamlit")

APP = Path(__file__).resolve().parent.parent / "src" / "poliscreen" / "ui" / "streamlit_app.py"
_TREE = ast.parse(APP.read_text(encoding="utf-8"))

# Keys reassigned on every pass so they survive changing stage. A widget cannot take one of these
# AND declare a default: Streamlit warns and drops the default, so the value asked for is ignored.
PREFIXES = ("pep_", "modo_", "cat_", "sec_", "rec_", "box_", "sites_", "rx_",
            "cx_", "cy_", "cz_", "sx_", "sy_", "sz_", "src_", "vis_", "cfg_", "tun_")

# A button has no value to carry, and Streamlit refuses to have one assigned: the reassignment
# loop that keeps the persistent widgets alive across stages raises on the first button whose key
# happens to start with one of those prefixes. It reaches the user as
# StreamlitValueAssignmentNotAllowedError, from a line that only draws a button.
PRESSED_NOT_STORED = ("button", "download_button", "form_submit_button", "link_button")


def _persistent_prefixes():
    """The prefixes the app itself reassigns, read from the app rather than copied.

    Copying them is how this test was first written and it was wrong the same day: the list here
    was stricter than the app's, so it flagged six keys of which only two could ever crash.
    """
    for node in ast.walk(_TREE):
        if (isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name)
                and node.targets[0].id == "_PERSISTENT_PREFIXES"):
            return tuple(ast.literal_eval(node.value))
    raise AssertionError("_PERSISTENT_PREFIXES is gone; does this rule still hold?")
DEFAULT_KWARGS = ("value", "index", "default")


def test_no_persistent_widget_declares_a_default():
    offenders = []
    for call in ast.walk(_TREE):
        if not isinstance(call, ast.Call):
            continue
        fn = call.func
        name = fn.attr if isinstance(fn, ast.Attribute) else ""
        key = next((k.value.value for k in call.keywords
                    if k.arg == "key" and isinstance(k.value, ast.Constant)), None)
        if not (isinstance(key, str) and key.startswith(PREFIXES)):
            continue
        if any(k.arg in DEFAULT_KWARGS for k in call.keywords):
            offenders.append(f"line {call.lineno}: st.{name}(key={key!r})")
    assert not offenders, f"seed these with S.setdefault instead: {offenders}"


def test_no_button_takes_a_persistent_key():
    """The loop that keeps widgets alive across stages assigns to every key with one of those
    prefixes. A button cannot be assigned to, so one named that way crashes the panel it is on --
    which is how `vis_tun_pml` took down the transport view."""
    offenders = []
    for call in ast.walk(_TREE):
        if not isinstance(call, ast.Call) or not isinstance(call.func, ast.Attribute):
            continue
        if call.func.attr not in PRESSED_NOT_STORED:
            continue
        key = next((k.value.value for k in call.keywords
                    if k.arg == "key" and isinstance(k.value, ast.Constant)), None)
        if isinstance(key, str) and key.startswith(_persistent_prefixes()):
            offenders.append(f"line {call.lineno}: st.{call.func.attr}(key={key!r})")
    assert not offenders, (
        "a button cannot hold a value, so its key must not start with a persistent prefix: "
        f"{offenders}")


def test_the_notice_uses_a_dialog():
    """Custom HTML did not survive Streamlit's sanitizing and a toast clips its own text:
    neither was ever visible. The dialog is the mechanism that works."""
    src = APP.read_text(encoding="utf-8")
    assert "def _notice_dialog" in src
    assert "st.toast" not in src


def test_the_notice_is_opened_after_the_stages():
    """It is set while the active stage runs; opened earlier, it would only appear one
    interaction later."""
    src = APP.read_text(encoding="utf-8")
    assert src.index('_STAGE_FN[S["stage"]]()') < src.index('S.pop("_notice")')
