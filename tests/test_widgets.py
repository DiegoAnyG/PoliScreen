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
            "cx_", "cy_", "cz_", "sx_", "sy_", "sz_", "src_", "vis_", "cfg_")
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
