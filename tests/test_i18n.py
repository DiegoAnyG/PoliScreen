"""Translation-catalog tests.

A wrong translation does not fail loudly: it either shows English where Spanish was
expected, or crashes at .format() time in a branch nobody exercised. These tests
check the two things that can actually break the interface.
"""
import ast
import re
from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from poliscreen.ui import i18n  # noqa: E402

APP = Path(__file__).resolve().parent.parent / "src" / "poliscreen" / "ui" / "streamlit_app.py"
_PLACEHOLDER = re.compile(r"\{([A-Za-z_][A-Za-z0-9_]*)")


def _placeholders(s):
    return sorted(_PLACEHOLDER.findall(s))


@pytest.fixture(autouse=True)
def _isolated_language():
    """The language lives in a global session_state; leaving it set would make the other
    interface tests render in another language."""
    previous = i18n.st.session_state.get("lang")
    yield
    if previous is None:
        i18n.st.session_state.pop("lang", None)
    else:
        i18n.st.session_state["lang"] = previous


def test_translations_keep_their_placeholders():
    """A translation that loses or renames a {placeholder} raises KeyError on .format()."""
    bad_entries = [en for en, es in i18n._ES.items() if _placeholders(en) != _placeholders(es)]
    assert not bad_entries, f"translations with mismatched placeholders: {bad_entries[:5]}"


def test_no_spanish_text_is_left_in_the_base():
    """English is the source language: a Spanish string wrapped in t() would never be
    translated back, and would show in Spanish to an English user."""
    markers = (" el ", " la ", " los ", " las ", " que ", "ó", "¿", "ñ", "á", "é", "ú")
    src = APP.read_text(encoding="utf-8")
    sospechosas = []
    for node in ast.walk(ast.parse(src)):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "t"
                and node.args and isinstance(node.args[0], ast.Constant)
                and isinstance(node.args[0].value, str)):
            s = node.args[0].value
            if any(m in s for m in markers):
                sospechosas.append(s[:60])
    assert not sospechosas, f"Spanish text in the English base: {sospechosas[:5]}"


def test_the_default_language_returns_the_source():
    i18n.st.session_state["lang"] = "en"
    assert i18n.t("Receptors") == "Receptors"


def test_spanish_translates_and_falls_back_to_english():
    i18n.st.session_state["lang"] = "es"
    assert i18n.t("Receptors") == "Receptores"
    assert i18n.t("QED") == "QED"           # not in the catalog: falls back, does not fail
