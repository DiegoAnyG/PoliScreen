"""admelab labels its results in Spanish; they are translated where they enter.

Missing a form leaves Spanish in an English table, which is how "desconocido" and
"fenolico" reached the products table.
"""
import pytest

from poliscreen.core.design import english_value as ev


@pytest.mark.parametrize("given,expected", [
    ("primario", "primary"),
    ("secundario", "secondary"),
    ("desconocido", "unknown"),          # masculine: only the feminine was mapped
    ("desconocida", "unknown"),
    ("fenolico", "phenolic"),
    ("buena", "good"),
])
def test_labels_are_translated(given, expected):
    assert ev(given) == expected


def test_a_verdict_with_an_explanation_is_translated():
    got = ev("dificil (fenol poco nucleofilo; usar cloruro de acilo/DMAP)")
    assert got.startswith("difficult (")
    assert "fenol poco" not in got


def test_the_ghs_category_keeps_its_code():
    assert ev("GHS-5 (Baja toxicidad)") == "GHS-5 (Low toxicity)"


def test_an_unknown_value_passes_through():
    """Guessing at a label we do not know would be worse than leaving it."""
    assert ev("una nota que no conocemos") == "una nota que no conocemos"
    assert ev(None) is None
