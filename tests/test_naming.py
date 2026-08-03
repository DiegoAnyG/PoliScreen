"""IUPAC name verification.

A name that OPSIN parses is not necessarily the name of the molecule being docked:
on a benzofuroxan the N-oxide sits on one of two nitrogens, and putting it on the
wrong one gives a valid name for a different compound with the same formula. These
tests pin the check that caught it.
"""
import pytest

from poliscreen.core import naming as nm

# Ester of the benzofuroxan-5-carboxylic acid core used by the project.
BUTYL_ESTER = "O=C(OCCCC)c1ccc2[n+]([O-])onc2c1"
RIGHT = "butyl 1-oxido-2,1,3-benzoxadiazol-1-ium-5-carboxylate"
WRONG = "butyl 3-oxido-2,1,3-benzoxadiazol-3-ium-5-carboxylate"

needs_opsin = pytest.mark.skipif(not nm.available(), reason="OPSIN not installed")


def test_the_oxide_locant_has_a_variant():
    assert RIGHT in nm.variants(WRONG)
    assert WRONG in nm.variants(RIGHT)


@needs_opsin
def test_the_wrong_oxide_locant_is_corrected():
    """It parses, so only comparing against the structure catches it."""
    (name, ok), = nm.verify([WRONG], [BUTYL_ESTER])
    assert ok and name == RIGHT


@needs_opsin
def test_a_correct_name_is_kept():
    (name, ok), = nm.verify([RIGHT], [BUTYL_ESTER])
    assert ok and name == RIGHT


@needs_opsin
def test_a_name_of_another_molecule_is_dropped():
    """Better no name than a wrong one in a table that ends up in a paper."""
    (name, ok), = nm.verify(["ethyl acetate"], [BUTYL_ESTER])
    assert not ok and name is None


# --------------------------------------------------------------- preferred radicals
@pytest.mark.parametrize("given,expected", [
    ("Furfuryl X", "furan-2-ylmethyl X"),
    ("Isoamyl X", "3-methylbutyl X"),
    ("isopropyl X", "propan-2-yl X"),
    ("1-butyl X", "butyl X"),
    ("1-pentyl X", "pentyl X"),
])
def test_trivial_radicals_become_the_preferred_form(given, expected):
    assert nm.preferred_radical(given) == expected


@pytest.mark.parametrize("kept", ["Benzyl X", "Cyclohexyl X", "tert-butyl X"])
def test_radicals_iupac_retains_are_left_alone(kept):
    assert nm.preferred_radical(kept) is None


def test_a_positional_locant_is_not_stripped():
    """'1-' is redundant on an unbranched chain but meaningful on a ring."""
    assert nm.preferred_radical("1-naphthyl X") is None


@needs_opsin
def test_the_preferred_radical_wins_when_both_verify():
    given = "1-butyl 3-oxido-2,1,3-benzoxadiazol-3-ium-5-carboxylate"
    (name, ok), = nm.verify([given], [BUTYL_ESTER])
    assert ok and name == RIGHT      # 'butyl', not '1-butyl'
