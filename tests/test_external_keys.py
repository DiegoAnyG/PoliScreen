"""Upstream renames a key, and the code goes quiet instead of loud.

Twice now, the same shape of fault. PubChem renamed IsomericSMILES to SMILES, and the reagent
discovery kept asking for the old spelling: the search still returned rows, every row was missing
the key, every row was skipped, and the interface reported that PubChem had nothing useful. The
design engine renamed its feasibility verdict the same way, and every product read `unknown` while
silently counting as synthesizable.

Neither raised. Both are only visible by reading a value that arrived under a name nobody asked
for. These pin the readers, not the network, so they cannot go flaky.
"""
from pathlib import Path

from poliscreen.core import reagents as rg

RUNNER = (Path(__file__).resolve().parent.parent / "src" / "poliscreen" / "core"
          / "_admelab_runner.py")


def test_pubchem_smiles_is_read_under_any_spelling():
    """The 2025 names first, then the ones they replaced."""
    for key in ("SMILES", "IsomericSMILES", "ConnectivitySMILES", "CanonicalSMILES"):
        assert rg.smiles_of({key: "CCO", "Title": "ethanol"}) == "CCO", f"{key} was not read"


def test_a_record_with_no_smiles_at_all_is_skipped_not_guessed():
    assert rg.smiles_of({"Title": "something", "CID": 1}) is None


def test_the_request_asks_for_the_name_pubchem_uses_now():
    """Asking for the retired name still works, which is exactly why this went unnoticed."""
    assert rg.SMILES_PROPERTY == "SMILES"


def test_the_feasibility_verdict_is_read_under_either_vocabulary():
    """The same fault in the design bridge: English now, Spanish before."""
    source = RUNNER.read_text(encoding="utf-8")
    for key in ("fischer_viability", "viabilidad_fischer"):
        assert key in source, f"{key} is no longer read; a rename would go silent again"
    for verdict in ("unfavorable", "difficult", "desfavorable", "dificil"):
        assert verdict in source, f"{verdict} is not excluded, so it would count as synthesizable"
