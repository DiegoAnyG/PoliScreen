"""Approved-drug screening: the filters, and the shape of what ChEMBL sends.

Nothing here touches the network. The parsing is pinned because the keys are somebody else's to
rename -- twice already this month an upstream did exactly that and the reader went quiet instead
of loud -- and the filters are pinned because a bound that silently admits everything is the same
class of fault as the synthesizability filter that stopped filtering.
"""
import pytest

from poliscreen.core import drugs as dg

pytest.importorskip("rdkit")

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
# Ciclosporin: an approved drug, and far outside the rule of five. The filters have to be able to
# say no to something real, not only to something invented for a test.
CICLOSPORIN = ("CC[C@H]1NC(=O)[C@H]([C@H](O)[C@H](C)C/C=C/C)N(C)C(=O)[C@H](C(C)C)N(C)C(=O)"
               "[C@H](CC(C)C)N(C)C(=O)[C@H](CC(C)C)N(C)C(=O)[C@@H](C)NC(=O)[C@H](C)NC(=O)"
               "[C@H](CC(C)C)N(C)C(=O)[C@@H](NC(=O)[C@H](CC(C)C)N(C)C(=O)CN(C)C1=O)C(C)C")


def test_descriptors_are_the_ones_the_filters_name():
    d = dg.descriptors(ASPIRIN)
    assert set(dg.PROPERTIES) <= set(d), "a property can be filtered on but never computed"
    assert 179 < d["MW"] < 181, d["MW"]
    assert d["HBD"] == 1 and d["HBA"] == 3


def test_an_unreadable_structure_yields_nothing_rather_than_zeros():
    """Zeros would pass every upper bound, which is how junk reaches the docking."""
    assert dg.descriptors("not a molecule") == {}
    assert dg.descriptors("") == {}


def test_a_molecule_with_no_value_fails_the_filter():
    assert not dg.passes({"MW": None}, {"MW": (0, 500)})
    assert not dg.passes({}, {"MW": (0, 500)})


def test_lipinski_admits_aspirin_and_refuses_ciclosporin():
    kept = dg.apply_filters([{"name": "aspirin", "smiles": ASPIRIN},
                             {"name": "ciclosporin", "smiles": CICLOSPORIN}], dg.LIPINSKI)
    assert [r["name"] for r in kept] == ["aspirin"]


def test_no_limits_keeps_everything_that_parses():
    rows = [{"name": "a", "smiles": ASPIRIN}, {"name": "junk", "smiles": "xyz"}]
    assert [r["name"] for r in dg.apply_filters(rows, {})] == ["a"]


def test_an_open_bound_only_constrains_one_side():
    assert dg.passes({"LogP": -3.0}, {"LogP": (None, 5.0)})
    assert not dg.passes({"LogP": -3.0}, {"LogP": (0.0, 5.0)})


def test_biotherapeutics_are_dropped_rather_than_failing_later():
    """Antibodies come back with no SMILES at all; 3D generation is the wrong place to find out."""
    page = {"molecules": [
        {"molecule_chembl_id": "CHEMBL1", "pref_name": "SMALL",
         "molecule_structures": {"canonical_smiles": ASPIRIN}},
        {"molecule_chembl_id": "CHEMBL2", "pref_name": "ANTIBODY", "molecule_structures": None},
        {"molecule_chembl_id": "CHEMBL3", "pref_name": "BLANK",
         "molecule_structures": {"canonical_smiles": "  "}},
    ]}
    assert [r["name"] for r in dg._rows_from_page(page)] == ["SMALL"]


def test_the_library_survives_a_round_trip_to_disk(tmp_path):
    """It is cached in the project so a run records which library it used."""
    rows = [{"name": "aspirin", "smiles": ASPIRIN, "chembl_id": "CHEMBL25", "source": "chembl"}]
    path = tmp_path / "chembl_approved.csv"
    dg.write_csv(path, rows)
    assert dg.read_csv(path) == rows


def test_a_cached_library_is_used_instead_of_the_network(tmp_path):
    path = tmp_path / "chembl_approved.csv"
    dg.write_csv(path, [{"name": "aspirin", "smiles": ASPIRIN, "chembl_id": "X", "source": "chembl"}])
    rows, notice = dg.fetch_approved(cache=path)
    assert notice is None and [r["name"] for r in rows] == ["aspirin"]
