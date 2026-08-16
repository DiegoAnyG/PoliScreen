"""The control's chemistry must not depend on who compiled the file reader.

A coordinate file has no field for bond orders, so they were inferred from geometry. Two builds of
the same converter read one fused aromatic ring differently, and the losing assignment left a ring
nitrogen free to be protonated at pH 7.4 -- which it is not. Since the control is what the ranking
is normalised against, two machines on the same commit were measuring against different molecules.
The PDB publishes the answer for every hetero component; these tests keep us asking for it.
"""
import io

import pytest

from poliscreen.core import receptor as rc

Chem = pytest.importorskip("rdkit.Chem")

PYRIDINE = "c1ccncc1"


def _block(smiles=PYRIDINE):
    return Chem.MolToMolBlock(Chem.MolFromSmiles(smiles))


def _serving(block, seen=None):
    def fake_urlopen(url, timeout=None):
        if seen is not None:
            seen.append(url)
        return io.BytesIO(block.encode())
    return fake_urlopen


def test_the_published_definition_is_what_sets_the_bond_orders(tmp_path, monkeypatch):
    seen = []
    monkeypatch.setattr(rc.urllib.request, "urlopen", _serving(_block(), seen))
    smiles = rc.ccd_template("zi9", cache=tmp_path)
    assert Chem.CanonSmiles(smiles) == Chem.CanonSmiles(PYRIDINE)
    assert seen and seen[0].endswith("ZI9_ideal.sdf"), "the code is upper-cased in the dictionary"


def test_it_is_cached_so_a_rerun_needs_no_network(tmp_path, monkeypatch):
    monkeypatch.setattr(rc.urllib.request, "urlopen", _serving(_block()))
    rc.ccd_template("ZI9", cache=tmp_path)
    assert (tmp_path / "ZI9_ccd.sdf").exists(), "the project must record the definition it used"

    def refuse(*a, **k):
        raise AssertionError("went to the network with a cached definition on disk")

    monkeypatch.setattr(rc.urllib.request, "urlopen", refuse)
    assert Chem.CanonSmiles(rc.ccd_template("ZI9", cache=tmp_path)) == Chem.CanonSmiles(PYRIDINE)


def test_no_network_is_not_a_failure(tmp_path, monkeypatch):
    """Offline, the old behaviour has to remain available -- degraded, not broken."""
    def down(*a, **k):
        raise OSError("no route to host")

    monkeypatch.setattr(rc.urllib.request, "urlopen", down)
    assert rc.ccd_template("ZI9", cache=tmp_path) is None


def test_a_meaningless_code_asks_nobody(monkeypatch):
    def refuse(*a, **k):
        raise AssertionError("looked up an empty component code")

    monkeypatch.setattr(rc.urllib.request, "urlopen", refuse)
    assert rc.ccd_template("") is None
    assert rc.ccd_template("   ") is None
