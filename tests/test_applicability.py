"""Applicability domain of the ADMET-AI predictions.

It qualifies the ADMET numbers rather than replacing them: outside the domain a prediction is
unsupported by the training data, not wrong. The columns therefore have to reach the ranking beside
the numbers they qualify — a hard-coded whitelist there used to drop anything it did not name.

The bridge itself needs admelab >= 0.3 and is skipped when that is not installed; the wiring around
it is checked either way, because the installer ships no admelab at all.
"""
import pandas as pd
import pytest

from poliscreen.core import pipeline as pl
from poliscreen.core.design import AdmelabBridge


class _Bridge:
    """Stands in for admelab: the pipeline must not care where the numbers came from."""

    def __init__(self, has_domain=True, rows=None, boom=None):
        self._has, self._rows, self._boom = has_domain, rows or [], boom

    def has_applicability(self):
        return self._has

    def applicability(self, smiles, **kw):
        if self._boom:
            raise self._boom
        return type("R", (), {"to_dataframe": lambda s: pd.DataFrame(self._rows)})()


def _analogs():
    return pd.DataFrame({"name": ["a1", "a2"],
                         "SMILES": ["CCO", "c1ccccc1"],
                         "LD50_mg_per_kg": [1000.0, 2000.0]})


def _cfg(use_ml=True):
    return pl.RunConfig(receptors=[], out_dir=".", use_ml=use_ml)


def test_the_domain_columns_are_attached_to_the_analogues(monkeypatch):
    rows = [{"SMILES": "CCO", "fraction_in_domain": 1.0, "median_similarity": 0.81},
            {"SMILES": "c1ccccc1", "fraction_in_domain": 0.5, "median_similarity": 0.22}]
    monkeypatch.setattr(pl, "AdmelabBridge", lambda: _Bridge(rows=rows))
    df = _analogs()
    pl._add_applicability(df, _cfg(), lambda *a: None)
    assert list(df["fraction_in_domain"]) == [1.0, 0.5]
    assert list(df["median_similarity"]) == [0.81, 0.22]


def test_they_survive_the_merge_into_the_ranking():
    """The whitelist in the ranking merge must name them, or they are silently dropped."""
    for c in pl.AD_COLUMNS:
        assert c in ("fraction_in_domain", "median_similarity")
    src = pl.__file__
    text = open(src, encoding="utf-8").read()
    merge = text[text.index("cols_admet = "):text.index("cols_admet = ") + 400]
    assert "AD_COLUMNS" in merge, "the ranking merge would drop the domain columns"


def test_nothing_is_computed_without_admet_ai(monkeypatch):
    """With --no-ml there are no ADMET-AI predictions, so there is nothing to qualify."""
    monkeypatch.setattr(pl, "AdmelabBridge", lambda: _Bridge(rows=[{"SMILES": "CCO"}]))
    df = _analogs()
    pl._add_applicability(df, _cfg(use_ml=False), lambda *a: None)
    assert "fraction_in_domain" not in df.columns


def test_an_older_admelab_is_reported_and_skipped(monkeypatch):
    monkeypatch.setattr(pl, "AdmelabBridge", lambda: _Bridge(has_domain=False))
    said = []
    df = _analogs()
    pl._add_applicability(df, _cfg(), lambda k, v: said.append(v))
    assert "fraction_in_domain" not in df.columns
    assert any("0.3" in s for s in said), "the user is not told why the columns are missing"


def test_a_broken_bridge_does_not_take_the_screening_down(monkeypatch):
    monkeypatch.setattr(pl, "AdmelabBridge",
                        lambda: _Bridge(boom=pl.AdmelabError("no admelab here")))
    df = _analogs()
    pl._add_applicability(df, _cfg(), lambda *a: None)
    assert "fraction_in_domain" not in df.columns
    assert list(df["LD50_mg_per_kg"]) == [1000.0, 2000.0]


@pytest.mark.skipif(not AdmelabBridge().has_applicability(),
                    reason="needs admelab >= 0.3 installed")
def test_the_real_engine_separates_familiar_from_exotic():
    """Aspirin is in the training sets; a polysilyne is not. If the metric cannot tell them
    apart it is measuring nothing."""
    b = AdmelabBridge()
    r = b.applicability(["CC(=O)Oc1ccccc1C(=O)O", "[Si](C)(C)(C)C#CC#CC#C[Si](C)(C)C"],
                        endpoints=["hERG", "AMES"])
    by = {row["SMILES"]: row for row in r.rows}
    known = by["CC(=O)Oc1ccccc1C(=O)O"]["median_similarity"]
    exotic = by["[Si](C)(C)(C)C#CC#CC#C[Si](C)(C)C"]["median_similarity"]
    assert known > exotic, f"aspirin {known} should sit closer than the silyne {exotic}"
