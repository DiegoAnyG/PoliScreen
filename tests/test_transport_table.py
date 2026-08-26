"""What the Results table shows, and what it stops showing.

Three things were wrong at once on a real run, and each was reported as something it was not:

- Both bounds ticked ran two jobs, so one calculation filled two rows of identical numbers, one of
  them marked as having no upper bound. CaverDock computes the lower bound on its way to the upper
  one; the separate run took 7 min 33 s per compound and added nothing.
- An upper bound that was calculated and refused to converge looked the same as one never asked
  for.
- A coverage check written for CaverWeb, where every combination is run in both directions,
  announced that combinations "were never calculated" when the user had simply not asked for them.
"""
import pandas as pd

from poliscreen.ui.streamlit_app import _one_row_per_route, _readable_transport

COLUMNS = ["receptor", "ligand", "tunnel", "direction", "Ea", "flags", "source"]


def frame(rows):
    return pd.DataFrame(rows, columns=COLUMNS)


def test_the_pair_of_bounds_collapses_to_the_richer_row():
    """Two folders, one calculation. The one that has the upper bound is the one to keep."""
    table = frame([
        ("8HTB", "benzo", 3, "in", 2.6, "lower_bound_only", "lb/analysis-lb.pdbqt"),
        ("8HTB", "benzo", 3, "in", 2.6, "", "ub/analysis-lb.pdbqt"),
    ])
    out = _one_row_per_route(table)
    assert len(out) == 1
    assert out.iloc[0]["flags"] == ""


def test_a_refusal_outranks_silence_about_it():
    """Between "no upper bound" and "the upper bound did not pass", the second is the result."""
    table = frame([
        ("8HTB", "pentanol", 3, "in", 2.5, "lower_bound_only", "lb/analysis-lb.pdbqt"),
        ("8HTB", "pentanol", 3, "in", 2.5, "upper_bound_failed", "ub/analysis-lb.pdbqt"),
    ])
    out = _one_row_per_route(table)
    assert len(out) == 1
    assert out.iloc[0]["flags"] == "upper_bound_failed"


def test_different_routes_are_never_collapsed():
    table = frame([
        ("8HTB", "benzo", 3, "in", 2.6, "", "a"),
        ("8HTB", "benzo", 3, "out", 9.1, "", "b"),
        ("8HTB", "benzo", 4, "in", 1.2, "", "c"),
        ("8HTB", "etanol", 3, "in", 3.1, "", "d"),
    ])
    assert len(_one_row_per_route(table)) == 4


def test_what_happened_is_said_in_words():
    table = frame([
        ("8HTB", "a", 3, "in", 2.6, "", "x"),
        ("8HTB", "b", 3, "in", 2.5, "upper_bound_failed", "x"),
        ("8HTB", "c", 3, "in", None, "failed", "x"),
    ])
    status = list(_readable_transport(table)["Status"])
    assert status[0] == ""
    assert "upper bound" in status[1]
    assert "did not finish" in status[2]


def test_the_machine_names_do_not_reach_the_reader():
    table = frame([("8HTB", "a", 3, "in", 2.6, "", "x")])
    out = _readable_transport(table)
    assert "flags" not in out and "source" not in out
    assert "Compound" in out and "Tunnel" in out
