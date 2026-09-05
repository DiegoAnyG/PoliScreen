"""Tests for Pareto Multi-Objective Optimization in PoliScreen.

Verifies that:
1. compute_pareto_ranks correctly identifies non-dominated candidates (Rank 1).
2. Multi-tier Pareto ranking (Rank 1, 2, 3...) accurately reflects hierarchical dominance.
3. Multi-objective trade-offs are preserved (neither dominates if one has better docking and the other better quality).
4. Edge cases (single row, empty, NaNs, duplicate values) are handled robustly.
5. _scatter_dock_inter renders a valid figure with the Pareto frontier line.
6. compute_ranking includes pareto_rank and is_pareto in its output.
"""
import numpy as np
import pandas as pd
import pytest

from poliscreen.core.screening import compute_pareto_ranks, compute_ranking
from poliscreen.ui.components.admet import _scatter_dock_inter


def test_pareto_clear_dominance():
    """A strictly superior compound dominates an inferior one across all objectives."""
    df = pd.DataFrame([
        {"compound": "Leader", "best_dock": -10.5, "inter_quality": 0.95, "confidence": 0.85},
        {"compound": "Inferior", "best_dock": -6.0, "inter_quality": 0.40, "confidence": 0.35},
    ])
    ranks, is_p = compute_pareto_ranks(df, objectives=["best_dock", "inter_quality", "confidence"],
                                       minimize_cols={"best_dock"})
    assert ranks.iloc[0] == 1
    assert is_p.iloc[0] == True
    assert ranks.iloc[1] == 2
    assert is_p.iloc[1] == False


def test_pareto_tradeoff_coexistence():
    """Compounds excelling in different objectives both belong to the Pareto frontier (Rank 1)."""
    df = pd.DataFrame([
        # Compound A has superior affinity (-12.0 kcal/mol) but lower interaction quality
        {"compound": "Affinity_Champ", "best_dock": -12.0, "inter_quality": 0.50, "confidence": 0.70},
        # Compound B has superior interaction quality (1.10) but lower affinity (-6.5 kcal/mol)
        {"compound": "Quality_Champ", "best_dock": -6.5, "inter_quality": 1.10, "confidence": 0.70},
        # Compound C is dominated by both (-5.0, 0.30, 0.40)
        {"compound": "Dominated", "best_dock": -5.0, "inter_quality": 0.30, "confidence": 0.40},
    ])
    ranks, is_p = compute_pareto_ranks(df, objectives=["best_dock", "inter_quality", "confidence"],
                                       minimize_cols={"best_dock"})
    assert ranks.iloc[0] == 1 and is_p.iloc[0] == True
    assert ranks.iloc[1] == 1 and is_p.iloc[1] == True
    assert ranks.iloc[2] == 2 and is_p.iloc[2] == False


def test_pareto_multitier_hierarchy():
    """Multi-tier sorting generates sequential Pareto fronts."""
    df = pd.DataFrame([
        {"compound": "Front1", "best_dock": -10.0, "inter_quality": 0.90, "confidence": 0.90},
        {"compound": "Front2", "best_dock": -8.0, "inter_quality": 0.70, "confidence": 0.70},
        {"compound": "Front3", "best_dock": -6.0, "inter_quality": 0.50, "confidence": 0.50},
    ])
    ranks, is_p = compute_pareto_ranks(df, objectives=["best_dock", "inter_quality", "confidence"],
                                       minimize_cols={"best_dock"})
    assert list(ranks) == [1, 2, 3]
    assert list(is_p) == [True, False, False]


def test_pareto_identical_duplicates():
    """Exact ties in all objectives belong to the same Pareto rank without error."""
    df = pd.DataFrame([
        {"compound": "TwinA", "best_dock": -9.0, "inter_quality": 0.80, "confidence": 0.75},
        {"compound": "TwinB", "best_dock": -9.0, "inter_quality": 0.80, "confidence": 0.75},
    ])
    ranks, is_p = compute_pareto_ranks(df, objectives=["best_dock", "inter_quality", "confidence"],
                                       minimize_cols={"best_dock"})
    assert list(ranks) == [1, 1]
    assert list(is_p) == [True, True]


def test_pareto_nan_handling():
    """Rows with NaN in objectives cannot falsely dominate valid rows."""
    df = pd.DataFrame([
        {"compound": "Valid", "best_dock": -8.0, "inter_quality": 0.75, "confidence": 0.60},
        {"compound": "HasNaN", "best_dock": -12.0, "inter_quality": np.nan, "confidence": 0.90},
    ])
    ranks, is_p = compute_pareto_ranks(df, objectives=["best_dock", "inter_quality", "confidence"],
                                       minimize_cols={"best_dock"})
    assert ranks.iloc[0] == 1
    assert is_p.iloc[0] == True
    # HasNaN must not be considered valid Pareto leader
    assert is_p.iloc[1] == False


def test_pareto_empty_and_single():
    """Handles empty dataframes and single candidate gracefully."""
    empty_df = pd.DataFrame()
    r_empty, p_empty = compute_pareto_ranks(empty_df, ["best_dock", "inter_quality"])
    assert r_empty.empty and p_empty.empty

    single = pd.DataFrame([{"best_dock": -7.5, "inter_quality": 0.65}])
    r_single, p_single = compute_pareto_ranks(single, ["best_dock", "inter_quality"], {"best_dock"})
    assert r_single.iloc[0] == 1
    assert p_single.iloc[0] == True


def test_scatter_dock_inter_renders_pareto_frontier():
    """_scatter_dock_inter generates a publication figure with the Pareto frontier line."""
    df = pd.DataFrame([
        {"compound": "Ctrl", "best_dock": -8.5, "best_inter": 1.0, "is_control": 1, "is_pareto": True},
        {"compound": "C1", "best_dock": -9.5, "best_inter": 0.85, "is_control": 0, "is_pareto": True},
        {"compound": "C2", "best_dock": -7.0, "best_inter": 0.95, "is_control": 0, "is_pareto": True},
        {"compound": "C3", "best_dock": -6.0, "best_inter": 0.40, "is_control": 0, "is_pareto": False},
    ])
    fig = _scatter_dock_inter(df)
    assert fig is not None
    ax = fig.axes[0]
    lines = [line.get_label() for line in ax.get_lines()]
    assert any("Pareto" in str(l) for l in lines)


def test_compute_ranking_includes_pareto_columns():
    """compute_ranking returns pareto_rank and is_pareto columns."""
    inter = pd.DataFrame([
        {"name": "Complex_C1_rec_8HTB-model1", "compound": "C1", "ckey": "c1", "receptor": "8HTB",
         "is_control": 0, "Tyr157_hbond": 1},
        {"name": "Complex_ctrl_rec_8HTB-model1", "compound": "ctrl", "ckey": "ctrl", "receptor": "8HTB",
         "is_control": 1, "Tyr157_hbond": 1},
    ])
    dc = pd.DataFrame([
        {"compound": "C1", "ckey": "c1", "receptor": "8HTB", "docking_score": -9.0},
        {"compound": "ctrl", "ckey": "ctrl", "receptor": "8HTB", "docking_score": -8.0},
    ])
    ref_info = {"8HTB": {"feats": ["Tyr157_hbond"], "src": "crystal", "ckey": "ctrl"}}
    weights = {"dock": 0.4, "inter": 0.6}
    icols = ["Tyr157_hbond"]
    dscore = {"c1_model1": -9.0, "ctrl_model1": -8.0}
    cat_map = {"8HTB": ["Tyr157"]}

    rk = compute_ranking(inter=inter, dc=dc, control_keys={"ctrl"}, control_assign={"ctrl": "8HTB"},
                         ref_info=ref_info, icols=icols, dscore=dscore, cat_map=cat_map, weights=weights)

    assert "pareto_rank" in rk.columns
    assert "is_pareto" in rk.columns
    assert (rk["is_pareto"] == True).any()
