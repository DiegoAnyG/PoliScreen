"""ADMET UI components: summary tables, scatter plots, and radar charts."""
from __future__ import annotations

import pandas as pd
import streamlit as st

from poliscreen.core import reagents as rg
from poliscreen.core import report as rp
from poliscreen.ui.i18n import t


def _shade(df, col, value="yours", color="rgba(255,205,60,0.20)"):
    """Highlight the rows whose column `col` equals `value` (to mark what the user contributed)."""
    if col not in df.columns:
        return df
    return df.style.apply(lambda r: [f"background-color: {color}" if str(r.get(col)) == value else ""
                                     for _ in r], axis=1)


def _scatter_dock_inter(sub):
    """Docking vs. interaction quality scatter with Pareto frontier envelope."""
    import matplotlib.pyplot as plt
    from poliscreen.core.screening import compute_pareto_ranks
    d = sub.copy()
    d["bd"] = pd.to_numeric(d.get("best_dock"), errors="coerce")
    d["bi"] = pd.to_numeric(d.get("best_inter"), errors="coerce")
    d = d.dropna(subset=["bd", "bi"])
    if d.empty:
        return None

    # Compute 2D Pareto frontier on (bd, bi) for the plot
    _, is_pareto_2d = compute_pareto_ranks(d, objectives=["bd", "bi"], minimize_cols={"bd"})
    d["_is_pareto_2d"] = is_pareto_2d

    fig, ax = plt.subplots(figsize=(6.4, 4.4), dpi=150)

    # Sort 2D Pareto optimal points by bd to draw the frontier line
    pareto_pts = d[d["_is_pareto_2d"]].sort_values("bd")
    if len(pareto_pts) > 1:
        ax.plot(pareto_pts["bd"], pareto_pts["bi"], color="#2563eb", linestyle="--",
                linewidth=1.8, alpha=0.85, label=t("Pareto frontier"), zorder=2)

    legend_handles = {}
    for _, r in d.iterrows():
        es_ctrl = r.get("is_control") == 1
        es_pareto = bool(r.get("_is_pareto_2d", False))

        if es_ctrl:
            lbl = t("Control")
            sc = ax.scatter(r["bd"], r["bi"], s=120, marker="D", c="#dc2626",
                            edgecolors="black", linewidths=0.8, zorder=5)
            if lbl not in legend_handles:
                legend_handles[lbl] = sc
            ax.annotate(str(r["compound"])[:14], (r["bd"], r["bi"]), fontsize=8,
                        fontweight="bold", color="#991b1b", xytext=(4, 4), textcoords="offset points")
        elif es_pareto:
            lbl = t("Pareto optimal")
            sc = ax.scatter(r["bd"], r["bi"], s=100, marker="o", c="#2563eb",
                            edgecolors="#f59e0b", linewidths=1.8, zorder=4)
            if lbl not in legend_handles:
                legend_handles[lbl] = sc
            ax.annotate(str(r["compound"])[:14], (r["bd"], r["bi"]), fontsize=8,
                        fontweight="bold", color="#1e3a8a", xytext=(4, 4), textcoords="offset points")
        else:
            lbl = t("Candidate")
            sc = ax.scatter(r["bd"], r["bi"], s=48, marker="o", c="#10b981",
                            edgecolors="#334155", linewidths=0.5, alpha=0.75, zorder=3)
            if lbl not in legend_handles:
                legend_handles[lbl] = sc
            ax.annotate(str(r["compound"])[:12], (r["bd"], r["bi"]), fontsize=6.5,
                        color="#64748b", alpha=0.7, xytext=(3, 3), textcoords="offset points")

    ax.set_xlabel(t("Docking (kcal/mol; more negative = better)"), fontsize=8.5)
    ax.set_ylabel(t("Interaction quality (0-1 vs. control)"), fontsize=8.5)
    ax.invert_xaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#94a3b8")
    ax.spines["bottom"].set_color("#94a3b8")
    ax.tick_params(colors="#475569", labelsize=8)
    ax.grid(True, linestyle=":", alpha=0.5, color="#cbd5e1")
    ax.legend(frameon=True, facecolor="white", edgecolor="#e2e8f0", fontsize=7.5, loc="lower left")
    ax.set_title(t("Docking vs. quality · Pareto frontier · ideal: top-right"), fontsize=9.5, fontweight="bold", pad=10)
    fig.tight_layout()
    return fig


def _render_adme(admet, items, keyp):
    """items: [(label, smiles)]. Shows a summary table of all + detail per compound."""
    rows_ = []
    for lb, smi in items:
        r = admet.get(rg.inchikey(smi)) or {}
        rows_.append({"compound": lb, "MW": r.get("MW"), "LogP": r.get("LogP"), "QED": r.get("QED"),
                      "LD50 (mg/kg)": r.get("LD50_mg_per_kg"), "GHS": r.get("GHS_category"),
                      "AMES": r.get("AMES"), "hERG": r.get("hERG"), "DILI": r.get("DILI")})
    st.markdown(t("**ADMET summary of all compounds**"))
    if not any(r.get(k) is not None for r in admet.values()
               for k in ("AMES", "hERG", "DILI", "LD50_mg_per_kg")):
        st.info(t("ADMET-AI is not installed on this machine: what you see are the properties "
                  "computed from the structure (MW, LogP, QED), not predicted endpoints. "
                  "docs/INSTALL.md explains how to add it."))
    st.dataframe(pd.DataFrame(rows_), width="stretch", height=min(320, 60 + 34 * len(rows_)))
    st.caption(t("AMES/hERG/DILI = toxicity probability (lower is better). LD50 in mg/kg (higher is better). Predicted on the WHOLE molecule (core + reagent), not the reagent alone."))
    labels = dict(items)
    sel = st.selectbox(t("View detail of"), list(labels), key=f"adme_det_{keyp}")
    row = admet.get(rg.inchikey(labels[sel]))
    if not row:
        return
    ca, cb = st.columns([1, 1])
    ca.pyplot(rp.radar_fig(row, title=sel))
    cb.metric(t("Oral LD50 (mg/kg)"), rp._f(row.get("LD50_mg_per_kg"), 0))
    cb.metric(t("GHS category"), str(row.get("GHS_category") or "-"))
    cb.metric(t("QED"), rp._f(row.get("QED")))
    cb.caption(t("Green = favorable · amber = intermediate · red = unfavorable."))
    cb.caption(t(rp.LD50_NOTICE))
    _col = {"good": "background-color:rgba(46,158,126,0.22)",
            "mid": "background-color:rgba(226,168,44,0.22)",
            "bad": "background-color:rgba(214,70,70,0.22)", "info": ""}
    for title_, fs in rp.sections(row):
        st.markdown(f"**{title_}**")
        dd = pd.DataFrame(fs, columns=["Property", "Value", "v"])
        sty = dd.style.apply(lambda r: [_col.get(r["v"], ""), _col.get(r["v"], ""), ""], axis=1)
        st.dataframe(sty, width="stretch", hide_index=True, column_config={"v": None})

