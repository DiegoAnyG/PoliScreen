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
    """Docking vs. interaction quality scatter. Shows the trade-off: top-right = good at both."""
    import matplotlib.pyplot as plt
    d = sub.copy()
    d["bd"] = pd.to_numeric(d.get("best_dock"), errors="coerce")
    d["bi"] = pd.to_numeric(d.get("best_inter"), errors="coerce")
    d = d.dropna(subset=["bd", "bi"])
    if d.empty:
        return None
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for _, r in d.iterrows():
        es_ctrl = r.get("is_control") == 1
        ax.scatter(r["bd"], r["bi"], s=110 if es_ctrl else 55,
                   c="#d62728" if es_ctrl else "#1b9e77", edgecolors="black", linewidths=0.6, zorder=3)
        ax.annotate(str(r["compound"])[:14], (r["bd"], r["bi"]), fontsize=7,
                    xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel(t("Docking (kcal/mol; more negative = better)"))
    ax.set_ylabel(t("Interaction quality (0-1 vs. control)"))
    ax.invert_xaxis()
    ax.grid(alpha=0.3)
    ax.set_title(t("Docking vs. quality · red = control · ideal: top-right"))
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

