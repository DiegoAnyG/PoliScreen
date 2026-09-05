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
    """Docking vs. interaction quality scatter with Pareto frontier, reference crosshairs,
    superiority shading, and dynamic collision-free label positioning."""
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

    fig, ax = plt.subplots(figsize=(6.8, 4.6), dpi=160)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # Detect control for reference crosshairs
    ctrl_rows = d[d["is_control"] == 1]
    has_ctrl = not ctrl_rows.empty
    if has_ctrl:
        ctrl_bd = float(ctrl_rows["bd"].iloc[0])
        ctrl_bi = float(ctrl_rows["bi"].iloc[0])
        ax.axhline(y=ctrl_bi, color="#ef4444", linestyle=":", linewidth=1.1, alpha=0.45, zorder=1)
        ax.axvline(x=ctrl_bd, color="#ef4444", linestyle=":", linewidth=1.1, alpha=0.45, zorder=1)

    # Sort 2D Pareto optimal points by bd to draw the frontier line
    pareto_pts = d[d["_is_pareto_2d"]].sort_values("bd")
    line_handle = None
    if len(pareto_pts) > 1:
        line_handle, = ax.plot(pareto_pts["bd"], pareto_pts["bi"], color="#2563eb", linestyle="--",
                               linewidth=2.0, alpha=0.9, label=t("Pareto frontier"), zorder=3)

    legend_handles = {}
    if line_handle:
        legend_handles[t("Pareto frontier")] = line_handle

    pts_to_annotate = []
    for _, r in d.iterrows():
        es_ctrl = r.get("is_control") == 1
        es_pareto = bool(r.get("_is_pareto_2d", False))
        bd_val = float(r["bd"])
        bi_val = float(r["bi"])
        name = str(r["compound"])[:15]

        if es_ctrl:
            lbl = t("Control")
            sc = ax.scatter(bd_val, bi_val, s=130, marker="D", c="#dc2626",
                            edgecolors="#7f1d1d", linewidths=1.2, zorder=6)
            if lbl not in legend_handles:
                legend_handles[lbl] = sc
            pts_to_annotate.append({"x": bd_val, "y": bi_val, "text": name, "color": "#991b1b",
                                    "fontsize": 8.0, "fontweight": "bold", "priority": 1,
                                    "badge_bg": "#fee2e2", "badge_ec": "#ef4444"})
        elif es_pareto:
            lbl = t("Pareto optimal")
            sc = ax.scatter(bd_val, bi_val, s=110, marker="o", c="#2563eb",
                            edgecolors="#f59e0b", linewidths=2.2, zorder=5)
            if lbl not in legend_handles:
                legend_handles[lbl] = sc
            pts_to_annotate.append({"x": bd_val, "y": bi_val, "text": name, "color": "#1e3a8a",
                                    "fontsize": 8.0, "fontweight": "bold", "priority": 2,
                                    "badge_bg": "#dbeafe", "badge_ec": "#3b82f6"})
        else:
            lbl = t("Candidate")
            sc = ax.scatter(bd_val, bi_val, s=55, marker="o", c="#10b981",
                            edgecolors="#0f766e", linewidths=0.6, alpha=0.85, zorder=4)
            if lbl not in legend_handles:
                legend_handles[lbl] = sc
            pts_to_annotate.append({"x": bd_val, "y": bi_val, "text": name, "color": "#334155",
                                    "fontsize": 7.0, "fontweight": "normal", "priority": 3,
                                    "badge_bg": "#ffffff", "badge_ec": "#cbd5e1"})

    ax.set_xlabel(t("Docking (kcal/mol; more negative = better)"), fontsize=8.8, fontweight="semibold")
    ax.set_ylabel(t("Interaction quality (0-1 vs. control)"), fontsize=8.8, fontweight="semibold")
    ax.invert_xaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#94a3b8")
    ax.spines["bottom"].set_color("#94a3b8")
    ax.tick_params(colors="#475569", labelsize=8)
    ax.grid(True, linestyle=":", alpha=0.5, color="#cbd5e1")

    # Shading the superiority quadrant (better docking AND better quality than control)
    if has_ctrl:
        x_lims = ax.get_xlim()
        y_lims = ax.get_ylim()
        more_neg_bound = min(x_lims)
        sup_x_start = min(ctrl_bd, more_neg_bound)
        sup_x_end = max(ctrl_bd, more_neg_bound)
        if y_lims[1] > y_lims[0] and ctrl_bi <= y_lims[1]:
            y_norm = max(0.0, min(1.0, (ctrl_bi - y_lims[0]) / (y_lims[1] - y_lims[0])))
            ax.axvspan(sup_x_start, sup_x_end, ymin=y_norm, ymax=1.0,
                       color="#10b981", alpha=0.06, zorder=0)

    # Pre-render canvas to compute exact pixel bounding boxes for dynamic collision avoidance
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    pts_to_annotate.sort(key=lambda p: (p["priority"], -p["y"]))
    placed_bboxes = []
    offsets_pool = [
        (8, 8), (8, -14), (-10, 8), (-10, -14),
        (0, 16), (0, -20), (14, 0), (-20, 0),
        (12, 22), (-14, 22), (12, -26), (-14, -26),
        (20, 10), (-24, 10), (20, -16), (-24, -16)
    ]

    for p in pts_to_annotate:
        pt_x, pt_y = p["x"], p["y"]
        name = p["text"]
        best_off = offsets_pool[0]
        min_overlap = float("inf")

        for off in offsets_pool:
            ann = ax.annotate(name, (pt_x, pt_y), xytext=off, textcoords="offset points",
                              fontsize=p["fontsize"],
                              bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#cbd5e1"))
            bbox = ann.get_window_extent(renderer)

            overlap = 0.0
            for pb in placed_bboxes:
                dx = min(bbox.x1, pb.x1) - max(bbox.x0, pb.x0)
                dy = min(bbox.y1, pb.y1) - max(bbox.y0, pb.y0)
                if dx > 0 and dy > 0:
                    overlap += dx * dy
            overlap += (abs(off[0]) + abs(off[1])) * 0.5

            if overlap < min_overlap:
                min_overlap = overlap
                best_off = off
            ann.remove()
            if overlap == (abs(off[0]) + abs(off[1])) * 0.5:
                break

        final_ann = ax.annotate(name, (pt_x, pt_y), xytext=best_off, textcoords="offset points",
                                fontsize=p["fontsize"], fontweight=p["fontweight"],
                                color=p["color"],
                                bbox=dict(boxstyle="round,pad=0.25,rounding_size=0.2",
                                          fc=p["badge_bg"], ec=p["badge_ec"],
                                          alpha=0.88, lw=0.6),
                                zorder=7)
        fig.canvas.draw()
        placed_bboxes.append(final_ann.get_window_extent(renderer))

    if legend_handles:
        ax.legend(handles=list(legend_handles.values()), labels=list(legend_handles.keys()),
                  frameon=True, facecolor="white", edgecolor="#e2e8f0", fontsize=7.5, loc="best")

    ax.set_title(t("Docking vs. quality · Pareto frontier · ideal: top-right"), fontsize=9.6, fontweight="bold", pad=10)
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

