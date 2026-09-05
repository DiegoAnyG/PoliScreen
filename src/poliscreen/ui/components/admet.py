"""ADMET UI components: summary tables, scatter plots, and radar charts."""
from __future__ import annotations

import base64
import json

import pandas as pd
import streamlit as st

from poliscreen.core import reagents as rg
from poliscreen.core import report as rp
from poliscreen.core import screening as sc
from poliscreen.core import viewer as vw
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

    n_total = len(d)
    if n_total > 10:
        cands = d[d["is_control"] != 1].copy()
        cands["_pareto_order"] = cands["_is_pareto_2d"].map({True: 0, False: 1})
        cands["_eff_order"] = -pd.to_numeric(cands.get("effectiveness_pct"), errors="coerce").fillna(0)
        top_candidates = set(cands.sort_values(["_pareto_order", "_eff_order"]).head(3)["compound"])
    else:
        top_candidates = set(d[d["is_control"] != 1]["compound"])

    pts_to_annotate = []
    for _, r in d.iterrows():
        name_raw = str(r["compound"])
        es_ctrl = r.get("is_control") == 1
        es_pareto = bool(r.get("_is_pareto_2d", False))
        bd_val = float(r["bd"])
        bi_val = float(r["bi"])
        eff_val = pd.to_numeric(pd.Series([r.get("effectiveness_pct")]), errors="coerce").iloc[0]
        eff_str = f" ({eff_val:.0f}%)" if pd.notna(eff_val) else ""
        name = f"{name_raw[:16]}{eff_str}"

        if es_ctrl:
            lbl = t("Control")
            sc_pt = ax.scatter(bd_val, bi_val, s=130, marker="D", c="#dc2626",
                               edgecolors="#7f1d1d", linewidths=1.2, zorder=6)
            if lbl not in legend_handles:
                legend_handles[lbl] = sc_pt
            pts_to_annotate.append({"x": bd_val, "y": bi_val, "text": name, "color": "#991b1b",
                                    "fontsize": 8.0, "fontweight": "bold", "priority": 1,
                                    "badge_bg": "#fee2e2", "badge_ec": "#ef4444"})
        elif es_pareto:
            lbl = t("Pareto optimal")
            sc_pt = ax.scatter(bd_val, bi_val, s=110, marker="o", c="#2563eb",
                               edgecolors="#f59e0b", linewidths=2.2, zorder=5)
            if lbl not in legend_handles:
                legend_handles[lbl] = sc_pt
            if name_raw in top_candidates:
                pts_to_annotate.append({"x": bd_val, "y": bi_val, "text": name, "color": "#1e3a8a",
                                        "fontsize": 8.0, "fontweight": "bold", "priority": 2,
                                        "badge_bg": "#dbeafe", "badge_ec": "#3b82f6"})
        else:
            lbl = t("Candidate")
            sc_pt = ax.scatter(bd_val, bi_val, s=55, marker="o", c="#10b981",
                               edgecolors="#0f766e", linewidths=0.6, alpha=0.85, zorder=4)
            if lbl not in legend_handles:
                legend_handles[lbl] = sc_pt
            if name_raw in top_candidates:
                pts_to_annotate.append({"x": bd_val, "y": bi_val, "text": name, "color": "#334155",
                                        "fontsize": 7.0, "fontweight": "normal", "priority": 3,
                                        "badge_bg": "#ffffff", "badge_ec": "#cbd5e1"})

    min_bd, max_bd = float(d["bd"].min()), float(d["bd"].max())
    min_bi, max_bi = float(d["bi"].min()), float(d["bi"].max())
    span_x = max(max_bd - min_bd, 1.2)
    span_y = max(max_bi - min_bi, 0.4)
    pad_x = max(span_x * 0.18, 0.9)
    pad_y = max(span_y * 0.20, 0.15)

    ax.set_xlabel(t("Docking (kcal/mol; more negative = better)"), fontsize=8.8, fontweight="semibold")
    ax.set_ylabel(t("Interaction quality (0-1 vs. control)"), fontsize=8.8, fontweight="semibold")
    ax.set_xlim(max_bd + pad_x, min_bd - pad_x)
    ax.set_ylim(max(0.0, min_bi - pad_y), max_bi + pad_y)
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


def _interactive_pareto_chart(sub, smap=None) -> str:
    """Generates an interactive, responsive HTML/SVG scatter plot with mouse zoom,

    pan, generous edge tolerance, and hover tooltips containing the 2D chemical
    structure image.
    """
    from poliscreen.core.screening import compute_pareto_ranks
    d = sub.copy()
    d["best_dock"] = pd.to_numeric(d.get("best_dock"), errors="coerce")
    d["best_inter"] = pd.to_numeric(d.get("best_inter"), errors="coerce")
    d = d.dropna(subset=["best_dock", "best_inter"])
    if d.empty:
        return ""

    # Compute 2D Pareto frontier on (best_dock, best_inter)
    _, is_pareto_2d = compute_pareto_ranks(d, objectives=["best_dock", "best_inter"], minimize_cols={"best_dock"})
    d["_is_pareto_2d"] = is_pareto_2d

    smap = smap or {}

    ctrl_rows = d[d["is_control"] == 1]
    has_ctrl = not ctrl_rows.empty
    ctrl_bd = float(ctrl_rows["best_dock"].iloc[0]) if has_ctrl else None
    ctrl_bi = float(ctrl_rows["best_inter"].iloc[0]) if has_ctrl else None

    # Calculate domain with generous tolerance/padding (at least 18-20% margin)
    min_bd, max_bd = float(d["best_dock"].min()), float(d["best_dock"].max())
    min_bi, max_bi = float(d["best_inter"].min()), float(d["best_inter"].max())

    span_x = max(max_bd - min_bd, 1.2)
    span_y = max(max_bi - min_bi, 0.4)
    pad_x = max(span_x * 0.18, 0.9)
    pad_y = max(span_y * 0.20, 0.15)

    domain_left = max_bd + pad_x      # less negative (left)
    domain_right = min_bd - pad_x     # more negative (right)
    domain_bottom = max(0.0, min_bi - pad_y)
    domain_top = max_bi + pad_y

    n_total = len(d)
    if n_total > 10:
        cands = d[d["is_control"] != 1].copy()
        cands["_pareto_order"] = cands["_is_pareto_2d"].map({True: 0, False: 1})
        cands["_eff_order"] = -pd.to_numeric(cands.get("effectiveness_pct"), errors="coerce").fillna(0)
        top_candidates = set(cands.sort_values(["_pareto_order", "_eff_order"]).head(3)["compound"])
    else:
        top_candidates = set(d[d["is_control"] != 1]["compound"])

    points_data = []
    for _, r in d.iterrows():
        name = str(r["compound"])
        es_ctrl = bool(r.get("is_control") == 1)
        es_pareto = bool(r.get("_is_pareto_2d", False))
        bd = float(r["best_dock"])
        bi = float(r["best_inter"])
        eff = float(r["effectiveness_pct"]) if pd.notna(r.get("effectiveness_pct")) else None
        conf = float(r["confidence"]) if pd.notna(r.get("confidence")) else None
        prank = int(r["pareto_rank"]) if pd.notna(r.get("pareto_rank")) else None

        show_label = es_ctrl or (name in top_candidates)
        label_text = f"{name[:14]} ({eff:.0f}%)" if (show_label and eff is not None) else (name[:14] if show_label else "")

        # 2D structure image in base64
        smi = smap.get(sc.normalize_key(name)) or ""
        img_b64 = ""
        if smi:
            try:
                png = vw.molecule_png(smi, size=140)
                if png:
                    img_b64 = f"data:image/png;base64,{base64.b64encode(png).decode('ascii')}"
            except Exception:
                img_b64 = ""

        points_data.append({
            "name": name,
            "bd": bd,
            "bi": bi,
            "eff": eff,
            "conf": conf,
            "prank": prank,
            "is_ctrl": es_ctrl,
            "is_pareto": es_pareto,
            "category": t("Control") if es_ctrl else (t("Pareto optimal") if es_pareto else t("Candidate")),
            "show_label": show_label,
            "label_text": label_text,
            "smiles": smi,
            "img": img_b64,
        })

    pareto_pts = [p for p in points_data if p["is_ctrl"] or p["is_pareto"]]
    pareto_pts.sort(key=lambda p: -p["bd"])

    chart_payload = json.dumps({
        "points": points_data,
        "pareto_line": [{"bd": p["bd"], "bi": p["bi"]} for p in pareto_pts],
        "ctrl": {"bd": ctrl_bd, "bi": ctrl_bi} if has_ctrl else None,
        "domain": {
            "left": domain_left,
            "right": domain_right,
            "bottom": domain_bottom,
            "top": domain_top
        },
        "labels": {
            "docking": t("Docking (kcal/mol; more negative = better)"),
            "quality": t("Interaction quality (0-1 vs. control)"),
            "effectiveness": t("Effectiveness (%)"),
            "confidence": t("Confidence"),
            "pareto_rank": t("Pareto rank"),
            "status": t("Status"),
            "reset": t("Reset view"),
        }
    })

    return f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<style>
  html, body {{
    margin: 0; padding: 0; width: 100%; height: 100%; overflow: hidden;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    background: transparent; color: #f1f5f9; user-select: none;
  }}
  #chart-container {{
    position: relative; width: 100%; height: 100%; min-height: 440px;
    background: rgba(15, 23, 42, 0.45); border-radius: 8px; border: 1px solid #334155;
    box-sizing: border-box; overflow: hidden;
  }}
  svg {{ width: 100%; height: 100%; display: block; }}
  .grid-line {{ stroke: #334155; stroke-dasharray: 2,3; stroke-width: 0.8; opacity: 0.55; }}
  .axis-line {{ stroke: #64748b; stroke-width: 1.2; }}
  .axis-text {{ fill: #94a3b8; font-size: 11px; font-weight: 500; }}
  .axis-title {{ fill: #cbd5e1; font-size: 11.5px; font-weight: 600; text-anchor: middle; }}
  .frontier-line {{ stroke: #2563eb; stroke-dasharray: 6,4; stroke-width: 2.2; fill: none; }}
  .ctrl-line {{ stroke: #ef4444; stroke-dasharray: 4,4; stroke-width: 1.2; opacity: 0.75; }}
  .sup-quad {{ fill: #10b981; fill-opacity: 0.08; }}
  
  .pt-marker {{ transition: transform 0.15s ease, r 0.15s ease; cursor: pointer; }}
  .pt-marker:hover {{ filter: drop-shadow(0 0 6px rgba(59, 130, 246, 0.9)); }}
  
  .badge-rect {{ fill: #1e293b; stroke: #475569; stroke-width: 0.8; opacity: 0.9; rx: 4; ry: 4; }}
  .badge-text {{ font-size: 10px; font-weight: 500; fill: #e2e8f0; pointer-events: none; }}
  
  /* Floating Tooltip */
  #tooltip {{
    position: absolute; display: none; z-index: 100; pointer-events: none;
    background: rgba(15, 23, 42, 0.97); border: 1px solid #475569;
    box-shadow: 0 12px 32px rgba(0,0,0,0.65); backdrop-filter: blur(12px);
    border-radius: 8px; padding: 10px 12px; font-size: 11px; color: #f8fafc;
    max-width: 250px;
  }}
  .tt-img {{
    width: 130px; height: 130px; background: #ffffff; border-radius: 6px;
    display: block; margin: 0 auto 8px auto; border: 1px solid #cbd5e1; object-fit: contain;
  }}
  .tt-title {{ font-size: 13px; font-weight: 700; color: #ffffff; margin-bottom: 4px; word-break: break-all; }}
  .tt-tag {{
    display: inline-block; padding: 2px 7px; border-radius: 4px; font-size: 10px;
    font-weight: 600; margin-bottom: 8px;
  }}
  .tag-ctrl {{ background: #fee2e2; color: #991b1b; border: 1px solid #ef4444; }}
  .tag-pareto {{ background: #dbeafe; color: #1e3a8a; border: 1px solid #3b82f6; }}
  .tag-cand {{ background: #d1fae5; color: #065f46; border: 1px solid #10b981; }}
  .tt-row {{ display: flex; justify-content: space-between; margin-bottom: 3px; gap: 12px; }}
  .tt-lbl {{ color: #94a3b8; font-weight: 500; }}
  .tt-val {{ color: #f8fafc; font-weight: 600; }}
  .tt-smi {{ font-size: 9px; color: #64748b; margin-top: 6px; word-break: break-all; font-family: monospace; }}

  /* Controls overlay */
  .chart-btn-bar {{
    position: absolute; top: 10px; right: 12px; display: flex; gap: 6px; z-index: 10;
  }}
  .chart-btn {{
    background: rgba(30, 41, 59, 0.85); border: 1px solid #475569; color: #cbd5e1;
    border-radius: 4px; padding: 4px 8px; font-size: 11px; cursor: pointer;
    backdrop-filter: blur(4px); transition: all 0.2s;
  }}
  .chart-btn:hover {{ background: #334155; color: #ffffff; }}
</style>
</head>
<body>
<div id="chart-container">
  <div class="chart-btn-bar">
    <button class="chart-btn" id="btn-reset">Reset</button>
  </div>
  <svg id="chart-svg"></svg>
  <div id="tooltip"></div>
</div>

<script>
const data = {chart_payload};
const container = document.getElementById("chart-container");
const svg = document.getElementById("chart-svg");
const tooltip = document.getElementById("tooltip");

let width = container.clientWidth;
let height = container.clientHeight || 440;
const margin = {{ top: 32, right: 35, bottom: 46, left: 62 }};

let viewLeft = data.domain.left;
let viewRight = data.domain.right;
let viewBottom = data.domain.bottom;
let viewTop = data.domain.top;

const initBounds = {{ left: viewLeft, right: viewRight, bottom: viewBottom, top: viewTop }};

function scaleX(val) {{
  return margin.left + ((val - viewLeft) / (viewRight - viewLeft)) * (width - margin.left - margin.right);
}}

function scaleY(val) {{
  return height - margin.bottom - ((val - viewBottom) / (viewTop - viewBottom)) * (height - margin.top - margin.bottom);
}}

function render() {{
  width = container.clientWidth;
  height = container.clientHeight || 440;
  svg.innerHTML = "";

  const plotW = width - margin.left - margin.right;
  const plotH = height - margin.top - margin.bottom;

  // Superiority quadrant (better docking AND better quality than control)
  if (data.ctrl) {{
    const cX = scaleX(data.ctrl.bd);
    const cY = scaleY(data.ctrl.bi);
    const qX = Math.min(cX, scaleX(viewRight));
    const qW = Math.abs(scaleX(viewRight) - cX);
    const qH = Math.max(0, cY - margin.top);
    if (qW > 0 && qH > 0) {{
      const rect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
      rect.setAttribute("x", cX);
      rect.setAttribute("y", margin.top);
      rect.setAttribute("width", qW);
      rect.setAttribute("height", qH);
      rect.setAttribute("class", "sup-quad");
      svg.appendChild(rect);
    }}
  }}

  // Grid and Ticks X
  const nTicksX = Math.max(4, Math.floor(plotW / 85));
  const stepX = (viewRight - viewLeft) / nTicksX;
  for (let i = 0; i <= nTicksX; i++) {{
    const v = viewLeft + i * stepX;
    const x = scaleX(v);
    if (x >= margin.left && x <= width - margin.right) {{
      const gl = document.createElementNS("http://www.w3.org/2000/svg", "line");
      gl.setAttribute("x1", x); gl.setAttribute("x2", x);
      gl.setAttribute("y1", margin.top); gl.setAttribute("y2", height - margin.bottom);
      gl.setAttribute("class", "grid-line");
      svg.appendChild(gl);

      const txt = document.createElementNS("http://www.w3.org/2000/svg", "text");
      txt.setAttribute("x", x); txt.setAttribute("y", height - margin.bottom + 16);
      txt.setAttribute("class", "axis-text"); txt.setAttribute("text-anchor", "middle");
      txt.textContent = v.toFixed(1);
      svg.appendChild(txt);
    }}
  }}

  // Grid and Ticks Y
  const nTicksY = Math.max(3, Math.floor(plotH / 60));
  const stepY = (viewTop - viewBottom) / nTicksY;
  for (let i = 0; i <= nTicksY; i++) {{
    const v = viewBottom + i * stepY;
    const y = scaleY(v);
    if (y <= height - margin.bottom && y >= margin.top) {{
      const gl = document.createElementNS("http://www.w3.org/2000/svg", "line");
      gl.setAttribute("x1", margin.left); gl.setAttribute("x2", width - margin.right);
      gl.setAttribute("y1", y); gl.setAttribute("y2", y);
      gl.setAttribute("class", "grid-line");
      svg.appendChild(gl);

      const txt = document.createElementNS("http://www.w3.org/2000/svg", "text");
      txt.setAttribute("x", margin.left - 8); txt.setAttribute("y", y + 4);
      txt.setAttribute("class", "axis-text"); txt.setAttribute("text-anchor", "end");
      txt.textContent = v.toFixed(2);
      svg.appendChild(txt);
    }}
  }}

  // Axis Lines
  const axX = document.createElementNS("http://www.w3.org/2000/svg", "line");
  axX.setAttribute("x1", margin.left); axX.setAttribute("x2", width - margin.right);
  axX.setAttribute("y1", height - margin.bottom); axX.setAttribute("y2", height - margin.bottom);
  axX.setAttribute("class", "axis-line");
  svg.appendChild(axX);

  const axY = document.createElementNS("http://www.w3.org/2000/svg", "line");
  axY.setAttribute("x1", margin.left); axY.setAttribute("x2", margin.left);
  axY.setAttribute("y1", margin.top); axY.setAttribute("y2", height - margin.bottom);
  axY.setAttribute("class", "axis-line");
  svg.appendChild(axY);

  // Axis Titles
  const titleX = document.createElementNS("http://www.w3.org/2000/svg", "text");
  titleX.setAttribute("x", margin.left + plotW / 2);
  titleX.setAttribute("y", height - 10);
  titleX.setAttribute("class", "axis-title");
  titleX.textContent = data.labels.docking;
  svg.appendChild(titleX);

  const titleY = document.createElementNS("http://www.w3.org/2000/svg", "text");
  titleY.setAttribute("x", -(margin.top + plotH / 2));
  titleY.setAttribute("y", 18);
  titleY.setAttribute("transform", "rotate(-90)");
  titleY.setAttribute("class", "axis-title");
  titleY.textContent = data.labels.quality;
  svg.appendChild(titleY);

  // Control reference crosshairs
  if (data.ctrl) {{
    const cx = scaleX(data.ctrl.bd);
    const cy = scaleY(data.ctrl.bi);
    const hL = document.createElementNS("http://www.w3.org/2000/svg", "line");
    hL.setAttribute("x1", margin.left); hL.setAttribute("x2", width - margin.right);
    hL.setAttribute("y1", cy); hL.setAttribute("y2", cy);
    hL.setAttribute("class", "ctrl-line");
    svg.appendChild(hL);

    const vL = document.createElementNS("http://www.w3.org/2000/svg", "line");
    vL.setAttribute("x1", cx); vL.setAttribute("x2", cx);
    vL.setAttribute("y1", margin.top); vL.setAttribute("y2", height - margin.bottom);
    vL.setAttribute("class", "ctrl-line");
    svg.appendChild(vL);
  }}

  // Pareto Frontier Line
  if (data.pareto_line && data.pareto_line.length > 1) {{
    let pathD = "";
    data.pareto_line.forEach((p, idx) => {{
      const px = scaleX(p.bd);
      const py = scaleY(p.bi);
      pathD += (idx === 0 ? "M " : "L ") + px + " " + py + " ";
    }});
    const pline = document.createElementNS("http://www.w3.org/2000/svg", "path");
    pline.setAttribute("d", pathD);
    pline.setAttribute("class", "frontier-line");
    svg.appendChild(pline);
  }}

  // Points and Labels
  data.points.forEach(pt => {{
    const px = scaleX(pt.bd);
    const py = scaleY(pt.bi);
    if (px < margin.left - 10 || px > width - margin.right + 10 || py < margin.top - 10 || py > height - margin.bottom + 10) return;

    let shape;
    if (pt.is_ctrl) {{
      shape = document.createElementNS("http://www.w3.org/2000/svg", "polygon");
      const s = 9;
      shape.setAttribute("points", `${{px}},${{py-s}} ${{px+s}},${{py}} ${{px}},${{py+s}} ${{px-s}},${{py}}`);
      shape.setAttribute("fill", "#dc2626");
      shape.setAttribute("stroke", "#7f1d1d");
      shape.setAttribute("stroke-width", "1.5");
    }} else if (pt.is_pareto) {{
      shape = document.createElementNS("http://www.w3.org/2000/svg", "polygon");
      const s = 8;
      shape.setAttribute("points", `${{px}},${{py-s}} ${{px+s}},${{py}} ${{px}},${{py+s}} ${{px-s}},${{py}}`);
      shape.setAttribute("fill", "#2563eb");
      shape.setAttribute("stroke", "#f59e0b");
      shape.setAttribute("stroke-width", "2.2");
    }} else {{
      shape = document.createElementNS("http://www.w3.org/2000/svg", "circle");
      shape.setAttribute("cx", px);
      shape.setAttribute("cy", py);
      shape.setAttribute("r", "6");
      shape.setAttribute("fill", "#10b981");
      shape.setAttribute("stroke", "#0f766e");
      shape.setAttribute("stroke-width", "1.0");
      shape.setAttribute("opacity", "0.88");
    }}
    shape.setAttribute("class", "pt-marker");

    shape.addEventListener("mouseenter", (e) => showTooltip(e, pt));
    shape.addEventListener("mousemove", (e) => positionTooltip(e));
    shape.addEventListener("mouseleave", hideTooltip);
    svg.appendChild(shape);

    // Badges only for top candidates & control
    if (pt.show_label && pt.label_text) {{
      const g = document.createElementNS("http://www.w3.org/2000/svg", "g");
      const offX = 10, offY = -8;
      const txt = document.createElementNS("http://www.w3.org/2000/svg", "text");
      txt.setAttribute("x", px + offX + 4);
      txt.setAttribute("y", py + offY + 11);
      txt.setAttribute("class", "badge-text");
      txt.textContent = pt.label_text;

      svg.appendChild(txt);
      const b = txt.getBBox();
      svg.removeChild(txt);

      const rect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
      rect.setAttribute("x", b.x - 4);
      rect.setAttribute("y", b.y - 2);
      rect.setAttribute("width", b.width + 8);
      rect.setAttribute("height", b.height + 4);
      rect.setAttribute("class", "badge-rect");
      if (pt.is_ctrl) {{ rect.setAttribute("stroke", "#ef4444"); rect.setAttribute("fill", "#450a0a"); }}
      else if (pt.is_pareto) {{ rect.setAttribute("stroke", "#3b82f6"); rect.setAttribute("fill", "#172554"); }}

      g.appendChild(rect);
      g.appendChild(txt);
      svg.appendChild(g);
    }}
  }});
}}

function showTooltip(e, pt) {{
  let tagClass = "tag-cand";
  if (pt.is_ctrl) tagClass = "tag-ctrl";
  else if (pt.is_pareto) tagClass = "tag-pareto";

  let html = "";
  if (pt.img) {{
    html += `<img class="tt-img" src="${{pt.img}}" alt="Structure" />`;
  }}
  html += `<div class="tt-title">${{pt.name}}</div>`;
  html += `<span class="tt-tag ${{tagClass}}">${{pt.category}}</span>`;
  if (pt.eff !== null) html += `<div class="tt-row"><span class="tt-lbl">${{data.labels.effectiveness}}</span><span class="tt-val">${{pt.eff.toFixed(1)}}%</span></div>`;
  html += `<div class="tt-row"><span class="tt-lbl">Docking</span><span class="tt-val">${{pt.bd.toFixed(2)}} kcal/mol</span></div>`;
  html += `<div class="tt-row"><span class="tt-lbl">Quality</span><span class="tt-val">${{pt.bi.toFixed(3)}}</span></div>`;
  if (pt.conf !== null) html += `<div class="tt-row"><span class="tt-lbl">${{data.labels.confidence}}</span><span class="tt-val">${{pt.conf.toFixed(2)}}</span></div>`;
  if (pt.prank !== null) html += `<div class="tt-row"><span class="tt-lbl">${{data.labels.pareto_rank}}</span><span class="tt-val">${{pt.prank}}</span></div>`;
  if (pt.smiles) html += `<div class="tt-smi">${{pt.smiles}}</div>`;

  tooltip.innerHTML = html;
  tooltip.style.display = "block";
  positionTooltip(e);
}}

function positionTooltip(e) {{
  const rect = container.getBoundingClientRect();
  const mouseX = e.clientX - rect.left;
  const mouseY = e.clientY - rect.top;
  const ttWidth = tooltip.offsetWidth || 240;
  const ttHeight = tooltip.offsetHeight || 300;

  let x = mouseX + 16;
  let y = mouseY - 20;

  if (x + ttWidth > width - 10) {{
    x = mouseX - ttWidth - 16;
  }}
  if (x < 10) x = 10;

  if (y + ttHeight > height - 10) {{
    y = height - ttHeight - 10;
  }}
  if (y < 10) y = 10;

  tooltip.style.left = x + "px";
  tooltip.style.top = y + "px";
}}

function hideTooltip() {{
  tooltip.style.display = "none";
}}

// Zoom & Pan
let isPanning = false;
let startX, startY;

container.addEventListener("wheel", (e) => {{
  e.preventDefault();
  const zoomFactor = e.deltaY > 0 ? 1.12 : 0.88;
  const rect = container.getBoundingClientRect();
  const mouseX = e.clientX - rect.left;
  const mouseY = e.clientY - rect.top;

  const curXVal = viewLeft + ((mouseX - margin.left) / (width - margin.left - margin.right)) * (viewRight - viewLeft);
  const curYVal = viewBottom + ((height - margin.bottom - mouseY) / (height - margin.top - margin.bottom)) * (viewTop - viewBottom);

  viewLeft = curXVal + (viewLeft - curXVal) * zoomFactor;
  viewRight = curXVal + (viewRight - curXVal) * zoomFactor;
  viewBottom = curYVal + (viewBottom - curYVal) * zoomFactor;
  viewTop = curYVal + (viewTop - curYVal) * zoomFactor;

  render();
}});

container.addEventListener("mousedown", (e) => {{
  if (e.target.closest(".chart-btn")) return;
  isPanning = true;
  startX = e.clientX;
  startY = e.clientY;
}});

window.addEventListener("mousemove", (e) => {{
  if (!isPanning) return;
  const dx = e.clientX - startX;
  const dy = e.clientY - startY;
  startX = e.clientX;
  startY = e.clientY;

  const xSpan = viewRight - viewLeft;
  const ySpan = viewTop - viewBottom;
  const plotW = width - margin.left - margin.right;
  const plotH = height - margin.top - margin.bottom;

  const dValX = (dx / plotW) * xSpan;
  const dValY = (dy / plotH) * ySpan;

  viewLeft -= dValX;
  viewRight -= dValX;
  viewBottom += dValY;
  viewTop += dValY;

  render();
}});

window.addEventListener("mouseup", () => {{ isPanning = false; }});

document.getElementById("btn-reset").addEventListener("click", () => {{
  viewLeft = initBounds.left;
  viewRight = initBounds.right;
  viewBottom = initBounds.bottom;
  viewTop = initBounds.top;
  render();
}});

window.addEventListener("resize", render);
render();
</script>
</body>
</html>"""


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

