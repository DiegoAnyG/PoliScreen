"""Geometric Interaction Footprint (Polygon) Modeling and Visualization.

Generates polygon-based pocket interaction diagrams where:
- Pocket amino acid residues form a fixed perimeter ring.
- Candidate compound binding is represented as a geometric polygon hull.
- Crystallographic control footprint is overlaid as a faint dashed reference silhouette.
- Interaction spokes (aristas) are color-coded by PLIP biophysical bond types.
- Catalytic and secondary residues are highlighted with distinctive badges.
- Multi-compound composite reports can be exported as PDF or high-res PNG.
"""
from __future__ import annotations

import io
import re
from typing import Any, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .screening import TYPE_STYLE


def residue_sort_key(res: str) -> tuple[int, str]:
    """Sort residues by their author sequence number, then residue name."""
    m = re.search(r"\d+", str(res))
    num = int(m.group()) if m else 9999
    return (num, str(res))


def extract_pocket_residues(
    inter: pd.DataFrame,
    receptor: str,
    extra_residues: Sequence[str] | None = None,
) -> list[str]:
    """Collect and naturally sort all residues contacted in the pocket of a receptor."""
    if inter is None or inter.empty:
        return sorted(list(set(extra_residues or [])), key=residue_sort_key)

    sub = inter[inter.get("receptor") == receptor] if "receptor" in inter.columns else inter
    if sub.empty:
        sub = inter

    icols = [
        c for c in sub.columns
        if "_" in c and c.rsplit("_", 1)[-1] in TYPE_STYLE
    ]

    contacted: set[str] = set()
    for col in icols:
        try:
            val = pd.to_numeric(sub[col], errors="coerce").fillna(0).sum()
            if val > 0:
                res = col.rsplit("_", 1)[0]
                contacted.add(res)
        except Exception:
            continue

    if extra_residues:
        for r in extra_residues:
            if r and str(r).strip():
                contacted.add(str(r).strip())

    return sorted(list(contacted), key=residue_sort_key)


def parse_contacts_from_row(row: pd.Series | dict) -> dict[str, list[tuple[str, int]]]:
    """Parse interaction contacts from an interacciones.csv row.

    Returns
    -------
    dict[str, list[tuple[str, int]]]
        {residue_name: [(interaction_type, count), ...]}
    """
    contacts: dict[str, list[tuple[str, int]]] = {}
    items = row.items() if hasattr(row, "items") else row.items()
    for col, val in items:
        if "_" not in str(col):
            continue
        parts = str(col).rsplit("_", 1)
        if len(parts) != 2:
            continue
        res, itype = parts[0], parts[1]
        if itype not in TYPE_STYLE:
            continue
        try:
            cnt = int(val)
            if cnt > 0:
                contacts.setdefault(res, []).append((itype, cnt))
        except (ValueError, TypeError):
            continue
    return contacts


def parse_contacts_from_features(features: Sequence[str]) -> dict[str, list[tuple[str, int]]]:
    """Parse interaction contacts from a list of feature strings ('Residue_type')."""
    contacts: dict[str, list[tuple[str, int]]] = {}
    for feat in features or []:
        if "_" not in str(feat):
            continue
        res, itype = str(feat).rsplit("_", 1)
        if itype in TYPE_STYLE:
            contacts.setdefault(res, []).append((itype, 1))
    return contacts


def draw_interaction_polygon(
    all_pocket_residues: Sequence[str],
    compound_contacts: dict[str, list[tuple[str, int]]],
    control_contacts: dict[str, list[tuple[str, int]]] | None = None,
    catalytic_residues: Sequence[str] | None = None,
    secondary_residues: Sequence[str] | None = None,
    title: str = "Compound",
    figsize: tuple[float, float] = (5.5, 5.5),
    dpi: int = 150,
    ax: plt.Axes | None = None,
    show_legend: bool = True,
) -> plt.Figure:
    """Render a single geometric interaction footprint polygon."""
    sorted_res = sorted(list(set(all_pocket_residues)), key=residue_sort_key)
    m = len(sorted_res)

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi, facecolor="#0f172a")
    else:
        fig = ax.figure  # type: ignore

    ax.set_facecolor("#0f172a")
    ax.set_aspect("equal")
    ax.axis("off")

    if m == 0:
        ax.text(0, 0, "No pocket interactions recorded", ha="center", va="center", color="#94a3b8")
        return fig

    catalytic = set(catalytic_residues or [])
    secondary = set(secondary_residues or [])
    ctrl_res = set(control_contacts.keys()) if control_contacts else set()
    cmp_res = set(compound_contacts.keys())

    # Angles clockwise from 90 deg (top)
    angles = np.linspace(90, 90 - 360, m, endpoint=False) * np.pi / 180.0
    res_coord: dict[str, tuple[float, float, float]] = {}
    r_ring = 1.0
    for res, ang in zip(sorted_res, angles):
        res_coord[res] = (r_ring * np.cos(ang), r_ring * np.sin(ang), ang)

    # Concentric guide rings
    for r_guide, alpha_g in [(0.4, 0.12), (0.7, 0.18), (1.0, 0.30)]:
        circle = plt.Circle(
            (0, 0), r_guide, color="#334155", fill=False, ls=":", lw=0.8, alpha=alpha_g, zorder=1
        )
        ax.add_patch(circle)

    # 1. Ghost Control Polygon (Reference silhouette)
    if ctrl_res:
        ctrl_ordered = [r for r in sorted_res if r in ctrl_res]
        if len(ctrl_ordered) >= 3:
            pts = np.array([[res_coord[r][0] * 0.86, res_coord[r][1] * 0.86] for r in ctrl_ordered])
            ctrl_poly = plt.Polygon(
                pts,
                closed=True,
                facecolor="#ef4444",
                edgecolor="#ef4444",
                alpha=0.10,
                ls="--",
                lw=1.5,
                zorder=2,
            )
            ax.add_patch(ctrl_poly)
            ax.plot(
                np.append(pts[:, 0], pts[0, 0]),
                np.append(pts[:, 1], pts[0, 1]),
                color="#ef4444",
                ls="--",
                lw=1.6,
                alpha=0.6,
                zorder=2,
            )
        elif len(ctrl_ordered) == 2:
            r1, r2 = ctrl_ordered
            p1 = [res_coord[r1][0] * 0.86, res_coord[r1][1] * 0.86]
            p2 = [res_coord[r2][0] * 0.86, res_coord[r2][1] * 0.86]
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="#ef4444", ls="--", lw=1.6, alpha=0.6, zorder=2)

    # 2. Candidate Polygon Hull
    cmp_ordered = [r for r in sorted_res if r in cmp_res]
    if len(cmp_ordered) >= 3:
        pts = np.array([[res_coord[r][0] * 0.86, res_coord[r][1] * 0.86] for r in cmp_ordered])
        cand_poly = plt.Polygon(
            pts,
            closed=True,
            facecolor="#2563eb",
            edgecolor="#3b82f6",
            alpha=0.22,
            ls="-",
            lw=2.2,
            zorder=3,
        )
        ax.add_patch(cand_poly)
        ax.plot(
            np.append(pts[:, 0], pts[0, 0]),
            np.append(pts[:, 1], pts[0, 1]),
            color="#60a5fa",
            ls="-",
            lw=2.0,
            alpha=0.9,
            zorder=3,
        )
    elif len(cmp_ordered) == 2:
        r1, r2 = cmp_ordered
        p1 = [res_coord[r1][0] * 0.86, res_coord[r1][1] * 0.86]
        p2 = [res_coord[r2][0] * 0.86, res_coord[r2][1] * 0.86]
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="#60a5fa", ls="-", lw=2.2, alpha=0.9, zorder=3)

    # 3. Draw Interaction Spokes (Aristas from Ligand center to Residues)
    used_types: set[str] = set()
    for res, bonds in compound_contacts.items():
        if res not in res_coord:
            continue
        rx, ry, ang = res_coord[res]
        k_bonds = len(bonds)
        for b_idx, (b_type, count) in enumerate(bonds):
            col, _lbl, ls = TYPE_STYLE.get(b_type, ("#cbd5e1", b_type, "-"))
            used_types.add(b_type)
            offset_ang = (b_idx - (k_bonds - 1) / 2.0) * 0.05
            target_x = 0.82 * (np.cos(ang + offset_ang))
            target_y = 0.82 * (np.sin(ang + offset_ang))
            lw = min(4.0, 1.6 + 0.8 * count)
            ax.plot([0, target_x], [0, target_y], color=col, ls=ls, lw=lw, alpha=0.9, zorder=4)

    # 4. Central Ligand Node
    ax.scatter([0], [0], s=800, c="#1e293b", edgecolors="#3b82f6", linewidths=2.2, zorder=6)
    ax.text(0, 0, "LIG", ha="center", va="center", color="#f8fafc", fontsize=9, fontweight="bold", zorder=7)

    # 5. Pocket Residue Nodes around the Perimeter
    for res in sorted_res:
        rx, ry, _ang = res_coord[res]
        is_hit = res in cmp_res
        is_cat = res in catalytic
        is_sec = res in secondary

        if is_cat:
            node_fc = "#854d0e" if is_hit else "#292524"
            node_ec = "#facc15"
            lw = 2.4
            txt_c = "#fef08a"
        elif is_sec:
            node_fc = "#0e7490" if is_hit else "#1e293b"
            node_ec = "#38bdf8"
            lw = 1.8
            txt_c = "#e0f2fe"
        elif is_hit:
            node_fc = "#1e3a8a"
            node_ec = "#3b82f6"
            lw = 1.6
            txt_c = "#ffffff"
        else:
            node_fc = "#1e293b"
            node_ec = "#475569"
            lw = 1.0
            txt_c = "#94a3b8"

        ax.scatter([rx], [ry], s=650, c=node_fc, edgecolors=node_ec, linewidths=lw, zorder=8)
        ax.text(rx, ry, res, ha="center", va="center", color=txt_c, fontsize=7.5, fontweight="bold", zorder=9)

        if is_cat:
            star_x, star_y = rx * 1.18, ry * 1.18
            ax.text(star_x, star_y, "★", ha="center", va="center", color="#facc15", fontsize=9, zorder=10)

    ax.set_xlim(-1.38, 1.38)
    ax.set_ylim(-1.42, 1.42)
    if title:
        ax.set_title(title, color="#f8fafc", fontsize=11, fontweight="bold", pad=12)

    # Chemical Legend
    if show_legend:
        handles = []
        if ctrl_res:
            handles.append(plt.Line2D([0], [0], color="#ef4444", ls="--", lw=1.8, label="Control ref"))
        if cmp_res:
            handles.append(plt.Line2D([0], [0], color="#3b82f6", ls="-", lw=2.0, label="Candidate hull"))

        for t_key in ["hbond", "saltbridge", "hydrophobic", "pistack", "pication", "halogen", "water"]:
            if t_key in used_types and t_key in TYPE_STYLE:
                col, lbl, ls = TYPE_STYLE[t_key]
                handles.append(plt.Line2D([0], [0], color=col, ls=ls, lw=2.0, label=lbl))

        if handles:
            ax.legend(
                handles=handles,
                loc="lower center",
                bbox_to_anchor=(0.5, -0.12),
                ncol=min(4, len(handles)),
                fontsize=7.5,
                frameon=False,
                labelcolor="#cbd5e1",
            )

    if standalone:
        plt.tight_layout()
    return fig


def draw_interaction_polygon_bytes(
    all_pocket_residues: Sequence[str],
    compound_contacts: dict[str, list[tuple[str, int]]],
    control_contacts: dict[str, list[tuple[str, int]]] | None = None,
    catalytic_residues: Sequence[str] | None = None,
    secondary_residues: Sequence[str] | None = None,
    title: str = "Compound",
    figsize: tuple[float, float] = (5.5, 5.5),
    dpi: int = 180,
    format: str = "png",
) -> bytes:
    """Generate image bytes for an interaction polygon figure."""
    fig = draw_interaction_polygon(
        all_pocket_residues=all_pocket_residues,
        compound_contacts=compound_contacts,
        control_contacts=control_contacts,
        catalytic_residues=catalytic_residues,
        secondary_residues=secondary_residues,
        title=title,
        figsize=figsize,
        dpi=dpi,
        show_legend=True,
    )
    buf = io.BytesIO()
    try:
        fig.savefig(buf, format=format, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
        return buf.getvalue()
    finally:
        plt.close(fig)


def generate_top_interactions_report(
    compounds_data: list[dict[str, Any]],
    all_pocket_residues: Sequence[str],
    control_contacts: dict[str, list[tuple[str, int]]] | None,
    catalytic_residues: Sequence[str] | None,
    secondary_residues: Sequence[str] | None,
    target_name: str = "Receptor",
    format: str = "pdf",
    dpi: int = 250,
) -> bytes:
    """Generate a multi-panel composite TOP interaction report (PDF or PNG).

    Parameters
    ----------
    compounds_data : list[dict]
        Each dict has: 'rank', 'name', 'eff', 'dock', 'quality', 'cat_cov', 'contacts'
    """
    n = len(compounds_data)
    if n == 0:
        fig, ax = plt.subplots(figsize=(6, 4), facecolor="#0f172a")
        ax.text(0.5, 0.5, "No compounds to display", ha="center", va="center", color="#f8fafc")
        ax.axis("off")
        buf = io.BytesIO()
        fig.savefig(buf, format=format, bbox_inches="tight")
        plt.close(fig)
        return buf.getvalue()

    cols = min(3, n)
    rows = (n + cols - 1) // cols
    fig = plt.figure(
        figsize=(5.2 * cols, 5.6 * rows + 0.9),
        facecolor="#0f172a",
        dpi=dpi,
    )
    fig.suptitle(
        f"PoliScreen · TOP Interaction Footprints · Target: {target_name}",
        color="#f8fafc",
        fontsize=16,
        fontweight="bold",
        y=0.98,
    )

    all_used_types: set[str] = set()

    for idx, item in enumerate(compounds_data):
        ax = fig.add_subplot(rows, cols, idx + 1)
        r_title = (
            f"#{item.get('rank', idx+1)} · {item.get('name', 'Compound')}\n"
            f"Eff: {item.get('eff', 0):.1f}% | Dock: {item.get('dock', 0):.2f} kcal/mol | Quality: {item.get('quality', 0):.3f}"
        )
        draw_interaction_polygon(
            all_pocket_residues=all_pocket_residues,
            compound_contacts=item.get("contacts", {}),
            control_contacts=control_contacts,
            catalytic_residues=catalytic_residues,
            secondary_residues=secondary_residues,
            title=r_title,
            ax=ax,
            show_legend=False,
        )
        for _res, bonds in item.get("contacts", {}).items():
            for b_type, _ in bonds:
                all_used_types.add(b_type)

    # Global bottom legend
    handles = [
        plt.Line2D([0], [0], color="#ef4444", ls="--", lw=2.0, label="Control ref"),
        plt.Line2D([0], [0], color="#3b82f6", ls="-", lw=2.2, label="Candidate hull"),
    ]
    for t_key in ["hbond", "saltbridge", "hydrophobic", "pistack", "pication", "halogen", "water"]:
        if t_key in all_used_types and t_key in TYPE_STYLE:
            col, lbl, ls = TYPE_STYLE[t_key]
            handles.append(plt.Line2D([0], [0], color=col, ls=ls, lw=2.0, label=lbl))

    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=min(6, len(handles)),
        fontsize=9,
        frameon=False,
        labelcolor="#cbd5e1",
        bbox_to_anchor=(0.5, 0.01),
    )

    plt.tight_layout(rect=[0, 0.04, 1, 0.96])
    buf = io.BytesIO()
    try:
        fig.savefig(buf, format=format, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
        return buf.getvalue()
    finally:
        plt.close(fig)

