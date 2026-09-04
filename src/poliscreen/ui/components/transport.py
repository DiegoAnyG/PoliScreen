"""Transport tunnel tables, selection, route preferences, and viewer groups."""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import streamlit as st

from poliscreen.core import caver as cv
from poliscreen.core import tunnels as tn
from poliscreen.core import viewer as vw
from poliscreen.ui.common import _download_table
from poliscreen.ui.i18n import t

# What each column is called on screen. The CSV keeps the machine names, which is what a script
# reading it expects; the table is for a person.
TRANSPORT_LABELS = {
    "receptor": "Receptor", "ligand": "Compound", "tunnel": "Tunnel", "direction": "Direction",
    "E_surface": "E_surface", "E_bound": "E_bound", "E_max": "E_max", "Ea": "Ea", "dE_BS": "dE_BS",
    "n_discs": "Discs", "span_A": "Span (Å)", "tunnel_length_A": "Length (Å)",
    "bottleneck_radius_A": "Bottleneck (Å)", "curvature": "Curvature", "priority": "Priority",
}


STATUS_TEXT = {
    "failed": "did not finish",
    "upper_bound_failed": "upper bound did not pass",
    "lower_bound_only": "lower bound only",
}


def _one_row_per_route(table):
    """One row per receptor, compound, tunnel and direction: the richest calculation of each."""
    if table.empty or "flags" not in table:
        return table
    keys = [c for c in ("receptor", "ligand", "tunnel", "direction") if c in table]
    if not keys:
        return table
    rank = table["flags"].fillna("").map(
        lambda f: 0 if "failed" in f.split() else 1 if "lower_bound_only" in f else 2)
    ordered = table.assign(_rank=rank).sort_values("_rank", ascending=False)
    return ordered.drop_duplicates(subset=keys, keep="first").drop(columns="_rank")


def _readable_transport(table):
    """The same rows, named for reading, with what happened said in words."""
    out = table.copy()

    def status(flags):
        for flag in str(flags or "").split():
            if flag in STATUS_TEXT:
                return t(STATUS_TEXT[flag])
        return ""

    out.insert(0, "Status", [status(f) for f in out["flags"]])
    out = out.drop(columns=[c for c in ("flags", "source") if c in out])
    return out.rename(columns={k: t(v) for k, v in TRANSPORT_LABELS.items()})


def route_preference(table):
    """Which route each compound takes, and by how much it beats its next one."""
    if table.empty or "Ea" not in table:
        return pd.DataFrame()
    usable = table.dropna(subset=["Ea"])
    rows = []
    for ligand, group in usable.groupby("ligand", dropna=False):
        ranked = group.sort_values("Ea")
        best = ranked.iloc[0]
        margin = (ranked.iloc[1]["Ea"] - best["Ea"]) if len(ranked) > 1 else None
        rows.append({
            "ligand": ligand,
            "tunnel": best.get("tunnel"),
            "direction": best.get("direction"),
            "Ea": best["Ea"],
            "dE_BS": best.get("dE_BS"),
            "margin": None if margin is None else round(float(margin), 2),
            "routes": len(ranked),
        })
    return pd.DataFrame(rows).sort_values("Ea") if rows else pd.DataFrame()


def _route_preferences(table):
    """The per-compound answer, above the full table."""
    preference = route_preference(table)
    if preference.empty or preference["routes"].max() < 2:
        return
    with st.expander(t("Which route each compound takes"), expanded=True):
        shown = preference.rename(columns={
            "ligand": t("Compound"), "tunnel": t("Tunnel"), "direction": t("Direction"),
            "margin": t("Beats the next by"), "routes": t("Routes tried")})
        st.dataframe(shown, width="stretch", hide_index=True)
        st.caption(t("Lowest Ea per compound. A margin under about 0.5 kcal/mol is a tie, not a "
                     "preference: that is the order of CaverDock's own repeatability."))


def _tunnel_number(cluster) -> int:
    """The cluster's own number, which is its rank: tun_cl_003 and tun_cl_003_1 are both 3."""
    digits = re.search(r"tun_cl_0*(\d+)", Path(cluster).stem)
    return int(digits.group(1)) if digits else 0


def _tunnel_groups(only=None, root=None):
    """The tunnels found for this project, as sphere groups the viewer can draw."""
    root = root if root is not None else st.session_state.get("tun_drawn")
    if not root or not Path(root).exists():
        return []
    groups_ = []
    for cluster in cv.clusters(root):
        n = _tunnel_number(cluster)
        if only is not None and n not in only:
            continue
        spheres = cv.tunnel_spheres(cluster)
        if spheres:
            groups_.append({"alpha": spheres, "name": cluster.stem, "number": n,
                            "color": vw.TUNNEL_PALETTE[(n - 1) % len(vw.TUNNEL_PALETTE)],
                            "opacity": 1.0})
    return groups_


def _with_tunnel_geometry(table, folder):
    """Join CAVER's own numbers for each tunnel onto the transport rows."""
    if not tn.available() or table.empty:
        return table
    from caver_translate.parse import parse_tunnels

    geometry = {}
    for summary in sorted(Path(folder).rglob("summary.txt")):
        for g in parse_tunnels(summary):
            geometry.setdefault(g.tunnel, g)
    if not geometry:
        return table

    out = table.copy()
    for column, attribute in (("tunnel_length_A", "length"),
                              ("bottleneck_radius_A", "bottleneck_radius"),
                              ("curvature", "curvature"), ("priority", "priority")):
        filled = [getattr(geometry[n], attribute) if n in geometry else existing
                  for n, existing in zip(out["tunnel"], out.get(column, [None] * len(out)))]
        out[column] = filled
    return out


def _tunnel_table(caver_out, found):
    """The tunnels, with a box to pick which ones are drawn and the colour they are drawn in."""
    S = st.session_state
    geometry = {}
    summary = Path(caver_out) / "summary.txt"
    if summary.exists() and tn.available():
        from caver_translate.parse import parse_tunnels
        geometry = {g.tunnel: g for g in parse_tunnels(summary)}

    numbers = [_tunnel_number(c) for c in found]
    widths = [0.6, 1.6, 1.2, 1.2, 1.2, 1.2]
    head = st.columns(widths, vertical_alignment="bottom")
    for col, label in zip(head, ["", t("Tunnel"), t("Bottleneck (Å)"), t("Length (Å)"),
                                 t("Curvature"), t("Priority")]):
        col.caption(label)

    shown = []
    for i, n in enumerate(numbers):
        g = geometry.get(n)
        dot = vw.emoji_for_color(vw.TUNNEL_PALETTE[(n - 1) % len(vw.TUNNEL_PALETTE)])
        S.setdefault(f"tun_draw_{n}", i == 0)
        row = st.columns(widths, vertical_alignment="center")
        if row[0].checkbox(t("Draw"), key=f"tun_draw_{n}", label_visibility="collapsed"):
            shown.append(n)
        row[1].markdown(f"{dot} **{n}**")
        for col, value in zip(row[2:], (g.bottleneck_radius if g else None,
                                        g.length if g else None,
                                        g.curvature if g else None,
                                        g.priority if g else None)):
            col.markdown(f"{value:g}" if isinstance(value, (int, float)) else "—")
    S["tun_shown"] = shown
    st.caption(t("Priority and length are read together: a tunnel with nothing to cross costs "
                 "nothing to cross."))

