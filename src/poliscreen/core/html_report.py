"""Standalone, Zero-Server Interactive HTML Report for PoliScreen.

Generates a completely self-contained .html file with:
- Embedded 3Dmol.js for real-time 3D receptor-ligand complex, pocket, and CAVER tunnel visualization.
- Embedded Plotly.js for interactive multi-objective Pareto frontier with rich 2D structure tooltips.
- Interactive geometric interaction polygon footprints (PLIP contact contours).
- Client-side searchable and sortable classification table with CSV and Excel export.
- CAVER transport tunnel cluster metrics and CaverDock kinetics (when present).
- ADMET drug-likeness profile and physicochemical metrics.
- Comprehensive in silico Methods and Software Reproducibility block.
- Target / Pocket switcher for multi-site and multi-target screening runs.

Requires no server, no Python, and no internet connection to open.
"""
from __future__ import annotations

import base64
import html
import io
import json
import os
import re
import zipfile
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np
import pandas as pd

from . import layout as lay
from . import polygon_interaction as poly
from . import report as rp
from . import screening as sc
from . import viewer as vw


def _load_3dmol_js() -> str:
    """Load bundled 3Dmol-min.js from package assets."""
    asset_path = Path(__file__).resolve().parent.parent / "assets" / "3Dmol-min.js"
    if asset_path.is_file():
        try:
            return asset_path.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            pass
    return ""


def _load_plotly_js() -> str:
    """Retrieve plotly.js from python package, bundled assets, or local workspace."""
    try:
        from plotly.offline import get_plotlyjs
        content = get_plotlyjs()
        if len(content) > 1000:
            return content
    except Exception:
        pass

    asset_path = Path(__file__).resolve().parent.parent / "assets" / "plotly.min.js"
    if asset_path.is_file():
        try:
            content = asset_path.read_text(encoding="utf-8", errors="ignore")
            if len(content) > 1000:
                return content
        except Exception:
            pass

    search_paths = [
        Path.home() / ".poliscreen" / "assets" / "plotly.min.js",
        Path.home() / "plotly.min.js",
    ]
    for p in search_paths:
        if p.is_file():
            try:
                content = p.read_text(encoding="utf-8", errors="ignore")
                if len(content) > 1000:
                    return content
            except Exception:
                pass

    return ""


def _find_receptor_pdb(proj: Path, preferred_name: Optional[str] = None) -> tuple[Optional[Path], str]:
    """Locate primary receptor PDB file in the project."""
    rec_dir = lay.artifact(proj, lay.RECEPTORS)
    if not rec_dir.is_dir():
        rec_dir = proj / "receptors"
        if not rec_dir.is_dir():
            rec_dir = proj / "receptores"

    if rec_dir.is_dir():
        if preferred_name:
            clean_pref = sc.base_of(preferred_name)
            for f in rec_dir.glob("*.pdb"):
                if sc.base_of(f.name) == clean_pref or clean_pref in f.name:
                    return f, f.read_text(errors="ignore")

        ready = list(rec_dir.glob("*_ready.pdb"))
        if ready:
            return ready[0], ready[0].read_text(errors="ignore")

        pdbs = list(rec_dir.glob("*.pdb"))
        if pdbs:
            return pdbs[0], pdbs[0].read_text(errors="ignore")

    for f in proj.glob("*.pdb"):
        if "complex" not in f.name.lower() and "pose" not in f.name.lower():
            return f, f.read_text(errors="ignore")

    return None, ""


def _find_pose_pdb(proj: Path, compound_name: str, model_idx: int = 1, receptor: Optional[str] = None) -> str:
    """Locate the ligand pose PDB text with strict pocket tag isolation."""
    ckey = sc.normalize_key(compound_name)
    poses_dir = proj / "poses"

    if poses_dir.is_dir():
        exact_tag = receptor if receptor else ""
        pocket_tag = f"~{receptor.split('~', 1)[1]}" if receptor and "~" in receptor else ""
        clean_rec = sc.base_of(receptor) if receptor else ""

        # Pass 1: Match EXACT target tag (e.g. "8HTB_ready~principal") AND compound AND model
        if exact_tag:
            for pf in sorted(poses_dir.glob("*.pdb")):
                pname = pf.stem
                if exact_tag in pname:
                    mod_num = sc.model_of(pname)
                    cmp_part = sc.compound_from_pose_name(pname)
                    if sc.normalize_key(cmp_part) == ckey and str(mod_num) == str(model_idx):
                        return pf.read_text(errors="ignore")

        # Pass 2: Match pocket tag (e.g. "~principal") if target had a specific pocket
        if pocket_tag:
            for pf in sorted(poses_dir.glob("*.pdb")):
                pname = pf.stem
                if pocket_tag in pname and (not clean_rec or clean_rec in pname):
                    mod_num = sc.model_of(pname)
                    cmp_part = sc.compound_from_pose_name(pname)
                    if sc.normalize_key(cmp_part) == ckey and str(mod_num) == str(model_idx):
                        return pf.read_text(errors="ignore")

        # Pass 3: If NO pocket tag was requested, avoid files that belong to specific pockets (~...)
        if not pocket_tag:
            for pf in sorted(poses_dir.glob("*.pdb")):
                pname = pf.stem
                if "~" in pname:
                    continue
                if clean_rec and clean_rec not in pname:
                    continue
                mod_num = sc.model_of(pname)
                cmp_part = sc.compound_from_pose_name(pname)
                if sc.normalize_key(cmp_part) == ckey and str(mod_num) == str(model_idx):
                    return pf.read_text(errors="ignore")

    cx_dir = lay.artifact(proj, lay.COMPLEXES)
    if cx_dir.is_dir():
        r_key = sc.base_of(receptor) if receptor else ""
        for cf in sorted(cx_dir.glob("*.pdb")):
            pname = lay.strip_complex_prefix(cf.stem)
            if r_key and r_key not in pname:
                continue
            cmp_part = sc.compound_from_pose_name(pname)
            if sc.normalize_key(cmp_part) == ckey:
                return cf.read_text(errors="ignore")

    return ""


def _extract_fpocket_spheres(proj: Path, max_spheres: int = 300) -> list[dict[str, Any]]:
    """Extract cavity alpha-spheres from project prep or pockets folder."""
    spheres: list[dict[str, Any]] = []
    pqr_candidates = list(proj.rglob("pocket*_vert.pqr")) + list(proj.rglob("*.pqr"))
    if not pqr_candidates:
        return []

    target_pqr = pqr_candidates[0]
    try:
        for line in target_pqr.read_text(errors="ignore").splitlines():
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                toks = line.split()
                r = float(toks[-1]) if toks and 0.5 <= float(toks[-1]) <= 5.0 else 1.6
                spheres.append({"x": round(x, 2), "y": round(y, 2), "z": round(z, 2), "r": round(r, 2)})
                if len(spheres) >= max_spheres:
                    break
            except (ValueError, IndexError):
                continue
    except Exception:
        pass
    return spheres


def _extract_fpocket_data_for_target(
    proj: Path,
    target_id: str,
    receptor_pdb: str,
    catalytic_res_nums: list[int],
    max_spheres: int = 400,
) -> tuple[list[dict[str, Any]], list[int], Optional[dict[str, Any]]]:
    """Extract authentic fpocket cavity detection data (alpha spheres, lining residues, and metadata).
    
    Checks:
    1. Any .poliscreen session archive (ui_state.json -> pockets)
    2. Loose fpocket output files (pocket*_vert.pqr, pocket*_atm.pdb)
    3. If fpocket is available on the machine, runs detection on receptor PDB
    """
    pockets_list: list[dict[str, Any]] = []

    # 1. Try reading from .poliscreen session archive
    zip_candidates = list(proj.glob("*.poliscreen")) + list(proj.parent.glob("*.poliscreen"))
    for zf_path in zip_candidates:
        try:
            with zipfile.ZipFile(zf_path, "r") as zf:
                if "ui_state.json" in zf.namelist():
                    st = json.loads(zf.read("ui_state.json").decode("utf-8"))
                    pks = st.get("pockets", {})
                    for rec_k, plist in pks.items():
                        if plist and isinstance(plist, list):
                            rec_base = Path(rec_k.replace("{project}/", "")).stem
                            if rec_base in target_id or target_id.startswith(rec_base) or len(pks) == 1:
                                pockets_list = plist
                                break
                if pockets_list:
                    break
        except Exception:
            pass

    # 2. Check loose fpocket output files
    if not pockets_list:
        pqr_candidates = sorted(proj.rglob("pocket*_vert.pqr"))
        if pqr_candidates:
            from .pockets import _residues_from_atm, _spheres_from_pqr
            for pqr_f in pqr_candidates:
                m_n = re.search(r"pocket(\d+)_vert\.pqr", pqr_f.name)
                n = int(m_n.group(1)) if m_n else 1
                atm_f = pqr_f.parent / f"pocket{n}_atm.pdb"
                sph = _spheres_from_pqr(pqr_f, max_n=max_spheres)
                res = _residues_from_atm(atm_f) if atm_f.exists() else []
                pockets_list.append({
                    "n": n,
                    "alpha_xyz": sph,
                    "residues": res,
                    "druggability": 0.5,
                    "volume": 0.0,
                })

    # 3. If still empty and fpocket is available, run detection on receptor PDB
    if not pockets_list and receptor_pdb:
        try:
            from . import pockets as pk
            if pk.fpocket_available():
                rec_stem = target_id.split("~")[0] if "~" in target_id else "receptor"
                local_pdb = proj / "receptors" / f"{rec_stem}.pdb"
                if not local_pdb.exists():
                    local_pdb = proj / "prep" / f"{rec_stem}.pdb"
                if not local_pdb.exists():
                    import tempfile
                    with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as tmp_pdb:
                        tmp_pdb.write(receptor_pdb)
                        local_pdb = Path(tmp_pdb.name)
                pockets_list = pk.detect(local_pdb, timeout=60)
        except Exception:
            pass

    # 4. Match target_id to the most appropriate pocket
    matched: Optional[dict[str, Any]] = None
    m_pk = re.search(r"(?:pk|pocket)[_\s~-]*(\d+)", target_id, re.IGNORECASE)
    if m_pk:
        target_n = int(m_pk.group(1))
        for p in pockets_list:
            if p.get("n") == target_n:
                matched = p
                break

    if not matched and catalytic_res_nums:
        best_p, best_cov = None, -1
        for p in pockets_list:
            p_nums = {int(re.search(r"\d+", str(r)).group()) for r in p.get("residues", []) if re.search(r"\d+", str(r))}
            cov = len(p_nums & set(catalytic_res_nums))
            if cov > best_cov:
                best_cov = cov
                best_p = p
        if best_p and best_cov > 0:
            matched = best_p

    if not matched and pockets_list:
        matched = pockets_list[0]

    if not matched:
        return [], [], None

    cavity_spheres: list[dict[str, Any]] = []
    for s in matched.get("alpha_xyz", [])[:max_spheres]:
        try:
            cavity_spheres.append({
                "x": round(float(s[0]), 2),
                "y": round(float(s[1]), 2),
                "z": round(float(s[2]), 2),
                "r": round(float(s[3]), 2),
            })
        except (IndexError, ValueError, TypeError):
            continue

    fp_nums: list[int] = []
    for r in matched.get("residues", []):
        m_r = re.search(r"\d+", str(r))
        if m_r:
            fp_nums.append(int(m_r.group()))
    fp_nums = sorted(list(set(fp_nums)))

    pocket_info = {
        "n": matched.get("n"),
        "druggability": matched.get("druggability"),
        "volume": matched.get("volume"),
        "spheres": len(cavity_spheres),
        "residues_count": len(fp_nums),
    }

    return cavity_spheres, fp_nums, pocket_info


def _compute_pocket_centroid(receptor_pdb: str, pocket_res_nums: list[int]) -> Optional[dict[str, float]]:
    """Calculate geometric centroid and bounding radius of pocket residues in the receptor PDB."""
    if not receptor_pdb or not pocket_res_nums:
        return None
    wanted = set(pocket_res_nums)
    coords = []
    for line in receptor_pdb.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                resi = int(line[22:26].strip())
                if resi in wanted:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
            except (ValueError, IndexError):
                continue
    if not coords:
        return None
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    cx = sum(xs) / len(coords)
    cy = sum(ys) / len(coords)
    cz = sum(zs) / len(coords)
    dists = [((c[0] - cx) ** 2 + (c[1] - cy) ** 2 + (c[2] - cz) ** 2) ** 0.5 for c in coords]
    r = max(4.0, min(14.0, float(np.percentile(dists, 75)) if len(dists) > 4 else 8.0))
    return {"x": round(cx, 2), "y": round(cy, 2), "z": round(cz, 2), "radius": round(r, 2)}


TUNNEL_COLORS = [
    "#ef4444",  # Bright Red (Tunnel 1 - matches CAVER)
    "#22c55e",  # Bright Green (Tunnel 2 - matches CAVER)
    "#3b82f6",  # Bright Electric Blue (Tunnel 3 - matches CAVER)
    "#eab308",  # Bright Amber / Yellow (Tunnel 4 - matches CAVER)
    "#f97316",  # Bright Orange (Tunnel 5 - matches CAVER)
    "#a855f7",  # Bright Purple (Tunnel 6 - matches CAVER)
    "#06b6d4",  # Bright Cyan (Tunnel 7)
    "#ec4899",  # Bright Pink (Tunnel 8)
]


def _extract_caver_tunnels(
    proj: Path, preferred_rec: Optional[str] = None
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Extract CAVER tunnel clusters, transport calculations, and full caver-translate report tables."""
    tunnels_dir = proj / lay.TUNNELS
    if not tunnels_dir.is_dir():
        tunnels_dir = proj / "tunnels"
    if not tunnels_dir.is_dir():
        return [], [], {"easiest": [], "routes": [], "flags": {}}

    # 1. Parse CAVER cluster summary
    sum_txt = next(tunnels_dir.rglob("summary.txt"), None)
    cluster_metrics: dict[int, dict[str, float]] = {}
    if sum_txt and sum_txt.is_file():
        try:
            in_tbl = False
            for line in sum_txt.read_text(errors="ignore").splitlines():
                s = line.strip()
                if s.startswith("ID") and "Avg_BR" in s:
                    in_tbl = True
                    continue
                if in_tbl and s.startswith("---"):
                    break
                if in_tbl and s:
                    parts = s.split()
                    if len(parts) >= 12 and parts[0].isdigit():
                        cid = int(parts[0])
                        cluster_metrics[cid] = {
                            "bottleneck": float(parts[3]),
                            "length": float(parts[6]),
                            "curvature": float(parts[8]),
                            "priority": float(parts[10]),
                            "throughput": float(parts[11]),
                        }
        except Exception:
            pass

    # 2. Parse tunnel sphere clusters
    clusters: list[dict[str, Any]] = []
    seen_ids = set()

    cluster_pdbs = sorted(tunnels_dir.rglob("tun_cl_*.pdb"))
    for pf in cluster_pdbs:
        m = re.search(r"tun_cl_(\d+)", pf.name)
        if not m:
            continue
        cid = int(m.group(1))
        if cid in seen_ids:
            continue
        seen_ids.add(cid)

        pts = []
        try:
            for line in pf.read_text(errors="ignore").splitlines():
                if line.startswith(("ATOM", "HETATM")):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    toks = line.split()
                    r = float(toks[-1]) if toks and 0.4 <= float(toks[-1]) <= 10.0 else 1.5
                    pts.append({"x": round(x, 2), "y": round(y, 2), "z": round(z, 2), "r": round(r, 2)})
        except Exception:
            continue

        if pts:
            cm = cluster_metrics.get(cid, {})
            clusters.append({
                "id": cid,
                "name": f"Tunnel {cid}",
                "color": TUNNEL_COLORS[(cid - 1) % len(TUNNEL_COLORS)],
                "spheres": pts,
                "bottleneck": cm.get("bottleneck"),
                "length": cm.get("length"),
                "curvature": cm.get("curvature"),
                "priority": cm.get("priority"),
                "throughput": cm.get("throughput"),
            })

    clusters.sort(key=lambda c: (-(c["priority"] or 0.0), c["id"]))

    # 3. Parse CaverDock transport calculations & full CaverTranslate tables
    transport_runs: list[dict[str, Any]] = []
    caver_report_data: dict[str, Any] = {"easiest": [], "routes": [], "flags": {}}

    try:
        from . import tunnels as tn
        if tn.available():
            tb, _ = tn.read(tunnels_dir)
            if not tb.empty:
                for _, r in tb.iterrows():
                    transport_runs.append({
                        "receptor": str(r.get("receptor", "")),
                        "ligand": str(r.get("ligand", "")),
                        "tunnel": int(r.get("tunnel", 0)),
                        "direction": str(r.get("direction", "in")),
                        "ea": float(r["Ea"]) if pd.notna(r.get("Ea")) else None,
                        "de_bs": float(r["dE_BS"]) if pd.notna(r.get("dE_BS")) else None,
                        "e_surface": float(r["E_surface"]) if pd.notna(r.get("E_surface")) else None,
                        "e_bound": float(r["E_bound"]) if pd.notna(r.get("E_bound")) else None,
                        "e_max": float(r["E_max"]) if pd.notna(r.get("E_max")) else None,
                    })
    except Exception:
        pass

    try:
        import dataclasses
        from caver_translate.parse import scan, parse_tunnels
        from caver_translate.report import rows, by_route, spark, FLAG_TEXT

        sum_candidates = list(tunnels_dir.rglob("*summary.txt"))
        parsed_t = []
        if sum_candidates:
            try:
                parsed_t = parse_tunnels(sum_candidates[0])
            except Exception:
                pass

        t_scan, jobs = scan(tunnels_dir)
        t_final = t_scan if t_scan else parsed_t
        if t_final and preferred_rec:
            clean_r = sc.base_of(preferred_rec)
            t_final = [dataclasses.replace(t, receptor=clean_r) for t in t_final]

        if jobs:
            records = rows(t_final, jobs)
            by_source = {j.source: j for j in jobs}
            ranked = sorted((r for r in records if r.get("Ea") is not None), key=lambda r: r["Ea"])

            easiest_list = []
            for r in ranked:
                j = by_source.get(r.get("source"))
                sp_svg = spark(j.profile, width=220, height=44) if (j and getattr(j, "profile", None)) else ""
                easiest_list.append({
                    "receptor": r.get("receptor", ""),
                    "ligand": r.get("ligand", ""),
                    "tunnel": r.get("tunnel", 0),
                    "direction": r.get("direction", "in"),
                    "ea": r.get("Ea"),
                    "de_bs": r.get("dE_BS"),
                    "e_surface": r.get("E_surface"),
                    "e_max": r.get("E_max"),
                    "e_bound": r.get("E_bound"),
                    "length": r.get("tunnel_length_A"),
                    "neck": r.get("bottleneck_radius_A"),
                    "flags": r["flags"].split() if r.get("flags") else [],
                    "spark": sp_svg,
                })
            caver_report_data["easiest"] = easiest_list

            routes = by_route(records)
            caver_report_data["routes"] = routes

            used_flags = {f for r in records for f in (r["flags"].split() if r.get("flags") else [])}
            caver_report_data["flags"] = {f: FLAG_TEXT.get(f, "") for f in used_flags}
    except Exception:
        pass

    return clusters, transport_runs, caver_report_data


def _build_caver_route_table_html(
    routes: list[dict[str, Any]], group_keys: tuple[str, str], ranked_by: str, headers: list[str]
) -> str:
    """Build formatted CaverTranslate route table with dark styling."""
    if not routes:
        return ""
    try:
        from caver_translate.report import _easiest_first
    except Exception:
        def _easiest_first(r): return r

    keyed: dict[tuple[Any, Any], list[dict[str, Any]]] = {}
    for r in routes:
        keyed.setdefault(tuple(r.get(g) for g in group_keys), []).append(r)

    rows_html = []
    for key in sorted(keyed, key=lambda k: tuple(str(part or "") for part in k)):
        for i, r in enumerate(_easiest_first(keyed[key])):
            head0 = html.escape(str(key[0])) if (i == 0 and key[0]) else ""
            head1 = html.escape(str(key[1])) if (i == 0 and key[1]) else ""
            val_ranked = html.escape(str(r.get(ranked_by) or ""))
            ea_in = f"{r['ea_in']:.2f}" if r.get("ea_in") is not None else "-"
            ea_out = f"{r['ea_out']:.2f}" if r.get("ea_out") is not None else "-"
            de_bs = f"{r['dE_BS']:.2f}" if r.get("dE_BS") is not None else "-"
            if "positive_surface" in r.get("flags", set()):
                de_bs = f"<span class='badge-sec'>{de_bs} (mouth clash)</span>"
            e_surf = f"{r['E_surface']:.2f}" if r.get("E_surface") is not None else "-"
            lng = f"{r['length']:.1f}" if r.get("length") is not None else "-"
            neck = f"{r['neck']:.2f}" if r.get("neck") is not None else "-"
            flags_txt = " / ".join(sorted(r.get("flags", set())))
            flags_html = f"<span class='tag-status tag-candidate'>{html.escape(flags_txt)}</span>" if flags_txt else "-"

            rows_html.append(
                f"<tr>"
                f"<td><strong>{head0}</strong></td>"
                f"<td><strong>{head1}</strong></td>"
                f"<td class='mono'>{val_ranked}</td>"
                f"<td class='mono' style='color:#38bdf8; font-weight:700;'>{ea_in}</td>"
                f"<td class='mono'>{ea_out}</td>"
                f"<td class='mono'>{de_bs}</td>"
                f"<td class='mono'>{e_surf}</td>"
                f"<td class='mono'>{lng}</td>"
                f"<td class='mono'>{neck}</td>"
                f"<td>{flags_html}</td>"
                f"</tr>"
            )

    return f"""
    <div class="table-wrapper" style="margin-top: 0.75rem; margin-bottom: 1.5rem;">
      <table class="data-table">
        <thead>
          <tr>
            <th>{headers[0]}</th>
            <th>{headers[1]}</th>
            <th>{headers[2]}</th>
            <th>Ea in (kcal/mol)</th>
            <th>Ea out (kcal/mol)</th>
            <th>dE_BS (kcal/mol)</th>
            <th>E_surface</th>
            <th>Length (Å)</th>
            <th>Neck (Å)</th>
            <th>Notes</th>
          </tr>
        </thead>
        <tbody>
          {''.join(rows_html)}
        </tbody>
      </table>
    </div>
    """


def _render_tunnels_section(
    caver_clusters: list[dict[str, Any]], caver_report_data: dict[str, Any], L: dict[str, str]
) -> str:
    """Compose semantic HTML section for CAVER transport tunnels and kinetics."""
    has_tunnels = len(caver_clusters) > 0
    caver_rep = caver_report_data or {}
    if not (has_tunnels or caver_rep.get("easiest")):
        return ""

    cl_rows = []
    for c in caver_clusters:
        prio = f"{c['priority']:.3f}" if c.get('priority') is not None else "-"
        btnk = f"{c['bottleneck']:.2f}" if c.get('bottleneck') is not None else "-"
        lng = f"{c['length']:.1f}" if c.get('length') is not None else "-"
        curv = f"{c['curvature']:.2f}" if c.get('curvature') is not None else "-"
        thrp = f"{c['throughput']:.3f}" if c.get('throughput') is not None else "-"
        cl_rows.append(
            f"<tr><td><span style='color:{c['color']}; font-weight:700;'>●</span> {c['name']}</td>"
            f"<td class='mono'>{prio}</td><td class='mono'>{btnk}</td><td class='mono'>{lng}</td>"
            f"<td class='mono'>{curv}</td><td class='mono'>{thrp}</td></tr>"
        )
    clusters_block = ""
    if cl_rows:
        clusters_block = f"""
        <div style="font-size: 0.92rem; font-weight: 700; color: #f8fafc; margin-top: 0.5rem; margin-bottom: 0.5rem;">CAVER Identified Tunnel Clusters</div>
        <div class="table-wrapper" style="margin-bottom: 1.5rem;">
          <table class="data-table">
            <thead>
              <tr>
                <th>{L['th_cluster']}</th>
                <th>{L['th_priority']}</th>
                <th>{L['th_bottleneck']}</th>
                <th>{L['th_length']}</th>
                <th>{L['th_curvature']}</th>
                <th>{L['th_throughput']}</th>
              </tr>
            </thead>
            <tbody>
              {''.join(cl_rows)}
            </tbody>
          </table>
        </div>
        """

    easiest_routes = caver_rep.get("easiest", [])
    easiest_block = ""
    if easiest_routes:
        e_rows = []
        for r in easiest_routes:
            ea_s = f"{r['ea']:.2f}" if r.get("ea") is not None else "-"
            de_s = f"{r['de_bs']:.2f}" if r.get("de_bs") is not None else "-"
            es_s = f"{r['e_surface']:.2f}" if r.get("e_surface") is not None else "-"
            em_s = f"{r['e_max']:.2f}" if r.get("e_max") is not None else "-"
            eb_s = f"{r['e_bound']:.2f}" if r.get("e_bound") is not None else "-"
            l_s = f"{r['length']:.1f}" if r.get("length") is not None else "-"
            n_s = f"{r['neck']:.2f}" if r.get("neck") is not None else "-"
            spark_html = r.get("spark", "")
            spark_cell = f"<div style='color:#38bdf8;'>{spark_html}</div>" if spark_html else "-"
            e_rows.append(
                f"<tr>"
                f"<td><strong>{r['receptor']}</strong></td>"
                f"<td><strong>{r['ligand']}</strong></td>"
                f"<td class='mono'>Tunnel {r['tunnel']}</td>"
                f"<td>{r['direction']}</td>"
                f"<td class='mono' style='color:#38bdf8; font-weight:700;'>{ea_s}</td>"
                f"<td class='mono'>{de_s}</td>"
                f"<td class='mono'>{es_s}</td>"
                f"<td class='mono'>{em_s}</td>"
                f"<td class='mono'>{eb_s}</td>"
                f"<td class='mono'>{l_s}</td>"
                f"<td class='mono'>{n_s}</td>"
                f"<td style='text-align:center;'>{spark_cell}</td>"
                f"</tr>"
            )
        easiest_block = f"""
        <div style="font-size: 0.92rem; font-weight: 700; color: #f8fafc; margin-top: 1rem; margin-bottom: 0.35rem;">Easiest Transport Routes (Ranked by Activation Barrier Ea)</div>
        <div style="font-size: 0.8rem; color: #94a3b8; margin-bottom: 0.65rem;">Ea is what entering costs and compares tunnels. Lower barrier indicates easier site access. Profiles show energy along the path with the barrier obstacle marked in red.</div>
        <div class="table-wrapper" style="margin-bottom: 1.5rem;">
          <table class="data-table">
            <thead>
              <tr>
                <th>Receptor</th>
                <th>Compound</th>
                <th>Tunnel</th>
                <th>Dir</th>
                <th>Ea (kcal/mol)</th>
                <th>dE_BS (kcal/mol)</th>
                <th>E_surface</th>
                <th>E_max</th>
                <th>E_bound</th>
                <th>Length (Å)</th>
                <th>Neck (Å)</th>
                <th>Energy Profile</th>
              </tr>
            </thead>
            <tbody>
              {''.join(e_rows)}
            </tbody>
          </table>
        </div>
        """

    routes = caver_rep.get("routes", [])
    pref_compound_block = ""
    pref_tunnel_block = ""
    if routes:
        t1 = _build_caver_route_table_html(routes, ("receptor", "ligand"), "tunnel", ["Receptor", "Compound", "Tunnel"])
        pref_compound_block = f"""
        <div style="font-size: 0.92rem; font-weight: 700; color: #f8fafc; margin-top: 1rem; margin-bottom: 0.35rem;">Which Tunnel Each Compound Prefers</div>
        <div style="font-size: 0.8rem; color: #94a3b8; margin-bottom: 0.65rem;">Every route one compound has, lowest barrier Ea first.</div>
        {t1}
        """
        t2 = _build_caver_route_table_html(routes, ("receptor", "tunnel"), "ligand", ["Receptor", "Tunnel", "Compound"])
        pref_tunnel_block = f"""
        <div style="font-size: 0.92rem; font-weight: 700; color: #f8fafc; margin-top: 1rem; margin-bottom: 0.35rem;">Which Compounds Prefer Each Tunnel</div>
        <div style="font-size: 0.8rem; color: #94a3b8; margin-bottom: 0.65rem;">The same routes read by tunnel: each tunnel and the compounds that navigate it most easily.</div>
        {t2}
        """

    flags_used = caver_rep.get("flags", {})
    flags_block = ""
    if flags_used:
        f_items = [f"<li><strong>{html.escape(k)}</strong>: {html.escape(v)}</li>" for k, v in sorted(flags_used.items())]
        flags_block = f"""
        <div style="font-size: 0.82rem; color: #94a3b8; background: #0f172a; border: 1px solid #1e293b; border-radius: 6px; padding: 0.75rem 1rem; margin-top: 1rem;">
          <strong style="color: #cbd5e1;">What the marks mean:</strong>
          <ul style="margin: 0.4rem 0 0 1.2rem; padding: 0;">
            {''.join(f_items)}
          </ul>
        </div>
        """

    return f"""
    <div class="panel-card" style="margin-bottom: 2.2rem;">
      <div class="section-head">
        <div class="section-title">{L['sec_tunnel_title']}</div>
        <div class="section-desc">{L['sec_tunnel_desc']}</div>
      </div>
      {clusters_block}
      {easiest_block}
      {pref_compound_block}
      {pref_tunnel_block}
      {flags_block}
    </div>
    """


def _render_admet_section(admet_table_rows: list[dict[str, Any]], L: dict[str, str]) -> str:
    """Compose semantic HTML section for ADMET properties and Lipinski Rule-of-5."""
    if not any(r.get("mw") is not None for r in admet_table_rows):
        return ""

    ad_rows = []
    for r in admet_table_rows:
        tag_class = "tag-control" if r["is_control"] else "tag-candidate"
        mw_s = f"{r['mw']:.1f}" if r.get('mw') is not None else "-"
        logp_s = f"{r['logp']:.2f}" if r.get('logp') is not None else "-"
        hbd_s = str(r['hbd']) if r.get('hbd') is not None else "-"
        hba_s = str(r['hba']) if r.get('hba') is not None else "-"
        sa_s = f"{r['sa_score']:.2f}" if r.get('sa_score') is not None else "-"
        pains_s = str(r['pains']) if r.get('pains') is not None else "-"
        if r.get('pass_ro5') is True:
            ro5_s = '<span style="color:#34d399; font-weight:600;">Pass</span>'
        elif r.get('pass_ro5') is False:
            ro5_s = '<span style="color:#f87171;">Fail</span>'
        else:
            ro5_s = '-'
        ad_rows.append(
            f"<tr><td><strong>{r['name']}</strong></td>"
            f"<td><span class='tag-status {tag_class}'>{r['status']}</span></td>"
            f"<td class='mono'>{mw_s}</td><td class='mono'>{logp_s}</td><td class='mono'>{hbd_s}</td>"
            f"<td class='mono'>{hba_s}</td><td class='mono'>{sa_s}</td><td class='mono'>{pains_s}</td>"
            f"<td>{ro5_s}</td></tr>"
        )
    return f"""
    <div class="panel-card" style="margin-bottom: 2.2rem;">
      <div class="section-head">
        <div class="section-title">{L['sec_admet_title']}</div>
        <div class="section-desc">{L['sec_admet_desc']}</div>
      </div>
      <div class="table-wrapper">
        <table class="data-table">
          <thead>
            <tr>
              <th>{L['th_compound']}</th>
              <th>{L['th_status']}</th>
              <th>{L['th_mw']}</th>
              <th>{L['th_logp']}</th>
              <th>{L['th_hbd']}</th>
              <th>{L['th_hba']}</th>
              <th>{L['th_sascore']}</th>
              <th>{L['th_pains']}</th>
              <th>{L['th_ro5']}</th>
            </tr>
          </thead>
          <tbody>
            {''.join(ad_rows)}
          </tbody>
        </table>
      </div>
    </div>
    """


def _prepare_target_dataset(
    proj_path: Path,
    target_id: str,
    rk_target: pd.DataFrame,
    inter: pd.DataFrame,
    meta: dict[str, Any],
    smap: dict[str, str],
    max_poses: int = 25,
    lang: str = "en",
    L: Optional[dict[str, str]] = None,
) -> dict[str, Any]:
    """Prepare and pre-render complete dataset for a specific target pocket."""
    clean_base = sc.base_of(target_id)
    site_part = target_id.split("~", 1)[1] if "~" in target_id else ""
    display_name = sc.display_name(clean_base) + (f" (site {site_part})" if site_part else "")

    # Ensure Pareto rankings are computed for this target
    if not rk_target.empty and "is_pareto" not in rk_target.columns:
        from .screening import compute_pareto_ranks
        obj_cols = [c for c in ["best_dock", "best_inter"] if c in rk_target.columns]
        if len(obj_cols) == 2:
            pranks, is_p = compute_pareto_ranks(rk_target, objectives=obj_cols, minimize_cols={"best_dock"})
            rk_target["pareto_rank"] = pranks
            rk_target["is_pareto"] = is_p

    rec_file, receptor_pdb = _find_receptor_pdb(proj_path, target_id)

    catalytic_map = meta.get("catalytic", {})
    secondary_map = meta.get("secondary", {})
    catalytic_res = catalytic_map.get(target_id, [])
    secondary_res = secondary_map.get(target_id, [])
    pocket_res_list = meta.get("pocket_residues", {}).get(target_id, [])

    pocket_res_nums = []
    for r_str in pocket_res_list:
        m = re.search(r"\d+", str(r_str))
        if m:
            pocket_res_nums.append(int(m.group()))
    pocket_res_nums = sorted(list(set(pocket_res_nums)))

    catalytic_res_nums = []
    for r_str in catalytic_res:
        m = re.search(r"\d+", str(r_str))
        if m:
            catalytic_res_nums.append(int(m.group()))
    catalytic_res_nums = sorted(list(set(catalytic_res_nums)))

    secondary_res_nums = []
    for r_str in secondary_res:
        m = re.search(r"\d+", str(r_str))
        if m:
            secondary_res_nums.append(int(m.group()))
    secondary_res_nums = sorted(list(set(secondary_res_nums)))

    control_contacts: dict[str, list[tuple[str, int]]] = {}
    crystal_feats = meta.get("crystal_feats", {}).get(target_id, [])
    if crystal_feats:
        control_contacts = poly.parse_contacts_from_features(crystal_feats)

    ctrl_rows = rk_target[rk_target.get("is_control") == 1].copy() if "is_control" in rk_target.columns else pd.DataFrame()
    cand_rows = rk_target[rk_target.get("is_control") != 1].copy() if "is_control" in rk_target.columns else rk_target.copy()

    ref_le = None
    if not ctrl_rows.empty and "LE" in ctrl_rows.columns and pd.notna(ctrl_rows["LE"].iloc[0]):
        try:
            ref_le = float(ctrl_rows["LE"].iloc[0])
        except Exception:
            ref_le = None

    w_dock = float(meta.get("weights", {}).get("dock", 0.5))
    w_inter = float(meta.get("weights", {}).get("inter", 0.5))
    tot_w = (w_dock + w_inter) if (w_dock + w_inter) > 0 else 1.0

    def calc_eff_le(le_val: Optional[float], inter_val: float, is_ctrl: bool) -> float:
        if is_ctrl:
            return 100.0
        if le_val is not None and ref_le is not None and ref_le > 0:
            return float(np.clip(100.0 * (w_dock * (le_val / ref_le) + w_inter * inter_val) / tot_w, 0.0, 150.0))
        return 0.0

    if not cand_rows.empty:
        cand_rows["_eff_val"] = pd.to_numeric(cand_rows.get("effectiveness_pct"), errors="coerce").fillna(0)
        cand_rows["_pareto_val"] = cand_rows.get("is_pareto", False).map({True: 0, False: 1})
        cand_sorted = cand_rows.sort_values(["_pareto_val", "_eff_val"], ascending=[True, False])
    else:
        cand_sorted = pd.DataFrame()

    active_items: list[dict[str, Any]] = []

    if not ctrl_rows.empty:
        crow = ctrl_rows.iloc[0]
        c_name = str(crow["compound"])
        active_items.append({
            "rank": "REF",
            "name": c_name,
            "is_control": True,
            "is_pareto": False,
            "dock": float(crow.get("best_dock") or 0.0),
            "quality": float(crow.get("best_inter") or 1.0),
            "eff": float(crow.get("effectiveness_pct") or 100.0),
            "eff_le": 100.0,
            "conf": float(crow.get("confidence") or 1.0),
            "pki": float(crow.get("pKi") or 0.0) if pd.notna(crow.get("pKi")) else None,
            "le": float(crow.get("LE") or 0.0) if pd.notna(crow.get("LE")) else None,
            "pareto_rank": int(crow.get("pareto_rank") or 1),
            "pose": int(crow.get("pose") or 1),
        })

    rank_idx = 1
    for _, crow in cand_sorted.head(max_poses).iterrows():
        c_name = str(crow["compound"])
        le_v = float(crow.get("LE") or 0.0) if pd.notna(crow.get("LE")) else None
        eff_le_v = calc_eff_le(le_v, float(crow.get("best_inter") or 0.0), is_ctrl=False)
        active_items.append({
            "rank": str(rank_idx),
            "name": c_name,
            "is_control": False,
            "is_pareto": bool(crow.get("is_pareto", False)),
            "dock": float(crow.get("best_dock") or 0.0),
            "quality": float(crow.get("best_inter") or 0.0),
            "eff": float(crow.get("effectiveness_pct") or 0.0),
            "eff_le": eff_le_v,
            "conf": float(crow.get("confidence") or 0.0) if pd.notna(crow.get("confidence")) else 0.0,
            "pki": float(crow.get("pKi") or 0.0) if pd.notna(crow.get("pKi")) else None,
            "le": le_v,
            "pareto_rank": int(crow.get("pareto_rank") or 99),
            "pose": int(crow.get("pose") or 1),
        })
        rank_idx += 1

    poses_pdb_map: dict[str, str] = {}
    polygon_svgs_map: dict[str, str] = {}
    structures_2d_map: dict[str, str] = {}
    contacts_detail_map: dict[str, list[dict[str, Any]]] = {}

    sr = inter[inter["receptor"] == target_id] if ("receptor" in inter.columns and target_id) else inter

    for item in active_items:
        cname = item["name"]
        ckey = sc.normalize_key(cname)

        smi = smap.get(ckey, "")
        if smi:
            try:
                png_bytes = vw.molecule_png(smi, size=180)
                if png_bytes:
                    structures_2d_map[cname] = f"data:image/png;base64,{base64.b64encode(png_bytes).decode('ascii')}"
            except Exception:
                pass

        pose_pdb = _find_pose_pdb(proj_path, cname, model_idx=item["pose"], receptor=target_id)
        if pose_pdb:
            poses_pdb_map[cname] = pose_pdb

        c_contacts = {}
        if not sr.empty:
            scmp = sr[sr["compound"] == cname]
            if scmp.empty:
                scmp = sr[sr["compound"].map(sc.normalize_key) == ckey]
            if not scmp.empty:
                pose_num = str(item.get("pose", 1))
                mod_pose = scmp[scmp["name"].apply(lambda n: sc.model_of(n) == pose_num)]
                row_p = mod_pose.iloc[0] if not mod_pose.empty else scmp.iloc[0]
                c_contacts = poly.parse_contacts_from_row(row_p)

        c_list = []
        for res, clist in c_contacts.items():
            m = re.search(r"\d+", str(res))
            r_num = int(m.group()) if m else 0
            is_c = r_num in catalytic_res_nums
            is_s = r_num in secondary_res_nums
            for ctype, cnt in clist:
                c_list.append({
                    "residue": res,
                    "type": ctype.capitalize(),
                    "count": cnt,
                    "is_cat": is_c,
                    "is_sec": is_s,
                    "color": sc.TYPE_STYLE.get(ctype.lower(), ("#94a3b8", ctype, "-"))[0],
                })
        c_list.sort(key=lambda x: (not x["is_cat"], not x["is_sec"], x["residue"]))
        contacts_detail_map[cname] = c_list

        item["contacts"] = c_contacts
        if item["is_control"] and not control_contacts and c_contacts:
            control_contacts = c_contacts

        try:
            fig = poly.draw_interaction_polygon(
                all_pocket_residues=None,
                compound_contacts=c_contacts,
                control_contacts=control_contacts if not item["is_control"] else None,
                catalytic_residues=catalytic_res,
                secondary_residues=secondary_res,
                title=f"{cname} (Pose {item['pose']})",
                figsize=(4.6, 4.6),
                dpi=120,
            )
            buf = io.StringIO()
            fig.savefig(buf, format="svg", bbox_inches="tight", facecolor=fig.get_facecolor(), edgecolor="none")
            import matplotlib.pyplot as plt
            plt.close(fig)
            polygon_svgs_map[cname] = buf.getvalue()
        except Exception:
            polygon_svgs_map[cname] = "<div class='no-svg'>Interaction polygon could not be rendered.</div>"

    # Extract authentic fpocket cavity detection (alpha spheres, lining residues, pocket metadata)
    cavity_spheres, fp_pocket_res_nums, pocket_info = _extract_fpocket_data_for_target(
        proj_path, target_id, receptor_pdb, catalytic_res_nums
    )
    if fp_pocket_res_nums:
        pocket_res_nums = fp_pocket_res_nums
    else:
        active_pocket_nums = set(catalytic_res_nums) | set(secondary_res_nums)
        for it in active_items:
            for res_name in it.get("contacts", {}).keys():
                m_res = re.search(r"\d+", str(res_name))
                if m_res:
                    active_pocket_nums.add(int(m_res.group()))
        if active_pocket_nums:
            pocket_res_nums = sorted(list(active_pocket_nums))

    pocket_centroid = _compute_pocket_centroid(receptor_pdb, pocket_res_nums)
    caver_clusters, caver_transport, caver_report_data = _extract_caver_tunnels(proj_path, target_id)

    n_total = len(cand_rows)
    best_dock = float(pd.to_numeric(cand_rows["best_dock"], errors="coerce").min()) if not cand_rows.empty and "best_dock" in cand_rows.columns else 0.0
    best_eff_raw = float(pd.to_numeric(cand_rows["effectiveness_pct"], errors="coerce").max()) if not cand_rows.empty and "effectiveness_pct" in cand_rows.columns else 0.0
    best_eff_le = float(max((it["eff_le"] for it in active_items if not it["is_control"]), default=0.0))
    n_pareto = int((cand_rows.get("is_pareto") == True).sum()) if "is_pareto" in cand_rows.columns else 0
    n_beat = int((pd.to_numeric(cand_rows["effectiveness_pct"], errors="coerce").fillna(0) >= 105).sum()) if not cand_rows.empty and "effectiveness_pct" in cand_rows.columns else 0
    n_reliable = int((pd.to_numeric(cand_rows["confidence"], errors="coerce").fillna(0) >= 0.5).sum()) if not cand_rows.empty and "confidence" in cand_rows.columns else 0

    kpis = {
        "total": n_total,
        "best_dock": best_dock,
        "best_eff_raw": best_eff_raw,
        "best_eff_le": best_eff_le,
        "pareto": n_pareto,
        "beat_control": n_beat,
        "reliable": n_reliable,
    }

    full_table_rows: list[dict[str, Any]] = []
    admet_table_rows: list[dict[str, Any]] = []

    for _, r in rk_target.iterrows():
        is_c = bool(r.get("is_control") == 1)
        is_p = bool(r.get("is_pareto", False))
        c_name = str(r.get("compound", ""))
        eff_v = float(r.get("effectiveness_pct")) if pd.notna(r.get("effectiveness_pct")) else 0.0
        dock_v = float(r.get("best_dock")) if pd.notna(r.get("best_dock")) else 0.0
        inter_v = float(r.get("best_inter")) if pd.notna(r.get("best_inter")) else 0.0
        conf_v = float(r.get("confidence")) if pd.notna(r.get("confidence")) else None
        pki_v = float(r.get("pKi")) if pd.notna(r.get("pKi")) else None
        le_v = float(r.get("LE")) if pd.notna(r.get("LE")) else None
        prank_v = int(r.get("pareto_rank")) if pd.notna(r.get("pareto_rank")) else None
        eff_le_v = calc_eff_le(le_v, inter_v, is_ctrl=is_c)

        status_str = "Control" if is_c else ("Pareto" if is_p else "Candidate")
        if lang == "es":
            status_str = "Control" if is_c else ("Pareto" if is_p else "Candidato")

        full_table_rows.append({
            "name": c_name,
            "is_control": is_c,
            "is_pareto": is_p,
            "status": status_str,
            "dock": dock_v,
            "quality": inter_v,
            "eff": eff_v,
            "eff_le": eff_le_v,
            "conf": conf_v,
            "pki": pki_v,
            "le": le_v,
            "pareto_rank": prank_v,
            "has_3d": c_name in poses_pdb_map,
        })

        admet_table_rows.append({
            "name": c_name,
            "is_control": is_c,
            "status": status_str,
            "mw": float(r.get("ro5_mw")) if pd.notna(r.get("ro5_mw")) else None,
            "logp": float(r.get("ro5_logp")) if pd.notna(r.get("ro5_logp")) else (float(r.get("logp")) if pd.notna(r.get("logp")) else None),
            "hbd": int(r.get("ro5_hbd")) if pd.notna(r.get("ro5_hbd")) else None,
            "hba": int(r.get("ro5_hba")) if pd.notna(r.get("ro5_hba")) else None,
            "sa_score": float(r.get("sa_score")) if pd.notna(r.get("sa_score")) else None,
            "pains": int(r.get("pains")) if pd.notna(r.get("pains")) else None,
            "pass_ro5": bool(r.get("ro5_pass") == 1 or r.get("ro5_pass") is True) if pd.notna(r.get("ro5_pass")) else None,
        })

    plotly_points = []
    for row in full_table_rows:
        ck = sc.normalize_key(row["name"])
        plotly_points.append({
            "name": row["name"],
            "x": row["dock"],
            "y": row["quality"],
            "dock": row["dock"],
            "quality": row["quality"],
            "eff": row["eff"],
            "eff_le": row["eff_le"],
            "le": row["le"],
            "conf": row["conf"] if row["conf"] is not None else 0.0,
            "prank": row["pareto_rank"] or 99,
            "is_control": row["is_control"],
            "is_pareto": row["is_pareto"],
            "status": row["status"],
            "img": structures_2d_map.get(row["name"], ""),
            "smiles": smap.get(ck, ""),
        })

    pareto_line = [
        {"x": p["x"], "y": p["y"]}
        for p in sorted(plotly_points, key=lambda pt: pt["x"])
        if p["is_pareto"] or p["is_control"]
    ]
    ctrl_pt = next((p for p in plotly_points if p["is_control"]), None)

    tunnels_html = _render_tunnels_section(caver_clusters, caver_report_data, L or {})
    admet_html = _render_admet_section(admet_table_rows, L or {})

    return {
        "target_id": target_id,
        "rec_name": display_name,
        "receptor_pdb": receptor_pdb,
        "pocket_res_nums": pocket_res_nums,
        "catalytic_res_nums": catalytic_res_nums,
        "secondary_res_nums": secondary_res_nums,
        "pocket_centroid": pocket_centroid,
        "pocket_info": pocket_info,
        "cavity_spheres": cavity_spheres,
        "caver_clusters": caver_clusters,
        "caver_transport": caver_transport,
        "has_tunnels": len(caver_clusters) > 0,
        "has_transport": len(caver_transport) > 0,
        "has_admet": any(r.get("mw") is not None for r in admet_table_rows),
        "tunnels_html": tunnels_html,
        "admet_html": admet_html,
        "kpis": kpis,
        "active_items": active_items,
        "poses_pdb": poses_pdb_map,
        "polygon_svgs": polygon_svgs_map,
        "contacts_detail": contacts_detail_map,
        "structures_2d": structures_2d_map,
        "full_table_rows": full_table_rows,
        "admet_table_rows": admet_table_rows,
        "plotly_data": {
            "points": plotly_points,
            "pareto_line": pareto_line,
            "ctrl": ctrl_pt,
        },
    }


def build_interactive_report(
    proj: str | Path,
    out_path: Optional[str | Path] = None,
    target: Optional[str] = None,
    max_poses: int = 25,
    lang: str = "en",
    all_targets: bool = True,
) -> str:
    """Build a complete, standalone, Zero-Server interactive HTML dossier for a PoliScreen project."""
    proj_path = Path(proj).resolve()
    if not proj_path.is_dir():
        raise FileNotFoundError(f"Project folder not found: {proj_path}")

    # 1. Load Metadata (run.json)
    meta_p = proj_path / "run.json"
    meta: dict[str, Any] = {}
    if meta_p.is_file():
        try:
            meta = json.loads(meta_p.read_text(encoding="utf-8"))
        except Exception:
            pass

    project_title = meta.get("project") or proj_path.name

    # 2. Load Ranking Data
    rk_p = proj_path / "ranking.csv"
    if not rk_p.is_file():
        rk_p = proj_path / "results.csv"
    if not rk_p.is_file():
        rk_p = lay.artifact(proj_path, lay.SUMMARY_CSV)

    if rk_p.is_file():
        rk_all = sc.normalize_columns(pd.read_csv(rk_p))
    else:
        rk_all = pd.DataFrame()

    # 3. Load Interactions Data
    inter_p = lay.artifact(proj_path, lay.INTERACTIONS_CSV)
    if inter_p.is_file():
        inter_raw = pd.read_csv(inter_p)
        ckeys = {sc.normalize_key(Path(c).stem) for c in meta.get("controls", [])}
        inter = sc.prepare_interactions(inter_raw, ckeys)
    else:
        inter = pd.DataFrame()

    # 4. Load Ligands Metadata (SMILES)
    smap: dict[str, str] = {}
    lig_meta_p = proj_path / "ligands_meta.csv"
    if lig_meta_p.is_file():
        try:
            lm_df = pd.read_csv(lig_meta_p)
            if {"name", "smiles"} <= set(lm_df.columns):
                for _, r in lm_df.iterrows():
                    smap[sc.normalize_key(r["name"])] = str(r["smiles"])
        except Exception:
            pass

    for k, v in sc.build_smiles_map(str(lay.artifact(proj_path, lay.RECEPTORS))).items():
        smap.setdefault(k, v)

    # 5. Localization dictionary setup
    i18n = {
        "en": {
            "title": "PoliScreen Virtual Screening Dossier",
            "zero_server": "Zero-Server Standalone",
            "project": "Project",
            "target": "Target",
            "target_pocket": "Target / Pocket",
            "kpi_total": "Evaluated Ligands",
            "kpi_total_desc": "Compounds docked in site",
            "kpi_dock": "Best Affinity (Vina)",
            "kpi_dock_desc": "Lowest binding energy",
            "kpi_eff": "Max Effectiveness",
            "kpi_eff_desc": "Normalized vs. crystal control",
            "kpi_pareto": "Pareto Leaders",
            "kpi_pareto_desc": "Multi-objective non-dominated",
            "kpi_beat": "Beat Control",
            "kpi_beat_desc": "Effectiveness >= 105%",
            "kpi_conf": "High Confidence",
            "kpi_conf_desc": "Consistency score >= 0.50",
            "sec_3d_title": "3D Structural Inspection & PLIP Geometric Footprint",
            "sec_3d_desc": "Interactive 3D molecular visualization (3Dmol.js) with real-time pose switching, active site pocket highlighting, transport tunnels, and PLIP geometric contact footprints.",
            "compound_label": "Compound",
            "btn_surface": "Surface",
            "btn_pocket": "Pocket",
            "btn_pocket_title": "Toggle fpocket cavity alpha spheres and zoom to binding site",
            "btn_tunnels": "Tunnels",
            "btn_spin": "Spin 360",
            "btn_center": "Center View",
            "loading_struct": "Loading structure...",
            "plip_title": "PLIP Interaction Footprint",
            "th_residue": "Residue",
            "th_type": "Interaction Type",
            "th_count": "Count",
            "no_contacts": "No direct PLIP contacts detected.",
            "cat_badge": "Catalytic",
            "sec_badge": "Secondary",
            "sec_pareto_title": "Multi-Objective Pareto Frontier Landscape",
            "sec_pareto_desc": "Hover over any compound to inspect its 2D chemical structure, affinity, and contact metrics. Click any point to load its 3D complex and PLIP polygon above.",
            "sec_table_title": "Ranked Candidates & Virtual Screening Table",
            "search_ph": "Search compound by name...",
            "filter_all": "All",
            "filter_pareto": "Pareto Leaders",
            "filter_beat": "Beat Control (>=105%)",
            "btn_csv": "Export CSV",
            "btn_excel": "Export Excel",
            "th_compound": "Compound",
            "th_status": "Status",
            "th_dock": "Docking (kcal/mol)",
            "th_quality": "Quality",
            "th_eff_raw": "Eff (Raw %)",
            "th_eff_le": "Eff (LE %)",
            "th_conf": "Confidence",
            "th_pki": "pKi",
            "th_le": "LE",
            "th_rank": "Pareto Rank",
            "sec_tunnel_title": "CAVER Transport Tunnels & CaverDock Kinetics",
            "sec_tunnel_desc": "Calculated transport pathways connecting the bulk solvent to the catalytic pocket, including bottleneck radii, curvature, and kinetic barriers.",
            "th_cluster": "Cluster",
            "th_priority": "Priority",
            "th_bottleneck": "Bottleneck Radius (A)",
            "th_length": "Length (A)",
            "th_curvature": "Curvature",
            "th_throughput": "Avg Throughput",
            "th_ligand": "Ligand",
            "th_route": "Tunnel Route",
            "th_dir": "Direction",
            "th_ea": "Ea (Barrier; kcal/mol)",
            "th_de": "dE_BS (Affinity delta; kcal/mol)",
            "sec_admet_title": "ADMET & Drug-Likeness Profile",
            "sec_admet_desc": "Physicochemical properties, Lipinski Rule-of-5 compliance, and synthetic accessibility metrics for candidate molecules.",
            "th_mw": "MW (Da)",
            "th_logp": "LogP",
            "th_hbd": "HBD",
            "th_hba": "HBA",
            "th_sascore": "SAscore",
            "th_pains": "PAINS",
            "th_ro5": "Ro5 Compliance",
            "sec_methods_title": "In Silico Methodology & Reproducibility",
            "sec_methods_desc": "Detailed dock parameters, active site grids, weighting schemas, and exact software versions recorded for scientific publication.",
            "how_to_cite": "How to cite this virtual screening pipeline:",
        },
        "es": {
            "title": "Dossier de Cribado Virtual PoliScreen",
            "zero_server": "Zero-Server Autónomo",
            "project": "Proyecto",
            "target": "Diana",
            "target_pocket": "Diana / Bolsillo",
            "kpi_total": "Ligandos Evaluados",
            "kpi_total_desc": "Moléculas acopladas en sitio",
            "kpi_dock": "Mejor Afinidad (Vina)",
            "kpi_dock_desc": "Mayor energía de unión",
            "kpi_eff": "Máx Efectividad",
            "kpi_eff_desc": "Normalizado vs. ligando control",
            "kpi_pareto": "Líderes de Pareto",
            "kpi_pareto_desc": "Óptimos multi-objetivo",
            "kpi_beat": "Superan al Control",
            "kpi_beat_desc": "Efectividad >= 105%",
            "kpi_conf": "Alta Confianza",
            "kpi_conf_desc": "Puntaje de consistencia >= 0.50",
            "sec_3d_title": "Inspección Estructural 3D y Huella Geométrica PLIP",
            "sec_3d_desc": "Visualización molecular interactiva (3Dmol.js) con conmutación de poses en tiempo real, resaltado de bolsillo catalítico, túneles de transporte y polígonos geométricos de contacto PLIP.",
            "compound_label": "Compuesto",
            "btn_surface": "Superficie",
            "btn_pocket": "Bolsillo",
            "btn_pocket_title": "Visualizar cavidad y esferas alfa de fpocket",
            "btn_tunnels": "Túneles",
            "btn_spin": "Girar 360",
            "btn_center": "Centrar",
            "loading_struct": "Cargando estructura...",
            "plip_title": "Huella Geométrica PLIP",
            "th_residue": "Residuo",
            "th_type": "Tipo de Interacción",
            "th_count": "Conteo",
            "no_contacts": "Sin contactos PLIP detectados.",
            "cat_badge": "Catalítico",
            "sec_badge": "Secundario",
            "sec_pareto_title": "Frontera de Pareto Multi-Objetivo",
            "sec_pareto_desc": "Pasa el cursor sobre cualquier punto para ver su estructura 2D, afinidad y contactos. Haz clic para cargar instantáneamente su complejo 3D arriba.",
            "sec_table_title": "Tabla de Clasificación y Resultados",
            "search_ph": "Buscar compuesto por nombre...",
            "filter_all": "Todos",
            "filter_pareto": "Líderes Pareto",
            "filter_beat": "Superan Control (>=105%)",
            "btn_csv": "Exportar CSV",
            "btn_excel": "Exportar Excel",
            "th_compound": "Compuesto",
            "th_status": "Estado",
            "th_dock": "Docking (kcal/mol)",
            "th_quality": "Calidad",
            "th_eff_raw": "Efect. Cruda (%)",
            "th_eff_le": "Efect. LE (%)",
            "th_conf": "Confianza",
            "th_pki": "pKi",
            "th_le": "LE",
            "th_rank": "Rango Pareto",
            "sec_tunnel_title": "Túneles de Transporte CAVER y Cinética CaverDock",
            "sec_tunnel_desc": "Rutas de transporte calculadas entre el solvente y el bolsillo catalítico, con radios de cuello de botella, curvatura y barreras cinéticas.",
            "th_cluster": "Clúster",
            "th_priority": "Prioridad",
            "th_bottleneck": "Radio Cuello (A)",
            "th_length": "Longitud (A)",
            "th_curvature": "Curvatura",
            "th_throughput": "Rendimiento Promedio",
            "th_ligand": "Ligando",
            "th_route": "Ruta de Túnel",
            "th_dir": "Dirección",
            "th_ea": "Ea (Barrera; kcal/mol)",
            "th_de": "dE_BS (Delta afinidad; kcal/mol)",
            "sec_admet_title": "Perfil ADMET y Propiedades Farmacológicas",
            "sec_admet_desc": "Propiedades fisicoquímicas, cumplimiento de la Regla de 5 de Lipinski y accesibilidad sintética.",
            "th_mw": "PM (Da)",
            "th_logp": "LogP",
            "th_hbd": "HBD",
            "th_hba": "HBA",
            "th_sascore": "SAscore",
            "th_pains": "PAINS",
            "th_ro5": "Cumple Ro5",
            "sec_methods_title": "Metodología In Silico y Reproducibilidad",
            "sec_methods_desc": "Parámetros de docking, cuadrículas de bolsillo, esquemas de ponderación y versiones de software registradas para publicación científica.",
            "how_to_cite": "Cómo citar este cribado virtual:",
        }
    }
    L = i18n.get(lang, i18n["en"])

    # 6. Discover targets
    if not rk_all.empty and "receptor" in rk_all.columns:
        all_targets_list = sorted(
            [str(x) for x in rk_all["receptor"].dropna().unique()],
            key=lambda r: (0 if ("principal" in str(r).lower() or "~" not in str(r)) else 1, str(r))
        )
    elif meta.get("receptors"):
        all_targets_list = [str(r) for r in meta["receptors"]]
    else:
        all_targets_list = []

    if not all_targets_list:
        all_targets_list = [target] if target else ["default"]

    active_target = target if (target and target in all_targets_list) else all_targets_list[0]
    target_keys = all_targets_list if all_targets else [active_target]

    # 7. Prepare target datasets
    targets_data: dict[str, Any] = {}
    for t_id in target_keys:
        if not rk_all.empty and "receptor" in rk_all.columns and t_id in rk_all["receptor"].values:
            sub_rk = rk_all[rk_all["receptor"] == t_id].copy()
        else:
            sub_rk = rk_all.copy()

        targets_data[t_id] = _prepare_target_dataset(
            proj_path=proj_path,
            target_id=t_id,
            rk_target=sub_rk,
            inter=inter,
            meta=meta,
            smap=smap,
            max_poses=max_poses,
            lang=lang,
            L=L,
        )

    # 8. Methods block
    methods_md = rp.methods_text(meta)
    methods_html = _markdown_to_html(methods_md)

    # 9. Load Embedded Scripts
    js_3dmol = _load_3dmol_js()
    js_plotly = _load_plotly_js()

    # 10. Assemble Master Document
    html_content = _render_master_html(
        project_title=project_title,
        active_target=active_target,
        all_targets_keys=target_keys,
        targets_data=targets_data,
        methods_html=methods_html,
        js_3dmol=js_3dmol,
        js_plotly=js_plotly,
        lang=lang,
        L=L,
    )

    if out_path:
        out_f = Path(out_path).resolve()
        out_f.parent.mkdir(parents=True, exist_ok=True)
        out_f.write_text(html_content, encoding="utf-8")

    return html_content


def _markdown_to_html(md: str) -> str:
    """Basic conversion of Markdown methods block to semantic HTML."""
    lines = md.splitlines()
    out: list[str] = []
    in_list = False

    for line in lines:
        s = line.strip()
        if not s:
            if in_list:
                out.append("</ul>")
                in_list = False
            continue

        if s.startswith("### "):
            if in_list:
                out.append("</ul>")
                in_list = False
            out.append(f"<h4>{html.escape(s[4:])}</h4>")
        elif s.startswith("## "):
            if in_list:
                out.append("</ul>")
                in_list = False
            out.append(f"<h3>{html.escape(s[3:])}</h3>")
        elif s.startswith("# "):
            if in_list:
                out.append("</ul>")
                in_list = False
            out.append(f"<h2>{html.escape(s[2:])}</h2>")
        elif s.startswith("- "):
            if not in_list:
                out.append("<ul class='methods-list'>")
                in_list = True
            text = _format_inline_markdown(s[2:])
            out.append(f"<li>{text}</li>")
        else:
            if in_list:
                out.append("</ul>")
                in_list = False
            text = _format_inline_markdown(s)
            out.append(f"<p>{text}</p>")

    if in_list:
        out.append("</ul>")

    return "\n".join(out)


def _format_inline_markdown(text: str) -> str:
    """Format bold and code in a string."""
    escaped = html.escape(text)
    escaped = re.sub(r"\*\*(.*?)\*\*", r"<strong>\1</strong>", escaped)
    escaped = re.sub(r"`(.*?)`", r"<code>\1</code>", escaped)
    escaped = re.sub(r"\*(.*?)\*", r"<em>\1</em>", escaped)
    return escaped


def _render_master_html(
    project_title: str,
    active_target: str,
    all_targets_keys: list[str],
    targets_data: dict[str, Any],
    methods_html: str,
    js_3dmol: str,
    js_plotly: str,
    lang: str = "en",
    L: Optional[dict[str, str]] = None,
) -> str:
    """Compose the complete zero-server standalone HTML document."""
    if not L:
        L = {}

    cur_d = targets_data[active_target]
    kpis = cur_d["kpis"]

    embed_3dmol = f"<script>{js_3dmol}</script>" if js_3dmol else "<script src='https://3Dmol.org/build/3Dmol-min.js'></script>"
    embed_plotly = f"<script>{js_plotly}</script>" if len(js_plotly) > 1000 else "<script src='https://cdn.plot.ly/plotly-2.35.2.min.js'></script>"

    target_options = [(tid, targets_data[tid]["rec_name"]) for tid in all_targets_keys]
    if len(all_targets_keys) > 1:
        target_opts_html = "".join(
            f'<option value="{tid}" {"selected" if tid == active_target else ""}>{html.escape(tname)}</option>'
            for tid, tname in target_options
        )
        target_nav_html = f"""
        <div class="nav-target-wrap">
          <label for="target-selector" class="nav-target-lbl">{L['target_pocket']}:</label>
          <select id="target-selector" class="select-target" onchange="switchTarget(this.value)">
            {target_opts_html}
          </select>
        </div>
        """
    else:
        target_nav_html = f'<span class="nav-target">{html.escape(cur_d["rec_name"])}</span>'

    json_targets_data = json.dumps(targets_data)
    json_all_targets = json.dumps(all_targets_keys)

    return f"""<!DOCTYPE html>
<html lang="{lang}">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=1280, initial-scale=0.35, minimum-scale=0.1, maximum-scale=3.0, user-scalable=yes">
<title>{L['title']} — {html.escape(project_title)}</title>
{embed_3dmol}
{embed_plotly}
<style>
  :root {{
    --bg-page: #080c14;
    --bg-card: #0f172a;
    --bg-card-subtle: #131f37;
    --bg-hover: #1e293b;
    --border: #1e293b;
    --border-light: #334155;
    --text: #f8fafc;
    --text-muted: #94a3b8;
    --text-dim: #64748b;
    --accent: #3b82f6;
    --accent-hover: #2563eb;
    --emerald: #10b981;
    --amber: #f59e0b;
    --rose: #ef4444;
    --cyan: #06b6d4;
    --font-mono: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", monospace;
  }}

  * {{ box-sizing: border-box; margin: 0; padding: 0; }}

  html, body {{
    min-width: 1280px;
    width: 100%;
    background-color: var(--bg-page);
    color: var(--text);
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Inter", "Helvetica Neue", Arial, sans-serif;
    line-height: 1.5;
    overflow-x: auto;
    -webkit-text-size-adjust: 100%;
  }}

  /* Top Navigation Bar */
  .navbar {{
    display: flex;
    justify-content: space-between;
    align-items: center;
    background: #060911;
    border-bottom: 1px solid var(--border);
    padding: 0.85rem 2rem;
    position: sticky;
    top: 0;
    z-index: 1000;
  }}
  .nav-brand {{
    display: flex;
    align-items: center;
    gap: 0.75rem;
  }}
  .nav-logo {{
    font-weight: 800;
    font-size: 1.2rem;
    letter-spacing: -0.03em;
    color: #ffffff;
    display: flex;
    align-items: center;
    gap: 0.45rem;
  }}
  .nav-tag {{
    font-size: 0.72rem;
    font-weight: 600;
    padding: 0.18rem 0.55rem;
    border-radius: 9999px;
    background: rgba(59, 130, 246, 0.15);
    color: #60a5fa;
    border: 1px solid rgba(59, 130, 246, 0.3);
  }}
  .nav-target {{
    font-size: 0.76rem;
    font-weight: 600;
    padding: 0.22rem 0.75rem;
    border-radius: 9999px;
    background: rgba(16, 185, 129, 0.15);
    color: #34d399;
    border: 1px solid rgba(16, 185, 129, 0.35);
  }}
  .nav-target-wrap {{
    display: flex;
    align-items: center;
    gap: 0.45rem;
    margin-left: 0.5rem;
  }}
  .nav-target-lbl {{
    font-size: 0.75rem;
    font-weight: 600;
    color: var(--text-muted);
    text-transform: uppercase;
    letter-spacing: 0.04em;
  }}
  .select-target {{
    background: #0f172a;
    color: #38bdf8;
    border: 1px solid #0284c7;
    padding: 0.3rem 0.75rem;
    border-radius: 6px;
    font-size: 0.82rem;
    font-weight: 700;
    cursor: pointer;
    outline: none;
    transition: all 0.15s ease;
  }}
  .select-target:hover, .select-target:focus {{
    border-color: #38bdf8;
    box-shadow: 0 0 0 2px rgba(56, 189, 248, 0.25);
  }}
  .nav-meta {{
    font-size: 0.84rem;
    color: var(--text-muted);
  }}

  /* Main Container */
  .container {{
    min-width: 1240px;
    max-width: 1560px;
    margin: 0 auto;
    padding: 1.8rem 2rem 4rem 2rem;
  }}

  /* Executive KPI Cards Grid */
  .kpi-grid {{
    display: grid;
    grid-template-columns: repeat(6, 1fr);
    gap: 1rem;
    margin-bottom: 2rem;
  }}
  .kpi-card {{
    background: var(--bg-card);
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 1.15rem 1.25rem;
    display: flex;
    flex-direction: column;
    justify-content: space-between;
    transition: transform 0.15s ease, border-color 0.15s ease;
  }}
  .kpi-card:hover {{
    border-color: var(--accent);
    transform: translateY(-2px);
  }}
  .kpi-label {{
    font-size: 0.76rem;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.05em;
    color: var(--text-muted);
    margin-bottom: 0.35rem;
  }}
  .kpi-val {{
    font-size: 1.55rem;
    font-weight: 700;
    letter-spacing: -0.02em;
    color: #ffffff;
    font-family: var(--font-mono);
  }}
  .kpi-desc {{
    font-size: 0.74rem;
    color: var(--text-dim);
    margin-top: 0.3rem;
  }}

  /* Section Headers */
  .section-head {{
    margin-bottom: 1.1rem;
  }}
  .section-title {{
    font-size: 1.22rem;
    font-weight: 700;
    letter-spacing: -0.02em;
    color: #ffffff;
    margin-bottom: 0.25rem;
  }}
  .section-desc {{
    font-size: 0.84rem;
    color: var(--text-muted);
  }}

  /* Dual Visual Grid (3D Complex & PLIP Polygon) */
  .visual-grid {{
    display: grid;
    grid-template-columns: 1.15fr 0.85fr;
    gap: 1.25rem;
    margin-bottom: 2.2rem;
    align-items: stretch;
  }}

  .panel-card {{
    background: var(--bg-card);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 1.2rem;
    display: flex;
    flex-direction: column;
    overflow: hidden;
  }}

  /* 3D Viewer Toolbar */
  .viewer-toolbar {{
    display: flex;
    flex-wrap: wrap;
    gap: 0.5rem;
    margin-bottom: 0.75rem;
    align-items: center;
  }}
  .select-custom {{
    background: #172554;
    color: #ffffff;
    border: 1px solid var(--border-light);
    padding: 0.45rem 0.85rem;
    border-radius: 6px;
    font-size: 0.82rem;
    font-weight: 600;
    cursor: pointer;
    outline: none;
    max-width: 280px;
  }}
  .btn {{
    background: #1e293b;
    color: #f1f5f9;
    border: 1px solid var(--border-light);
    padding: 0.45rem 0.85rem;
    border-radius: 6px;
    font-size: 0.8rem;
    font-weight: 500;
    cursor: pointer;
    transition: all 0.15s ease;
  }}
  .btn:hover:not(:disabled) {{
    background: #334155;
    border-color: var(--accent);
    color: #ffffff;
  }}
  .btn.active {{
    background: var(--accent);
    border-color: #60a5fa;
    color: #ffffff;
  }}
  .btn:disabled {{
    opacity: 0.45;
    cursor: not-allowed;
  }}

  #viewer-container {{
    width: 100%;
    height: 570px;
    min-height: 540px;
    flex: 1;
    background: #060911;
    border-radius: 6px;
    position: relative;
    border: 1px solid var(--border);
  }}
  .viewer-caption {{
    font-size: 0.78rem;
    color: var(--text-muted);
    margin-top: 0.6rem;
    font-family: var(--font-mono);
  }}

  /* Polygon Container */
  .polygon-box {{
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    min-height: 420px;
    flex: 1;
    background: #060911;
    border-radius: 6px;
    padding: 0.8rem;
    border: 1px solid var(--border);
  }}
  .polygon-box svg {{
    max-width: 100%;
    height: auto;
    display: block;
  }}

  /* PLIP Contacts Table */
  .contacts-scroll {{
    margin-top: 0.75rem;
    max-height: 140px;
    overflow-y: auto;
    border: 1px solid var(--border);
    border-radius: 6px;
  }}
  .contacts-table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 0.78rem;
  }}
  .contacts-table th, .contacts-table td {{
    padding: 0.35rem 0.65rem;
    text-align: left;
    border-bottom: 1px solid var(--border);
  }}
  .contacts-table th {{
    background: #0b1329;
    color: var(--text-muted);
    font-weight: 600;
    position: sticky;
    top: 0;
  }}

  .badge-cat {{
    background: rgba(245, 158, 11, 0.2);
    color: #fbbf24;
    padding: 0.1rem 0.35rem;
    border-radius: 4px;
    font-size: 0.68rem;
    border: 1px solid rgba(245, 158, 11, 0.4);
  }}
  .badge-sec {{
    background: rgba(6, 182, 212, 0.2);
    color: #38bdf8;
    padding: 0.1rem 0.35rem;
    border-radius: 4px;
    font-size: 0.68rem;
    border: 1px solid rgba(6, 182, 212, 0.4);
  }}

  /* Plotly Chart Card */
  .chart-card {{
    margin-bottom: 2.2rem;
    position: relative;
  }}
  #plotly-pareto-chart {{
    width: 100%;
    height: 500px;
    background: var(--bg-card);
    border-radius: 8px;
  }}

  /* Rich Pareto Floating Tooltip */
  #pareto-tooltip {{
    position: absolute;
    display: none;
    z-index: 999;
    pointer-events: none;
    background: rgba(15, 23, 42, 0.98);
    border: 1px solid #475569;
    box-shadow: 0 14px 35px rgba(0,0,0,0.75);
    backdrop-filter: blur(14px);
    border-radius: 8px;
    padding: 10px 12px;
    font-size: 11px;
    color: #f8fafc;
    max-width: 250px;
    transition: opacity 0.12s ease;
  }}
  .tt-img {{
    width: 130px;
    height: 130px;
    background: #ffffff;
    border-radius: 6px;
    display: block;
    margin: 0 auto 8px auto;
    border: 1px solid #cbd5e1;
    object-fit: contain;
  }}
  .tt-title {{
    font-size: 13px;
    font-weight: 700;
    color: #ffffff;
    margin-bottom: 4px;
    word-break: break-all;
  }}
  .tt-tag {{
    display: inline-block;
    padding: 2px 7px;
    border-radius: 4px;
    font-size: 10px;
    font-weight: 600;
    margin-bottom: 8px;
  }}
  .tag-control {{ background: #fee2e2; color: #991b1b; border: 1px solid #ef4444; }}
  .tag-pareto {{ background: #dbeafe; color: #1e3a8a; border: 1px solid #3b82f6; }}
  .tag-candidate {{ background: #d1fae5; color: #065f46; border: 1px solid #10b981; }}
  .tt-row {{
    display: flex;
    justify-content: space-between;
    margin-bottom: 3px;
    gap: 12px;
  }}
  .tt-lbl {{ color: #94a3b8; font-weight: 500; }}
  .tt-val {{ color: #f8fafc; font-weight: 600; font-family: var(--font-mono); }}
  .tt-smi {{
    font-size: 9px;
    color: #64748b;
    margin-top: 6px;
    word-break: break-all;
    font-family: var(--font-mono);
  }}

  /* Data Table Section */
  .table-card {{
    background: var(--bg-card);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 1.4rem;
    margin-bottom: 2.2rem;
  }}
  .table-toolbar {{
    display: flex;
    flex-wrap: wrap;
    justify-content: space-between;
    align-items: center;
    gap: 1rem;
    margin-bottom: 1.1rem;
  }}
  .toolbar-left {{
    display: flex;
    flex-wrap: wrap;
    gap: 0.75rem;
    align-items: center;
  }}
  .search-input {{
    background: #060911;
    color: #f1f5f9;
    border: 1px solid var(--border-light);
    padding: 0.5rem 0.9rem;
    border-radius: 6px;
    font-size: 0.84rem;
    min-width: 270px;
  }}
  .filter-pills {{
    display: flex;
    gap: 0.4rem;
  }}
  .pill {{
    background: #1e293b;
    color: var(--text-muted);
    border: 1px solid var(--border-light);
    padding: 0.35rem 0.75rem;
    border-radius: 9999px;
    font-size: 0.76rem;
    font-weight: 500;
    cursor: pointer;
  }}
  .pill.active {{
    background: var(--accent);
    color: #ffffff;
    border-color: #60a5fa;
  }}
  .export-btns {{
    display: flex;
    gap: 0.5rem;
  }}
  .btn-export {{
    background: #0f172a;
    color: #cbd5e1;
    border: 1px solid var(--border-light);
    padding: 0.4rem 0.8rem;
    border-radius: 6px;
    font-size: 0.78rem;
    font-weight: 600;
    cursor: pointer;
  }}
  .btn-export:hover {{
    background: #1e293b;
    border-color: var(--accent);
    color: #ffffff;
  }}

  .table-wrapper {{
    overflow-x: auto;
    border: 1px solid var(--border);
    border-radius: 6px;
  }}
  .data-table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 0.82rem;
    text-align: left;
  }}
  .data-table th, .data-table td {{
    padding: 0.65rem 0.85rem;
    border-bottom: 1px solid var(--border);
  }}
  .data-table th {{
    background: #0b1329;
    color: #cbd5e1;
    font-weight: 600;
    cursor: pointer;
    user-select: none;
    white-space: nowrap;
  }}
  .data-table th:hover {{
    color: #60a5fa;
  }}
  .data-table tr:hover {{
    background: var(--bg-hover);
    cursor: pointer;
  }}
  .data-table tr.selected-row {{
    background: rgba(59, 130, 246, 0.16) !important;
    border-left: 3px solid var(--accent);
  }}
  .mono {{ font-family: var(--font-mono); }}

  .tag-status {{
    display: inline-block;
    padding: 0.18rem 0.55rem;
    border-radius: 4px;
    font-size: 0.72rem;
    font-weight: 600;
  }}

  /* Methods Block */
  .methods-card {{
    background: var(--bg-card);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 1.4rem;
  }}
  .methods-list {{
    margin-left: 1.3rem;
    margin-top: 0.4rem;
    margin-bottom: 0.9rem;
    color: #cbd5e1;
    font-size: 0.85rem;
  }}
  .methods-list li {{
    margin-bottom: 0.35rem;
  }}
  code {{
    background: #060911;
    color: #38bdf8;
    padding: 0.15rem 0.35rem;
    border-radius: 4px;
    font-size: 0.8rem;
    font-family: var(--font-mono);
  }}
  .cite-card {{
    margin-top: 1.5rem;
    background: #060911;
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 1rem;
    font-size: 0.82rem;
    color: var(--text-muted);
  }}
</style>
</head>
<body>

<!-- Navigation Bar -->
<div class="navbar">
  <div class="nav-brand">
    <span class="nav-logo">PoliScreen Dossier</span>
    <span class="nav-tag">{L['zero_server']}</span>
    {target_nav_html}
  </div>
  <div class="nav-meta">
    {L['project']}: <strong>{html.escape(project_title)}</strong>
  </div>
</div>

<div class="container">

  <!-- KPI Executive Metrics Cards -->
  <div class="kpi-grid">
    <div class="kpi-card">
      <div class="kpi-label">{L['kpi_total']}</div>
      <div class="kpi-val" id="kpi-total">{kpis['total']}</div>
      <div class="kpi-desc">{L['kpi_total_desc']}</div>
    </div>
    <div class="kpi-card">
      <div class="kpi-label">{L['kpi_dock']}</div>
      <div class="kpi-val" id="kpi-dock" style="color: #60a5fa;">{kpis['best_dock']:.2f} <span style="font-size: 0.85rem;">kcal/mol</span></div>
      <div class="kpi-desc">{L['kpi_dock_desc']}</div>
    </div>
    <div class="kpi-card">
      <div class="kpi-label">{L['kpi_eff']}</div>
      <div class="kpi-val" id="kpi-eff" style="color: #34d399; font-size: 1.35rem;">
        {kpis['best_eff_le']:.0f}% <span style="font-size: 0.78rem; color: #94a3b8; font-weight: 500;">(LE)</span> / {kpis['best_eff_raw']:.0f}% <span style="font-size: 0.78rem; color: #94a3b8; font-weight: 500;">(Raw)</span>
      </div>
      <div class="kpi-desc">{L['kpi_eff_desc']}</div>
    </div>
    <div class="kpi-card">
      <div class="kpi-label">{L['kpi_pareto']}</div>
      <div class="kpi-val" id="kpi-pareto" style="color: #fbbf24;">{kpis['pareto']}</div>
      <div class="kpi-desc">{L['kpi_pareto_desc']}</div>
    </div>
    <div class="kpi-card">
      <div class="kpi-label">{L['kpi_beat']}</div>
      <div class="kpi-val" id="kpi-beat">{kpis['beat_control']}</div>
      <div class="kpi-desc">{L['kpi_beat_desc']}</div>
    </div>
    <div class="kpi-card">
      <div class="kpi-label">{L['kpi_conf']}</div>
      <div class="kpi-val" id="kpi-conf">{kpis['reliable']}</div>
      <div class="kpi-desc">{L['kpi_conf_desc']}</div>
    </div>
  </div>

  <!-- Interactive 3D Viewer & PLIP Geometric Footprint Panel -->
  <div class="section-head">
    <div class="section-title">{L['sec_3d_title']}</div>
    <div class="section-desc">{L['sec_3d_desc']}</div>
  </div>

  <div class="visual-grid">
    <!-- 3D Molecular Complex Viewer -->
    <div class="panel-card">
      <div class="viewer-toolbar">
        <label for="pose-selector" style="font-size: 0.8rem; font-weight: 600; color: #cbd5e1;">{L['compound_label']}:</label>
        <select id="pose-selector" class="select-custom" onchange="onPoseSelected(this.value)">
        </select>
        <button id="btn-surface" class="btn" onclick="toggleSurface()">{L['btn_surface']}</button>
        <button id="btn-pocket" class="btn" onclick="togglePocket()" title="{L['btn_pocket_title']}">{L['btn_pocket']}</button>
        <button id="btn-tunnels" class="btn" onclick="toggleTunnels()" {'disabled' if not cur_d['has_tunnels'] else ''}>{L['btn_tunnels']}{' (N/A)' if not cur_d['has_tunnels'] else ''}</button>
        <button id="btn-spin" class="btn" onclick="toggleSpin()">{L['btn_spin']}</button>
        <button class="btn" onclick="centerView()">{L['btn_center']}</button>
      </div>
      <div id="viewer-container"></div>
      <div id="viewer-caption" class="viewer-caption">
        {L['loading_struct']}
      </div>
    </div>

    <!-- PLIP Polygon & Contacts Breakdown -->
    <div class="panel-card">
      <div style="font-size: 0.95rem; font-weight: 700; color: #f8fafc; margin-bottom: 0.5rem; display: flex; justify-content: space-between; align-items: center;">
        <span>{L['plip_title']}</span>
        <span id="polygon-badge" class="tag-status tag-pareto">Pareto Lead</span>
      </div>
      <div id="polygon-container" class="polygon-box">
      </div>
      <div class="contacts-scroll">
        <table class="contacts-table">
          <thead>
            <tr>
              <th>{L['th_residue']}</th>
              <th>{L['th_type']}</th>
              <th>{L['th_count']}</th>
            </tr>
          </thead>
          <tbody id="contacts-tbody">
          </tbody>
        </table>
      </div>
    </div>
  </div>

  <!-- Multi-Objective Pareto Frontier Landscape (Plotly.js) -->
  <div class="panel-card chart-card">
    <div class="section-head">
      <div class="section-title">{L['sec_pareto_title']}</div>
      <div class="section-desc">{L['sec_pareto_desc']}</div>
    </div>
    <div id="plotly-pareto-chart"></div>
    <div id="pareto-tooltip"></div>
  </div>

  <!-- CAVER Transport Tunnels Section -->
  <div id="tunnels-section-container">
    {cur_d['tunnels_html']}
  </div>

  <!-- ADMET Drug-Likeness Section -->
  <div id="admet-section-container">
    {cur_d['admet_html']}
  </div>

  <!-- Searchable & Sortable Classification Table -->
  <div class="table-card">
    <div class="section-head">
      <div class="section-title">{L['sec_table_title']}</div>
    </div>
    <div class="table-toolbar">
      <div class="toolbar-left">
        <input type="text" id="table-search" class="search-input" placeholder="{L['search_ph']}" oninput="onSearchInput(this.value)">
        <div class="filter-pills">
          <button class="pill active" onclick="setFilter('all', this)">{L['filter_all']}</button>
          <button class="pill" onclick="setFilter('pareto', this)">{L['filter_pareto']}</button>
          <button class="pill" onclick="setFilter('beat', this)">{L['filter_beat']}</button>
        </div>
      </div>
      <div class="export-btns">
        <button class="btn-export" onclick="exportTableCSV()">{L['btn_csv']}</button>
        <button class="btn-export" onclick="exportTableExcel()">{L['btn_excel']}</button>
      </div>
    </div>

    <div class="table-wrapper">
      <table class="data-table" id="main-results-table">
        <thead>
          <tr>
            <th onclick="sortTable('name')">{L['th_compound']}</th>
            <th onclick="sortTable('status')">{L['th_status']}</th>
            <th onclick="sortTable('dock')">{L['th_dock']}</th>
            <th onclick="sortTable('quality')">{L['th_quality']}</th>
            <th onclick="sortTable('eff')">{L['th_eff_raw']}</th>
            <th onclick="sortTable('eff_le')">{L['th_eff_le']}</th>
            <th onclick="sortTable('conf')">{L['th_conf']}</th>
            <th onclick="sortTable('pki')">{L['th_pki']}</th>
            <th onclick="sortTable('le')">{L['th_le']}</th>
            <th onclick="sortTable('pareto_rank')">{L['th_rank']}</th>
          </tr>
        </thead>
        <tbody id="table-tbody">
        </tbody>
      </table>
    </div>
  </div>

  <!-- In Silico Methods & Reproducibility Block -->
  <div class="methods-card">
    <div class="section-head">
      <div class="section-title">{L['sec_methods_title']}</div>
      <div class="section-desc">{L['sec_methods_desc']}</div>
    </div>
    <div class="methods-body">
      {methods_html}
    </div>

    <div class="cite-card">
      <strong>{L['how_to_cite']}</strong><br>
      Anaya-Guerrero, D. C. <em>PoliScreen: Reproducible Virtual Screening and Multi-Objective Analogue Optimization</em> (2026).
    </div>
  </div>

</div>

<!-- Embedded Data and Interactivity Scripts -->
<script>
// --- Multi-Target Data Payloads ---
const PROJECT_TITLE = "{html.escape(project_title)}";
const TARGETS_DATA = {json_targets_data};
const ALL_TARGETS = {json_all_targets};
let currentTarget = "{active_target}";

// Pointers to active target's datasets
let KPIS = TARGETS_DATA[currentTarget].kpis;
let ACTIVE_ITEMS = TARGETS_DATA[currentTarget].active_items;
let POSES_PDB = TARGETS_DATA[currentTarget].poses_pdb;
let RECEPTOR_PDB = TARGETS_DATA[currentTarget].receptor_pdb;
let POCKET_RESIDUES_NUMS = TARGETS_DATA[currentTarget].pocket_res_nums;
let CATALYTIC_RESIDUES_NUMS = TARGETS_DATA[currentTarget].catalytic_res_nums;
let SECONDARY_RESIDUES_NUMS = TARGETS_DATA[currentTarget].secondary_res_nums;
let CAVITY_SPHERES = TARGETS_DATA[currentTarget].cavity_spheres;
let CAVER_TUNNELS = TARGETS_DATA[currentTarget].caver_clusters;
let POLYGON_SVGS = TARGETS_DATA[currentTarget].polygon_svgs;
let CONTACTS_DETAIL = TARGETS_DATA[currentTarget].contacts_detail;
let PLOTLY_DATA = TARGETS_DATA[currentTarget].plotly_data;
let STRUCTURES_2D = TARGETS_DATA[currentTarget].structures_2d;
let FULL_TABLE_ROWS = TARGETS_DATA[currentTarget].full_table_rows;

let currentCompound = ACTIVE_ITEMS.length > 0 ? ACTIVE_ITEMS[0].name : "";
let currentFilter = "all";
let currentSearch = "";
let sortColumn = "dock";
let sortAscending = true;

// 3D Viewer variables
let glviewer = null;
let receptorModel = null;
let ligandModel = null;
let surfaceObject = null;
let pocketShapes = [];
let tunnelShapes = [];
let isSpinning = false;
let showSurface = false;
let showPocket = false;
let showTunnels = false;

// Utility number formatter
function fmtNum(val, decimals = 2, fallback = "-") {{
  if (val === null || val === undefined || isNaN(Number(val))) return fallback;
  return Number(val).toFixed(decimals);
}}

// Initialize when DOM is ready
function initializeDossier() {{
  try {{ initPoseSelector(); }} catch(e) {{ console.error("Error in initPoseSelector:", e); }}
  try {{ init3DViewer(); }} catch(e) {{ console.error("Error in init3DViewer:", e); }}
  try {{ initPlotlyChart(); }} catch(e) {{ console.error("Error in initPlotlyChart:", e); }}
  try {{ renderTable(); }} catch(e) {{ console.error("Error in renderTable:", e); }}
  try {{ selectCompound(currentCompound); }} catch(e) {{ console.error("Error in selectCompound:", e); }}
}}

if (document.readyState === "loading") {{
  document.addEventListener("DOMContentLoaded", initializeDossier);
}} else {{
  initializeDossier();
}}

// --- Target / Pocket Switcher ---
function switchTarget(targetName) {{
  if (!TARGETS_DATA[targetName]) return;
  currentTarget = targetName;
  const d = TARGETS_DATA[targetName];

  KPIS = d.kpis;
  ACTIVE_ITEMS = d.active_items;
  POSES_PDB = d.poses_pdb;
  POCKET_RESIDUES_NUMS = d.pocket_res_nums;
  CATALYTIC_RESIDUES_NUMS = d.catalytic_res_nums;
  SECONDARY_RESIDUES_NUMS = d.secondary_res_nums;
  CAVITY_SPHERES = d.cavity_spheres;
  CAVER_TUNNELS = d.caver_clusters;
  POLYGON_SVGS = d.polygon_svgs;
  CONTACTS_DETAIL = d.contacts_detail;
  PLOTLY_DATA = d.plotly_data;
  STRUCTURES_2D = d.structures_2d;
  FULL_TABLE_ROWS = d.full_table_rows;

  const selTarget = document.getElementById("target-selector");
  if (selTarget && selTarget.value !== targetName) selTarget.value = targetName;

  // Update Top KPIs
  const elTotal = document.getElementById("kpi-total"); if (elTotal) elTotal.innerText = KPIS.total;
  const elDock = document.getElementById("kpi-dock"); if (elDock) elDock.innerHTML = fmtNum(KPIS.best_dock, 2) + ' <span style="font-size: 0.85rem;">kcal/mol</span>';
  const elEff = document.getElementById("kpi-eff"); if (elEff) elEff.innerHTML = fmtNum(KPIS.best_eff_le, 0) + '% <span style="font-size: 0.78rem; color: #94a3b8; font-weight: 500;">(LE)</span> / ' + fmtNum(KPIS.best_eff_raw, 0) + '% <span style="font-size: 0.78rem; color: #94a3b8; font-weight: 500;">(Raw)</span>';
  const elPareto = document.getElementById("kpi-pareto"); if (elPareto) elPareto.innerText = KPIS.pareto;
  const elBeat = document.getElementById("kpi-beat"); if (elBeat) elBeat.innerText = KPIS.beat_control;
  const elConf = document.getElementById("kpi-conf"); if (elConf) elConf.innerText = KPIS.reliable;

  // Update Tunnels and ADMET blocks in DOM
  const tunContainer = document.getElementById("tunnels-section-container");
  if (tunContainer) tunContainer.innerHTML = d.tunnels_html || "";
  const admetContainer = document.getElementById("admet-section-container");
  if (admetContainer) admetContainer.innerHTML = d.admet_html || "";

  // Update Tunnels button
  const btnTun = document.getElementById("btn-tunnels");
  if (btnTun) {{
    btnTun.disabled = !d.has_tunnels;
    btnTun.innerText = d.has_tunnels ? "{L['btn_tunnels']}" : "{L['btn_tunnels']} (N/A)";
    if (!d.has_tunnels && showTunnels) {{
      toggleTunnels();
    }}
  }}

  // Reload receptor if it changed
  if (d.receptor_pdb && d.receptor_pdb !== RECEPTOR_PDB) {{
    RECEPTOR_PDB = d.receptor_pdb;
    if (glviewer && receptorModel) {{
      glviewer.removeModel(receptorModel);
      receptorModel = glviewer.addModel(RECEPTOR_PDB, "pdb");
      glviewer.setStyle({{ model: receptorModel }}, {{
        cartoon: {{ color: "spectrum", opacity: 0.95 }}
      }});
    }}
  }}

  // Re-populate Pose Selector
  initPoseSelector();

  // Reset compound to first candidate of new pocket
  currentCompound = ACTIVE_ITEMS.length > 0 ? ACTIVE_ITEMS[0].name : "";
  selectCompound(currentCompound);

  // If Pocket is currently shown, refresh highlights
  if (showPocket) {{
    showPocket = false;
    togglePocket();
  }}

  // If Tunnels are shown, re-trigger
  if (showTunnels && d.has_tunnels) {{
    showTunnels = false;
    toggleTunnels();
  }}

  // Re-draw Plotly Pareto Chart
  initPlotlyChart();

  // Re-render Table
  renderTable();
}}

// --- Pose Selector Initialization ---
function initPoseSelector() {{
  const sel = document.getElementById("pose-selector");
  if (!sel) return;
  sel.innerHTML = "";
  ACTIVE_ITEMS.forEach(it => {{
    const opt = document.createElement("option");
    opt.value = it.name;
    const prefix = it.is_control ? "[REF] " : (it.is_pareto ? "[Pareto] " : "");
    opt.textContent = prefix + it.name + " (" + fmtNum(it.dock, 2) + " kcal/mol)";
    sel.appendChild(opt);
  }});
}}

// --- 3Dmol Viewer Setup ---
function init3DViewer() {{
  const container = document.getElementById("viewer-container");
  if (!container) return;
  const molLib = window.$3Dmol || window["3Dmol"] || (typeof $3Dmol !== "undefined" ? $3Dmol : null);
  if (!molLib) {{
    container.innerHTML = "<div style='color:#ef4444; padding:2rem; text-align:center;'>3Dmol.js is not available offline.</div>";
    return;
  }}

  glviewer = molLib.createViewer(container, {{ backgroundColor: "#060911" }});

  // Load Receptor
  if (RECEPTOR_PDB && RECEPTOR_PDB.length > 50) {{
    receptorModel = glviewer.addModel(RECEPTOR_PDB, "pdb");
    glviewer.setStyle({{ model: receptorModel }}, {{
      cartoon: {{ color: "spectrum", opacity: 0.95 }}
    }});
  }}

  // Load Initial Ligand
  loadLigand(currentCompound);

  glviewer.resize();
  glviewer.zoomTo();
  glviewer.render();

  window.addEventListener("resize", function() {{
    if (glviewer) {{
      glviewer.resize();
      glviewer.render();
    }}
  }});
}}

function loadLigand(cmpName) {{
  if (!glviewer) return;

  if (ligandModel) {{
    glviewer.removeModel(ligandModel);
    ligandModel = null;
  }}

  const pdbStr = POSES_PDB[cmpName];
  if (pdbStr && pdbStr.length > 30) {{
    ligandModel = glviewer.addModel(pdbStr, "pdb");
    glviewer.setStyle({{ model: ligandModel }}, {{
      stick: {{ radius: 0.26, colorscheme: "default" }}
    }});
    glviewer.zoomTo({{ model: ligandModel }});
  }} else {{
    if (receptorModel) glviewer.zoomTo({{ model: receptorModel }});
  }}

  // Maintain receptor cartoon ribbon and current compound interacting residue highlights
  reapplyHighlights(cmpName);

  glviewer.render();

  // Caption update
  const cap = document.getElementById("viewer-caption");
  const it = ACTIVE_ITEMS.find(x => x.name === cmpName);
  if (it && cap) {{
    cap.innerHTML = "<strong>" + it.name + "</strong> (Pose " + it.pose + ") · Affinity: <strong>" +
                    fmtNum(it.dock, 2) + " kcal/mol</strong> · Eff (Raw): <strong>" +
                    fmtNum(it.eff, 0) + "%</strong> · Eff (LE): <strong>" +
                    fmtNum(it.eff_le, 0) + "%</strong>";
  }}
}}

function toggleSurface() {{
  showSurface = !showSurface;
  const btn = document.getElementById("btn-surface");
  if (btn) btn.classList.toggle("active", showSurface);

  if (!glviewer || !receptorModel) return;
  const molLib = window.$3Dmol || window["3Dmol"] || (typeof $3Dmol !== "undefined" ? $3Dmol : null);

  if (!showSurface) {{
    if (surfaceObject !== null) {{
      try {{ glviewer.removeSurface(surfaceObject); }} catch(e) {{}}
      surfaceObject = null;
    }}
  }} else if (molLib) {{
    glviewer.addSurface(molLib.SurfaceType.VDW, {{
      opacity: 0.62,
      color: "#b9c6d6"
    }}, {{ model: receptorModel }}).then(function(sId) {{
      surfaceObject = sId;
      glviewer.render();
    }});
  }}

  reapplyHighlights(currentCompound);
  if (ligandModel) {{
    glviewer.setStyle({{ model: ligandModel }}, {{
      stick: {{ radius: 0.26, colorscheme: "default" }}
    }});
  }}
  glviewer.render();
}}

function reapplyHighlights(cmpName) {{
  if (!glviewer || !receptorModel) return;

  glviewer.setStyle({{ model: receptorModel }}, {{
    cartoon: {{ color: "spectrum", opacity: 0.95 }}
  }});

  const contacts = (CONTACTS_DETAIL && cmpName) ? (CONTACTS_DETAIL[cmpName] || []) : [];
  const activeCatRes = [];
  const activeNonCatRes = [];

  contacts.forEach(c => {{
    const m = (c.residue || "").match(/\\d+/);
    if (!m) return;
    const num = parseInt(m[0], 10);
    if (c.is_cat || (CATALYTIC_RESIDUES_NUMS && CATALYTIC_RESIDUES_NUMS.includes(num))) {{
      if (!activeCatRes.includes(num)) activeCatRes.push(num);
    }} else {{
      if (!activeNonCatRes.includes(num)) activeNonCatRes.push(num);
    }}
  }});

  // When pocket is displayed or if compound has no catalytic contacts, ensure catalytic site residues are visible
  if (showPocket || activeCatRes.length === 0) {{
    if (CATALYTIC_RESIDUES_NUMS) {{
      CATALYTIC_RESIDUES_NUMS.forEach(num => {{
        if (!activeCatRes.includes(num)) activeCatRes.push(num);
      }});
    }}
  }}

  // 1. Additional interacting residues in sleek cyan sticks
  if (activeNonCatRes.length > 0) {{
    glviewer.addStyle({{ model: receptorModel, resi: activeNonCatRes }}, {{
      stick: {{ radius: 0.18, colorscheme: "cyanCarbon" }}
    }});
  }}

  // 2. Catalytic residues in radiant golden/amber sticks
  if (activeCatRes.length > 0) {{
    glviewer.addStyle({{ model: receptorModel, resi: activeCatRes }}, {{
      stick: {{ radius: 0.22, colorscheme: "yellowCarbon" }}
    }});
  }}
}}

function togglePocket() {{
  showPocket = !showPocket;
  const btn = document.getElementById("btn-pocket");
  if (btn) btn.classList.toggle("active", showPocket);

  if (!glviewer || !receptorModel) return;

  // Clear previous alpha spheres
  pocketShapes.forEach(s => {{
    try {{ glviewer.removeShape(s); }} catch(e) {{}}
  }});
  pocketShapes = [];

  if (showPocket) {{
    // 1. Alpha spheres from fpocket cavity detection
    if (CAVITY_SPHERES && CAVITY_SPHERES.length > 0) {{
      CAVITY_SPHERES.forEach(sp => {{
        const shape = glviewer.addSphere({{
          center: {{ x: sp.x, y: sp.y, z: sp.z }},
          radius: sp.r,
          color: "#f4d35e",
          opacity: 0.50
        }});
        pocketShapes.push(shape);
      }});
    }}

    // 2. Re-apply highlights for current compound + catalytic site
    reapplyHighlights(currentCompound);

    // 3. Zoom to pocket cavity
    if (POCKET_RESIDUES_NUMS && POCKET_RESIDUES_NUMS.length > 0) {{
      glviewer.zoomTo({{ model: receptorModel, resi: POCKET_RESIDUES_NUMS }});
    }}
  }} else {{
    reapplyHighlights(currentCompound);
    centerView();
  }}

  if (ligandModel) {{
    glviewer.setStyle({{ model: ligandModel }}, {{
      stick: {{ radius: 0.26, colorscheme: "default" }}
    }});
  }}

  glviewer.render();
}}

function toggleTunnels() {{
  if (!CAVER_TUNNELS || CAVER_TUNNELS.length === 0) return;
  showTunnels = !showTunnels;
  const btn = document.getElementById("btn-tunnels");
  if (btn) btn.classList.toggle("active", showTunnels);

  if (!glviewer) return;

  tunnelShapes.forEach(s => {{
    try {{ glviewer.removeShape(s); }} catch(e) {{}}
  }});
  tunnelShapes = [];

  if (showTunnels) {{
    CAVER_TUNNELS.forEach(tun => {{
      tun.spheres.forEach(sp => {{
        const shape = glviewer.addSphere({{
          center: {{ x: sp.x, y: sp.y, z: sp.z }},
          radius: sp.r,
          color: tun.color,
          opacity: 0.65
        }});
        tunnelShapes.push(shape);
      }});
    }});
  }}

  glviewer.render();
}}

function toggleSpin() {{
  isSpinning = !isSpinning;
  const btn = document.getElementById("btn-spin");
  if (btn) btn.classList.toggle("active", isSpinning);
  if (glviewer) glviewer.spin(isSpinning ? "y" : false);
}}

function centerView() {{
  if (!glviewer) return;
  if (ligandModel) {{
    glviewer.zoomTo({{ model: ligandModel }});
  }} else if (receptorModel) {{
    glviewer.zoomTo({{ model: receptorModel }});
  }}
  glviewer.render();
}}

function onPoseSelected(cmpName) {{
  selectCompound(cmpName);
}}

// --- Selection Synchronizer (NO FORCED SCROLL) ---
function selectCompound(cmpName) {{
  currentCompound = cmpName;

  // 1. Sync dropdown
  const sel = document.getElementById("pose-selector");
  if (sel) sel.value = cmpName;

  // 2. Load in 3D viewer
  loadLigand(cmpName);

  // 3. Update Polygon SVG
  const polyBox = document.getElementById("polygon-container");
  const svgStr = POLYGON_SVGS[cmpName];
  if (polyBox) {{
    if (svgStr) {{
      polyBox.innerHTML = svgStr;
    }} else {{
      polyBox.innerHTML = "<div style='color:var(--text-muted);'>No interaction polygon for this compound.</div>";
    }}
  }}

  // 4. Update Polygon badge
  const it = ACTIVE_ITEMS.find(x => x.name === cmpName);
  const badge = document.getElementById("polygon-badge");
  if (it && badge) {{
    if (it.is_control) {{
      badge.className = "tag-status tag-control";
      badge.textContent = "Crystallographic Control";
    }} else if (it.is_pareto) {{
      badge.className = "tag-status tag-pareto";
      badge.textContent = "Pareto Leader";
    }} else {{
      badge.className = "tag-status tag-candidate";
      badge.textContent = "Candidate";
    }}
  }}

  // 5. Update Contacts table
  const tbody = document.getElementById("contacts-tbody");
  if (tbody) {{
    tbody.innerHTML = "";
    const contacts = CONTACTS_DETAIL[cmpName] || [];
    if (contacts.length === 0) {{
      tbody.innerHTML = "<tr><td colspan='3' style='text-align:center; color:var(--text-muted);'>No direct contacts detected.</td></tr>";
    }} else {{
      contacts.forEach(c => {{
        const tr = document.createElement("tr");
        let resBadge = "";
        if (c.is_cat) resBadge = " <span class='badge-cat'>Catalytic</span>";
        if (c.is_sec) resBadge = " <span class='badge-sec'>Secondary</span>";
        tr.innerHTML = "<td><strong>" + c.residue + "</strong>" + resBadge + "</td>" +
                       "<td><span style='color:" + c.color + "; font-weight:700;'>●</span> " + c.type + "</td>" +
                       "<td class='mono'>" + c.count + "</td>";
        tbody.appendChild(tr);
      }});
    }}
  }}

  // 6. Highlight row in Table (without auto-scroll)
  const rows = document.querySelectorAll("#table-tbody tr");
  rows.forEach(r => {{
    if (r.dataset.name === cmpName) {{
      r.classList.add("selected-row");
    }} else {{
      r.classList.remove("selected-row");
    }}
  }});
}}

// --- Plotly Chart Setup & Rich Floating Tooltip ---
function initPlotlyChart() {{
  const chartDiv = document.getElementById("plotly-pareto-chart");
  if (!chartDiv) return;
  const plotLib = window.Plotly || (typeof Plotly !== "undefined" ? Plotly : null);
  if (!plotLib) {{
    chartDiv.innerHTML = "<div style='color:#ef4444; padding:2rem; text-align:center;'>Plotly.js is not available offline.</div>";
    return;
  }}

  const pts = PLOTLY_DATA.points || [];
  const paretoLine = PLOTLY_DATA.pareto_line || [];
  const ctrl = PLOTLY_DATA.ctrl;

  const lineTrace = {{
    x: paretoLine.map(p => p.x),
    y: paretoLine.map(p => p.y),
    mode: "lines",
    name: "Pareto Frontier",
    line: {{ color: "#3b82f6", dash: "dash", width: 2.2 }},
    hoverinfo: "none"
  }};

  const cands = pts.filter(p => !p.is_control && !p.is_pareto);
  const candTrace = {{
    x: cands.map(p => p.x),
    y: cands.map(p => p.y),
    mode: "markers",
    name: "Candidates",
    text: cands.map(p => p.name),
    customdata: cands,
    hoverinfo: "none",
    marker: {{ color: "#10b981", size: 9, opacity: 0.85, line: {{ color: "#065f46", width: 1 }} }}
  }};

  const paretoPts = pts.filter(p => p.is_pareto);
  const paretoTrace = {{
    x: paretoPts.map(p => p.x),
    y: paretoPts.map(p => p.y),
    mode: "markers+text",
    name: "Pareto Leaders",
    text: paretoPts.map(p => p.name),
    textposition: "top center",
    textfont: {{ size: 11, color: "#f8fafc" }},
    customdata: paretoPts,
    hoverinfo: "none",
    marker: {{ color: "#2563eb", size: 13, line: {{ color: "#fbbf24", width: 2.2 }} }}
  }};

  const ctrlPts = pts.filter(p => p.is_control);
  const ctrlTrace = {{
    x: ctrlPts.map(p => p.x),
    y: ctrlPts.map(p => p.y),
    mode: "markers+text",
    name: "Control",
    text: ctrlPts.map(p => p.name),
    textposition: "top center",
    textfont: {{ size: 11, color: "#f87171" }},
    customdata: ctrlPts,
    hoverinfo: "none",
    marker: {{ color: "#ef4444", size: 14, symbol: "diamond", line: {{ color: "#7f1d1d", width: 1.5 }} }}
  }};

  const shapes = [];
  if (ctrl) {{
    shapes.push({{
      type: "line",
      x0: ctrl.x, x1: ctrl.x,
      yref: "paper",
      y0: 0, y1: 1,
      line: {{ color: "#ef4444", width: 1.2, dash: "dot" }}
    }});
    shapes.push({{
      type: "line",
      xref: "paper",
      x0: 0, x1: 1,
      y0: ctrl.y, y1: ctrl.y,
      line: {{ color: "#ef4444", width: 1.2, dash: "dot" }}
    }});
  }}

  const layout = {{
    paper_bgcolor: "transparent",
    plot_bgcolor: "#060911",
    margin: {{ l: 65, r: 30, t: 30, b: 60 }},
    xaxis: {{
      title: {{ text: "Vina Binding Affinity (kcal/mol; more negative = tighter binding)", font: {{ color: "#cbd5e1", size: 12 }} }},
      gridcolor: "#1e293b",
      color: "#94a3b8",
      zerolinecolor: "#334155",
      autorange: "reversed"
    }},
    yaxis: {{
      title: {{ text: "PLIP Interaction Quality (0-1 vs. Reference)", font: {{ color: "#cbd5e1", size: 12 }} }},
      gridcolor: "#1e293b",
      color: "#94a3b8",
      zerolinecolor: "#334155"
    }},
    legend: {{
      font: {{ color: "#f8fafc", size: 11 }},
      bgcolor: "rgba(15, 23, 42, 0.85)",
      bordercolor: "#334155",
      borderwidth: 1
    }},
    shapes: shapes
  }};

  const config = {{
    responsive: true,
    displayModeBar: true,
    displaylogo: false,
    modeBarButtonsToRemove: ["lasso2d", "select2d"]
  }};

  const pPromise = plotLib.newPlot(chartDiv, [lineTrace, candTrace, paretoTrace, ctrlTrace], layout, config);

  const setupEvents = function() {{
    if (typeof chartDiv.on === "function") {{
      chartDiv.on("plotly_hover", function(data) {{
        if (data.points && data.points.length > 0) {{
          const pt = data.points[0];
          if (pt.customdata) {{
            showTooltip(data.event, pt.customdata);
          }}
        }}
      }});

      chartDiv.on("plotly_unhover", function() {{
        hideTooltip();
      }});

      chartDiv.on("plotly_click", function(data) {{
        if (data.points && data.points.length > 0) {{
          const pt = data.points[0];
          if (pt.customdata && pt.customdata.name) {{
            selectCompound(pt.customdata.name);
          }}
        }}
      }});
    }}
  }};

  if (pPromise && typeof pPromise.then === "function") {{
    pPromise.then(setupEvents);
  }} else {{
    setupEvents();
  }}
}}

// --- Rich Custom Tooltip Functions ---
function showTooltip(e, pt) {{
  const tooltip = document.getElementById("pareto-tooltip");
  if (!tooltip) return;

  let tagClass = "tag-candidate";
  if (pt.is_control) tagClass = "tag-control";
  else if (pt.is_pareto) tagClass = "tag-pareto";

  let htmlStr = "";
  if (pt.img) {{
    htmlStr += '<img class="tt-img" src="' + pt.img + '" alt="Structure" />';
  }}
  htmlStr += '<div class="tt-title">' + pt.name + '</div>';
  htmlStr += '<span class="tt-tag ' + tagClass + '">' + pt.status + '</span>';
  htmlStr += '<div class="tt-row"><span class="tt-lbl">Effectiveness (Raw)</span><span class="tt-val">' + fmtNum(pt.eff, 1) + '%</span></div>';
  htmlStr += '<div class="tt-row"><span class="tt-lbl">Effectiveness (LE)</span><span class="tt-val">' + fmtNum(pt.eff_le, 1) + '%</span></div>';
  htmlStr += '<div class="tt-row"><span class="tt-lbl">Docking</span><span class="tt-val">' + fmtNum(pt.dock, 2) + ' kcal/mol</span></div>';
  htmlStr += '<div class="tt-row"><span class="tt-lbl">Quality</span><span class="tt-val">' + fmtNum(pt.quality, 3) + '</span></div>';
  if (pt.conf !== null && pt.conf !== undefined) htmlStr += '<div class="tt-row"><span class="tt-lbl">Confidence</span><span class="tt-val">' + fmtNum(pt.conf, 2) + '</span></div>';
  if (pt.le !== null && pt.le !== undefined) htmlStr += '<div class="tt-row"><span class="tt-lbl">LE</span><span class="tt-val">' + fmtNum(pt.le, 3) + '</span></div>';
  if (pt.prank !== null && pt.prank !== undefined) htmlStr += '<div class="tt-row"><span class="tt-lbl">Pareto Rank</span><span class="tt-val">#' + pt.prank + '</span></div>';
  if (pt.smiles) htmlStr += '<div class="tt-smi">' + pt.smiles + '</div>';

  tooltip.innerHTML = htmlStr;
  tooltip.style.display = "block";
  positionTooltip(e);
}}

function positionTooltip(e) {{
  const tooltip = document.getElementById("pareto-tooltip");
  const container = document.querySelector(".chart-card");
  if (!container || !tooltip) return;
  const rect = container.getBoundingClientRect();
  const mouseX = e.clientX - rect.left;
  const mouseY = e.clientY - rect.top;
  const ttWidth = tooltip.offsetWidth || 240;
  const ttHeight = tooltip.offsetHeight || 300;

  let x = mouseX + 16;
  let y = mouseY - 20;

  if (x + ttWidth > rect.width - 15) {{
    x = mouseX - ttWidth - 16;
  }}
  if (x < 10) x = 10;

  if (y + ttHeight > rect.height - 15) {{
    y = rect.height - ttHeight - 15;
  }}
  if (y < 10) y = 10;

  tooltip.style.left = x + "px";
  tooltip.style.top = y + "px";
}}

function hideTooltip() {{
  const tooltip = document.getElementById("pareto-tooltip");
  if (tooltip) tooltip.style.display = "none";
}}

// --- Dynamic Table Logic ---
function getFilteredRows() {{
  return FULL_TABLE_ROWS.filter(r => {{
    if (currentSearch && !r.name.toLowerCase().includes(currentSearch.toLowerCase())) {{
      return false;
    }}
    if (currentFilter === "pareto" && !r.is_pareto && !r.is_control) {{
      return false;
    }}
    if (currentFilter === "beat" && r.eff < 105 && !r.is_control) {{
      return false;
    }}
    return true;
  }});
}}

function renderTable() {{
  const tbody = document.getElementById("table-tbody");
  if (!tbody) return;
  tbody.innerHTML = "";

  let filtered = getFilteredRows();

  filtered.sort((a, b) => {{
    let va = a[sortColumn];
    let vb = b[sortColumn];
    if (va === null || va === undefined) return 1;
    if (vb === null || vb === undefined) return -1;
    if (typeof va === "string") {{
      return sortAscending ? va.localeCompare(vb) : vb.localeCompare(va);
    }}
    return sortAscending ? (va - vb) : (vb - va);
  }});

  filtered.forEach(row => {{
    const tr = document.createElement("tr");
    tr.dataset.name = row.name;
    if (row.name === currentCompound) tr.classList.add("selected-row");

    let tagClass = "tag-candidate";
    if (row.is_control) tagClass = "tag-control";
    else if (row.is_pareto) tagClass = "tag-pareto";

    tr.innerHTML =
      "<td><strong>" + row.name + "</strong>" + (row.has_3d ? " <span style='color:#3b82f6; font-size:0.75rem;'>[3D]</span>" : "") + "</td>" +
      "<td><span class='tag-status " + tagClass + "'>" + row.status + "</span></td>" +
      "<td class='mono'>" + fmtNum(row.dock, 2) + "</td>" +
      "<td class='mono'>" + fmtNum(row.quality, 3) + "</td>" +
      "<td class='mono'>" + fmtNum(row.eff, 0) + "%</td>" +
      "<td class='mono'>" + fmtNum(row.eff_le, 0) + "%</td>" +
      "<td class='mono'>" + fmtNum(row.conf, 2) + "</td>" +
      "<td class='mono'>" + fmtNum(row.pki, 2) + "</td>" +
      "<td class='mono'>" + fmtNum(row.le, 3) + "</td>" +
      "<td class='mono'>" + (row.pareto_rank !== null && row.pareto_rank !== undefined ? "#" + row.pareto_rank : "-") + "</td>";

    tr.onclick = function() {{
      selectCompound(row.name);
    }};

    tbody.appendChild(tr);
  }});
}}

function onSearchInput(val) {{
  currentSearch = val;
  renderTable();
}}

function setFilter(filterType, btn) {{
  currentFilter = filterType;
  document.querySelectorAll(".filter-pills .pill").forEach(p => p.classList.remove("active"));
  btn.classList.add("active");
  renderTable();
}}

function sortTable(col) {{
  if (sortColumn === col) {{
    sortAscending = !sortAscending;
  }} else {{
    sortColumn = col;
    sortAscending = (col === "dock" || col === "pareto_rank");
  }}
  renderTable();
}}

// --- Table Export Functions ---
function exportTableCSV() {{
  const headers = ["Compound", "Status", "Docking_kcal_mol", "Interaction_Quality", "Effectiveness_Raw_pct", "Effectiveness_LE_pct", "Confidence", "pKi", "LE", "Pareto_Rank"];
  const rows = getFilteredRows().map(r => [
    '"' + r.name.replace(/"/g, '""') + '"',
    '"' + r.status + '"',
    fmtNum(r.dock, 2),
    fmtNum(r.quality, 3),
    fmtNum(r.eff, 1),
    fmtNum(r.eff_le, 1),
    fmtNum(r.conf, 2, ""),
    fmtNum(r.pki, 2, ""),
    fmtNum(r.le, 3, ""),
    r.pareto_rank !== null && r.pareto_rank !== undefined ? r.pareto_rank : ""
  ]);
  const csvContent = "\\uFEFF" + [headers.join(","), ...rows.map(e => e.join(","))].join("\\r\\n");
  const activeRecName = (TARGETS_DATA[currentTarget] && TARGETS_DATA[currentTarget].rec_name) ? TARGETS_DATA[currentTarget].rec_name.replace(/[^a-zA-Z0-9_-]/g, "_") : "target";
  downloadBlob(csvContent, "text/csv;charset=utf-8;", PROJECT_TITLE + "_" + activeRecName + "_ranking.csv");
}}

function exportTableExcel() {{
  const headers = ["Compound", "Status", "Docking (kcal/mol)", "Interaction Quality", "Effectiveness Raw (%)", "Effectiveness LE (%)", "Confidence", "pKi", "LE", "Pareto Rank"];
  const rows = getFilteredRows();
  let tableHtml = '<html xmlns:o="urn:schemas-microsoft-com:office:office" xmlns:x="urn:schemas-microsoft-com:office:excel" xmlns="http://www.w3.org/TR/REC-html40">';
  tableHtml += '<head><meta charset="utf-8"><!--[if gte mso 9]><xml><x:ExcelWorkbook><x:ExcelWorksheets><x:ExcelWorksheet><x:Name>Ranking</x:Name><x:WorksheetOptions><x:DisplayGridlines/></x:WorksheetOptions></x:ExcelWorksheet></x:ExcelWorksheets></x:ExcelWorkbook></xml><![endif]--></head><body>';
  tableHtml += '<table border="1"><thead><tr style="background:#0f172a;color:#ffffff;font-weight:bold;">';
  headers.forEach(h => {{ tableHtml += `<th>${{h}}</th>`; }});
  tableHtml += '</tr></thead><tbody>';
  rows.forEach(r => {{
    tableHtml += `<tr><td>${{r.name}}</td><td>${{r.status}}</td><td>${{fmtNum(r.dock, 2)}}</td><td>${{fmtNum(r.quality, 3)}}</td><td>${{fmtNum(r.eff, 1)}}</td><td>${{fmtNum(r.eff_le, 1)}}</td><td>${{fmtNum(r.conf, 2, "")}}</td><td>${{fmtNum(r.pki, 2, "")}}</td><td>${{fmtNum(r.le, 3, "")}}</td><td>${{r.pareto_rank !== null && r.pareto_rank !== undefined ? r.pareto_rank : ""}}</td></tr>`;
  }});
  tableHtml += '</tbody></table></body></html>';
  const activeRecName = (TARGETS_DATA[currentTarget] && TARGETS_DATA[currentTarget].rec_name) ? TARGETS_DATA[currentTarget].rec_name.replace(/[^a-zA-Z0-9_-]/g, "_") : "target";
  downloadBlob(tableHtml, "application/vnd.ms-excel", PROJECT_TITLE + "_" + activeRecName + "_ranking.xls");
}}

function downloadBlob(content, mimeType, filename) {{
  const blob = new Blob([content], {{ type: mimeType }});
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}}
</script>
</body>
</html>
"""


def build_all_target_reports(
    proj: str | Path,
    out_dir: Optional[str | Path] = None,
    max_poses: int = 25,
    lang: str = "en",
) -> dict[str, Path]:
    """Build interactive HTML dossiers for every target receptor evaluated in the project."""
    proj_path = Path(proj).resolve()
    target_out_dir = Path(out_dir).resolve() if out_dir else proj_path
    target_out_dir.mkdir(parents=True, exist_ok=True)

    rk_p = proj_path / "ranking.csv"
    if not rk_p.is_file():
        rk_p = proj_path / "results.csv"

    receptors: list[str] = []
    if rk_p.is_file():
        try:
            df = sc.normalize_columns(pd.read_csv(rk_p))
            if "receptor" in df.columns:
                receptors = sorted(df["receptor"].dropna().unique().tolist())
        except Exception:
            pass

    if not receptors:
        meta_p = proj_path / "run.json"
        if meta_p.is_file():
            try:
                meta = json.loads(meta_p.read_text(encoding="utf-8"))
                receptors = [sc.base_of(r) for r in meta.get("receptors", [])]
            except Exception:
                pass

    if not receptors:
        out_f = target_out_dir / f"{proj_path.name}_Dossier.html"
        build_interactive_report(proj_path, out_path=out_f, max_poses=max_poses, lang=lang)
        return {"default": out_f}

    generated: dict[str, Path] = {}
    for rec in receptors:
        site_suffix = f"_{rec.split('~', 1)[1]}" if "~" in rec else ""
        clean_name = sc.display_name(sc.base_of(rec)).replace(" ", "_") + site_suffix
        out_f = target_out_dir / f"{proj_path.name}_{clean_name}_Dossier.html"
        build_interactive_report(
            proj_path,
            out_path=out_f,
            target=rec,
            max_poses=max_poses,
            lang=lang,
            all_targets=True,
        )
        generated[rec] = out_f

    return generated
