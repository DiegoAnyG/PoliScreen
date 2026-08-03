"""Cavity (pocket) detection with fpocket and their druggability properties.

Like MolModa: lists the pockets with their Druggability Score, volume and center, to center the
docking box on the most promising site without guessing coordinates. fpocket is free (conda-forge).
"""
from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import numpy as np


def fpocket_available() -> bool:
    return shutil.which("fpocket") is not None


def _spheres_from_pqr(pqr, max_n: int = 500) -> list:
    """Alpha-spheres of a pocket: [(x, y, z, radius)]. They define the real VOLUME of the cavity;
    rendered translucent they give an isosurface (like CaverWeb/MolModa), not a box."""
    out = []
    for l in Path(pqr).read_text(errors="ignore").splitlines():
        if not l.startswith(("ATOM", "HETATM")):
            continue
        try:
            x, y, z = float(l[30:38]), float(l[38:46]), float(l[46:54])
            toks = l.split()
            r = float(toks[-1]) if toks and 0.5 <= float(toks[-1]) <= 6.0 else 1.6
        except (ValueError, IndexError):
            continue
        out.append((round(x, 2), round(y, 2), round(z, 2), round(r, 2)))
    return out[:max_n]


def _residues_from_atm(atm) -> list:
    """Residues lining the cavity, read from fpocket's pocket{n}_atm.pdb.
    They let us say whether a cavity contains catalytic residues, not just where it is."""
    res, seen_items = [], set()
    for l in Path(atm).read_text(errors="ignore").splitlines():
        if not l.startswith(("ATOM", "HETATM")):
            continue
        rt, rn = l[17:20].strip(), l[22:26].strip()
        if rt and rn and (rt, rn) not in seen_items:
            seen_items.add((rt, rn))
            res.append(f"{rt.capitalize()}{rn}")
    return sorted(res, key=lambda r: (int("".join(c for c in r if c.isdigit()) or 0), r))


def _dims_from_pqr(pqr, pad: float = 6.0, lo: float = 14.0, hi: float = 30.0) -> Optional[tuple]:
    """Center and PER-AXIS size of a pocket from its vertices (alpha spheres). Returns
    (cx,cy,cz,sx,sy,sz): an anisotropic box that follows the pocket shape (like MolModa), not a cube.
    The center is that of the bounding box; each side = axis extent + margin, clamped to [lo,hi]."""
    pts = []
    for l in Path(pqr).read_text(errors="ignore").splitlines():
        if l.startswith(("ATOM", "HETATM")):
            try:
                pts.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
            except ValueError:
                continue
    if not pts:
        return None
    a = np.array(pts)
    mn, mx = a.min(0), a.max(0)
    c = (mn + mx) / 2.0
    ext = mx - mn
    dims = np.clip(ext + pad, lo, hi)
    return (round(float(c[0]), 2), round(float(c[1]), 2), round(float(c[2]), 2),
            round(float(dims[0]), 1), round(float(dims[1]), 1), round(float(dims[2]), 1),
            round(float(ext[0]), 1), round(float(ext[1]), 1), round(float(ext[2]), 1))


def _parse_info(info_txt) -> dict:
    props, cur = {}, None
    for line in Path(info_txt).read_text(errors="ignore").splitlines():
        m = re.match(r"Pocket\s+(\d+)", line.strip())
        if m:
            cur = int(m.group(1)); props[cur] = {}; continue
        if cur is not None and ":" in line:
            k, _, v = line.strip().partition(":")
            try:
                props[cur][k.strip()] = float(v.strip())
            except ValueError:
                pass
    return props


def detect(pdb, timeout: int = 300) -> list:
    """Returns the pockets ordered by druggability (desc).

    Each one: {n, druggability, score, volume, spheres, cx, cy, cz, size, label}.
    Empty list if fpocket is absent or finds nothing.
    """
    if not fpocket_available():
        return []
    pdb = Path(pdb)
    tmp = Path(tempfile.mkdtemp())
    try:
        local = tmp / f"{pdb.stem}.pdb"
        shutil.copy(pdb, local)
        subprocess.run(["fpocket", "-f", str(local)], capture_output=True, text=True, timeout=timeout)
        outd = tmp / f"{pdb.stem}_out"
        info = outd / f"{pdb.stem}_info.txt"
        if not info.exists():
            return []
        props = _parse_info(info)
        pockets = []
        for n in sorted(props):
            vert = outd / "pockets" / f"pocket{n}_vert.pqr"
            ctr = _dims_from_pqr(vert) if vert.exists() else None
            if ctr is None:
                continue
            cx, cy, cz, sx, sy, sz, ex, ey, ez = ctr
            alpha = _spheres_from_pqr(vert) if vert.exists() else []
            p = props[n]
            dr = p.get("Druggability Score")
            atm = outd / "pockets" / f"pocket{n}_atm.pdb"
            pockets.append({
                "n": n, "druggability": dr, "score": p.get("Score"), "volume": p.get("Volume"),
                "spheres": int(p.get("Number of Alpha Spheres", 0)),
                "cx": cx, "cy": cy, "cz": cz, "sx": sx, "sy": sy, "sz": sz,
                "ex": ex, "ey": ey, "ez": ez, "minimo_aplicado": bool(max(ex, ey, ez) + 6.0 < 14.0),
                "alpha_xyz": alpha,
                "residues": _residues_from_atm(atm) if atm.exists() else [],
                "props": dict(p),
                "size": round(max(sx, sy, sz), 1),
                "label": f"Pocket {n} · druggability {dr:.2f} · vol {p.get('Volume', 0):.0f} · "
                         f"box {sx:.0f}x{sy:.0f}x{sz:.0f}" if dr is not None else f"Pocket {n}",
            })
        pockets.sort(key=lambda x: (x["druggability"] if x["druggability"] is not None else -1), reverse=True)
        return pockets
    except subprocess.TimeoutExpired:
        return []
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def size_from_volume(volume, factor: float = 1.8, pad: float = 8.0,
                     lo: float = 16.0, hi: float = 30.0) -> float:
    """Box edge from the pocket volume (A^3): side of the equivalent cube x factor + margin.
    The generous factor leaves room for exploration; clamped to [lo, hi] to stay within the site."""
    if not volume or volume <= 0:
        return 22.0
    return round(float(np.clip(volume ** (1.0 / 3.0) * factor + pad, lo, hi)), 1)


def pocket_box(pocket: dict, sizing: str = "shape"):
    """Turns a pocket into a docking box (Box).

    sizing='shape': ANISOTROPIC box that follows the pocket shape (sx,sy,sz per axis; like MolModa).
    sizing='cube': cube of the largest side. sizing='volume': cube derived from the fpocket volume.
    """
    from .docking import Box
    if sizing == "volume":
        s = size_from_volume(pocket.get("volume"))
        return Box(pocket["cx"], pocket["cy"], pocket["cz"], s, s, s)
    if sizing == "cube":
        s = pocket.get("size", 22.0)
        return Box(pocket["cx"], pocket["cy"], pocket["cz"], s, s, s)
    sx = pocket.get("sx", pocket.get("size", 22.0))
    sy = pocket.get("sy", pocket.get("size", 22.0))
    sz = pocket.get("sz", pocket.get("size", 22.0))
    return Box(pocket["cx"], pocket["cy"], pocket["cz"], sx, sy, sz)
