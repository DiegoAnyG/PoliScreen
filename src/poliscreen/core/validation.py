"""Redocking validation.

If the co-crystallized control does not recover its pose when re-docked, that target's setup is
wrong (box off the site, poorly prepared receptor, ligand from another crystal) and no downstream
number is reliable. It runs on every screening to warn early, not late.
"""
from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from .screening import (base_of, compound_from_pose_name, display_name, model_of, normalize_key,
                        receptor_from_name, redocking_rmsd)

RMSD_ACCEPTABLE = 2.0


def _formula(pdb) -> dict:
    """Heavy-atom count per element, to compare two structures."""
    out = {}
    try:
        for l in Path(pdb).read_text(errors="ignore").splitlines():
            if l.startswith(("ATOM", "HETATM")):
                e = (l[76:78].strip() or l[12:16].strip()[:1]).upper()
                if e and e != "H":
                    out[e] = out.get(e, 0) + 1
    except Exception:
        pass
    return out


def _reason_no_rmsd(ref, pose) -> str:
    """Why there is no RMSD. obrms returns infinity when it cannot match the atoms, a sign that two
    different molecules are being compared."""
    fr, fp = _formula(ref), _formula(pose)
    if fr and fp and fr != fp:
        sobran = {e: n - fp.get(e, 0) for e, n in fr.items() if n > fp.get(e, 0)}
        detalle = ", ".join(f"{e}×{n}" for e, n in sorted(sobran.items())) or "different composition"
        return (f"reference and pose are not the same molecule "
                f"({sum(fr.values())} vs {sum(fp.values())} heavy atoms; {detalle} extra "
                f"in the reference)")
    return "not computable"


def redock_validation(controls: Sequence, control_assign: dict, poses_dir,
                      rmsd_ok: float = RMSD_ACCEPTABLE, max_poses: int = 10) -> pd.DataFrame:
    """RMSD of each control against its co-crystallized file.

    rmsd_pose1 is the best-energy pose; rmsd_min the closest to the crystal among the first ones.
    Valid if the RMSD drops below rmsd_ok (2 A by convention).
    """
    poses_dir = Path(poses_dir)
    rows_ = []
    for ref in controls:
        ref = Path(ref)
        ck = normalize_key(ref.stem)
        target_ = control_assign.get(ck)
        all_of = [p for p in poses_dir.glob("*-model*.pdb")
                 if normalize_key(compound_from_pose_name(p.stem)) == ck
                 and (not target_ or base_of(receptor_from_name(p.stem)) == base_of(target_))]
        by_site = {}
        for p in all_of:
            by_site.setdefault(receptor_from_name(p.stem), []).append(p)
        if not by_site:
            rows_.append({"control": ref.stem, "target": target_ or "unassigned", "rmsd_pose1_A": np.nan,
                          "rmsd_min_A": np.nan, "validated": "missing poses"})
            continue
        best_ = None
        for site, ps in by_site.items():
            ps = sorted(ps, key=lambda p: model_of(p.stem))[:max_poses]
            rs = [redocking_rmsd(ref, p) for p in ps]
            v = [r for r in rs if pd.notna(r)]
            key_ = min(v) if v else float("inf")
            if best_ is None or key_ < best_[0]:
                best_ = (key_, site, rs)
        target_ = best_[1]
        rmsds = best_[2]
        validos = [r for r in rmsds if pd.notna(r)]
        r1 = rmsds[0] if rmsds else np.nan
        rmin = min(validos) if validos else np.nan
        if pd.notna(r1) and r1 <= rmsd_ok:
            verdict = "YES"
        elif pd.notna(rmin) and rmin <= rmsd_ok:
            verdict = "YES (in another pose)"
        elif validos:
            verdict = "NO"
        else:
            verdict = _reason_no_rmsd(ref, by_site[target_][0])
        rows_.append({"control": ref.stem, "target": target_ or "unassigned",
                      "rmsd_pose1_A": round(r1, 3) if pd.notna(r1) else np.nan,
                      "rmsd_min_A": round(rmin, 3) if pd.notna(rmin) else np.nan,
                      "validated": verdict})
    return pd.DataFrame(rows_)


# Both spellings are accepted: tables written by earlier versions say SI.
_VALID_PREFIXES = ("YES", "SI")
LEGACY_COLUMNS = {"diana": "target", "validado": "validated"}


def is_valid(verdict) -> bool:
    """True if the verdict means the control recovered its pose, in any of the stored spellings."""
    return str(verdict).strip().upper().startswith(_VALID_PREFIXES)


def normalize(val: pd.DataFrame) -> pd.DataFrame:
    """Renames the columns of a table written by an older version so it can be read today."""
    if val is None or val.empty:
        return val
    faltan = {k: v for k, v in LEGACY_COLUMNS.items() if k in val.columns and v not in val.columns}
    return val.rename(columns=faltan) if faltan else val


def summary(val: pd.DataFrame) -> dict:
    """The verdict as data: {'ok', 'n', 'n_failing', 'targets'}.

    The sentence is built by whoever displays it. Assembling it here would hand the
    interface a finished English string with no way to translate it.
    """
    if val is None or val.empty:
        return {"ok": None, "n": 0, "n_failing": 0, "targets": []}
    val = normalize(val)
    failing = val[~val["validated"].map(is_valid)]
    return {"ok": failing.empty, "n": len(val), "n_failing": len(failing),
            "targets": [display_name(x) for x in failing["target"].astype(str)]}


def summary_text(val: pd.DataFrame) -> str:
    """The same verdict as an English sentence, for the console."""
    s = summary(val)
    if s["ok"] is None:
        return "No controls: the setup cannot be validated."
    n, m = s["n"], s["n_failing"]
    if s["ok"]:
        subject_ = "The control recovers" if n == 1 else f"The {n} controls recover"
        return f"{subject_} the crystallographic pose: the setup is reliable."
    how_many = "The control does not recover" if n == 1 else f"{m} of {n} controls do NOT recover"
    return (f"WARNING: {how_many} the pose ({', '.join(s['targets'])}). "
            "Check the box or the preparation of that target before trusting the ranking.")
