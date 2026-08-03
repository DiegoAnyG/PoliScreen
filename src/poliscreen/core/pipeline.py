"""Orchestrator: chains design, docking, interactions and scoring.

Writes everything to folders, no database. Each stage is resumable: if the result is already on disk
it is not recomputed, so re-launching a run is cheap.
"""
from __future__ import annotations

import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path

from . import layout as lay
from typing import Callable, Optional, Sequence

import pandas as pd

from . import adcp
from . import docking as dk
from . import interactions as ix
from . import ligands as lg
from . import peptides as pep
from . import screening as sc
from . import validation as vl
from .design import AdmelabBridge

DEFAULT_WEIGHTS = dict(sc.DEFAULT_WEIGHTS)


@dataclass
class RunConfig:
    receptors: Sequence[Path]
    out_dir: Path
    lead: Optional[str] = None
    ligands: Sequence[Path] = field(default_factory=list)
    controls: Sequence[Path] = field(default_factory=list)
    control_map: dict = field(default_factory=dict)
    boxes: dict = field(default_factory=dict)
    site_boxes: dict = field(default_factory=dict)
    catalytic: dict = field(default_factory=dict)
    weights: dict = field(default_factory=lambda: dict(DEFAULT_WEIGHTS))
    n_analogs: int = 20
    n_substitutions: Sequence[int] = (1,)
    use_ml: bool = True
    seed: int = 42
    exhaustiveness: int = 24
    n_poses: int = 10
    energy_range: float = 3.0
    ph: float = 7.4
    cpu: int = 1
    workers: int = 0
    reuse: bool = False     # everything is recomputed by default: reusing poses from another box falsifies the result
    rescoring_cnn: bool = False
    adcp_steps: int = 250_000
    adcp_replicas: int = 20


DERIVED = ("poses", lay.COMPLEXES, "xml_plip", "san", "prep", "ligands", "xtal")
DERIVED_FILES = (lay.DOCKING_CSV, lay.INTERACTIONS_CSV, "ranking.csv",
                      lay.SUMMARY_CSV, lay.ANALOGUES_CSV, "run.json")


def clean(out_dir) -> None:
    """Deletes what was computed, not what the user provided (receptors and controls are kept)."""
    out = Path(out_dir)
    for d in DERIVED:
        shutil.rmtree(out / d, ignore_errors=True)
    for f in DERIVED_FILES:
        (out / f).unlink(missing_ok=True)


@dataclass
class RunResult:
    out_dir: Path
    analogs: Optional[pd.DataFrame] = None
    docking: Optional[pd.DataFrame] = None
    interactions: Optional[pd.DataFrame] = None
    ranking: Optional[pd.DataFrame] = None
    validation: Optional[pd.DataFrame] = None
    ref_info: dict = field(default_factory=dict)
    docking_errors: list = field(default_factory=list)


def _assign_controls(controls: Sequence[Path], receptors: Sequence[Path], given: dict) -> dict:
    """Matches each control with its receptor. Priority: manually specified > geometry > name.

    Geometry is the reliable criterion: a co-crystallized control shares the coordinate system of ITS
    receptor, so its atoms overlap with that structure's and stay far from the others. The name does
    not help, because the control is usually named after its ligand (control_ZI9) and the receptor
    after its PDB (8HTB), with no word in common.
    """
    out = dict(given)
    stems = [Path(r).stem for r in receptors]
    for c in controls:
        ck = sc.normalize_key(Path(c).stem)
        if ck in out:
            continue
        if len(stems) == 1:
            out[ck] = stems[0]
            continue
        hit = _receptor_by_geometry(c, receptors)
        if hit is None:
            hit = next((Path(r).stem for r in receptors
                        if any(len(t) >= 4 and t.lower() in Path(c).stem.lower()
                               for t in re.split(r"[_\-.]", Path(r).stem))), None)
        if hit:
            out[ck] = hit
    return out


def _receptor_by_geometry(control, receptors) -> Optional[str]:
    """Receptor in whose space the control falls, or None. Compares the control centroid with each
    receptor's atoms: it belongs to the nearest, and only if it is inside (below a threshold), so as
    not to force the assignment of a loose control."""
    try:
        cc = dk.coords_from_file(control)
    except Exception:
        cc = []
    if not cc:
        return None
    n = len(cc)
    centro = [sum(p[i] for p in cc) / n for i in range(3)]
    best_, best_d = None, 1e18
    for r in receptors:
        try:
            pts = dk._coords(r)
        except Exception:
            continue
        if not pts:
            continue
        d = min((centro[0] - p[0]) ** 2 + (centro[1] - p[1]) ** 2 + (centro[2] - p[2]) ** 2
                for p in pts)
        if d < best_d:
            best_, best_d = Path(r).stem, d
    return best_ if best_d <= 64.0 else None


def _split_peptides(proj, ligands_) -> tuple:
    """({name: (sequence, cyclic)} for ADCP, [paths] for Vina).

    The provenance is read from ligands_meta.csv (the interface writes it when building the library):
    it survives restarts and tells a peptide from any ligand without deducing it from the structure.
    Only those within its range go to ADCP; below the minimum, Vina is practicable.
    """
    meta_p = Path(proj) / "ligands_meta.csv"
    if not (adcp.available() and meta_p.exists()):
        return {}, list(ligands_)
    try:
        meta = sc.normalize_columns(pd.read_csv(meta_p))
    except Exception:
        return {}, list(ligands_)
    if "source" not in meta.columns or "product" not in meta.columns:
        return {}, list(ligands_)

    seqs = {}
    for _, r in meta[meta["source"].astype(str).str.startswith("pep")].iterrows():
        s = "".join(ch for ch in str(r.get("product") or "").upper() if ch.isalpha())
        if adcp.MIN_RESIDUES <= len(s) <= adcp.MAX_RESIDUES:
            seqs[str(r["name"])] = (s, str(r["name"]).lower().startswith("ciclo"))

    a_adcp, a_vina = {}, []
    for l in ligands_:
        name_ = Path(l).stem
        if name_ in seqs:
            a_adcp[name_] = seqs[name_]
            continue
        # A peptide control must also go to ADCP, or its energy would not be comparable with the rest.
        est = pep.sequence_from_structure(l)
        if est and adcp.MIN_RESIDUES <= len(est[0]) <= adcp.MAX_RESIDUES:
            a_adcp[name_] = est
        else:
            a_vina.append(l)
    return a_adcp, a_vina


def run(cfg: RunConfig, on_step: Optional[Callable[[str, str], None]] = None) -> RunResult:
    """Runs the full cycle. on_step(stage, detail) reports progress."""
    def step(name, detail=""):
        if on_step:
            on_step(name, detail)

    out = Path(cfg.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if not cfg.reuse:
        clean(out)
        step("cleanup", "run from scratch; nothing previous is reused")
    res = RunResult(out_dir=out)

    lig_files = [Path(p) for p in cfg.ligands]
    if cfg.lead:
        step("design", "generating analogues and predicting ADMET")
        d = AdmelabBridge().design(cfg.lead, use_ml=cfg.use_ml,
                                   n_substitutions=list(cfg.n_substitutions),
                                   max_rows=cfg.n_analogs)
        res.analogs = d.to_dataframe()
        smis = d.smiles()[: cfg.n_analogs]
        step("design", f"{d.n_generated} generated, {len(smis)} selected")
        step("3d", "generating conformers")
        made = lg.materialize(smis, out / "ligands", seed=cfg.seed)
        if res.analogs is not None and not res.analogs.empty:
            names = {s: n for n, _, s in made}
            res.analogs.insert(0, "name", res.analogs["SMILES"].map(names) if "SMILES" in res.analogs else None)
            res.analogs.to_csv(lay.artifact(out, lay.ANALOGUES_CSV), index=False)
        lig_files += [p for _, p, _ in made]
        step("3d", f"{len(made)} 3D structures")

    dock_set = list(lig_files) + [Path(c) for c in cfg.controls]
    if not dock_set:
        raise ValueError("Nothing to dock: give a lead molecule or ligands.")
    if not cfg.receptors:
        raise ValueError("No receptors.")

    # Box priority: specified > co-crystallized control > automatic. The control marks the real site.
    control_keys = {sc.normalize_key(Path(c).stem) for c in cfg.controls}
    control_assign = _assign_controls(cfg.controls, cfg.receptors, cfg.control_map)
    boxes = dict(cfg.boxes)
    for r in cfg.receptors:
        if str(r) in boxes:
            step("box", f"{Path(r).name}: specified by the user")
            continue
        ctrl = next((c for c in cfg.controls
                     if control_assign.get(sc.normalize_key(Path(c).stem)) == Path(r).stem), None)
        if ctrl:
            boxes[str(r)] = dk.box_from_file(ctrl)
            step("box", f"{Path(r).name}: centered on the control {Path(ctrl).name}")
        else:
            boxes[str(r)] = dk.auto_box(r)
            step("box", f"{Path(r).name}: automatic, no control assigned")

    import re as _re
    def _sid(stem, label):
        return f"{stem}~{_re.sub(r'[^A-Za-z0-9]+', '', str(label))[:24] or 'site'}"
    def _dedup(pares):
        seen_items, out = set(), []
        for lab, bx in pares:
            k = tuple(round(v, 1) for v in (bx.cx, bx.cy, bx.cz, bx.sx, bx.sy, bx.sz))
            if k not in seen_items:
                seen_items.add(k); out.append((lab, bx))
        return out
    targets, stem_to_file, site_box, sites_of = [], {}, {}, {}
    for r in cfg.receptors:
        sb = cfg.site_boxes.get(str(r)) if cfg.site_boxes else None
        pares = _dedup(sb) if sb else []
        if len(pares) <= 1:
            lst = [(Path(r).stem, pares[0][1] if pares else boxes[str(r)])]
        else:
            lst = [(_sid(Path(r).stem, lab), bx) for lab, bx in pares]
        sites_of[str(r)] = [sid for sid, _ in lst]
        for sid, bx in lst:
            targets.append((r, sid, bx)); stem_to_file[sid] = str(r); site_box[sid] = bx
    hibrido = any(len(v) > 1 for v in sites_of.values())

    step("docking", f"{len(dock_set)} ligands x {len(targets)} site(s)"
                    + (" (hybrid)" if hibrido else ""))
    engine = dk.VinaEngine(cpu=cfg.cpu, seed=cfg.seed, exhaustiveness=cfg.exhaustiveness,
                           n_poses=cfg.n_poses, energy_range=cfg.energy_range)
    # The ligand decides the engine, not the user: peptides go to ADCP, Vina does not cover their flexibility.
    peptides, remaining = _split_peptides(out, dock_set)
    rows, errors = [], []
    if peptides:
        step("docking", f"{len(peptides)} peptide(s) -> ADCP; "
                        f"{len(remaining)} ligand(s) -> Vina")
        # ADCP parallelizes its replicas and stays reproducible; in Vina more than one thread breaks determinism.
        import os as _os
        _nuc = max(1, (_os.cpu_count() or 2) - 2)
        f_ad, e_ad = adcp.dock_sites(targets, peptides, out,
                                      {sid: r for r, sid, _b in targets},
                                      n_poses=cfg.n_poses, n_cores=_nuc, seed_=cfg.seed,
                                      n_steps=cfg.adcp_steps,
                                      n_replicas=max(cfg.adcp_replicas, cfg.n_poses))
        rows += f_ad
        errors += e_ad
    if remaining:
        def _dock_progress(done, total, base, _err):
            if done == 0:
                step("docking", f"{total} job(s), {base.split('=')[-1]} at a time · "
                                f"exhaustiveness={cfg.exhaustiveness}, {cfg.n_poses} pose(s), "
                                f"{cfg.cpu} thread(s) per docking")
            elif done == total or done % max(1, total // 20) == 0:
                step("docking-progress", f"{done}/{total}")

        f_vi, e_vi = dk.dock_batch(cfg.receptors, remaining, boxes, out, engine=engine,
                                   workers=cfg.workers, ph=cfg.ph, targets=targets,
                                   on_progress=_dock_progress)
        for f in f_vi:
            f.setdefault("engine", "vina")
        rows += f_vi
        errors += e_vi
    res.docking = pd.DataFrame(rows)
    if not res.docking.empty:
        res.docking.to_csv(lay.artifact(out, lay.DOCKING_CSV), index=False)
    step("docking", f"{len(rows)} poses" + (f", {len(errors)} failed" if errors else ""))
    if errors:
        res.docking_errors = [(b, e) for b, e in errors]
        for base, err in errors[:8]:
            step("docking: failure", f"{base}: {err}")
    if rows:
        import collections
        by_site = collections.Counter(r["receptor"] for r in rows)
        expected_items = len(dock_set)
        for sid, _n in by_site.items():
            distinct = len({r["compound_name"] for r in rows if r["receptor"] == sid})
            if distinct < expected_items:
                step("docking: warning",
                     f"{sid}: {distinct} of {expected_items} compounds with poses")

    step("complexes", "fusing receptor and pose")
    ix.fuse_batch(cfg.receptors, out / "poses", lay.artifact(out, lay.COMPLEXES), cache_dir=out / "prep",
                  stem_to_file=stem_to_file)
    complexes = sorted(lay.artifact(out, lay.COMPLEXES).glob("*.pdb"))
    step("plip", f"{len(complexes)} complexes")
    res.interactions = ix.plip_batch(complexes, out, workers=cfg.workers, force=not cfg.reuse)

    step("scoring", "crystallographic reference + pocket interactions")
    base_crystal = ix.crystal_fingerprints(cfg.receptors, cfg.controls, control_assign,
                                           out / "xtal", cache_dir=out / "prep") if cfg.controls else {}
    # Each site takes its reference from its own control, its cofactor or its catalytic residues.
    crystal_feats, ref_src = {}, {}
    for r in cfg.receptors:
        stem = Path(r).stem
        bf = base_crystal.get(stem)
        ctrl = next((c for c in cfg.controls
                     if control_assign.get(sc.normalize_key(Path(c).stem)) == stem), None)
        cc = None
        if ctrl:
            try:
                pts = dk.coords_from_file(ctrl)
                cc = [sum(a) / len(a) for a in zip(*pts)]
            except Exception:
                cc = None
        for sid in sites_of[str(r)]:
            bx = site_box[sid]
            in_box = cc is not None and (abs(cc[0] - bx.cx) <= bx.sx / 2
                                          and abs(cc[1] - bx.cy) <= bx.sy / 2
                                          and abs(cc[2] - bx.cz) <= bx.sz / 2)
            if bf and in_box:
                crystal_feats[sid] = bf; ref_src[sid] = f"crystallographic ({Path(ctrl).stem})"
            else:
                feats, etq = ix.hetero_fingerprint(r, bx, out / "xtal")
                if feats:
                    crystal_feats[sid] = feats; ref_src[sid] = f"cofactor {etq}"
                else:
                    ref_src[sid] = "pocket residues"
    if crystal_feats:
        step("reference", f"ligand fingerprint in {len(crystal_feats)} of {len(site_box)} site(s)")
    pocket_res = {}
    for sid, bx in site_box.items():
        pocket_res[sid] = sorted(dk.residues_in_box(stem_to_file[sid], bx))
    inter = sc.prepare_interactions(res.interactions, control_keys)
    dc = res.docking.copy()
    if not dc.empty:
        dc["ckey"] = dc["compound_name"].apply(sc.normalize_key)
    step("confidence", "geometric stability of poses (obrms)")
    pose_stab = sc.pose_stability_map(dc, out / "poses")

    cnn_map = {}
    if cfg.rescoring_cnn and not dc.empty and dk.gnina_available():
        targets_to_score, pose_receptors = [], {}
        for (sid, ck), sub in dc.groupby(["receptor", "ckey"]):
            best_ = sub.sort_values("docking_score").iloc[0]["pose_name"]
            p = out / "poses" / f"{best_}.pdb"
            if p.exists():
                targets_to_score.append((ck, p))
                pose_receptors[ck] = stem_to_file.get(sid, cfg.receptors[0])
        step("rescoring", f"neural network (gnina) over {len(targets_to_score)} poses")
        for ck, p in targets_to_score:
            r = dk.rescore_poses(pose_receptors[ck], [p])
            for _name, vals in r.items():
                cnn_map[ck] = vals
        step("rescoring", f"{len(cnn_map)} compounds re-scored")
    elif cfg.rescoring_cnn and not dk.gnina_available():
        step("rescoring", "gnina not available; the second scoring is skipped")
    ref_info, icols, dscore = sc.build_ref_info(inter, dc, control_keys, control_assign,
                                                crystal_feats=crystal_feats)
    cat = {r: (cfg.catalytic.get(r) or cfg.catalytic.get(sc.base_of(r))
               or ref_info.get(r, {}).get("autocat", [])) for r in ref_info}
    res.ref_info = ref_info

    smiles_map = {}
    if res.analogs is not None and not res.analogs.empty and {"name", "SMILES"} <= set(res.analogs.columns):
        for _, r in res.analogs.iterrows():
            if r.get("name"):
                smiles_map[sc.normalize_key(r["name"])] = r["SMILES"]
    for folder in {Path(p).parent for p in list(lig_files) + [Path(c) for c in cfg.controls]}:
        for k, v in sc.build_smiles_map(folder).items():
            smiles_map.setdefault(k, v)

    res.ranking = sc.compute_ranking(inter, dc, control_keys, control_assign,
                                     ref_info, icols, dscore, cat, cfg.weights,
                                     smiles_map=smiles_map, pocket_res_map=pocket_res,
                                     pose_stability=pose_stab, cnn_map=cnn_map)

    if res.analogs is not None and not res.analogs.empty and "name" in res.analogs.columns:
        adm = res.analogs.copy()
        adm["ckey"] = adm["name"].map(lambda n: sc.normalize_key(n) if n else None)
        cols_admet = [c for c in ("MW", "LogP", "TPSA", "QED", "Lipinski_violations",
                                  "LD50_mg_per_kg", "GHS_category") if c in adm.columns]
        if cols_admet:
            res.ranking = res.ranking.merge(adm[["ckey"] + cols_admet].drop_duplicates("ckey"),
                                            on="ckey", how="left")
    res.ranking.to_csv(out / "ranking.csv", index=False)

    import json
    (out / "run.json").write_text(json.dumps({
        "lead": cfg.lead,
        # File names only: run.json travels inside sessions and export packages, and an absolute
        # path would carry the user's folder layout with it.
        "receptors": [Path(r).name for r in cfg.receptors],
        "controls": [Path(c).name for c in cfg.controls],
        "control_assign": control_assign,
        "control_keys": sorted(control_keys),
        "catalytic": cat,
        "reference": {r: ref_info.get(r, {}).get("src") for r in ref_info},
        "crystal_feats": {r: sorted(f) for r, f in crystal_feats.items()},
        "pocket_residues": pocket_res,
        "sites": {Path(r).name: sids for r, sids in sites_of.items()},
        "site_reference": ref_src,
        "hibrido": hibrido,
        "pose_stability": pose_stab,
        "weights": cfg.weights,
        "seed": cfg.seed,
        "exhaustiveness": cfg.exhaustiveness,
        "n_poses": cfg.n_poses,
        "energy_range": cfg.energy_range,
        "ph": cfg.ph,
        "rescoring_cnn": bool(cnn_map),
    }, indent=2))

    cols = [c for c in ("receptor", "compound", "best_dock", "pred_ki_M", "best_inter",
                        "cat_coverage", "effectiveness_pct", "key_interaction", "type")
            if c in res.ranking.columns]
    res.ranking[cols].to_csv(lay.artifact(out, lay.SUMMARY_CSV), index=False)

    if cfg.controls:
        step("validation", "redocking the controls")
        res.validation = vl.redock_validation(cfg.controls, control_assign, out / "poses")
        res.validation.to_csv(lay.artifact(out, lay.VALIDATION_CSV), index=False)
        step("validation", vl.summary_text(res.validation))

    step("done", str(out))
    return res
