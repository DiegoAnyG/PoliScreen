"""Stage 4: Results view (ranking, PLIP interactions, 3D complex viewer, methods export, and tunnels)."""
from __future__ import annotations

import io
import json
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

from poliscreen.core import caver as cv
from poliscreen.core import layout as lay
from poliscreen.core import pipeline as pl
from poliscreen.core import reagents as rg
from poliscreen.core import report as rp
from poliscreen.core import screening as sc
from poliscreen.core import tunnels as tn
from poliscreen.core import validation as vl
from poliscreen.core import viewer as vw
from poliscreen.core.design import AdmelabBridge
from poliscreen.ui.common import (
    _download_image,
    _download_table,
    _empty_state,
    _fmt_ki,
    _how_to_cite,
    _rname,
    _scene_height,
    _viewer_height,
)
from poliscreen.ui.components.admet import _render_adme, _scatter_dock_inter, _shade
from poliscreen.ui.components.transport import (
    _one_row_per_route,
    _readable_transport,
    _route_preferences,
    _with_tunnel_geometry,
    route_preference,
)
from poliscreen.ui.i18n import t


def _pleiotropic_summary(rk, targets_):
    """Which compound binds well in SEVERAL targets at once."""
    sub = rk[rk.get("is_control", 0) != 1].copy() if "is_control" in rk.columns else rk.copy()
    sub = sub[pd.notna(sub.get("effectiveness_pct"))]
    if sub.empty:
        return
    sub["_target"] = sub["receptor"].map(lambda r: sc.display_name(sc.base_of(r)))
    best_ = (sub.groupby(["compound", "_target"])["effectiveness_pct"].max()
             .reset_index())
    piv = best_.pivot(index="compound", columns="_target", values="effectiveness_pct")
    presente_en_todas = piv.dropna()
    st.markdown(t("### Pleiotropic summary — activity across several targets"))
    if presente_en_todas.empty:
        st.caption(t("No compound docked in all targets; there is no broad-spectrum comparison."))
        return
    presente_en_todas = presente_en_todas.assign(
        **{"minimum": presente_en_todas.min(axis=1).round(1),
           "mean": presente_en_todas.mean(axis=1).round(1)})
    table_ = presente_en_todas.sort_values("minimum", ascending=False).reset_index()
    table_.columns = ["compound"] + [f"{c} (%)" if c in targets_ else c for c in table_.columns[1:]]
    st.caption(t("Effectiveness (%) of each compound in each target, ordered by the **minimum** across targets: broad-spectrum ones on top. Only those docked in all."))
    st.dataframe(table_.round(1), width="stretch", hide_index=True,
                 height=min(340, 60 + 34 * len(table_)))
    _download_table(table_, "pleiotropico", key="pleio")
    best_broad = table_.iloc[0]
    st.success(t('Best broad-spectrum: **{v1}** (minimum {v3:.0f} % across {v5} targets).').format(v1=best_broad['compound'], v3=best_broad['minimum'], v5=len(targets_)))
    st.divider()


def _results_tunnels(proj: Path):
    """CAVER and CaverDock output, read into a results table."""
    S = st.session_state
    st.subheader(t("Transport tunnels"))

    if not tn.available():
        st.info(t("Tunnel reading needs **caver-translate**, a separate package with no "
                  "dependencies of its own. The rest of PoliScreen works without it."))
        st.code(tn.INSTALL_HINT, language="bash")
        return

    here = proj / lay.TUNNELS
    S.setdefault("tun_folder", str(here) if here.is_dir() else "")
    folder_str = st.text_input(t("Results folder"), key="tun_folder",
                               help=t("This project's runs by default. Point it elsewhere to read "
                                      "a CaverWeb download."))
    if not folder_str.strip():
        st.info(t("No transport calculated yet. Run one in step 3."))
        return

    folder = Path(folder_str.strip().strip('"'))
    if not folder.is_dir():
        st.error(t("Not a folder: {p}").format(p=folder))
        return

    try:
        table, _cov = tn.read(folder)
    except Exception as e:
        st.error(t("Could not read that folder: {e}").format(e=e))
        return

    if table.empty:
        st.info(t("Nothing to read there yet."))
        return

    table = _one_row_per_route(_with_tunnel_geometry(table, folder))
    failed = int(table["flags"].fillna("").str.contains("failed").sum())

    c = st.columns(3)
    c[0].metric(t("Calculations"), len(table))
    c[1].metric(t("Routes"), int(table["tunnel"].nunique(dropna=True)))
    c[2].metric(t("Did not finish"), failed)

    _route_preferences(table)

    shown = table.sort_values("Ea", na_position="last")
    st.dataframe(_readable_transport(shown), width="stretch", hide_index=True)
    st.caption(t("Ea is what entering costs and compares tunnels. dE_BS is how much better the "
                 "site is than the outside. Detail in Help › Transport tunnels."))
    _download_table(shown, "tuneles", key="tunnels")

    if st.button(t("Write the full report"), key="tun_export"):
        out = tn.export(folder, proj / lay.TUNNELS / "report")
        st.success(t("Written to {p}").format(p=out))
        page = out / "report.html"
        if page.exists():
            st.download_button(t("report.html"), page.read_bytes(), file_name="report.html",
                               mime="text/html", key="tun_html")


def _results_screening(proj: Path):
    """Main results screening panel: metrics, ranking, and 2D interaction diagrams."""
    S = st.session_state
    st.subheader(t("Results"))
    meta_p = proj / "run.json"
    inter_p = lay.artifact(proj, lay.INTERACTIONS_CSV)
    dock_p = lay.artifact(proj, lay.DOCKING_CSV)
    if not (meta_p.exists() and inter_p.exists()):
        st.info(t("No results in this folder yet. Run step 3."))
        return

    meta = json.loads(meta_p.read_text())
    inter_raw = pd.read_csv(inter_p)
    dc = pd.read_csv(dock_p) if dock_p.exists() else pd.DataFrame()
    ckeys = {sc.normalize_key(Path(c).stem) for c in meta.get("controls", [])}
    cassign = meta.get("control_assign", {})
    inter = sc.prepare_interactions(inter_raw, ckeys)
    if not dc.empty:
        dc["ckey"] = dc["compound_name"].apply(sc.normalize_key)
    pocket_res_map = meta.get("pocket_residues", {})
    ref_info, icols, dscore = sc.build_ref_info(inter, dc, ckeys, cassign,
                                                crystal_feats=meta.get("crystal_feats"))

    st.markdown(t("**Catalytic / anchor residues**. The score also rewards the quality of the pocket's other interactions."))
    st.caption(t("Auto-suggested from the directional interactions of the crystallographic ligand. Edit them if you know your target's real catalytic site."))
    cat, sec = {}, {}
    cols = st.columns(max(1, len(ref_info)))
    for i, R in enumerate(sorted(ref_info)):
        options = sorted(set(ref_info[R].get("residues", [])) | set(pocket_res_map.get(R, [])),
                          key=lambda r: (sc.resnum(r), r))
        suggested = ref_info[R].get("autocat", [])
        prev = ([x for x in meta.get("catalytic", {}).get(R, []) if x in options]
                or [x for x in suggested if x in options])
        cat[R] = cols[i].multiselect(f"{R}  ·  ref: {ref_info[R].get('src', '?')}", options,
                                     default=prev, key=f"cat_{R}")
        free_slots = [x for x in options if x not in cat[R]]
        k_sec = f"sec_{R}"
        if k_sec in S:
            S[k_sec] = [x for x in S[k_sec] if x in free_slots]
        prev_s = [x for x in meta.get("secondary", {}).get(R, []) if x in free_slots]
        sec[R] = cols[i].multiselect(t('{v0} · secondary (bonus, not required)').format(v0=R), free_slots,
                                     default=prev_s, key=f"sec_{R}")

    val_p = lay.artifact(proj, lay.VALIDATION_CSV)
    if val_p.exists():
        val = vl.normalize(pd.read_csv(val_p))
        _v = vl.summary(val)
        if _v["ok"] is None:
            msg = t("No controls: the setup cannot be validated.")
        elif _v["ok"]:
            msg = t("The control recovers the crystallographic pose: the setup is reliable.") \
                if _v["n"] == 1 else \
                t("The {n} controls recover the crystallographic pose: the setup is reliable.").format(n=_v["n"])
        else:
            _who = (t("The control does not recover") if _v["n"] == 1
                    else t("{m} of {n} controls do NOT recover").format(m=_v["n_failing"], n=_v["n"]))
            msg = t("WARNING: {who} the pose ({targets}). Check the box or the preparation of "
                    "that target before trusting the ranking.").format(
                        who=_who, targets=", ".join(_v["targets"]))
        (st.success if _v["ok"] is not False else st.error)(msg)
        with st.expander(t("Redocking validation detail")):
            st.dataframe(val.assign(target=val["target"].map(_rname)),
                         width="stretch", hide_index=True)
            st.caption(t("RMSD against the co-crystallized ligand. Valid below 2 Å."))

    crystal = meta.get("crystal_feats", {})
    if crystal:
        rows_ = []
        for R in sorted(crystal):
            cks = {ck for ck, rc in cassign.items() if rc == R} or ckeys
            sub = inter[(inter["receptor"] == R) & (inter["ckey"].isin(cks))]
            if sub.empty:
                continue
            s = sub["name"].map(lambda n: dscore.get(sc.pose_key(n), float("nan")))
            best = sub.loc[s.idxmin()] if s.notna().any() else sub.iloc[0]
            pose_feats = [c for c in icols if best[c] > 0]
            rec = sc.fp_recovery(crystal[R], pose_feats)
            rows_.append({"receptor": R, "recovery": rec["recovery"], "Tanimoto": rec["tanimoto"],
                          "reproduced": f"{rec['shared']}/{rec['ref_n']}", "extra (non-crystal)": rec["extra"]})
        if rows_:
            st.markdown(t("**Interaction validation** — docked control vs. crystallographic ligand."))
            st.dataframe(pd.DataFrame(rows_), width="stretch", hide_index=True)
            st.caption(t("`recovery` = fraction of the crystallographic interactions reproduced by the control's docked pose; `Tanimoto` also includes the extra contacts docking adds. "))

    st.markdown(t("**Weighting**"))
    mw = meta.get("weights", {})
    metric_afin = st.radio(
        t("Affinity-axis metric"), ["dock", "le"], horizontal=True,
        index=1 if str(mw.get("dock_metric", "dock")).lower() == "le" else 0,
        format_func=lambda m: t("Raw score (kcal/mol)") if m == "dock" else t("Ligand efficiency (LE)"),
        help=t("Vina's raw score favors large molecules (size bias). LE = -ΔG/heavy atoms corrects it. Recommended if your library varies widely in size; both columns are reported."))
    c1, c2, c3, c4 = st.columns(4)
    w_dock = c1.slider(t("Docking weight"), 0.0, 1.0, float(mw.get("dock", 0.5)), 0.05)
    w_inter = c2.slider(t("Interactions weight"), 0.0, 1.0, float(mw.get("inter", 0.5)), 0.05)
    w_adme = c3.slider(t("ADME weight"), 0.0, 1.0, float(mw.get("adme", 0.0)), 0.05,
                       help=t("Physicochemical (drug-likeness) quality of the compound. Guards against rewarding only large/greasy molecules."))
    w_tox = c4.slider(t("Toxicity weight"), 0.0, 1.0, float(mw.get("tox", 0.0)), 0.05,
                      help=t("Requires ADMET predicted (Ligands tab); otherwise this axis is ignored."))
    c5, c6, c7 = st.columns(3)
    w_cat = c5.slider(t("Catalytic-residue weight"), 1.0, 6.0, float(mw.get("w_cat", 3.0)), 0.5,
                      help=t("How much an interaction with a catalytic (gate) residue is worth vs. an ordinary pocket one."))
    w_sec = c6.slider(t("Secondary-residue weight"), 1.0, 3.0, float(mw.get("w_sec", 1.5)), 0.25,
                      help=t("How much an interaction with a SECONDARY anchor is worth vs. an ordinary pocket contact (×1)."))
    cat_gate = c7.slider(t("Catalytic strictness"), 0.0, 1.0, float(mw.get("cat_gate", 0.5)), 0.05,
                         help=t("0 = missing a catalytic residue is not penalized; 1 = missing all nullifies the score."))
    _axw = {"docking": w_dock, "interaction": w_inter, "ADME": w_adme, "tox": w_tox}
    _tot = sum(_axw.values())
    if _tot > 0:
        st.caption("Real contribution of each axis: "
                   + " · ".join(f"{k} {v / _tot * 100:.0f}%" for k, v in _axw.items() if v > 0))
    else:
        st.warning(t("All axis weights are 0: there will be no score. Raise at least one."))
    with st.expander(t("Weights by interaction type (advanced)")):
        st.caption(t("Merit value per type (0-1). Default: salt bridge > H-bond > π > halogen > hydrophobic. Literature-guided; adjust to your judgment."))
        st.caption(t("`water` (water-mediated bridges) only matters if you keep water molecules when preparing the receptor. In the usual flow they are removed, so this weight has no effect."))
        tw = {}; tcols = st.columns(4)
        for j, (tk, tv) in enumerate(sc.TYPE_WEIGHTS.items()):
            tw[tk] = tcols[j % 4].number_input(
                tk, 0.0, 1.0, float((mw.get("type_weights") or {}).get(tk, tv)), 0.05, key=f"tw_{tk}")
    w = dict(pl.DEFAULT_WEIGHTS)
    w.update(dock=w_dock, inter=w_inter, adme=w_adme, tox=w_tox, dock_metric=metric_afin,
             w_cat=w_cat, w_sec=w_sec, cat_gate=cat_gate, type_weights=tw)
    S["_ui_weights"] = w

    smap = sc.build_smiles_map(str(lay.artifact(proj, lay.INPUT_LIGANDS)))
    for _k, _v in sc.build_smiles_map(str(lay.artifact(proj, lay.RECEPTORS))).items():
        smap.setdefault(_k, _v)
    ml0 = proj / "ligands_meta.csv"
    if ml0.exists():
        _m0 = pd.read_csv(ml0)
        for _n, _s in zip(_m0.get("name", []), _m0.get("smiles", [])):
            k0 = sc.normalize_key(_n)
            if k0 and pd.notna(_s) and k0 not in smap:
                smap[k0] = str(_s)
    pose_stab = meta.get("pose_stability", {})
    reliable_map = {}
    _vp = lay.artifact(proj, lay.VALIDATION_CSV)
    if _vp.exists():
        _v = vl.normalize(pd.read_csv(_vp))
        for _, _r in _v.iterrows():
            _rm = pd.to_numeric(pd.Series([_r.get("rmsd_min_A")]), errors="coerce").iloc[0]
            if pd.notna(_r.get("target")):
                reliable_map[str(_r["target"])] = bool(pd.notna(_rm) and _rm < 2.0)
    rk = sc.compute_ranking(inter, dc, ckeys, cassign, ref_info, icols, dscore, cat, w,
                            smiles_map=smap, pocket_res_map=pocket_res_map, sec_map=sec,
                            pose_stability=pose_stab, reliable_map=reliable_map)
    S["_rk_live"] = rk.copy()
    rk["Ki"] = rk["pred_ki_M"].map(_fmt_ki) if "pred_ki_M" in rk.columns else None
    faltan = []
    if w_adme > 0 and pd.to_numeric(rk.get("adme"), errors="coerce").isna().all():
        faltan.append("ADME")
    if w_tox > 0 and pd.to_numeric(rk.get("ld50_mgkg"), errors="coerce").isna().all():
        faltan.append("toxicidad")
    if faltan:
        st.warning(t('You weight **{v1}** but there is no data for that axis in this run: it is ignored in the score. Predict ADMET first, or lower its weight to 0 so Methods does not declare it.').format(v1=' and '.join(faltan)))

    with st.expander(t("Export Methods (for the paper)")):
        st.caption(t("Parameters, box, weights, reference and exact software versions. Reproducibility ready to paste into the Methods section."))
        methods_text_ = rp.methods_text(meta, weights=w, catalytic=cat, secondary=sec)
        st.download_button(t("Download Methods.md"), methods_text_, file_name="PoliScreen_Methods.md",
                           mime="text/markdown")
        st.code(methods_text_, language="markdown")

    engines_ = {}
    _dock_p = lay.artifact(proj, lay.DOCKING_CSV)
    if _dock_p.exists():
        try:
            _dd = pd.read_csv(_dock_p)
            if "engine" in _dd.columns:
                engines_ = {sc.normalize_key(c): m for c, m in
                           zip(_dd["compound_name"], _dd["engine"]) if pd.notna(m)}
        except Exception:
            engines_ = {}
    if engines_:
        rk["engine"] = rk["compound"].map(lambda c: engines_.get(sc.normalize_key(c), ""))
        _distintos = sorted({m for m in rk["engine"] if m})
        if len(_distintos) > 1:
            st.warning(
                t('This table mixes **{v1}**. Their energies come from different scoring functions and are not comparable: `best_dock`, `pKi` and `LE` only make sense within each engine. To compare across them use `effectiveness_pct`, computed from the contacts and independent of the engine.').format(v1=' and '.join(_distintos)))

    chosen_items = [c for c in ("compound", "IUPAC", "engine", "pose", "best_dock", "pKi", "LE",
                            "best_inter", "cat_coverage", "effectiveness_pct", "percentile",
                            "confidence", "consensus", "key_interaction", "sa_score", "pains",
                            "type")
                if c in rk.columns]

    meta_lig = proj / "ligands_meta.csv"
    tuyos, iupac_map, real_name = set(), {}, {}
    if meta_lig.exists():
        m = pd.read_csv(meta_lig)
        tuyos = {sc.normalize_key(n) for n, f in zip(m.get("name", []), m.get("source", []))
                 if str(f) == "yours"}
        if "iupac_name" in m.columns:
            iupac_map = {sc.normalize_key(n): v for n, v in zip(m.get("name", []), m.get("iupac_name", []))
                         if pd.notna(v) and str(v).strip()}
        real_name = {sc.normalize_key(n): n for n in m.get("name", []) if pd.notna(n)}
    if real_name:
        rk["compound"] = rk["compound"].map(lambda c: real_name.get(sc.normalize_key(c), c))
    if iupac_map:
        rk["IUPAC"] = rk["compound"].map(lambda c: iupac_map.get(sc.normalize_key(c), ""))
    if tuyos:
        rk["source"] = rk["compound"].map(lambda c: "yours" if sc.normalize_key(c) in tuyos else "")
        st.caption(t("Highlighted rows = compounds made with reagents you provided."))

    if any("~" in str(x) for x in rk["receptor"].unique()):
        st.info(t("**Hybrid docking**: each block is a different pocket of the same receptor. Compare a compound's effectiveness across sites to see where it prefers to bind."))
    _dianas = sorted({sc.display_name(sc.base_of(x)) for x in rk["receptor"].unique()})
    if len(_dianas) > 1:
        _pleiotropic_summary(rk, _dianas)
    for R in sorted(rk["receptor"].unique()):
        sub = rk[rk["receptor"] == R].copy()
        _rn = _rname(R)
        _et = (f"{_rn.split('~')[0]} · " + t("site") + f" **{_rn.split('~', 1)[1]}**"
               if "~" in _rn else f"**{_rn}**")
        _refsrc = meta.get("site_reference", {}).get(R) or ref_info.get(R, {}).get("src", "?")
        st.markdown(t('{v0} · interaction reference: `{v2}`').format(v0=_et, v2=_refsrc))
        noc = sub[sub["is_control"] != 1]
        if not noc.empty:
            m1, m2, m3, m4 = st.columns(4)
            try:
                bd = noc.loc[pd.to_numeric(noc["best_dock"], errors="coerce").idxmin()]
                m1.metric(t("Best docking"), str(bd["compound"])[:18], f"{bd['best_dock']:.2f} kcal/mol",
                          delta_color="inverse")
            except Exception:
                pass
            try:
                bi = noc.loc[pd.to_numeric(noc["best_inter"], errors="coerce").idxmax()]
                m2.metric(t("Best interaction quality"), str(bi["compound"])[:18], f"{bi['best_inter']:.2f}")
            except Exception:
                pass
            try:
                be = noc.loc[pd.to_numeric(noc["effectiveness_pct"], errors="coerce").idxmax()]
                m3.metric(t("Best effectiveness"), str(be["compound"])[:18], f"{be['effectiveness_pct']:.0f} %")
            except Exception:
                pass
            try:
                bc = noc.loc[pd.to_numeric(noc["confidence"], errors="coerce").idxmax()]
                m4.metric(t("Highest confidence"), str(bc["compound"])[:18], f"{bc['confidence']:.2f}")
            except Exception:
                pass
        view_ = sub[chosen_items]
        st.dataframe(_shade(view_.assign(source=sub.get("source", "")), "source") if tuyos else view_,
                     width="stretch", height=min(400, 60 + 34 * len(sub)))
        _download_table(view_, f"ranking_{R}", key=f"rk_{R}")
        g1, g2 = st.columns(2)
        ch = sub.dropna(subset=["effectiveness_pct"]).set_index("compound")["effectiveness_pct"]
        if not ch.empty:
            g1.bar_chart(ch, height=260)
        fig = _scatter_dock_inter(sub)
        if fig:
            g2.pyplot(fig)
    st.download_button(t("Download ranking (CSV)"), rk.to_csv(index=False).encode(), "ranking.csv")

    items_all = [(c, smap[sc.normalize_key(c)]) for c in rk["compound"].unique()
                 if sc.normalize_key(c) in smap and pd.notna(smap.get(sc.normalize_key(c)))]
    if items_all:
        with st.expander(t("ADMET report (compounds + core + control, those you choose)")):
            all_names = [c for c, _ in items_all]
            chosen_adme = st.multiselect(t("Which ligands to predict ADMET for?"), all_names,
                                           default=all_names, key="adme_sel_res")
            items = [(c, s) for c, s in items_all if c in chosen_adme]
            if st.button(t("Predict ADMET"), key="pred_res") and items:
                with st.spinner(t('Predicting with ADMET-AI for {v1} ligand(s)...').format(v1=len(items))):
                    pr = AdmelabBridge().predict([s for _, s in items], use_ml=True)
                S["admet"] = {**(S.get("admet") or {}), **{rg.inchikey(r.get("SMILES")): r for r in pr.rows}}
            if S.get("admet") and items:
                _render_adme(S["admet"], items, keyp="res")

    st.markdown("---")
    st.markdown(t("**Interaction diagram** of a specific pose."))
    d1, d2, d3 = st.columns(3)
    R = d1.selectbox(t("Receptor"), sorted(inter["receptor"].unique()))
    sr = inter[inter["receptor"] == R]
    compounds = sorted(sr["compound"].unique())
    first = 0
    if "is_control" in sr.columns:
        controls = [c for c in compounds if bool(sr[sr["compound"] == c]["is_control"].any())]
        if controls:
            first = compounds.index(controls[0])
    cmp_ = d2.selectbox(t("Compound"), compounds, index=first,
                        help=t("Opens on the control, which is the reference the other "
                               "diagrams are judged against."))
    scmp = sr[sr["compound"] == cmp_]
    mods = sorted({sc.model_of(n) for n in scmp["name"]})
    mod = d3.selectbox(t("Pose"), mods)
    row = scmp[scmp["name"].apply(lambda n: sc.model_of(n) == mod)]
    if not row.empty:
        reference_ = ref_info.get(R, {}).get("feats", [])
        fig_int = sc.draw_2d(row.iloc[0], f"{R} · {cmp_} · pose {mod}", reference=reference_)
        st.pyplot(fig_int, width="content")
        try:
            _b = io.BytesIO(); fig_int.savefig(_b, format="png", dpi=160, bbox_inches="tight")
            _download_image(_b.getvalue(), f"interaccion_{cmp_}_pose{mod}", key=f"int_{R}_{cmp_}_{mod}")
        except Exception:
            pass
        st.caption(t("Green = reproduces a control interaction (same residue and same bond). Gray = extra contact or the same residue with a different bond type."))

    st.markdown("---")
    _how_to_cite()


def render_results_tools(proj: Path):
    """Tool panel for Results stage: Screening ranking & transport tabs."""
    tab_screening, tab_tunnels = st.tabs([t("Screening"), t("Transport tunnels")])
    with tab_tunnels:
        _results_tunnels(proj)
    with tab_screening:
        _results_screening(proj)


def _complex_viewer(proj: Path):
    """3D Complex viewer with surface and heteroatom toggles."""
    S = st.session_state
    proj_p = proj
    inter_p = lay.artifact(proj_p, lay.INTERACTIONS_CSV)
    if not inter_p.exists():
        _empty_state("Run a screening and you can browse the 3D complexes here.")
        return
    inter = pd.read_csv(inter_p)
    if "name" not in inter.columns or inter.empty:
        _empty_state("The interactions table is empty.")
        return
    inter = inter.assign(_rec=inter["name"].map(sc.receptor_from_name),
                         _cmp=inter["name"].map(sc.compound_from_pose_name),
                         _mod=inter["name"].map(sc.model_of))

    _asig, _ctrlk = {}, set()
    _rj = proj_p / "run.json"
    if _rj.exists():
        try:
            _d = json.loads(_rj.read_text())
            _asig = _d.get("control_assign") or {}
            _ctrlk = set(_d.get("control_keys") or _asig.keys())
        except Exception:
            pass

    def _visible(rec, cmp):
        ck = sc.normalize_key(cmp)
        if ck not in _ctrlk:
            return True
        return _asig.get(ck) == sc.base_of(rec)

    c1, c2, c3 = st.columns([1.3, 1.7, 0.9])
    R = c1.selectbox(t("Receptor"), sorted(inter["_rec"].unique()), format_func=_rname,
                     key="vis_cx_rec")
    sr = inter[inter["_rec"] == R]
    _cmps = sorted(c for c in sr["_cmp"].unique() if _visible(R, c))
    C = c2.selectbox(t("Compound"), _cmps or sorted(sr["_cmp"].unique()), key="vis_cx_cmp")
    scmp = sr[sr["_cmp"] == C]
    M = c3.selectbox(t("Pose"), sorted(scmp["_mod"].unique()), key="vis_cx_pose")
    o1, o2 = st.columns(2)
    S.setdefault("vis_cx_surface", False)
    sup = o1.checkbox(t("Show the surface"), key="vis_cx_surface",
                      help=t("Translucent molecular surface of the receptor. With the ribbon alone you cannot tell whether the ligand is inside the cavity or resting outside."))
    S.setdefault("vis_cx_het", True)
    het = o2.checkbox(t("Cofactors and hetero"), key="vis_cx_het")
    row_ = scmp[scmp["_mod"] == M]
    if row_.empty:
        _empty_state("There is no pose with that combination.")
        return

    name_ = str(row_.iloc[0]["name"])
    pose_f = proj_p / "poses" / f"{lay.strip_complex_prefix(name_)}.pdb"
    rec_f = next((p for p in (lay.artifact(proj_p, lay.RECEPTORS)).glob(f"{R}.*")
                  if p.suffix.lower() in (".pdb", ".pdbqt")), None)

    _h = _viewer_height(190)
    try:
        if pose_f.exists() and rec_f is not None:
            html = vw.view_html(receptor=rec_f, ligand_=pose_f, show_waters=False,
                                show_hetero=het, surface=sup, height_=_h)
        else:
            fus = lay.artifact(proj_p, lay.COMPLEXES) / f"{name_}.pdb"
            if not fus.exists():
                _empty_state("Cannot find this pose's files in the project folder.")
                return
            html = vw.view_html(receptor=fus, show_waters=False, show_hetero=het,
                                surface=sup, height_=_h)
        st.iframe(html, height=_h + 12)
    except Exception as e:
        st.error(t('Could not draw the complex: {v1}').format(v1=e))
        return

    r0 = row_.iloc[0]
    feats = [c for c in row_.columns
             if "_" in str(c) and str(c).rsplit("_", 1)[-1] in sc.TYPE_STYLE
             and pd.notna(r0[c]) and r0[c] > 0]
    _cab = f"{C} · pose {M} · {R}"
    if feats:
        _cab += f" — **{len(feats)} contactos**"
    st.caption(_cab)


def _transport_viewer(proj: Path):
    """Transport poses and energy profile."""
    S = st.session_state
    runs = tn.runs_in(proj / lay.TUNNELS)
    if not runs:
        _empty_state("Run a transport and its poses will appear here.")
        return

    run = st.selectbox(t("Calculation"), runs, format_func=tn.short_name,
                       key="vis_tun_run", label_visibility="collapsed")

    profile = tn.profile_of(run)
    if not profile:
        st.warning(t("That calculation produced no profile."))
        return

    bound = tn.orientation_of(run)
    most = tn.most_extra(profile)
    if most:
        S.setdefault("vis_tun_poses", 3 + tn.suggested_extra(profile))
        S["vis_tun_poses"] = max(3, min(int(S["vis_tun_poses"]), 3 + most))
        extra = st.slider(t("Poses"), min_value=3, max_value=3 + most, key="vis_tun_poses",
                          help=t("Three are always there: the mouth, the barrier and the site. "
                                 "The rest are context, spaced between them.")) - 3
    else:
        extra = 0
    poses = tn.chosen_poses(profile, bound=bound, extra=int(extra))

    trajectory = next(iter(sorted(Path(run).glob("*-lb.pdbqt"))), None)
    receptor = next(iter(sorted(Path(run).glob("*_ready.pdb"))), None) \
        or next(iter(sorted(Path(run).glob("*.pdb"))), None)
    tunnel = next(iter(sorted(Path(run).glob("tun_cl_*.pdb"))), None)

    try:
        blocks = tn.pose_blocks(trajectory, [s for s, _t, _l, _r in poses]) if trajectory else []
        colours, marks_3d, context = [], [], 0
        for (_s, tag, label, _r), block in zip(poses, blocks):
            if tag in tn.TAG_COLOR:
                colours.append(tn.pose_color(tag))
            else:
                colours.append(tn.pose_color(tag, context))
                context += 1
            marks_3d.append((tn.centroid_of(block), label.split("  ")[0].strip() or tag))
        spheres = [{"alpha": cv.tunnel_spheres(tunnel), "color": "#CED4DA",
                    "opacity": 0.55}] if tunnel else None
        _h = _scene_height()
        st.iframe(vw.view_html(receptor=receptor, ligand_=None, cavities=spheres,
                               show_waters=False, show_hetero=True, height_=_h,
                               extra_models=[(b, "pdb") for b in blocks if b],
                               model_colors=colours, callouts=marks_3d),
                  height=_h + 12)
    except Exception as e:
        st.error(t('Could not draw: {v1}').format(v1=e))

    marks = tn.landmarks(profile, bound=bound)
    fig = tn.draw_profile(profile, bound=bound, title=tn.short_name(run))
    st.pyplot(fig, width="stretch")
    if marks:
        e_surface = marks["surface"][1]
        cols = st.columns(4)
        cols[0].metric("E_surface", f"{e_surface:.1f}")
        cols[1].metric("E_max", f"{marks['barrier'][1]:.1f}")
        cols[2].metric("E_bound", f"{marks['site'][1]:.1f}")
        cols[3].metric("Ea", f"{marks['barrier'][1] - e_surface:.1f}",
                       delta=f"dE_BS {marks['site'][1] - e_surface:.1f}", delta_color="off")
        st.caption(t("The four numbers are lower-bound quantities, which is why they are marked "
                     "on that line."))

    out = Path(run) / "figure.pml"
    tn.pymol_script(run, receptor, tunnel, poses, out)
    st.download_button(t("Download as a PyMOL script"), out.read_bytes(),
                       file_name=f"{tn.short_name(run).replace(' · ', '_')}.pml",
                       mime="text/plain", key="pml_download",
                       help=t("The receptor, the tunnel and exactly these poses. It loads them "
                              "itself: `pymol figure.pml`."))


def _visual_summary(proj: Path):
    """Results summary for the viewer panel."""
    S = st.session_state
    live = S.get("_rk_live")
    if live is not None and not live.empty:
        rk = live.copy()
    else:
        rk_p = proj / "ranking.csv"
        if not rk_p.exists():
            _empty_state("Run a screening and the results summary will appear here.")
            return
        rk = sc.normalize_columns(pd.read_csv(rk_p))
    sites = sorted(rk["receptor"].unique()) if "receptor" in rk.columns else []
    if len(sites) > 1:
        sel = st.selectbox(t("Site"), sites, key="vis_res_site",
                           format_func=lambda s: _rname(s).split("~", 1)[-1])
        rk = rk[rk["receptor"] == sel]
        st.caption(t("Summary of site **{site}**. Switch site to compare where each compound binds best.")
                   .format(site=_rname(sel)))
    noc = (rk[rk["is_control"] != 1] if "is_control" in rk.columns else rk).copy()
    if noc.empty or "effectiveness_pct" not in noc.columns:
        st.info(t("No compounds to summarize yet."))
        return

    ef = pd.to_numeric(noc["effectiveness_pct"], errors="coerce")
    conf = pd.to_numeric(noc.get("confidence"), errors="coerce") if "confidence" in noc else None
    superan = int((ef >= 105).sum())
    fiables = int((conf >= 0.5).sum()) if conf is not None else 0

    best_ = noc.loc[ef.idxmax()]
    enc, img = st.columns([2, 1])
    enc.markdown(f"### {str(best_['compound'])[:38]}")
    enc.caption(t("Highest-effectiveness compound"))
    smi_top = None
    ml = proj / "ligands_meta.csv"
    if ml.exists():
        _m = pd.read_csv(ml)
        if {"name", "smiles"} <= set(_m.columns):
            _k = sc.normalize_key(best_["compound"])
            _hit = _m[_m["name"].map(lambda n: sc.normalize_key(n) == _k)]
            if not _hit.empty:
                smi_top = str(_hit.iloc[0]["smiles"])
    if smi_top:
        png = vw.molecule_png(smi_top, size=170)
        if png:
            img.image(png)
    k1, k2, k3 = st.columns(3)
    k1.metric(t("Effectiveness"), f"{ef.max():.0f} %")
    if "best_dock" in noc.columns:
        k2.metric(t("Affinity"), f"{float(best_['best_dock']):.1f}", "kcal/mol", delta_color="off")
    if conf is not None and pd.notna(best_.get("confidence")):
        k3.metric(t("Confidence"), f"{float(best_['confidence']):.2f}")

    st.divider()
    c1, c2, c3 = st.columns(3)
    c1.metric(t("Compounds"), len(noc))
    c2.metric(t("Beat the control"), superan, help=t("Effectiveness ≥ 105 % vs. the crystallographic ligand."))
    if conf is not None:
        c3.metric(t("Confidence ≥ 0.5"), fiables, help=t("Concordant evidence: reliable result."))

    st.markdown(t("**Top five**"))
    top = noc.assign(_ef=ef).nlargest(5, "_ef")
    cap = max(float(top["_ef"].max()), 1e-9)
    for i, (_, row_) in enumerate(top.iterrows(), 1):
        name_ = str(row_["compound"])[:34]
        val = float(row_["_ef"])
        cf = float(row_["confidence"]) if (conf is not None and pd.notna(row_.get("confidence"))) else None
        label_for = f"**{i}. {name_}** — {val:.0f} %"
        if cf is not None:
            label_for += f" · confidence {cf:.2f}"
        st.progress(max(0.0, min(val / cap, 1.0)), text=label_for)

    if conf is not None and fiables == 0:
        st.warning(t("No compound reaches confidence 0.5. With few poses the metric loses resolution: raise \"Poses per ligand\" in step 3 and run again."))

    cols_res = [c for c in ("compound", "best_dock", "pKi", "LE", "best_inter", "effectiveness_pct",
                            "percentile", "confidence", "cnn_score", "consensus") if c in noc.columns]
    _download_table(noc[cols_res].sort_values("effectiveness_pct", ascending=False),
                     "resumen_" + str(noc["receptor"].iloc[0]).replace("~", "_"), key="summary_view")


def render_results_viewer(proj: Path):
    """Viewer panel for Results stage: Summary, 3D complex, and Transport."""
    views = ["Summary", "3D complex"]
    if tn.available() and tn.runs_in(proj / lay.TUNNELS):
        views.append("Transport")
    view_ = st.radio(t("View"), views, horizontal=True,
                     format_func=t, key="vis_res_view", label_visibility="collapsed")
    if view_ == "3D complex":
        _complex_viewer(proj)
    elif view_ == "Transport":
        _transport_viewer(proj)
    else:
        _visual_summary(proj)

