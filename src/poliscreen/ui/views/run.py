"""Stage 3: Run view (box configuration, hybrid docking, parameters, and execution)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from poliscreen.core import adcp
from poliscreen.core import caver as cv
from poliscreen.core import docking as dk
from poliscreen.core import layout as lay
from poliscreen.core import peptides as pp
from poliscreen.core import pipeline as pl
from poliscreen.core import pockets as pk
from poliscreen.core import screening as sc
from poliscreen.core import viewer as vw
from poliscreen.ui.common import (
    _already_done,
    _download_table,
    _empty_state,
    _mark_done,
    _notify,
    _rname,
    _run_summary,
    _viewer_height,
)
from poliscreen.ui.components.transport import _tunnel_groups, _tunnel_table
from poliscreen.ui.i18n import t


def _batch_chemotypes():
    """(has_vina, has_peptides) present in the batch."""
    S = st.session_state
    peps = S.get("ligand_mode") == "Generate peptides" or bool(S.get("pep_seqs"))
    vina = False
    for c in S.get("controls", []):
        cl = str(c).lower()
        if cl.endswith(".pdb") and pp.sequence_from_structure(c):
            peps = True
        elif cl.endswith((".sdf", ".mol2", ".mol")):
            vina = True
    if S.get("products") or (S.get("ligand_mode") in ("Build by reaction", "Upload ready ligands")
                             and S.get("ligands")):
        vina = True
    if not peps and not vina:
        vina = True
    return vina, peps


def _docking_params():
    """Docking parameters form and dictionary."""
    has_vina, has_peptides = _batch_chemotypes()
    adcp_ok = adcp.available()
    uses_adcp = has_peptides and adcp_ok
    uses_vina = has_vina or (has_peptides and not adcp_ok)

    exhaust, energy_range, ph, cpu, workers = 24, 3.0, 7.4, 1, 0
    adcp_steps, adcp_reps = 250_000, 20

    with st.expander(t("Advanced docking settings")):
        st.caption(t("The defaults are fine for a first exploration. Raise them for a definitive screen."))
        if uses_vina and uses_adcp:
            st.info(t("**Mixed** run: small molecules dock with **Vina** and peptides with **ADCP**. Each engine's settings appear separately below."))
        elif uses_adcp:
            st.info(t("**Peptide** screening with **ADCP**. Only its settings are shown; Vina's do not apply."))
        elif has_peptides and not adcp_ok:
            st.warning(t("There are peptides and ADCP is not installed: they will dock with **Vina**, whose sampling does not cover that flexibility. Install it with scripts/get_adcp.sh."))

        seed = st.number_input(t("Seed"), value=42, step=1,
                               help=t("Fixes the randomness: the same seed gives the same result in both engines."))
        n_poses = st.slider(t("Poses per ligand"), 1, 20, 5,
                            help=t("Below 3 the confidence metric loses resolution."))

        if uses_vina:
            st.markdown(t("**Vina** — small molecules"))
            exhaust = st.slider(t("Exhaustiveness"), 8, 64, 8, 8,
                                help=t("Higher = finer and slower search. 8 is for exploring; "
                                       "raise it to 24 for a run whose ranking you will report, "
                                       "where a single unlucky seed should not decide the order."))
            energy_range = st.slider(t("Energy range (kcal/mol)"), 1.0, 8.0, 3.0, 0.5,
                                     help=t("Energy window relative to the best pose for reporting alternative modes."))
            ph = st.slider(t("Protonation pH"), 5.0, 9.0, 7.4, 0.1,
                           help=t("pH at which OpenBabel protonates before docking (physiological ≈ 7.4)."))
            cpu = st.number_input(t("Threads per docking"), 1, 16, 1,
                                  help=t("1 keeps the result reproducible. Raise it only if you do not mind."))
            workers = st.number_input(t("Dockings in parallel (0 = automatic)"), 0, 32, 0)
            if cpu > 1:
                st.warning(t("With more than one thread per docking, Vina stops being deterministic."))

        if uses_adcp:
            st.markdown(t("**ADCP** — peptides"))
            st.caption(t("It uses the machine's cores automatically and is reproducible with the seed; Vina's thread settings do not affect it."))
            adcp_steps = st.select_slider(
                t("Steps per replica"), [50_000, 100_000, 250_000, 500_000, 1_000_000, 2_500_000],
                value=250_000, format_func=lambda v: f"{v // 1000} k",
                help=t("Length of each search. Raise it if the control does not recover its pose or if the energy keeps improving as you increase it."))
            adcp_reps = st.slider(t("Independent replicas"), 4, 100, 20, 2,
                                  help=t("Parallel searches from different starting points. More replicas lower the chance of getting stuck in a local minimum."))

        st.markdown(t("**Second opinion (neural network)**"))
        has_gnina = dk.gnina_available()
        rescnn = st.checkbox(t("Re-score the poses with gnina (CNN, GPU)"), value=False,
                             disabled=not has_gnina,
                             help=t("It does not re-dock: it keeps the poses and evaluates them with a neural network trained on crystallographic complexes. Adds independent evidence to the confidence metric."))
        if not has_gnina:
            st.caption(t("gnina is not installed. It is optional: without it, confidence uses the other evidence."))
        elif rescnn and has_peptides:
            st.warning(t("You are screening peptides. gnina's network was trained on small-molecule complexes, so here it scores outside its domain: low values do not necessarily mean the pose is bad. Use it to compare, not as a criterion, and declare it in Methods."))
        elif rescnn:
            st.caption(t("The best pose of each compound is re-scored (~2 s per compound)."))
    return dict(seed=int(seed), exhaustiveness=int(exhaust), n_poses=int(n_poses),
                energy_range=float(energy_range), ph=float(ph), cpu=int(cpu), workers=int(workers),
                rescoring_cnn=bool(rescnn), adcp_steps=int(adcp_steps),
                adcp_replicas=int(adcp_reps))


def _engine_status():
    """What is installed, said once, with the one line that turns each on."""
    rows = [
        (t("CAVER (tunnels)"), cv.caver_available(), "POLISCREEN_CAVER",
         t("GPL-3 and cross-platform, but a Java program: the runtime is larger than the rest of "
           "PoliScreen. Point POLISCREEN_CAVER at caver.jar, with java on PATH.")),
        (t("CaverDock (transport)"), cv.caverdock_available(), "POLISCREEN_CAVERDOCK",
         t("A Linux Apptainer image under an academic licence, so it cannot be redistributed. "
           "Point POLISCREEN_CAVERDOCK at the .sif, with apptainer on PATH.")),
    ]
    for label, ok, _var, how in rows:
        if ok:
            st.success(f"{label} — {t('ready')}")
        else:
            st.info(f"**{label}** — {t('not installed')}. {how}")
    return all(ok for _l, ok, _v, _h in rows)


def _run_tunnels(proj: Path):
    """CAVER and CaverDock search tab in Run stage."""
    S = st.session_state
    st.subheader(t("Transport tunnels"))
    st.caption(t("A docking score says how well a compound sits in the site. This asks whether it "
                 "can reach it: CAVER finds the routes through the protein, CaverDock costs one."))

    has_caver = cv.caver_available()
    has_dock = cv.caverdock_available()
    if not (has_caver and has_dock):
        _engine_status()

    existing = S.get("tun_drawn")
    if existing and cv.clusters(existing):
        _tunnel_table(existing, cv.clusters(existing))

    recs = [Path(p) for p in S.get("receptors", []) if Path(p).exists()]
    if not recs:
        st.info(t("Stage a receptor in step 1 first."))
        return

    out_root = proj / lay.TUNNELS
    rec = st.selectbox(t("Receptor"), recs, format_func=_rname, key="tun_rec")

    st.markdown(t("**1 · What CAVER should look at**"))
    st.caption(t("Hydrogens make every atom effectively larger and close the narrow routes, so "
                 "they are removed. A receptor prepared for docking is not the right input here."))

    original = next((p for p in sorted(Path(rec).parent.glob(f"{Path(rec).stem.split('_ready')[0]}*.pdb"))
                     if "_ready" not in p.name and "control" not in p.name), None)
    sources = [p for p in (original, Path(rec)) if p is not None]
    source = st.selectbox(t("Structure"), sources, format_func=lambda p: p.name, key="tun_src",
                          help=t("The original download, not the docking-ready file."))

    present = cv.hetero_groups(source)
    keep = []
    if present:
        S.setdefault("tun_het", [])
        keep = st.multiselect(
            t("Heterogroups to keep"), [name for name, _n in present],
            format_func=lambda n: f"{n} ({dict(present)[n]} atoms)", key="tun_het",
            help=t("A cofactor sitting in a channel closes it. Keeping one says the route is "
                   "blocked in the physiological state; removing it says it is not."))
        st.caption(t("What is kept decides how many routes exist. Nothing is kept by default."))
    S.setdefault("tun_wat", False)
    waters = st.checkbox(t("Keep waters"), key="tun_wat")

    st.markdown(t("**2 · Where to measure from**"))
    ctrls = [Path(p) for p in S.get("controls", []) if Path(p).exists()]
    stored = (S.get("_boxes") or {}).get(str(rec))
    options = ([t("Centre of the search box")] if stored else []) \
        + ([t("Centre of a control")] if ctrls else []) + [t("Centre of chosen residues")]
    how = st.radio(t("Starting point"), options, key="tun_start", horizontal=False)

    start = None
    if how == t("Centre of the search box"):
        start = dk.Box(**stored)
    elif how == t("Centre of a control"):
        which = st.selectbox(t("Control"), ctrls, format_func=lambda p: p.name, key="tun_ctrl")
        start = cv.ligand_atoms(which)
    else:
        labels = sorted({f"{ln[17:20].strip()}{ln[22:26].strip()}"
                         for ln in Path(source).read_text(errors="ignore").splitlines()
                         if ln.startswith("ATOM")},
                        key=lambda r: (sc.resnum(r), r))
        chosen = st.multiselect(t("Residues"), labels, key="tun_res",
                                help=t("The catalytic or anchor residues, if you know them. The "
                                       "starting point is the middle of the ones you pick."))
        start = cv.atoms_of(source, chosen) if chosen else None

    if start is None:
        st.info(t("Choose the residues to measure from."))
        return
    try:
        x, y, z = cv.start_point(start)
    except cv.CaverError as e:
        st.error(str(e))
        return
    st.caption(t("Starting point: {x}, {y}, {z}.").format(x=x, y=y, z=z))

    if not cv.inside_structure((x, y, z), source):
        st.error(t("That point has no protein around it. CAVER would measure outwards into open "
                   "space and report it as one very wide tunnel. Check the structure and the "
                   "control you picked belong to the same chain."))
        return

    caver_out = out_root / f"caver_{Path(rec).stem}" / "out"
    found = cv.clusters(caver_out)
    if st.button(t("Find tunnels"), key="tun_find", disabled=not has_caver):
        with st.spinner(t("Running CAVER...")):
            try:
                prepared = cv.prepare_for_caver(
                    source, caver_out.parent / "prepared" / f"{Path(source).stem}.pdb",
                    keep_hetero=keep, keep_waters=waters)
                cv.find_tunnels(prepared, (x, y, z), caver_out.parent)
            except cv.CaverError as e:
                st.error(str(e))
                return
        found = cv.clusters(caver_out)
        S["tun_drawn"] = str(caver_out)
        st.success(t("{n} tunnels found.").format(n=len(found)))

    if not found:
        return
    S["tun_drawn"] = str(caver_out)

    st.markdown(t("**2 · Push a compound through one**"))
    if not has_dock:
        return
    st.caption(t("Route measured on `{caver}`. Docking against `{dock}`, which is the "
                 "docking-ready receptor: the energies need its hydrogens and charges.").format(
                     caver=Path(source).name, dock=Path(rec).name))

    ligs = [Path(p) for p in S.get("ligands", []) if Path(p).exists()]
    if not ligs:
        st.info(t("No compounds staged. Add them in step 2."))
        return

    c = st.columns(2)
    S.setdefault("tun_pick", found[:1])
    S.setdefault("tun_lig", ligs[:1])
    tunnels_pick = c[0].multiselect(t("Tunnels"), found,
                                    format_func=lambda p: p.stem, key="tun_pick")
    ligands_pick = c[1].multiselect(t("Compounds"), ligs,
                                    format_func=lambda p: p.stem, key="tun_lig")

    c2 = st.columns(5)
    for key_, default_ in (("tun_in", True), ("tun_out", True),
                           ("tun_lb", True), ("tun_ub", False)):
        S.setdefault(key_, default_)
    want_in = c2[0].checkbox(t("Entering"), key="tun_in")
    want_out = c2[1].checkbox(t("Leaving"), key="tun_out")
    want_lb = c2[2].checkbox(t("Lower bound"), key="tun_lb")
    want_ub = c2[3].checkbox(t("Upper bound"), key="tun_ub")
    cpus = c2[4].number_input(t("MPI processes"), min_value=2, max_value=16, step=1, key="tun_cpus")

    directions = [d for d, on in (("in", want_in), ("out", want_out)) if on]
    bounds = ["ub"] if want_ub else (["lb"] if want_lb else [])
    jobs = [(cmp_, tun, d, b) for cmp_ in ligands_pick for tun in tunnels_pick
            for d in directions for b in bounds]

    if not jobs:
        st.info(t("Pick at least one tunnel, one compound, one direction and one bound."))
        return
    if want_ub and want_lb:
        st.caption(t("Both bounds is one calculation: the upper bound produces the lower one on "
                     "its way."))
    minutes = sum(28 if b == "ub" else 8 for _c, _t, _d, b in jobs)
    st.caption(t("{n} calculations, roughly {m} minutes in total.").format(n=len(jobs), m=minutes))

    if not cv.reproducible(int(cpus)):
        st.warning(t("More than two processes: CaverDock warns that the seed no longer makes the "
                     "run repeatable. Faster, and not reproducible."))

    if st.button(t("Run transport"), key="tun_run", type="primary"):
        done, failed = [], []
        with st.status(t("Running CaverDock..."), expanded=True) as status:
            bar = st.progress(0.0)
            for i, (cmp_, tun, d, b) in enumerate(jobs):
                label = f"{Path(cmp_).stem} · {Path(tun).stem} · {d} · {b}"
                st.write(t("{n} of {total}: {what}").format(n=i + 1, total=len(jobs), what=label))
                try:
                    run = cv.transport(rec, cmp_, tun, out_root, direction=d, bound=b,
                                       cpus=int(cpus))
                    done.append(run)
                except cv.CaverError as e:
                    failed.append((label, str(e)))
                    st.write(f"  {t('failed')}: {str(e)[:160]}")
                bar.progress((i + 1) / len(jobs))
            state = "complete" if not failed else "error"
            status.update(label=t("{n} of {total} finished").format(n=len(done), total=len(jobs)),
                          state=state)
        if done:
            st.success(t("{n} written to {p}").format(n=len(done), p=out_root))
            _notify(t("Transport finished. Read it in Results, Transport tunnels."), str(out_root))
        for label, message in failed:
            st.error(f"{label}: {message}")

    st.caption(t("Every run goes into `{p}`. Point the Results tab at that folder and they come "
                 "out as one table, with the combinations not yet run counted as missing.").format(
                     p=out_root))


def _run_log_panel():
    """Log of the last run, preserved across tabs."""
    S = st.session_state
    state_ = S.get("_log_state", "completo")
    with st.status("Screening completed" if state_ == "completo" else "The run failed",
                   state="complete" if state_ == "completo" else "error", expanded=False):
        for n, d in S.get("_log_run", []):
            st.write(f"**{n}** · {d}")
        _download_table(pd.DataFrame(S["_log_run"], columns=["stage", "detail"]),
                         "registro_corrida", key="log_run")


def _run_screening(proj: Path):
    """Screening configuration and execution form in Run stage."""
    S = st.session_state
    st.subheader(t("Run the screening"))
    recs = [Path(p) for p in S["receptors"]]
    ctrls = [Path(p) for p in S["controls"]]
    ligs = [Path(p) for p in S["ligands"]]
    st.write(t("Receptors: **{r}** · Controls: **{c}** · ").format(r=len(recs), c=len(ctrls))
             + (t("Lead: `{lead}`").format(lead=S.get("lead")) if S.get("lead")
                else t("Ligands: **{n}**").format(n=len(ligs))))
    if recs and not ctrls:
        st.warning(t("No control loaded. The control is docked alongside the ligands and defines the reference; without it there is no baseline or validation. Extract the co-crystallized one in step 1 (or upload it below). If you already extracted it, check it is in the project's `receptores/` folder."))

    boxes = {}
    site_boxes = {}
    if recs:
        st.markdown(t("**Search box** — where to search inside the protein."))
        st.caption(t("Most reliable centered on the co-crystallized ligand: it marks the real site. The geometric center or a cofactor point elsewhere."))
        manual = {}
        if len(recs) > 1 and ctrls:
            auto = pl._assign_controls(ctrls, recs, {})
            unplaced = [c for c in ctrls if sc.normalize_key(c.stem) not in auto]
            if unplaced:
                with st.expander(t("Assign the missing controls to their receptor"), expanded=True):
                    st.caption(t("These controls could not be placed by geometry; indicate which receptor each belongs to."))
                    _rec_labels = ["(none)"] + [r.stem for r in recs]
                    for c in unplaced:
                        sel = st.selectbox(f"Control \"{c.stem}\"", _rec_labels,
                                           key=f"ctrlrec_{sc.normalize_key(c.stem)}")
                        if sel != "(none)":
                            manual[sc.normalize_key(c.stem)] = sel
            else:
                _pares = ", ".join(f"{c.stem} → {auto[sc.normalize_key(c.stem)]}" for c in ctrls)
                st.caption(t('Controls assigned automatically by geometry: {v1}.').format(v1=_pares))
        assignment = pl._assign_controls(ctrls, recs, manual)
        S["_control_map"] = manual
        S.setdefault("pockets", {})
        for _i_r, r in enumerate(recs):
            st.divider()
            _th, _tv = st.columns([3, 1], vertical_alignment="center")
            _th.markdown(f"### ▸ {_rname(r)}")
            if _tv.button(t("Show in the viewer"), key=f"verrec_{r.name}", width="stretch"):
                S["vis_box_rec"] = str(r)
                st.rerun()
            groups_ = dk.hetero_groups(r)
            ctrl = next((c for c in ctrls if assignment.get(sc.normalize_key(c.stem)) == r.stem), None)
            b1, b2 = st.columns([1, 3])
            ya_pk = bool(S["pockets"].get(str(r)))
            if b1.button(t("Detect pockets"), key=f"pk_{r.name}", type="primary",
                         disabled=not pk.fpocket_available() or ya_pk,
                         help="Cavities already detected for this receptor." if ya_pk else None):
                S["vis_box_rec"] = str(r)
                _why = []
                with st.spinner(t("Searching cavities with fpocket...")):
                    S["pockets"][str(r)] = pk.detect(r, on_notice=_why.append)
                S.setdefault("pockets_why", {})[str(r)] = _why[0] if _why else None
                S["vis_show_cav"] = True
                st.rerun()
            pkts = S["pockets"].get(str(r), [])
            if not pkts and not pk.fpocket_available():
                b2.caption(t("fpocket not installed: `conda install -n cribado -c conda-forge fpocket`."))
            elif not pkts and S.get("pockets_why", {}).get(str(r)):
                st.warning(S["pockets_why"][str(r)])
            pk_opts = {p["label"]: p for p in pkts}
            opts = ([f"Center on the control ({ctrl.name})"] if ctrl else []) \
                + list(pk_opts.keys()) + ["Automatic"] + list(groups_.keys())
            pick = st.selectbox(t("Box source"), opts, key=f"box_{r.name}",
                                on_change=lambda rr=str(r): S.__setitem__("vis_box_rec", rr))

            if pkts:
                cats_r = {x.lower() for x in S.get(f"cat_{r.stem}", [])}
                chosen = pk_opts.get(pick, {}).get("n")
                hib_labels = set(S.get(f"sites_{r.name}") or [])
                used_set = {chosen} | {p["n"] for p in pkts if p["label"] in hib_labels}
                used_set.discard(None)
                rows_, cav = [], []
                for i, p in enumerate(pkts[:8]):
                    in_use = p["n"] in used_set
                    color = (vw.CHOSEN_COLOR if p["n"] == chosen
                             else vw.CAVITY_PALETTE[i % len(vw.CAVITY_PALETTE)])
                    cav.append({"alpha": p.get("alpha_xyz"), "color": color, "chosen": in_use})
                    pr = p.get("props", {})
                    resid = p.get("residues") or []
                    has_catalytic = bool(cats_r and {x.lower() for x in resid} & cats_r)
                    row_ = {
                        "Color": vw.emoji_for_color(color), "Pocket": p["n"],
                        "Used": ("main" if p["n"] == chosen else "hybrid") if in_use else "",
                        "Druggability": p["druggability"], "Score": p["score"],
                        "Volume (Å³)": round(p.get("volume") or 0),
                        "Cavity (Å)": "%.0f×%.0f×%.0f" % (p.get("ex", 0), p.get("ey", 0), p.get("ez", 0)),
                        "Box (Å)": ("%.0f×%.0f×%.0f" % (p["sx"], p["sy"], p["sz"])
                                     + (" *" if p.get("minimo_aplicado") else "")),
                        "α-spheres": p.get("spheres"),
                        "Hydrophobicity": pr.get("Hydrophobicity score"),
                        "Polarity": pr.get("Polarity score"),
                        "Charge": pr.get("Charge score"),
                        "Apolar SASA": pr.get("Apolar SASA"),
                        "Residues": ", ".join(resid[:14]) + ("…" if len(resid) > 14 else ""),
                    }
                    if cats_r:
                        row_["Catalytic"] = "yes" if has_catalytic else "no"
                    rows_.append(row_)
                S.setdefault("_cavities", {})[str(r)] = cav
                dfp = pd.DataFrame(rows_)
                st.dataframe(
                    dfp, width="stretch", hide_index=True, height=230,
                    column_config={"Color": st.column_config.TextColumn(
                        "Color", width="small",
                        help=t("Color the cavity is drawn with in the viewer"))})
                _download_table(pd.DataFrame([{"Pocket": p["n"], **p.get("props", {}),
                                                "Residues": ", ".join(p.get("residues") or [])}
                                               for p in pkts]),
                                 f"cavidades_{r.stem}", key=f"cav_{r.name}")
                st.caption(t("All cavities are drawn at once in the right panel. The one **used for docking** is highlighted and more opaque. `Cavity` is its real extent; `Box` is the search region, with a 14 Å minimum because below that a ligand would not fit (marked `*` when that minimum was applied)."))
                with st.expander(t("All properties fpocket computes")):
                    st.dataframe(pd.DataFrame([{"pocket": p["n"], **p.get("props", {})} for p in pkts]),
                                 width="stretch", hide_index=True)
            else:
                S.get("_cavities", {}).pop(str(r), None)

            if ctrl and pick.startswith("Center on the control"):
                base = dk.box_from_file(ctrl)
            elif pick in pk_opts:
                base = pk.pocket_box(pk_opts[pick])
            elif pick in groups_:
                base = groups_[pick]
            else:
                base = dk.auto_box(r)
            if S.get(f"src_{r.name}") != pick or f"cx_{r.name}" not in S:
                S[f"src_{r.name}"] = pick
                for k, v in (("cx", base.cx), ("cy", base.cy), ("cz", base.cz),
                             ("sx", base.sx), ("sy", base.sy), ("sz", base.sz)):
                    S[f"{k}_{r.name}"] = float(v)
            gc, gs = st.columns(2)
            with gc:
                st.markdown(t("**Center** — where the box sits (Å)"))
                st.caption(t("Moves the box through space. The axes are shown in the viewer on the right."))
                cc = st.columns(3)
                cx = cc[0].number_input(t("← X →"), step=1.0, key=f"cx_{r.name}", format="%.1f",
                                        help=t("Left / right (red axis)."))
                cy = cc[1].number_input(t("↓ Y ↑"), step=1.0, key=f"cy_{r.name}", format="%.1f",
                                        help=t("Down / up (green axis)."))
                cz = cc[2].number_input(t("⊙ Z ⊗"), step=1.0, key=f"cz_{r.name}", format="%.1f",
                                        help=t("Into / out of the screen (blue axis)."))
            with gs:
                st.markdown(t("**Size** — how much the box spans (Å)"))
                st.caption(t("Grows or shrinks each side. If the ligand does not fit, Vina fails."))
                cs = st.columns(3)
                sx = cs[0].number_input(t("width X"), min_value=6.0, step=1.0, key=f"sx_{r.name}", format="%.1f")
                sy = cs[1].number_input(t("height Y"), min_value=6.0, step=1.0, key=f"sy_{r.name}", format="%.1f")
                sz = cs[2].number_input(t("depth Z"), min_value=6.0, step=1.0, key=f"sz_{r.name}", format="%.1f")
            boxes[str(r)] = dk.Box(cx, cy, cz, sx, sy, sz)
            S.setdefault("_boxes", {})[str(r)] = boxes[str(r)].as_dict()
            st.caption(t("The box is drawn over the receptor in the right panel."))
            _lig_este = list(S["ligands"]) + ([str(ctrl)] if ctrl else [])
            _minimo = dk.min_box(_lig_este)
            if _minimo and min(sx, sy, sz) < _minimo:
                st.warning(
                    t('The largest ligand is **{v1:.0f} Å** on its major axis and the box is {v3:.0f} Å on its shortest side. It fits, but cannot reorient: the search is restricted to the orientations that fit. Raise all three sides to at least **{v5:.0f} Å**.').format(v1=_minimo - 4, v3=min(sx, sy, sz), v5=_minimo))

            if pkts:
                hybrid_available = [p["label"] for p in pkts if p["label"] != pick]
                extra = st.multiselect(
                    t("Also dock in other pockets (hybrid docking)"),
                    hybrid_available, key=f"sites_{r.name}",
                    help=t("Each chosen pocket is docked separately and gets its own ranking. Reveals whether a compound prefers the catalytic site or slips into an allosteric one."))
                if extra:
                    lst = [("principal", boxes[str(r)])]
                    for lab in extra:
                        pdd = next((p for p in pkts if p["label"] == lab), None)
                        if pdd:
                            lst.append((f"Pk{pdd['n']}", pk.pocket_box(pdd)))
                    site_boxes[str(r)] = lst
                    st.caption(t('Hybrid docking: {v1} sites (main + {v3} pocket(s)).').format(v1=len(lst), v3=len(extra)))
                S[f"_hib_sel_{r.name}"] = set(extra)

    params = _docking_params()
    c1, c2 = st.columns([2, 1])
    reuse = c1.checkbox(t("Reuse previous calculations from this folder"), value=False,
                        help=t("Off, each run recomputes everything. Enable only if nothing has changed: reusing poses made with another box gives false results."))
    from poliscreen.ui.common import _confirm_delete
    if c2.button(t("Delete this folder's results")):
        _confirm_delete(proj)
    st.caption(t("Everything is saved in `{proj}` — poses, complexes, PLIP XML and the CSV tables.").format(proj=proj))

    if not recs:
        st.info(t("Prepare at least one receptor in step 1."))
    elif not (S.get("lead") or ligs):
        st.info(t("Choose compounds in step 2."))
    else:
        firma = (tuple(sorted(str(x) for x in recs)), tuple(sorted(str(x) for x in ligs)),
                 tuple(sorted(str(x) for x in ctrls)),
                 tuple(sorted((k, tuple(v.as_dict().values())) for k, v in boxes.items())),
                 tuple(sorted((k, len(v)) for k, v in site_boxes.items())),
                 tuple(sorted(params.items())), reuse, str(proj))
        done_ = _already_done("run", firma)
        if st.button(t("Run"), type="primary", disabled=done_,
                     help="Already run with this configuration. Change something to launch again."
                          if done_ else None):
            cfg = pl.RunConfig(receptors=recs, out_dir=proj, lead=S.get("lead") or None, ligands=ligs,
                               controls=ctrls, boxes=boxes, site_boxes=site_boxes,
                               control_map=S.get("_control_map") or {},
                               n_analogs=int(S.get("n_analogs", 20)),
                               n_substitutions=S.get("n_sub", [1]) or [1], use_ml=bool(S.get("use_ml", True)),
                               reuse=reuse, **params)
            S["_log_run"] = []

            with st.status(t("Running..."), expanded=True) as status:
                _bar = st.empty()

                def _paso(n, d):
                    if n == "docking-progress":
                        done_n, total_n = (int(x) for x in d.split("/"))
                        _bar.progress(done_n / max(total_n, 1),
                                      text=t("Docking {done} of {total}").format(done=done_n, total=total_n))
                        return
                    S["_log_run"].append((n, d))
                    st.write(f"**{n}** · {d}")

                try:
                    pl.run(cfg, on_step=_paso)
                    _mark_done("run", firma)
                    status.update(label=t("Screening completed"), state="complete")
                    S["_log_state"] = "completo"
                    _notify(t("Screening completed. Go to the results tab."), _run_summary(S.get("_log_run")))
                    st.success(t("Done. Go to the results tab."))
                except Exception as e:
                    status.update(label=t("Failed"), state="error")
                    S["_log_state"] = "error"
                    S["_log_run"].append(("error", str(e)))
                    _notify(t("The screening failed. Check the message."), str(e)[:200])
                    st.error(str(e))
        elif S.get("_log_run"):
            _run_log_panel()
        if done_:
            st.caption(t("Screening completed with this configuration. Change a parameter to enable the button again."))


def render_run_tools(proj: Path):
    """Tool panel for Run stage."""
    tab_screening, tab_tunnels = st.tabs([t("Screening"), t("Transport tunnels")])
    with tab_tunnels:
        _run_tunnels(proj)
    with tab_screening:
        _run_screening(proj)


def render_run_viewer(proj: Path):
    """Viewer panel for Run stage: 3D search box, cavities, and tunnels."""
    S = st.session_state
    boxes_ = S.get("_boxes") or {}
    cav_map = S.get("_cavities") or {}
    if boxes_:
        c1, c2, c3, c4 = st.columns([2, 1, 1, 1])
        rsel = c1.selectbox(t("Receptor"), list(boxes_), format_func=_rname,
                            key="vis_box_rec", label_visibility="collapsed")
        groups_by_receptor = cav_map.get(rsel)
        S.setdefault("vis_show_cav", bool(groups_by_receptor))
        ver_cav = c2.checkbox(t("Cavities"), key="vis_show_cav")
        routes = _tunnel_groups(set(S.get("tun_shown") or []))
        any_found = bool(S.get("tun_drawn") and cv.clusters(S["tun_drawn"]))
        S.setdefault("vis_show_tun", True)
        ver_tun = c3.checkbox(t("Tunnels"), key="vis_show_tun", disabled=not any_found,
                              help=None if any_found else t("Find them in Run first."))
        S.setdefault("vis_axes_box", True)
        axes_ = c4.checkbox(t("XYZ axes"), key="vis_axes_box")
        groups_ = list(groups_by_receptor) if (ver_cav and groups_by_receptor) else []
        drawn = routes if (ver_tun and routes) else []
        try:
            _h = _viewer_height(150)
            st.iframe(vw.view_html(receptor=rsel, box_=boxes_[rsel],
                                   cavities=(groups_ + drawn) or None,
                                   show_waters=False, axes_=axes_, height_=_h), height=_h + 12)
            b = boxes_[rsel]
            st.caption(f"Box (mauve): center ({b['cx']}, {b['cy']}, {b['cz']}) · "
                       f"{b['sx']} × {b['sy']} × {b['sz']} Å"
                       + (f" · {len(groups_)} cavities; the one used is highlighted." if groups_ else "")
                       + (t(" · {n} tunnels").format(n=len(drawn)) if drawn else ""))
        except Exception as e:
            st.error(t('Could not draw: {v1}').format(v1=e))
    else:
        _empty_state("Define the search box and it will be drawn here over the receptor.")

