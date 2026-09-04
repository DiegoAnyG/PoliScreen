"""Stage 1: Receptors view (preparation, controls, and 3D visualization)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from poliscreen.core import adcp
from poliscreen.core import layout as lay
from poliscreen.core import peptides as pp
from poliscreen.core import pipeline as pl
from poliscreen.core import receptor as rc
from poliscreen.core import screening as sc
from poliscreen.core import viewer as vw
from poliscreen.ui.common import (
    _already_done,
    _controls_of,
    _empty_state,
    _forget_all_receptors,
    _forget_receptor,
    _mark_done,
    _rname,
    _viewer_height,
)
from poliscreen.ui.i18n import t


def render_receptors_tools(proj: Path):
    """Tool panel for Receptors stage: PDB fetch/upload, inspect, clean, extract control."""
    S = st.session_state
    st.subheader(t("Prepare a receptor"))
    st.caption(t("Type a PDB identifier or upload your own file. Waters are removed, hydrogens added, and the original residue numbering is kept."))
    c1, c2 = st.columns([1, 2])
    pdb_id = c1.text_input(t("PDB identifier"), placeholder=t("4D44"))
    up = c2.file_uploader(t("...or upload a .pdb file"), type=["pdb"])

    src = None
    if up is not None:
        src = lay.artifact(proj, lay.RECEPTORS) / up.name
        src.parent.mkdir(parents=True, exist_ok=True)
        src.write_bytes(up.getvalue())
    elif pdb_id.strip() and st.button(t("Download from the PDB")):
        try:
            src = rc.fetch_pdb(pdb_id, lay.artifact(proj, lay.RECEPTORS))
            S["src_pdb"] = str(src)
        except rc.ReceptorError as e:
            st.error(str(e))
    if src is None and S.get("src_pdb"):
        src = Path(S["src_pdb"])

    if src and src.exists():
        st.success(t('Structure loaded: {v1}').format(v1=src.name))
        info = rc.inspect(src)
        st.write(t('**{v1}** atoms · chains **{v3}** · **{v5}** waters').format(v1=info.n_atoms, v3=', '.join(info.chains), v5=info.n_waters))
        if info.het:
            st.dataframe(pd.DataFrame([{"group": h.resname, "chain": h.chain, "number": h.resseq,
                                        "atoms": h.n_atoms, "key": h.key} for h in info.het]),
                         width="stretch", height=220)
        kb = sc.normalize_key(src.stem)
        c1, c2, c3 = st.columns(3)
        chains = c1.multiselect(t("Chains to keep"), info.chains, default=info.chains[:1],
                                key=f"rec_chains_{kb}")
        mods = rc.modified_residues(src)
        modified_keys = {(m.chain, m.resseq) for m in mods}
        keys = [h.key for h in info.het if (h.chain, str(h.resseq).strip()) not in modified_keys]
        keep = c2.multiselect(t("Keep (cofactors)"), keys, key=f"rec_keep_{kb}",
                              help=t("A site cofactor, e.g. NADP."))
        _cad_libres = [c for c in info.chains if c not in chains]
        _opc_ctrl = keys + [f"chain:{c}" for c in _cad_libres]

        def _fmt_ctrl(o):
            if o.startswith("chain:"):
                c = o.split(":", 1)[1]
                n = sum(1 for _l in src.read_text(errors="ignore").splitlines()
                        if _l.startswith("ATOM") and _l[21] == c and _l[12:16].strip() == "CA")
                return f"Chain {c} · {n} residues"
            return f"{o} (hetero group)"

        _sel_ctrl = c3.multiselect(t("Extract as control"), _opc_ctrl, key=f"rec_extract_{kb}",
                                   format_func=_fmt_ctrl,
                                   help=t("The co-crystallized ligand that defines the reference fingerprint. It can be a hetero group or a peptide chain; both appear here."))
        extract = [o for o in _sel_ctrl if not o.startswith("chain:")]
        control_chain = [o.split(":", 1)[1] for o in _sel_ctrl if o.startswith("chain:")]
        smiles = st.text_input(t("SMILES of the extracted ligand (optional)"), key=f"rec_smiles_{kb}",
                               help=t("Fixes bond orders, which the PDB does not store."))
        keep_mod = []
        if mods:
            st.markdown(t("**Modified residues of the chain**"))
            st.caption(t("Detected in the structure. Checked, they are kept with their modification; unchecked, they are replaced by the amino acid they derive from and **the modification is lost** — which is often the function, as in a phosphorylated activation loop."))
            keep_mod = st.multiselect(t("Keep with its modification"),
                                      [m.key for m in mods],
                                      default=[m.key for m in mods],
                                      key=f"rec_mod_{kb}",
                                      format_func=lambda k: next(
                                          (m.label for m in mods if m.key == k), k))
        for _c in control_chain:
            _n = sum(1 for _l in src.read_text(errors="ignore").splitlines()
                     if _l.startswith("ATOM") and _l[21] == _c and _l[12:16].strip() == "CA")
            if _n > pp.MAX_CHAIN_LENGTH:
                st.warning(t('Chain {v1} has {v3} residues: too long to treat it as a reference ligand.').format(v1=_c, v3=_n))
        firma_prep = (str(src), tuple(chains), tuple(keep), tuple(extract), tuple(control_chain),
                      tuple(keep_mod), smiles)
        prep_hecho = _already_done("prep_" + kb, firma_prep)
        if prep_hecho:
            st.caption(t("Receptor already prepared with this selection. Change something to prepare it again."))
        if st.button(t("Prepare receptor"), type="primary", disabled=prep_hecho):
            with st.spinner(t("Preparing...")):
                dest = lay.artifact(proj, lay.RECEPTORS) / f"{src.stem}{sc.READY_SUFFIX}.pdb"
                rc.prepare(src, dest, keep_chains=chains or None, keep_het=keep, ph=7.4,
                           keep_modified=keep_mod, on_notice=st.warning)
                if str(dest) not in S["receptors"]:
                    S["receptors"].append(str(dest))
                S["last_prepared"] = str(dest)
                S["last_original"] = str(src)
                for k in extract:
                    het = info.find(k)
                    p = rc.extract_ligand(src, het, lay.artifact(proj, lay.RECEPTORS) / f"control_{het.resname}.sdf",
                                          smiles=smiles or None, on_notice=st.info)
                    if str(p) not in S["controls"]:
                        S["controls"].append(str(p))
                for c in control_chain:
                    p = rc.extract_chain(src, c, lay.artifact(proj, lay.RECEPTORS) / f"control_Chain{c}.pdb",
                                         on_notice=st.warning)
                    if str(p) not in S["controls"]:
                        S["controls"].append(str(p))
                    _seq = pp.sequence_from_structure(p)
                    if _seq:
                        st.info(f"Chain {c} extracted as control: `{_seq[0]}` "
                                f"({len(_seq[0])} residues)"
                                + (", it will be docked with ADCP." if adcp.available()
                                   and adcp.MIN_RESIDUES <= len(_seq[0]) <= adcp.MAX_RESIDUES
                                   else "."))
            _mark_done("prep_" + kb, firma_prep)
            st.success(f"Done: {dest.name}")
            st.rerun()

    if S.get("last_prepared") and Path(S["last_prepared"]).exists():
        st.markdown("---")
        st.subheader(t("Preparation check"))
        antes = vw.structure_summary(S["last_original"])
        despues = vw.structure_summary(S["last_prepared"])
        comp = pd.DataFrame([
            {"": "Atoms", "before": str(antes["atomos"]), "after": str(despues["atomos"])},
            {"": "Hydrogens", "before": str(antes["hidrogenos"]), "after": str(despues["hidrogenos"])},
            {"": "Waters", "before": str(antes["aguas"]), "after": str(despues["aguas"])},
            {"": "Chains", "before": ", ".join(antes["chains"]), "after": ", ".join(despues["chains"])},
            {"": "Hetero groups", "before": ", ".join(sorted(antes["heterogrupos"])) or "-",
             "after": ", ".join(sorted(despues["heterogrupos"])) or "-"},
        ])
        st.dataframe(comp, width="stretch", hide_index=True)
        if despues["aguas"] == 0 and despues["hidrogenos"] > 0:
            st.success(t("No waters and hydrogens added."))
        else:
            st.warning(t("Check: there should be 0 waters and hydrogens present."))

        st.caption(t("The structure is shown in the right panel; there you can change the view and style."))

    if S["receptors"] or S["controls"]:
        st.markdown("---")
        _t1, _t2 = st.columns([3, 1], vertical_alignment="center")
        _t1.markdown(t("**Prepared receptors**"))
        if _t2.button(t("Remove all"), key="wipe_receptors", width="stretch",
                      help=t("Deletes every prepared receptor and control from the project folder.")):
            _forget_all_receptors()
            st.rerun()
        for _p in list(S["receptors"]):
            _c1, _c2 = st.columns([5, 1], vertical_alignment="center")
            _ctrl_p = _controls_of(Path(_p), S["receptors"], S["controls"])
            _c1.write(_rname(_p) + (f" · {', '.join(_rname(c) for c in _ctrl_p)}" if _ctrl_p else ""))
            if _c2.button("🗑", key=f"drop_receptor_{sc.normalize_key(_rname(_p))}", width="stretch",
                          help=t("Removes this receptor and the controls extracted from it.")):
                _forget_receptor(_p)
                st.rerun()
        _huerfanos = [c for c in S["controls"]
                      if not any(c in _controls_of(Path(r), S["receptors"], S["controls"])
                                 for r in S["receptors"])]
        if _huerfanos:
            st.caption(t("Controls with no receptor: {v1}").format(
                v1=", ".join(_rname(c) for c in _huerfanos)))


def render_receptors_viewer(proj: Path):
    """Viewer panel for Receptors stage: 3D structure with controls and cofactors."""
    S = st.session_state
    prepared = [p for p in S.get("receptors", []) if Path(p).exists()]
    if prepared:
        _ult = S.get("last_prepared")
        if len(prepared) > 1:
            _idx = prepared.index(_ult) if _ult in prepared else len(prepared) - 1
            S.setdefault("vis_rec_sel", prepared[_idx])
            rsel = st.selectbox(t("Receptor"), prepared,
                                format_func=_rname, key="vis_rec_sel")
        else:
            rsel = prepared[0]

        _stem = Path(rsel).stem
        _orig = next((str(o) for o in [Path(rsel).with_name(_stem[:-6] + ".pdb")
                                       if _stem.endswith(("_ready", "_listo")) else None,
                                       Path(S["last_original"]) if rsel == _ult
                                       and S.get("last_original") else None]
                      if o and o.exists()), None)
        _ctrls = [Path(c) for c in S["controls"]]
        _asig = pl._assign_controls(_ctrls, [Path(p) for p in prepared], S.get("_control_map") or {})
        _ctrl = next((str(c) for c in _ctrls
                      if _asig.get(sc.normalize_key(c.stem)) == _stem), None)

        c1, c2 = st.columns([3, 1])
        options = ["Prepared"] + (["Original"] if _orig else []) \
            + (["With its control"] if _ctrl else [])
        cual = c1.radio(t("View"), options, horizontal=True, format_func=t,
                        key="vis_show_rec", label_visibility="collapsed")
        S.setdefault("vis_axes_rec", True)
        axes_ = c2.checkbox(t("XYZ axes"), key="vis_axes_rec")
        try:
            receptor = _orig if (cual == "Original" and _orig) else rsel
            ligand_ = _ctrl if cual == "With its control" else None
            _h = _viewer_height(120)
            st.iframe(vw.view_html(receptor=receptor, ligand_=ligand_,
                                         show_waters=False, axes_=axes_, height_=_h), height=_h + 12)
            if cual == "With its control" and _ctrl:
                st.caption(t('Control of this receptor: `{v1}`.').format(v1=Path(_ctrl).stem))
        except Exception as e:
            st.error(t('Could not draw the structure: {v1}').format(v1=e))
    else:
        _empty_state("Prepare a receptor and it will appear here in 3D.")

