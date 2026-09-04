"""Stage 2: Ligands view (reaction builder, approved drugs, peptides, upload, and viewer)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from poliscreen.core import adcp
from poliscreen.core import drugs as dg
from poliscreen.core import layout as lay
from poliscreen.core import ligands as lig
from poliscreen.core import naming as nm
from poliscreen.core import peptides as pp
from poliscreen.core import reactions as rx
from poliscreen.core import reagents as rg
from poliscreen.core import screening as sc
from poliscreen.core import viewer as vw
from poliscreen.core.design import AdmelabBridge
from poliscreen.ui.common import (
    _already_done,
    _download_table,
    _empty_state,
    _mark_done,
    _notify,
    _to_smiles,
)
from poliscreen.ui.components.admet import _render_adme, _shade
from poliscreen.ui.i18n import t


def _merge_ligand_meta(proj: Path, rows):
    """Adds these rows to the project's ligand table, keeping what other sources already wrote."""
    path = proj / "ligands_meta.csv"
    previous = []
    if path.exists():
        try:
            previous = pd.read_csv(path).to_dict("records")
        except Exception:
            previous = []
    fresh = {r["name"] for r in rows}
    merged = [r for r in previous if r.get("name") not in fresh] + list(rows)
    path.write_text(pd.DataFrame(merged).to_csv(index=False))


def _drug_mode(proj: Path):
    """Approved drugs as a ligand source, filtered by property."""
    S = st.session_state
    st.caption(t("Compounds already approved as medicines, from ChEMBL. Nothing is designed here: "
                 "they exist, so there is no synthesis feasibility to judge — only whether they "
                 "fit the properties you want."))
    cache = proj / "chembl_approved.csv"
    library = dg.read_csv(cache)

    top = st.columns([2, 1])
    wanted = top[0].number_input(t("How many to bring from the library"), min_value=100,
                                 max_value=5000, value=2000, step=100,
                                 help=t("ChEMBL holds about 4200 approved compounds. They are "
                                        "saved inside this project, so the run records exactly "
                                        "which library it used and the next one is instant."))
    if top[1].button(t("Download library"), disabled=bool(library)):
        with st.spinner(t("Asking ChEMBL...")):
            library, notice = dg.fetch_approved(cache=cache, max_records=int(wanted))
        if notice:
            st.warning(notice)
    if not library:
        st.info(t("Download the library to start. It needs the internet once; after that this "
                  "project works offline."))
        return
    st.caption(t("{v0} compounds in the library of this project. Delete "
                 "`chembl_approved.csv` from the project folder to refresh it.").format(
                     v0=len(library)))

    st.markdown(t("##### Filters"))
    cols = st.columns(2)
    limits = {}
    if cols[0].checkbox(t("Lipinski (rule of five)"), value=True,
                        help=t("MW under 500, LogP under 5, at most 5 donors and 10 acceptors.")):
        limits.update(dg.LIPINSKI)
    if cols[1].checkbox(t("Veber (oral bioavailability)"), value=False,
                        help=t("At most 10 rotatable bonds and 140 A^2 of polar surface. Looks at "
                               "flexibility, which Lipinski does not.")):
        limits.update(dg.VEBER)

    with st.expander(t("Adjust the ranges by hand")):
        st.caption(t("Anything set here replaces the preset for that property. The values are the "
                     "same ones the ranking table reports later, computed the same way."))
        ranges = {"MW": (0.0, 1000.0, (0.0, 500.0)), "LogP": (-5.0, 10.0, (-5.0, 5.0)),
                  "TPSA": (0.0, 300.0, (0.0, 140.0)), "HBD": (0.0, 20.0, (0.0, 5.0)),
                  "HBA": (0.0, 30.0, (0.0, 10.0)), "RotB": (0.0, 30.0, (0.0, 10.0))}
        picked = st.multiselect(t("Properties to bound"), list(ranges), default=[],
                                format_func=lambda k: k)
        for prop in picked:
            lo, hi, default = ranges[prop]
            limits[prop] = st.slider(prop, min_value=lo, max_value=hi, value=default,
                                     key=f"drug_range_{prop}")

    with st.spinner(t("Applying the filters...")):
        kept = dg.apply_filters(library, limits)
    if not kept:
        st.warning(t("No compound passes these filters. Loosen one of them."))
        return
    st.success(t("{v0} of {v2} compounds pass.").format(v0=len(kept), v2=len(library)))

    shown = pd.DataFrame(kept)[["name", "chembl_id"] + list(dg.PROPERTIES)]
    st.dataframe(shown.head(200), width="stretch", hide_index=True)
    if len(shown) > 200:
        st.caption(t("Showing the first 200. All {v0} are docked if you continue.").format(
            v0=len(kept)))

    how_many = st.number_input(t("How many to dock"), min_value=1, max_value=len(kept),
                               value=min(25, len(kept)), step=1,
                               help=t("Taken from the top of the filtered table. Docking is the "
                                      "slow step: start small, widen once the box is right."))
    chosen = kept[:int(how_many)]
    signature = (tuple(c["chembl_id"] for c in chosen),)
    already = _already_done("use_drugs", signature)
    if already:
        st.caption(t("These drugs are already loaded for the screening. Change the selection to "
                     "regenerate them."))
    if st.button(t("Add these drugs to the screening"), type="primary", disabled=already):
        d = lay.artifact(proj, lay.INPUT_LIGANDS)
        names_ = [lig.safe_name(c.get("name") or c.get("chembl_id")) for c in chosen]
        with st.spinner(t('Generating 3D of {v1} compounds...').format(v1=len(chosen))):
            made = lig.materialize([c["smiles"] for c in chosen], d, names=names_)
        done_set = {nm_ for nm_, _, _ in made}
        S["ligands"] = [str(p) for p in sorted(d.iterdir()) if p.is_file()]
        rows_ = [{"name": nm_, "smiles": c.get("smiles"), "source": "chembl",
                  "product": c.get("name"), "iupac_name": None, "feasibility": None}
                 for nm_, c in zip(names_, chosen) if nm_ in done_set]
        _merge_ligand_meta(proj, rows_)
        _mark_done("use_drugs", signature)
        _notify(t('{v0} compounds built and ready for step 3.').format(v0=len(made)))
        st.success(t("{v0} drugs added. The screening now has {v2} compounds in total, from "
                     "every source you have used.").format(v0=len(made), v2=len(S["ligands"])))


def _peptide_mode(proj: Path):
    """Peptide design: a path independent of reaction synthesis."""
    S = st.session_state
    st.caption(t("Peptides undergo no chemical reactions: they are built directly from the sequence. Between 1 and 20 residues."))
    entrada = st.radio(t("How to obtain the sequences"), ["Generate library", "Write sequences"],
                       horizontal=True, format_func=t, key="pep_input")

    sequences, notice, problems = [], "", []
    if entrada == "Write sequences":
        txt = st.text_area(t("One sequence per line, in one-letter code"),
                           placeholder=t("KWKLFKKI\nGIGKFLHSAK\nRRWWRF"), height=130, key="pep_text")
        raw_rows = [s.strip().upper() for s in txt.splitlines() if s.strip()]
        bad_entries = []
        for s in raw_rows:
            outside = set(s) - set(pp.AMINO_ACIDS)
            if outside:
                bad_entries.append(f"{s} (invalid symbols: {', '.join(sorted(outside))})")
            elif not (pp.MIN_LENGTH <= len(s) <= pp.MAX_LENGTH):
                bad_entries.append(f"{s} (length {len(s)}; the maximum is {pp.MAX_LENGTH})")
            else:
                sequences.append(s)
        if bad_entries:
            st.warning("These lines are ignored: " + " · ".join(bad_entries[:6]))
    else:
        c1, c2, c3 = st.columns(3)
        S.setdefault("pep_len", 7)
        largo = c1.number_input(t("Residues per peptide"), pp.MIN_LENGTH, pp.MAX_LENGTH, key="pep_len")
        S.setdefault("pep_n", 50)
        how_many = c2.number_input(t("How many peptides"), 1, 2000, key="pep_n")
        S.setdefault("pep_seed", 42)
        seed_ = c3.number_input(t("Seed"), step=1, key="pep_seed",
                                  help=t("Same seed and same rules = same library."))
        with st.expander(t("Composition: which amino acids it may use"), expanded=True):
            classes_ = st.multiselect(t("Allowed classes (empty = all 20)"), list(pp.CLASSES),
                                    format_func=lambda k: pp.CLASSES[k], key="pep_classes")
            excl = st.multiselect(t("Exclude specific residues"), sorted(pp.AMINO_ACIDS),
                                  format_func=lambda a: f"{a} · {pp.AMINO_ACIDS[a][0]}", key="pep_exclude")
            alf = pp.alphabet(include=classes_, exclude_residues=excl)
            st.caption(t('Resulting alphabet ({v1}): {v3}').format(v1=len(alf), v3=', '.join(alf) if alf else 'empty'))
        with st.expander(t("Sequence rules")):
            r1, r2 = st.columns(2)
            no_rep = r1.checkbox(t("No repeated residues"), key="pep_no_repeat")
            S.setdefault("pep_consecutive", 0)
            maxcons = r1.number_input(t("Max identical in a row (0 = no limit)"), 0, 10, key="pep_consecutive")
            S.setdefault("pep_maxres", 0)
            maxres = r2.number_input(t("Max times per residue (0 = no limit)"), 0, 20, key="pep_maxres")
            pre = r2.text_input(t("Starts with"), key="pep_prefix", placeholder=t("e.g. KK")).upper()
            suf = r1.text_input(t("Ends with"), key="pep_suffix", placeholder=t("e.g. GG")).upper()
        with st.expander(t("Physicochemical filters")):
            st.caption(t("In antimicrobial peptides, positive net charge and moderate hydrophobicity are the traits most associated with activity."))
            f1, f2 = st.columns(2)
            use_quick = f1.checkbox(t("Filter by net charge"), key="pep_use_q")
            S.setdefault("pep_q", (2.0, 9.0))
            q_rng = f1.slider(t("Net charge at pH 7.4"), -10.0, 10.0, step=0.5,
                              key="pep_q", disabled=not use_quick)
            use_gnina = f2.checkbox(t("Filter by hydropathy (GRAVY)"), key="pep_use_g")
            S.setdefault("pep_g", (-1.0, 1.0))
            g_rng = f2.slider(t("GRAVY"), -4.5, 4.5, step=0.1, key="pep_g", disabled=not use_gnina)

        rules = pp.Rules(length_=int(largo), alphabet=alf, no_repeats=no_rep,
                           max_consecutive=int(maxcons), max_per_residue=int(maxres),
                           prefix_=pre, suffix_=suf,
                           charge_min=q_rng[0] if use_quick else None,
                           charge_max=q_rng[1] if use_quick else None,
                           gravy_min=g_rng[0] if use_gnina else None,
                           gravy_max=g_rng[1] if use_gnina else None)
        problems = rules.validate()
        for p in problems:
            st.error(p)

    st.markdown("---")
    st.markdown(t("**Terminus chemistry**"))
    e1, e2, e3 = st.columns(3)
    n_ac = e1.checkbox(t("Acetylate N-terminus"), key="pep_n_acetyl",
                       help=t("Protects against aminopeptidases."))
    c_am = e2.checkbox(t("Amidate C-terminus"), key="pep_cam",
                       help=t("Removes the terminal negative charge: +1 net charge, which usually increases antimicrobial activity."))
    ciclo = e3.checkbox(t("Cyclize head-to-tail"), key="pep_cyclic",
                        help=t("Rigidifies the peptide and greatly reduces degrees of freedom, which also makes docking more reliable."))

    if entrada == "Generate library" and not problems:
        st.caption(t('Available combinatorial space: ~{v1:.0f} sequences.').format(v1=rules.space()))
        firma = (int(largo), int(how_many), int(seed_), tuple(alf), no_rep, int(maxcons),
                 int(maxres), pre, suf, use_quick, q_rng, use_gnina, g_rng)
        done_ = S.get("_pep_signature") == firma and S.get("pep_seqs")
        if st.button(t("Generate library"), type="primary"):
            with st.spinner(t("Generating sequences...")):
                sequences, notice = pp.generate(rules, int(how_many), seed=int(seed_))
            S["pep_seqs"] = sequences
            S["pep_notice"] = notice
            S["_pep_signature"] = firma
        if done_:
            st.caption(t('Library generated with these parameters ({v1} sequences).').format(v1=len(S['pep_seqs'])))
        sequences = sequences or S.get("pep_seqs", [])
        notice = notice or S.get("pep_notice", "")

    if notice:
        st.warning(notice)
    if not sequences:
        return

    rows_ = [pp.properties(s, c_amida=c_am, n_acetil=n_ac, cyclic=ciclo) for s in sequences]
    df = pd.DataFrame(rows_)[["name", "sequence", "length", "net_charge", "gravy",
                              "hydrophobic_moment", "hydrophobic_fraction", "boman_index"]]
    st.dataframe(df, width="stretch", hide_index=True, height=260)
    _download_table(df, "peptides", key="pep_table")
    st.caption(t("`momento_hidrofobico` measures amphipathicity (hydrophobic vs. polar face when folded into a helix); `indice_boman` estimates the tendency to bind other proteins: above 2.5 kcal/mol is considered promiscuous."))

    level, msg = pp.docking_feasibility(int(df["length"].max()), n_peptides=len(df),
                                       has_adcp=adcp.available())
    (st.success if level == "good" else st.warning if level == "mid" else st.error)(
        f"**Docking of {len(df)} peptides of {int(df['length'].max())} residues:** {msg}")

    S["_pep_preview"] = [(f["name"], f["sequence"]) for f in rows_[:24]]
    S["_pep_chemistry"] = (bool(n_ac), bool(c_am), bool(ciclo))
    _firma_pep = (tuple(sequences), n_ac, c_am, ciclo)
    _pep_listo = _already_done("use_peptides", _firma_pep)
    if st.button(t("Use these peptides in the screening"), type="primary", disabled=_pep_listo,
                 help="Change the sequences or terminus chemistry to build them again."
                      if _pep_listo else None):
        with st.spinner(t('Building the 3D structure of {v1} peptides...').format(v1=len(sequences))):
            smiles, names_, failures = [], [], 0
            for f in rows_:
                smi = pp.to_smiles(f["sequence"], n_acetil=n_ac, c_amida=c_am, cyclic=ciclo)
                if smi:
                    smiles.append(smi); names_.append(f["name"])
                else:
                    failures += 1
            made = lig.materialize(smiles, lay.artifact(proj, lay.INPUT_LIGANDS), names=names_)
        done_set = {nm_ for nm_, _, _ in made}
        S["ligands"] = [str(p) for _, p, _ in made]
        meta = pd.DataFrame([{"name": f["name"], "smiles": smi, "source": "peptide",
                              "product": f["sequence"], "iupac_name": None,
                              "feasibility": f"{f['length']} residues · charge {f['net_charge']}"}
                             for f, smi in zip(rows_, smiles) if f["name"] in done_set])
        (proj / "ligands_meta.csv").write_text(meta.to_csv(index=False))
        _mark_done("use_peptides", _firma_pep)
        lost = [n for n in names_ if n not in done_set]
        _notify(t('{v0} peptides built and ready for step 3.').format(v0=len(made)))
        st.success(f"{len(made)} peptides ready for step 3."
                   + (f" {failures} could not be built." if failures else ""))
        if lost:
            st.warning("Could not generate the 3D structure of: " + ", ".join(lost)
                       + ". They are long, flexible chains; try cyclizing them to rigidify them.")
    if _pep_listo:
        st.caption(t('{v0} peptides built with these parameters.').format(v0=len(S['ligands'])))


def render_ligands_tools(proj: Path):
    """Tool panel for Ligands stage."""
    S = st.session_state
    st.subheader(t("What do you want to dock?"))
    modo = st.radio(t("Source of the compounds"),
                    ["Build by reaction", "Screen approved drugs", "Generate peptides",
                     "Upload ready ligands"],
                    horizontal=True, format_func=t, key="ligand_mode")
    st.caption(t("Sources add up. Compounds from one source stay when you switch to another, so a "
                 "run can hold five products from the builder, five approved drugs and anything "
                 "you uploaded — and the results table says which is which."))

    if modo == "Screen approved drugs":
        S["lead"] = None
        _drug_mode(proj)
    elif modo == "Generate peptides":
        _peptide_mode(proj)
    elif modo == "Build by reaction":
        S["lead"] = None
        izq, der = st.columns(2)
        with izq:
            st.markdown(t("#### Core (your starting molecule)"))
            nuc = st.text_input(t("Core SMILES"), value=S.get("core_smiles", ""),
                                placeholder=t("O=C(O)c1ccc2[n+]([O-])onc2c1"))
            fnuc = st.file_uploader(t("...or core file"), type=["sdf", "mol2", "mol"])
            if fnuc is not None:
                d = proj / "nucleo"
                d.mkdir(parents=True, exist_ok=True)
                fp = d / fnuc.name
                fp.write_bytes(fnuc.getvalue())
                nuc = _to_smiles(fp) or nuc
            S["core_smiles"] = nuc
            rxkey = st.selectbox(t("Reaction"), list(rx.REACTIONS), key="rx_reaction",
                                 format_func=lambda k: t(rx.get(k).name_))
            reaction = rx.get(rxkey)
            st.caption(reaction.description_)
            if nuc:
                aplica = any(r.key == rxkey for r in rx.applicable(nuc))
                if not aplica:
                    st.warning(t('The core has no {v1}; this reaction does not apply.').format(v1=t(reaction.lead_grupo)))
                else:
                    sites = rx.lead_sites(nuc, reaction)
                    st.success(t('The core can undergo {v1}: {v3} reactive site(s).').format(v1=t(reaction.name_), v3=len(sites)))
                    idx = 0
                    if len(sites) > 1:
                        idx = st.selectbox(t("Growth point"), range(len(sites)),
                                           format_func=lambda i: f"atoms {sites[i]['atomos']}")
                    hl = sites[idx]["atomos"] if sites else []
                    S["_core_png"] = vw.molecule_png_indexed(nuc, highlight=hl, size=420)
        with der:
            if reaction.kind == "coupling":
                st.markdown(t("#### Reagents that couple"))
                use_internal_ = st.checkbox(t("Internal library"), key="use_internal_library",
                                       help=t('{v0} curated reagents.').format(v0=len(rg.load_internal(reaction)) if nuc else 0))
                ups = st.file_uploader(t("Your reagents (csv/xlsx with columns name and smiles · sdf · mol2 · smi)"),
                                       type=["csv", "xlsx", "xls", "sdf", "mol2", "mol", "smi"],
                                       accept_multiple_files=True)
                with st.expander(t("What columns must my Excel/CSV have?")):
                    st.markdown(
                        t("Two columns: a **name** and a **SMILES**. Accepted headers are:\n- Name: `name`, `nombre`, `compound`, `compuesto`, `Alcohol origen`, `Nombre clave`\n- SMILES: `smiles`, `smile`, `SMILES alcohol`\n\nDeduplicated by structure (InChIKey) and filtered to those bearing the reaction group (for esterification, an alcohol/phenol OH; acids and amines are discarded)."))
                    st.dataframe(pd.DataFrame({"name": ["Benzyl", "Menthol", "Cyclohexanol"],
                                                "smiles": ["OCc1ccccc1", "CC(C)C1CCC(C)CC1O", "OC1CCCCC1"]}),
                                 width="stretch", hide_index=True)
                use_pubchem_ = st.checkbox(t("Supplement with PubChem (experimental, needs internet)"), value=False)
                pc_max = st.number_input(t("PubChem maximum"), 5, 100, 25) if use_pubchem_ else 25
                upaths = []
                if ups:
                    d = proj / "reactivos"; d.mkdir(parents=True, exist_ok=True)
                    for u in ups:
                        (d / u.name).write_bytes(u.getvalue()); upaths.append(str(d / u.name))
                if st.button(t("Gather reagents"), type="primary"):
                    with st.spinner(t("Gathering and deduplicating...")):
                        reags, info = rg.build(reaction, use_internal=use_internal_, user_paths=upaths,
                                               use_pubchem=use_pubchem_, pubchem_max=int(pc_max))
                    S["reagents"] = [(r.name, r.smiles, r.source, r.inchikey) for r in reags]
                    S["reag_info"] = info
                    if info.get("aviso_pubchem"):
                        st.warning(info["aviso_pubchem"])
                if S.get("reag_info"):
                    info = S["reag_info"]
                    st.write(f"**{info['total']} reagents** — " +
                             " · ".join(f"{k}: {v}" for k, v in info["por_fuente"].items()))
                    dfa = pd.DataFrame([{"name": n, "SMILES": s, "source": src} for n, s, src, ik in S["reagents"]])
                    st.dataframe(_shade(dfa, "source"), width="stretch", height=240)
                    st.caption(t("Highlighted = reagents you provided."))
            else:
                st.markdown(t("#### Substituents"))
                st.caption(t("Decoration uses small internal groups (F, Cl, CN, OMe...); you upload no reagents."))
                c1, c2 = st.columns(2)
                S["n_analogs"] = c1.number_input(t("How many analogues"), 1, 200, S.get("n_analogs", 20))
                S["n_sub"] = c2.multiselect(t("Number of substitutions"), [1, 2, 3], default=S.get("n_sub", [1]))
                S["use_ml"] = st.checkbox(t("Predict ADMET with AI (slower the first time)"), value=S.get("use_ml", True))
                b = AdmelabBridge()
                if not b.available():
                    st.error(t("Cannot find the design engine (admelab)."))
                elif nuc and st.button(t("Generate analogues")):
                    with st.spinner(t("Generating and predicting properties...")):
                        d = b.design(nuc, use_ml=bool(S["use_ml"]),
                                     n_substitutions=S.get("n_sub", [1]) or [1], max_rows=int(S["n_analogs"]))
                    S["products"] = [dict(producto=(r.get("name") or f"analogo{i + 1:03d}"), smiles=r["SMILES"],
                                          fuente="internal", synthesizable=True, viabilidad="decoration")
                                     for i, r in enumerate(d.rows) if r.get("SMILES")]
                    if d.n_generated < int(S["n_analogs"]):
                        st.warning(t('{v0} analogues generated: with {v2} substitution(s) the chemical space runs out there. Try 2.').format(v0=d.n_generated, v2=S.get('n_sub', [1])))

        if reaction.kind == "coupling" and nuc and S.get("reagents"):
            st.markdown("---")
            b = AdmelabBridge()
            if not b.available():
                st.error(t("Cannot find the reaction engine (admelab)."))
            else:
                firma_p = (nuc, rxkey, tuple(sorted(ik for _n, _s, _src, ik in S["reagents"])))
                hecho_p = _already_done("productos", firma_p) and S.get("products")
                if hecho_p:
                    st.caption(t('{v0} products already built with this core and these reagents. Change one to rebuild.').format(v0=len(S['products'])))
                if st.button(t("Build products"), type="primary", disabled=bool(hecho_p)):
                    alcs = [{"name": n, "smiles": s} for n, s, src, ik in S["reagents"]]
                    with st.spinner(t("Building the series...")):
                        prods = b.esterify(nuc, alcs, policy="preferred")
                    src_by_ik = {ik: src for n, s, src, ik in S["reagents"]}
                    for p in prods:
                        p["source"] = src_by_ik.get(rg.inchikey(p.get("alcohol_smiles", "") or ""), "?")
                        p["product"] = p.get("alcohol")
                    _mark_done("productos", firma_p)
                    S["products"] = prods

        prods = S.get("products")
        if prods:
            st.markdown("---")
            n_ok = sum(1 for p in prods if p.get("synthesizable"))
            st.info(t('{v0} of {v2} products are synthesizable by this reaction.').format(v0=n_ok, v2=len(prods)))
            if not any(p.get("alcohol_smiles") for p in prods):
                st.caption(t("The IUPAC name is only generated for coupling reactions, where it is composed from the two fragments and verified with OPSIN. In decoration, products are identified by their SMILES."))
            if any(p.get("alcohol_smiles") for p in prods):
                if not nm.available():
                    st.button(t("Name (IUPAC, verified with OPSIN)"), disabled=True)
                    st.caption(t("IUPAC naming is off in this build. Products are identified by "
                                 "their SMILES, which is what gets docked. See Help."))
                elif st.button(t("Name (IUPAC, verified with OPSIN)")):
                    with st.spinner(t("Naming and verifying by round-trip...")):
                        named = AdmelabBridge().name_esters(
                            [p["smiles"] for p in prods], [p.get("alcohol_smiles") or "" for p in prods],
                            acid_smiles=nuc, alcohol_names=[p.get("product") for p in prods], use_web=True)
                    by = {n["smiles"]: n for n in named}
                    _names = [by.get(p["smiles"], {}).get("iupac_name") for p in prods]
                    checked = nm.verify(_names, [p["smiles"] for p in prods])
                    for p, (_nm, _ok) in zip(prods, checked):
                        p["iupac_name"], p["iupac_verif"] = _nm, _ok
                    S["products"] = prods
                    nver = sum(1 for p in prods if p.get("iupac_verif"))
                    st.success(t('{v0} of {v2} with a verified IUPAC name. The rest keep their label (the alcohol name); they are niche and OPSIN does not always cover them offline.').format(v0=nver, v2=len(prods)))
            dfp = pd.DataFrame(prods)
            cols = [c for c in ("product", "iupac_name", "source", "oh_type", "feasibility", "synthesizable", "smiles")
                    if c in dfp.columns]
            st.dataframe(_shade(dfp[cols].rename(columns={"product": "product", "source": "source",
                                                          "feasibility": "feasibility",
                                                          "synthesizable": "synthesizable"}), "source"),
                         width="stretch", height=320)
            st.caption(t("Highlighted = products with YOUR reagents. `synthesizable`=False are infeasible by this reaction."))
            st.caption(t("The 2D structures of the products are shown in the right panel so you can check the bond and stereochemistry."))

            with st.expander(t("ADMET report (predicts ~40 endpoints with AI for ALL at once)")):
                if st.button(t("Predict ADMET")):
                    with st.spinner(t("Predicting with ADMET-AI for all (the model downloads the first time)...")):
                        pr = AdmelabBridge().predict([p["smiles"] for p in prods], use_ml=True)
                    S["admet"] = {rg.inchikey(r.get("SMILES")): r for r in pr.rows}
                if S.get("admet"):
                    _render_adme(S["admet"], [(p.get("product") or f"prod{i}", p["smiles"])
                                              for i, p in enumerate(prods)], keyp="lig")

            c_sel1, c_sel2 = st.columns(2)
            solo_ok = c_sel1.checkbox(t("Dock only the synthesizable ones"), value=True)
            include_core = c_sel2.checkbox(t("Add the bare core (reference)"), value=True,
                                          help=t("Docks the unesterified core as a baseline: reveals how much activity the scaffold contributes on its own, apart from the tail."))
            use_signature = (tuple(p.get("smiles") for p in prods), solo_ok, bool(include_core and nuc))
            used_ = _already_done("use_products", use_signature)
            if used_:
                st.caption(t('These products are already loaded for the screening ({v1} compounds). Change the selection to regenerate them.').format(v1=len(S['ligands'])))
            if st.button(t("Use these products in the screening"), type="primary", disabled=used_):
                chosen_ones = [p for p in prods if (p.get("synthesizable") or not solo_ok)]
                if include_core and nuc:
                    chosen_ones = [dict(producto="free_core", smiles=nuc, fuente="core",
                                     iupac_name=None, viabilidad="reference (not esterified)",
                                     synthesizable=True)] + chosen_ones
                names_ = [lig.safe_name(p.get("product") or f"prod{i}") for i, p in enumerate(chosen_ones)]
                with st.spinner(t('Generating 3D of {v1} compounds...').format(v1=len(chosen_ones))):
                    made = lig.materialize([p["smiles"] for p in chosen_ones], lay.artifact(proj, lay.INPUT_LIGANDS), names=names_)
                done_set = {nm_ for nm_, _, _ in made}
                S["ligands"] = [str(p) for _, p, _ in made]
                meta = pd.DataFrame([{"name": nm_, "smiles": p.get("smiles"), "source": p.get("source", "?"),
                                      "product": p.get("product"), "iupac_name": p.get("iupac_name"),
                                      "feasibility": p.get("feasibility")}
                                     for (nm_, p) in zip(names_, chosen_ones) if nm_ in done_set])
                (proj / "ligands_meta.csv").write_text(meta.to_csv(index=False))
                _mark_done("use_products", use_signature)
                extra = " (includes the bare core as reference)" if (include_core and nuc) else ""
                _notify(t('{v0} compounds built and ready for step 3.').format(v0=len(made)))
                st.success(t('{v0} compounds ready for step 3{v2}.').format(v0=len(made), v2=extra))

    else:
        S["lead"] = None
        ups = st.file_uploader(t("Upload ligands"), type=["mol2", "sdf", "mol", "smi"], accept_multiple_files=True)
        if ups:
            d = lay.artifact(proj, lay.INPUT_LIGANDS)
            d.mkdir(parents=True, exist_ok=True)
            for u in ups:
                (d / u.name).write_bytes(u.getvalue())
            S["ligands"] = [str(p) for p in sorted(d.iterdir())]
            smap = sc.build_smiles_map(str(d))
            rows_ = []
            for p in S["ligands"]:
                name_ = Path(p).stem
                smi = smap.get(sc.normalize_key(name_))
                rows_.append({"name": name_, "smiles": smi, "source": "subido",
                              "product": name_, "iupac_name": None, "feasibility": None})
            without_smiles = sum(1 for f in rows_ if not f["smiles"])
            (proj / "ligands_meta.csv").write_text(pd.DataFrame(rows_).to_csv(index=False))
            if without_smiles:
                st.warning(t('{v0} of {v2} ligands gave no readable structure; there will be no ADMET or descriptors for them.').format(v0=without_smiles, v2=len(rows_)))
        if S["ligands"]:
            st.write(t("**{n} ligands:** ").format(n=len(S["ligands"]))
                     + ", ".join(Path(p).name for p in S["ligands"][:8]))
            ml = proj / "ligands_meta.csv"
            if ml.exists():
                mdf = pd.read_csv(ml)
                con = mdf["smiles"].notna().sum() if "smiles" in mdf.columns else 0
                st.caption(t('Structure read from {v1} of {v3}: allows computing ADMET, ligand efficiency, SAscore and PAINS alerts.').format(v1=con, v3=len(mdf)))
                items = [(r["name"], r["smiles"]) for _, r in mdf.iterrows() if pd.notna(r.get("smiles"))]
                if items:
                    with st.expander(t("ADMET report of the uploaded ligands")):
                        if st.button(t("Predict ADMET"), key="adme_uploaded"):
                            with st.spinner(t("Predicting with ADMET-AI...")):
                                pr = AdmelabBridge().predict([s for _, s in items], use_ml=True)
                            S["admet"] = {**(S.get("admet") or {}),
                                          **{rg.inchikey(r.get("SMILES")): r for r in pr.rows}}
                        if S.get("admet"):
                            _render_adme(S["admet"], items, keyp="sub")

    cup = st.file_uploader(t("Controls (co-crystallized ligand)"), type=["mol2", "sdf", "mol"],
                           accept_multiple_files=True,
                           help=t("If you already extracted it in step 1, no need to upload anything."))
    if cup:
        d = lay.artifact(proj, lay.RECEPTORS)
        d.mkdir(parents=True, exist_ok=True)
        new_items = []
        for u in cup:
            p = d / u.name
            try:
                p.write_bytes(u.getvalue())
            except Exception as e:
                st.error(t('Could not save {v1}: {v3}').format(v1=u.name, v3=e))
                continue
            if str(p) not in S["controls"]:
                S["controls"].append(str(p))
                new_items.append(u.name)
        if new_items:
            st.success(t('Control(s) loaded: {v1}. Total controls: {v3}.').format(v1=', '.join(new_items), v3=len(S['controls'])))
        st.caption(t("Loaded controls and those extracted in step 1 are docked alongside the ligands and define the reference fingerprint. With several receptors, assign each to its receptor in step 3."))


def render_ligands_viewer(proj: Path):
    """Viewer panel for Ligands stage: 2D grids of peptides, products, and uploaded ligands."""
    S = st.session_state
    prods = S.get("products")
    nuc_png = S.get("_core_png")
    peps = S.get("_pep_preview")
    if peps and S.get("ligand_mode") == "Generate peptides":
        st.markdown(t("**Generated sequences**"))
        legend = [("#3d7ea6", "hydrophobic"), ("#b5453c", "charge +"),
                   ("#3f7d4e", "charge -"), ("#7a6ba8", "polar"), ("#8a8a8a", "G/P")]
        def _color(a):
            cls = pp.AMINO_ACIDS[a][1]
            if "cargado_pos" in cls: return "#b5453c"
            if "cargado_neg" in cls: return "#3f7d4e"
            if "hidrofobico" in cls: return "#3d7ea6"
            if "especial" in cls:    return "#8a8a8a"
            return "#7a6ba8"
        html = []
        for nom, seq in peps:
            letters = "".join(
                f"<span style='display:inline-block;width:1.35em;text-align:center;"
                f"background:{_color(a)}22;color:{_color(a)};border-radius:3px;margin:1px;"
                f"font-weight:600'>{a}</span>" for a in seq)
            html.append(f"<div style='margin:.35rem 0'><code style='opacity:.6'>{nom}</code> "
                        f"<span style='font-family:monospace;font-size:1.05rem'>{letters}</span></div>")
        st.markdown("".join(html), unsafe_allow_html=True)
        st.markdown(" ".join(f"<span style='color:{c};font-size:.8rem'>■ {n}</span>"
                             for c, n in legend), unsafe_allow_html=True)
        _nac, _cam, _cic = S.get("_pep_chemistry", (False, False, False))
        _firma_grid = (tuple(s for _, s in peps), _nac, _cam, _cic)
        if S.get("_pep_grid_signature") != _firma_grid:
            smis, etiquetas = [], []
            for nom, seq in peps:
                s = pp.to_smiles(seq, n_acetil=_nac, c_amida=_cam, cyclic=_cic)
                if s:
                    smis.append(s); etiquetas.append(nom)
            S["_pep_grid"] = vw.grid_png(smis, legends=etiquetas, cols=3, sub=250) if smis else None
            S["_pep_grid_n"] = len(smis)
            S["_pep_grid_signature"] = _firma_grid
        if S.get("_pep_grid"):
            chemotypes = ("head-to-tail cyclized" if _cic else
                    ", ".join(filter(None, ["N-acetylated" if _nac else "",
                                            "C-amidated" if _cam else ""])) or "free termini")
            st.image(S["_pep_grid"],
                     caption=f"Structure of {S.get('_pep_grid_n', 0)} peptides · {chemotypes}.")
    elif prods:
        png = vw.grid_png([p.get("smiles") for p in prods],
                          legends=[str(p.get("product") or "") for p in prods])
        if png:
            st.image(png, caption=f"{len(prods)} products built. "
                                  "Check the ester bond and stereochemistry.")
    elif nuc_png:
        st.image(nuc_png, caption=t("Core with atom indices; in color, the reactive site."))
    elif S["ligands"]:
        st.success(t('{v0} ligands ready to dock.').format(v0=len(S['ligands'])))
        _ml = proj / "ligands_meta.csv"
        _sig = str(_ml.stat().st_mtime) if _ml.exists() else ""
        if _sig and S.get("_lig_grid_signature") != _sig:
            try:
                _m = sc.normalize_columns(pd.read_csv(_ml))
                _pairs = [(str(n), str(s)) for n, s in zip(_m.get("name", []), _m.get("smiles", []))
                          if isinstance(s, str) and s and s.lower() != "nan"]
            except Exception:
                _pairs = []
            S["_lig_grid"] = (vw.grid_png([s for _n, s in _pairs],
                                          legends=[_n for _n, s in _pairs]) if _pairs else None)
            S["_lig_grid_signature"] = _sig
        if S.get("_lig_grid"):
            st.image(S["_lig_grid"], caption=t("Structures of the ligands to be docked."))
        else:
            st.caption(", ".join(Path(p).stem for p in S["ligands"][:20]))
    else:
        _empty_state("Build or upload ligands and you will see their structures here.")

