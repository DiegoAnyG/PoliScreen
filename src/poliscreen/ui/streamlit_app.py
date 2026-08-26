"""PoliScreen interface. Wraps the core; contains no science of its own.

Launch:  poliscreen ui      (or: streamlit run .../ui/streamlit_app.py)
"""
from __future__ import annotations

import json
import os
import re
import shutil
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import pandas as pd
import streamlit as st

from poliscreen import __version__
from poliscreen.core import adcp
from poliscreen.core import caver as cv
from poliscreen.core import docking as dk
from poliscreen.core import drugs as dg
from poliscreen.core import ligands as lig
from poliscreen.core import pipeline as pl
from poliscreen.core import peptides as pp
from poliscreen.core import pockets as pk
from poliscreen.core import reactions as rx
from poliscreen.core import reagents as rg
from poliscreen.core import report as rp
from poliscreen.core import receptor as rc
from poliscreen.core import layout as lay
from poliscreen.core import naming as nm
from poliscreen.core import screening as sc
from poliscreen.core import tunnels as tn
from poliscreen.core import session as ss
from poliscreen.core import validation as vl
from poliscreen.core import viewer as vw
from poliscreen.core.design import AdmelabBridge
from poliscreen.ui import ayuda
from poliscreen.ui.i18n import LANGUAGES, get_lang, t


def _shade(df, col, value="yours", color="rgba(255,205,60,0.20)"):
    """Highlight the rows whose column `col` equals `value` (to mark what the user contributed)."""
    if col not in df.columns:
        return df
    return df.style.apply(lambda r: [f"background-color: {color}" if str(r.get(col)) == value else ""
                                     for _ in r], axis=1)


def _scatter_dock_inter(sub):
    """Docking vs. interaction quality scatter. Shows the trade-off: top-right = good at both."""
    import matplotlib.pyplot as plt
    d = sub.copy()
    d["bd"] = pd.to_numeric(d.get("best_dock"), errors="coerce")
    d["bi"] = pd.to_numeric(d.get("best_inter"), errors="coerce")
    d = d.dropna(subset=["bd", "bi"])
    if d.empty:
        return None
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for _, r in d.iterrows():
        es_ctrl = r.get("is_control") == 1
        ax.scatter(r["bd"], r["bi"], s=110 if es_ctrl else 55,
                   c="#d62728" if es_ctrl else "#1b9e77", edgecolors="black", linewidths=0.6, zorder=3)
        ax.annotate(str(r["compound"])[:14], (r["bd"], r["bi"]), fontsize=7,
                    xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel(t("Docking (kcal/mol; more negative = better)"))
    ax.set_ylabel(t("Interaction quality (0-1 vs. control)"))
    ax.invert_xaxis()
    ax.grid(alpha=0.3)
    ax.set_title(t("Docking vs. quality · red = control · ideal: top-right"))
    fig.tight_layout()
    return fig


def _render_adme(admet, items, keyp):
    """items: [(label, smiles)]. Shows a summary table of all + detail per compound."""
    rows_ = []
    for lb, smi in items:
        r = admet.get(rg.inchikey(smi)) or {}
        rows_.append({"compound": lb, "MW": r.get("MW"), "LogP": r.get("LogP"), "QED": r.get("QED"),
                      "LD50 (mg/kg)": r.get("LD50_mg_per_kg"), "GHS": r.get("GHS_category"),
                      "AMES": r.get("AMES"), "hERG": r.get("hERG"), "DILI": r.get("DILI")})
    st.markdown(t("**ADMET summary of all compounds**"))
    # Everything ADMET-AI predicts comes back empty when it is not installed, and the one-click
    # installer does not ship it: the descriptors are computed from the structure, the endpoints
    # are not. An empty column says nothing about why it is empty.
    if not any(r.get(k) is not None for r in admet.values()
               for k in ("AMES", "hERG", "DILI", "LD50_mg_per_kg")):
        st.info(t("ADMET-AI is not installed on this machine: what you see are the properties "
                  "computed from the structure (MW, LogP, QED), not predicted endpoints. "
                  "docs/INSTALL.md explains how to add it."))
    st.dataframe(pd.DataFrame(rows_), width="stretch", height=min(320, 60 + 34 * len(rows_)))
    st.caption(t("AMES/hERG/DILI = toxicity probability (lower is better). LD50 in mg/kg (higher is better). Predicted on the WHOLE molecule (core + reagent), not the reagent alone."))
    labels = dict(items)
    sel = st.selectbox(t("View detail of"), list(labels), key=f"adme_det_{keyp}")
    row = admet.get(rg.inchikey(labels[sel]))
    if not row:
        return
    ca, cb = st.columns([1, 1])
    ca.pyplot(rp.radar_fig(row, title=sel))
    cb.metric(t("Oral LD50 (mg/kg)"), rp._f(row.get("LD50_mg_per_kg"), 0))
    cb.metric(t("GHS category"), str(row.get("GHS_category") or "-"))
    cb.metric(t("QED"), rp._f(row.get("QED")))
    cb.caption(t("Green = favorable · amber = intermediate · red = unfavorable."))
    cb.caption(t(rp.LD50_NOTICE))
    _col = {"good": "background-color:rgba(46,158,126,0.22)",
            "mid": "background-color:rgba(226,168,44,0.22)",
            "bad": "background-color:rgba(214,70,70,0.22)", "info": ""}
    for title_, fs in rp.sections(row):
        st.markdown(f"**{title_}**")
        dd = pd.DataFrame(fs, columns=["Property", "Value", "v"])
        sty = dd.style.apply(lambda r: [_col.get(r["v"], ""), _col.get(r["v"], ""), ""], axis=1)
        st.dataframe(sty, width="stretch", hide_index=True, column_config={"v": None})


def _to_smiles(path):
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    suf = Path(path).suffix.lower()
    m = (Chem.MolFromMol2File(str(path)) if suf == ".mol2" else
         next(iter(Chem.SDMolSupplier(str(path))), None) if suf == ".sdf" else
         Chem.MolFromMolFile(str(path)) if suf == ".mol" else None)
    return Chem.MolToSmiles(m) if m is not None else None


def _fmt_ki(ki):
    """Ki in readable units. Informative only: it does not take part in the ranking."""
    if ki is None or pd.isna(ki):
        return ""
    if ki < 1e-9:
        return f"{ki * 1e12:.1f} pM"
    if ki < 1e-6:
        return f"{ki * 1e9:.1f} nM"
    if ki < 1e-3:
        return f"{ki * 1e6:.2f} uM"
    return f"{ki * 1e3:.2f} mM"


def _rname(p) -> str:
    """Receptor display name (8HTB_ready.pdb -> 8HTB). The file on disk does not change; the
    prepared-file tag is internal and must never reach the screen."""
    return sc.display_name(p)


CITATIONS = f"""**PoliScreen v{__version__}** — Anaya Guerrero DC.

**Cite also the tools PoliScreen runs.** Failing to do so is a common omission in
peer review:

- **AutoDock Vina 1.2** — Eberhardt J, Santos-Martins D, Tillack AF, Forli S. *AutoDock Vina 1.2.0:
  New Docking Methods, Expanded Force Field, and Python Bindings.* J Chem Inf Model. 2021;61(8):3891–3898.
- **AutoDock Vina (original)** — Trott O, Olson AJ. J Comput Chem. 2010;31(2):455–461.
- **AutoDock CrankPep (ADCP)** (if you docked peptides) — Zhang Y, Sanner MF. *Docking Flexible
  Cyclic Peptides with AutoDock CrankPep.* J Chem Theory Comput. 2019;15(10):5161–5168.
- **AGFR / AutoGridFR** (target preparation for ADCP) — Zhang Y, Forli S, Omelchenko A,
  Sanner MF. *AutoGridFR.* J Comput Chem. 2019;40(32):2882–2891.
- **GNINA** (if you re-scored with the neural network) — McNutt AT, et al. *GNINA 1.0: molecular docking
  with deep learning.* J Cheminform. 2021;13(1):43.
- **PLIP** — Adasme MF, et al. *PLIP 2021: expanding the scope of the protein–ligand interaction
  profiler.* Nucleic Acids Res. 2021;49(W1):W530–W534.
- **Open Babel** — O'Boyle NM, et al. *Open Babel: An open chemical toolbox.* J Cheminform. 2011;3:33.
- **RDKit** — RDKit: Open-source cheminformatics. https://www.rdkit.org
- **fpocket** — Le Guilloux V, Schmidtke P, Tuffery P. *Fpocket: An open source platform for ligand
  pocket detection.* BMC Bioinformatics. 2009;10:168.
- **OpenMM / PDBFixer** — Eastman P, et al. *OpenMM 7.* PLoS Comput Biol. 2017;13(7):e1005659.
- **SAscore** — Ertl P, Schuffenhauer A. J Cheminform. 2009;1:8.
- **PAINS** — Baell JB, Holloway GA. J Med Chem. 2010;53(7):2719–2740.
- **ADMET-AI** (if you used ADMET prediction) — Swanson K, et al. Bioinformatics. 2024;40(7):btae416.

The exact version of each tool is exported with **File → Methods**, so the article's
Methods section is reproducible."""

ACKNOWLEDGMENTS = """PoliScreen builds on these open-source projects:

- **Streamlit** (Apache-2.0) — web interface
- **3Dmol.js / py3Dmol** (BSD-3-Clause) — 3D molecular viewer
- **pandas** (BSD-3-Clause) and **NumPy** (BSD-3-Clause) — data handling
- **Matplotlib** (Matplotlib/PSF license) — interaction diagrams
- **OpenPyXL** (MIT) — XLSX export
- **OPSIN** (MIT) — IUPAC name verification

The scientific tools the engine runs (Vina, ADCP, gnina, PLIP, RDKit, Open Babel,
fpocket, OpenMM) have their full citation in «How to cite»."""


def _how_to_cite(expanded_: bool = False):
    with st.expander(t("How to cite"), expanded=expanded_):
        st.markdown(CITATIONS)
    with st.expander(t("Acknowledgments")):
        st.markdown(t(ACKNOWLEDGMENTS))


def _download_table(df, name_: str, key: str):
    """Offers downloading a DataFrame as CSV or XLSX below the table."""
    import io
    c = st.columns([2, 1, 1])
    c[0].caption(t("Download as:"))
    c[1].download_button(t("CSV"), df.to_csv(index=False).encode("utf-8"),
                         file_name=f"{name_}.csv", mime="text/csv", key=f"dl_csv_{key}",
                         width="stretch")
    buf = io.BytesIO()
    try:
        with pd.ExcelWriter(buf, engine="openpyxl") as w:
            df.to_excel(w, index=False, sheet_name="PoliScreen")
        c[2].download_button(t("XLSX"), buf.getvalue(), file_name=f"{name_}.xlsx", key=f"dl_xlsx_{key}",
                             mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                             width="stretch")
    except Exception:
        pass


def _download_image(png: bytes, name_: str, key: str):
    """Offers downloading an image (PNG) below the viewer. jpeg adds nothing here: they are diagrams
    with a flat background, PNG keeps sharpness without artifacts."""
    if not png:
        return
    c = st.columns([3, 1])
    c[0].caption(t("Download as:"))
    c[1].download_button(t("PNG"), png, file_name=f"{name_}.png", mime="image/png",
                         key=f"dl_png_{key}", width="stretch")


def _already_done(key_: str, firma) -> bool:
    """True if the action already ran with THESE parameters. Used to disable the button until
    something changes: avoids repeating a long computation by mistake and double clicks."""
    return S.get("_signature_" + key_) == firma


def _mark_done(key_: str, firma):
    S["_signature_" + key_] = firma


def _controls_of(rec: Path, receptors: list, controls: list) -> list:
    """Controls belonging to a receptor, by the criterion the pipeline itself uses: geometry.
    The file name cannot say it — a control is named after its ligand (control_ZI9), never after
    the structure it came from."""
    assign = pl._assign_controls([Path(c) for c in controls],
                                 [Path(r) for r in receptors], {})
    return [c for c in controls
            if assign.get(sc.normalize_key(Path(c).stem)) == rec.stem]


def _forget_receptor(path_str: str) -> None:
    """Deletes a prepared receptor together with its controls and the choices made for it.

    The controls go with it: extracted from that structure, they share its coordinate system, and
    left behind they would be assigned by geometry to whichever receptor remains.
    """
    rec = Path(path_str)
    doomed = [path_str] + _controls_of(rec, S["receptors"], S["controls"])
    for f in doomed:
        Path(f).unlink(missing_ok=True)
    S["receptors"] = [p for p in S["receptors"] if p != path_str]
    S["controls"] = [p for p in S["controls"] if p not in doomed]
    S.setdefault("_forget_prep", []).append(sc.normalize_key(_rname(rec)))
    if S.get("last_prepared") == path_str:
        S.pop("last_prepared", None)
        S.pop("last_original", None)


def _forget_all_receptors() -> None:
    """Empties the project's receptor folder: prepared structures and controls alike. The way out
    when a control was left orphaned because its receptor was already gone."""
    for f in list(S["receptors"]) + list(S["controls"]):
        Path(f).unlink(missing_ok=True)
    S["_forget_prep"] = [sc.normalize_key(_rname(p)) for p in S["receptors"]]
    S["receptors"], S["controls"] = [], []
    S.pop("last_prepared", None)
    S.pop("last_original", None)


def _notify(headline: str, detail: str = "") -> None:
    """Report the end of a long operation. Shown as a modal, the same mechanism as the deletion
    prompt: a toast clips anything longer than a line, and custom HTML does not survive
    Streamlit's sanitizing, so neither was ever seen."""
    S["_notice"] = (headline, detail)


# Steps whose detail says how much was processed. The rest of the log (box source, per-ligand
# failures) stays in the run log, which is kept below the Run button.
_SUMMARY_STEPS = ("3d", "docking", "plip", "rescoring", "validation")


def _run_summary(log) -> str:
    """One line per phase, with what each one processed."""
    seen, out = set(), []
    for name_, detail in reversed(list(log or [])):
        if name_ in _SUMMARY_STEPS and name_ not in seen and detail:
            seen.add(name_)
            out.append(f"{name_}: {detail}")
    return " · ".join(reversed(out))


def _empty_state(message: str):
    """Empty state of the viewer: the logo as a watermark and a line with what is missing.
    If no logo is installed in assets/ a monochrome fallback glyph is used."""
    message = t(message)
    logo = vw.logo_path()
    if logo and logo.suffix.lower() != ".svg":
        c = st.columns([1, 2, 1])[1]
        c.image(str(logo), width="stretch")
        c.markdown(f"<div style='text-align:center;opacity:.75;font-size:.92rem'>{message}</div>",
                   unsafe_allow_html=True)
        return
    brand = logo.read_text(errors="ignore") if logo else vw.logo_svg()
    st.markdown(
        f"<div style='display:flex;flex-direction:column;align-items:center;justify-content:center;"
        f"padding:3rem 1rem;opacity:.75;text-align:center'>{brand}"
        f"<div style='margin-top:.9rem;font-size:.92rem;opacity:.8'>{message}</div></div>",
        unsafe_allow_html=True)


@st.dialog("PoliScreen")
def _notice_dialog(headline: str, detail: str):
    st.markdown(f"### {headline}")
    if detail:
        st.caption(detail)
    if st.button(t("Close"), type="primary", width="stretch"):
        st.rerun()


@st.dialog("Confirm deletion")
def _confirm_delete(folder):
    """Modal window before a destructive action: deleting results cannot be undone."""
    st.warning(t("This deletes the poses, complexes, PLIP XML and all result tables in this folder. **This action cannot be undone.**"))
    st.caption(f"Folder: `{folder}`")
    st.caption(t("Prepared receptors, controls and input ligands are kept."))
    st.info(t("To keep this analysis, cancel and use File -> Save session first."))
    c1, c2 = st.columns(2)
    if c1.button(t("Yes, delete"), type="primary", width="stretch"):
        pl.clean(folder)
        _notify(t("Results deleted. Receptors, controls and ligands are kept."))
        st.rerun()
    if c2.button(t("Cancel"), width="stretch"):
        st.rerun()


def _human_size(b: int) -> str:
    return f"{b / 1e6:.1f} MB" if b >= 1e6 else (f"{b / 1e3:.0f} kB" if b else "—")


@st.dialog("Download results", width="large")
def _downloads_dialog(folder):
    """Export selector. Builds a single ZIP in memory with whatever the user checks.

    The project folder already contains almost all of these files; the package exists to take the
    analysis to another machine, attach it to a paper or archive it without the heavy intermediates.
    That is why each item says whether it is worth downloading and how much it takes, instead of
    offering them blindly.
    """
    cat = ss.catalog(folder)
    has_ = {k: v for k, v in cat.items() if v["has"]}
    if not has_:
        st.info(t("Nothing to export in this folder yet."))
        return
    st.caption(t('Project folder: `{v1}`').format(v1=folder))
    st.caption(t("Everything here already lives in that folder, except the Methods section, written on export. The package is for moving the analysis elsewhere or attaching it to a manuscript without the heavy intermediates."))

    c1, c2 = st.columns(2)
    mark = c1.button(t("Select the recommended"), width="stretch")
    clean_ = c2.button(t("Clear all"), width="stretch")
    for k, v in has_.items(): 
        if mark or clean_:
            S[f"dl_{k}"] = bool(mark and v["reason"])
        else:
            S.setdefault(f"dl_{k}", bool(v["reason"]))

    rec = {k: v for k, v in has_.items() if v["reason"]}
    otros = {k: v for k, v in has_.items() if not v["reason"]}
    if rec:
        st.markdown(t("**Recommended**"))
        for k, v in rec.items():
            st.checkbox(f"{v['desc']} · {_human_size(v['bytes'])}", key=f"dl_{k}")
            st.caption(v["reason"])
    if otros:
        st.markdown(t("**Optional**"))
        for k, v in otros.items():
            st.checkbox(f"{v['desc']} · {_human_size(v['bytes'])}", key=f"dl_{k}")
            if v["regenerable"]:
                st.caption(t("Regenerated by running again; takes space in the package."))

    chosen_items = [k for k in has_ if S.get(f"dl_{k}")]
    total = sum(has_[k]["bytes"] for k in chosen_items)
    st.divider()
    if not chosen_items:
        st.caption(t("You have not selected anything."))
        return
    st.caption(t('{v0} item(s) · {v2} uncompressed').format(v0=len(chosen_items), v2=_human_size(total)))
    if st.button(t("Prepare package"), type="primary", width="stretch"):
        try:
            meta = folder / "run.json"
            mtxt = rp.methods_text(json.loads(meta.read_text())) if meta.exists() else None
            datos, included_items = ss.package_bytes(folder, chosen_items, methods_text=mtxt)
            S["_zip"] = (f"{Path(folder).name}_PoliScreen.zip", datos, included_items)
        except Exception as e:
            st.error(t('Could not build the package: {v1}').format(v1=e))
    if S.get("_zip"):
        name_, datos, included_items = S["_zip"]
        st.download_button(t('Download {v1} ({v3})').format(v1=name_, v3=_human_size(len(datos))), datos,
                           file_name=name_, mime="application/zip",
                           type="primary", width="stretch")
        st.caption("Contains: " + ", ".join(included_items))


st.set_page_config(page_title="PoliScreen", layout="wide", initial_sidebar_state="collapsed")

st.markdown(
    """
    <style>
      /* Streamlit's header is a FIXED bar covering the whole top width. Even when it looks
         empty, it intercepted the menu bar's clicks. It is made transparent to events and only
         the real buttons (the three-dot one) receive them again; the content is also pushed down
         just enough not to overlap it. */
      [data-testid="stAppDeployButton"] { display: none; }
      header[data-testid="stHeader"] { background: transparent; pointer-events: none; }
      header[data-testid="stHeader"] button { pointer-events: auto; }
      .block-container { padding-top: 2.4rem; padding-bottom: 0.3rem;
                         padding-left: 1.6rem; padding-right: 1.6rem; }
      /* Streamlit adornments that add nothing: the fullscreen button over each image and the
         anchor that appears next to titles and only links to the page itself. */
      [data-testid="StyledFullScreenButton"] { display: none; }
      [data-testid="stHeaderActionElements"] { display: none; }
      h1 a.anchor-link, h2 a.anchor-link, h3 a.anchor-link,
      h4 a.anchor-link, h5 a.anchor-link { display: none; }
      /* Dropdown lists open inside a scrolling panel, which puts them on their own compositing
         layer and left the text looking blurred and washed out. */
      div[data-baseweb="popover"], div[data-baseweb="popover"] * {
          opacity: 1 !important; filter: none !important; backdrop-filter: none !important; }
      div[data-baseweb="popover"] ul, div[data-baseweb="popover"] li {
          transform: none !important; text-rendering: geometricPrecision; }
      /* Header: faint separator line */
      .st-key-barra_menu { border-bottom: 1px solid rgba(128,128,128,.28); padding-bottom: .35rem; }
      /* Bottom stage bar: text only, flat and centered */
      .st-key-barra_etapas { border-top: 1px solid rgba(128,128,128,.28); padding-top: .35rem; }
      .st-key-barra_etapas button { border: none !important; background: transparent !important;
                                    font-weight: 600; letter-spacing: .01em; }
      .st-key-barra_etapas button:hover { background: rgba(128,128,128,.14) !important; }
      /* Active stage: accent color and underline. A faint gray background was indistinguishable. */
      .st-key-barra_etapas button[kind="primary"] {
          background: rgba(46,158,126,.16) !important;
          color: #1f7a63 !important;
          box-shadow: inset 0 -3px 0 0 #2E9E7E !important; }
      .st-key-barra_etapas button[kind="primary"]:hover {
          background: rgba(46,158,126,.24) !important; }
    </style>
    """,
    unsafe_allow_html=True,
)
S = st.session_state
S.setdefault("receptors", [])
S.setdefault("controls", [])
S.setdefault("ligands", [])

S["lang"] = LANGUAGES.get(S.get("_lang_pick", "English"), "en")

STAGES = ["Receptors", "Ligands", "Run", "Results"]
S.setdefault("stage", STAGES[0])

if S.get("_proj_pending"):
    S["proj_dir"] = S.pop("_proj_pending")
    S["_proj_loaded"] = None
for _k, _v in (S.pop("_pending_widgets", None) or {}).items():
    S[_k] = _v

# Removing a receptor also drops the choices made for it, but a widget key cannot be deleted in
# the run that drew the widget, so the removal is queued and settled here, before anything is drawn.
for _kb in S.pop("_forget_prep", []):
    for _k in (f"rec_chains_{_kb}", f"rec_keep_{_kb}", f"rec_extract_{_kb}",
               f"rec_smiles_{_kb}", f"rec_mod_{_kb}", f"_signature_prep_{_kb}"):
        S.pop(_k, None)

_PERSISTENT_PREFIXES = ("pep_", "modo_", "cat_", "sec_", "rec_", "box_", "sites_", "rx_",
                          "cx_", "cy_", "cz_", "sx_", "sy_", "sz_", "src_", "vis_", "cfg_")
for _k in [k for k in S.keys() if isinstance(k, str) and k.startswith(_PERSISTENT_PREFIXES)]:
    S[_k] = S[_k]

_barra = st.container(key="menu_bar")
_logo_f = vw.logo_path()
_wm = vw.wordmark_path()
_izq, _hueco, _m1, _m2, _m3, _m4 = _barra.columns([3.6, 1.4, 1.1, 1.1, 1.6, 1.1],
                                                  vertical_alignment="center")


def _img_inline(path_, height_, cls_name="", estilo=""):
    import base64
    mime = "svg+xml" if Path(path_).suffix.lower() == ".svg" else Path(path_).suffix.lower().lstrip(".")
    b64 = base64.b64encode(Path(path_).read_bytes()).decode()
    return f"<img class='{cls_name}' src='data:image/{mime};base64,{b64}' style='height:{height_}px;{estilo}'>"


_INV = "filter:invert(1) hue-rotate(180deg) brightness(1.1);"


def _png_size(path_):
    """Pixel size straight from the PNG header, so the mask below can be given a width."""
    b = Path(path_).read_bytes()[16:24]
    return int.from_bytes(b[:4], "big"), int.from_bytes(b[4:], "big")


def _img_tinted(path_, height_):
    """A single-colour brand image painted with the theme's own text colour.

    The colour used to be decided in Python from `st.context.theme` and written into the tag.
    That is stale the moment the theme changes: switching it does not re-run the script, so the
    wordmark kept the previous theme's colour until something else forced a rerun -- changing the
    language was what usually did it. Painted as a mask filled with `currentColor`, the browser
    resolves it on every repaint and no rerun is involved.
    """
    import base64
    w, h = _png_size(path_)
    url = f"url('data:image/png;base64,{base64.b64encode(Path(path_).read_bytes()).decode()}')"
    return (f"<span style=\"display:inline-block;height:{height_}px;width:{height_ * w / h:.0f}px;"
            f"background-color:currentColor;-webkit-mask:{url} no-repeat center/contain;"
            f"mask:{url} no-repeat center/contain\"></span>")


# The theme is pinned in .streamlit/config.toml, so nothing here has to guess it. Detecting it at
# runtime was the old approach and it was wrong twice over: the value followed the system
# preference rather than the theme actually rendered, and it was read once per script run, so
# switching the theme left the answer stale until something forced a rerun.
#
# The wordmark needs none of this -- it is painted with currentColor. The logo is a single-tone
# drawing with internal shading, which a mask would flatten, so it is inverted instead. If the
# theme is ever overridden to light, this inversion is the one line to drop with it.
_est_marca = _INV

_marca = []
if _logo_f and _logo_f.suffix.lower() != ".svg":
    _marca.append(_img_inline(_logo_f, 72, "", _est_marca))
if not _wm:
    _marca.append("<span style='font-size:2.1rem;font-weight:700'>PoliScreen</span>")
elif _wm.suffix.lower() == ".png":
    _marca.append(_img_tinted(_wm, 56))
else:
    _marca.append(_img_inline(_wm, 56, "", _est_marca))
_marca.append("<span style='color:currentColor;opacity:.65;font-size:.85rem;"
              f"align-self:flex-end;padding-bottom:.5rem'>v{__version__}</span>")
_izq.markdown("<div style='display:flex;align-items:center;gap:.7rem;margin:-.3rem 0 .4rem'>"
              + "".join(_marca) + "</div>", unsafe_allow_html=True)
_menu_archivo = _m1.popover(t("File"), width="stretch")
_menu_datos = _m2.popover(t("Data"), width="stretch")
_menu_cfg = _m3.popover(t("Settings"), width="stretch")
_menu_ayuda = _m4.popover(t("Help"), width="stretch")

with _menu_archivo:
    st.markdown(t("**Project**"))
    # Through default_root(), which honours POLISCREEN_PROJECTS. Building the path from home()
    # here ignored it, and inside Docker that wrote the results to the container's own home
    # instead of the mounted volume, where they were lost when the container went away.
    S.setdefault("proj_dir", str(ss.default_project()))
    _escrito = st.text_input(
        t("Project folder"), key="proj_dir",
        help=t("Any folder you can write to; `C:\\...` is used as it is.") if os.name == "nt"
        else t("Path inside Linux (WSL). If you paste a Windows path (`\\\\wsl.localhost\\...` or `C:\\...`) it is translated automatically."))
    proj, _path_notice = ss.normalize_path(_escrito)
    if _path_notice:
        st.warning(_path_notice)
        S["_proj_pending"] = str(proj)
    try:
        proj.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        st.error(t('Cannot create that folder: {v1}').format(v1=e))
        st.stop()
    st.caption(t('Results in `{v1}`').format(v1=proj))
    # A project is named after the day it starts. Reopen the interface after midnight and this
    # resolves to a new, empty folder: the analysis finished last night is one directory over, and
    # without this line nothing on screen says where it went.
    try:
        _empty = not any(proj.iterdir())
    except OSError:
        _empty = False
    if _empty:
        _prev = ss.previous_project(proj)
        if _prev:
            st.info(t('This project is empty — it is a new day. The previous one, with your '
                      'analysis in it, is `{v1}`: paste that path above to go back to it.'
                      ).format(v1=_prev))
    elif lay.artifact(proj, lay.DOCKING_CSV).exists():
        # The other half of the same edge. The day's folder already holds a finished analysis, and
        # a second one started here would mix with it -- the reason results were being deleted by
        # hand to make room. Offered, not automatic: reopening this folder is the normal case.
        _next = ss.next_project(proj)
        st.caption(t('This project already holds an analysis. To start another without mixing '
                     'them, use `{v1}`.').format(v1=_next))
        if st.button(t('Start a new analysis'), key='new_proj'):
            S['_proj_pending'] = str(_next)
            st.rerun()
    if S.get("_proj_loaded") != str(proj):
        S["_proj_loaded"] = str(proj)
        rec_dir, lig_dir = lay.artifact(proj, lay.RECEPTORS), lay.artifact(proj, lay.INPUT_LIGANDS)
        S["receptors"] = (sorted(str(p) for suf in ("_ready", "_listo")
                                 for p in rec_dir.glob(f"*{suf}.pdb"))
                          if rec_dir.exists() else [])
        # Not every .sdf here is a control: ccd_template caches the dictionary entry beside them,
        # and its idealised coordinates sit nowhere near the crystallographic pose. Loaded as a
        # control it becomes a second, phantom reference, and centring the box on it searches
        # empty space -- 9 A off the real site on 8HTB, with nothing saying so.
        S["controls"] = (sorted(str(p) for p in rec_dir.iterdir()
                                if not rc.is_ccd_cache(p)
                                and (p.suffix.lower() in (".sdf", ".mol2", ".mol")
                                     or (p.suffix.lower() == ".pdb" and p.name.startswith("control_"))))
                         if rec_dir.exists() else [])
        S["ligands"] = (sorted(str(p) for p in lig_dir.iterdir()
                               if p.suffix.lower() in (".sdf", ".mol2", ".mol", ".smi", ".pdbqt"))
                        if lig_dir.exists() else [])
    if S["receptors"] or S["controls"] or S["ligands"]:
        st.caption(t('Recovered from disk: {v1} receptor(s), {v3} control(s), {v5} ligand(s).').format(v1=len(S['receptors']), v3=len(S['controls']), v5=len(S['ligands'])))

    with st.expander(t("Session and export")):
        sub = st.file_uploader(t("Open session (.poliscreen)"), type=["poliscreen", "zip"],
                               help=t("Restores a previous analysis: tables, receptors and ligands. You can change the weighting without repeating the docking."))
        if sub is not None and st.button(t("Restore this session")):
            dest_ = ss.default_root() / Path(sub.name).stem
            tmp_s = dest_.parent / f"_{Path(sub.name).stem}.poliscreen"
            tmp_s.parent.mkdir(parents=True, exist_ok=True)
            tmp_s.write_bytes(sub.getvalue())
            try:
                ss.load_session(tmp_s, dest_)
                info = ss.session_info(tmp_s)
                est = ss.restore_paths(ss.read_state(tmp_s), dest_)
                for key_ in ("products", "core_smiles", "reagents", "reag_info", "pockets",
                              "pep_seqs", "pep_notice", "_log_run", "_log_state"):
                    if est.get(key_) is not None:
                        S[key_] = est[key_]
                S["_pending_widgets"] = est.get("widgets") or {}
                S["_proj_pending"] = str(dest_)
                st.success(t("Session '{v1}' restored in {v3}").format(
                    v1=info.get("project") or info.get("proyecto") or "?", v3=dest_))
                st.rerun()
            except Exception as e:
                st.error(t('Could not restore the session: {v1}').format(v1=e))
            finally:
                tmp_s.unlink(missing_ok=True)

        full_ = st.checkbox(t("Include poses and complexes (heavy session)"), value=False,
                               help=t("Unchecked, the session is a few MB and enough to reopen and re-score. Checked, it also lets you re-examine 3D structures."))
        if st.button(t("Save session")):
            try:
                _rj = proj / "run.json"
                if S.get("_ui_weights") and _rj.exists():
                    try:
                        _d = json.loads(_rj.read_text())
                        _d["weights"] = S["_ui_weights"]
                        _rj.write_text(json.dumps(_d, indent=2))
                    except Exception:
                        pass
                ui_keys = ("products", "core_smiles", "reagents", "reag_info", "pockets",
                             "pep_seqs", "pep_notice", "_log_run", "_log_state")
                estado_ui = {k: S.get(k) for k in ui_keys if S.get(k) is not None}
                estado_ui["widgets"] = {k: S[k] for k in S.keys()
                                        if isinstance(k, str) and k.startswith(_PERSISTENT_PREFIXES)
                                        and isinstance(S[k], (str, int, float, bool, list))}
                _tmp = Path(tempfile.mkdtemp())
                sfile = ss.save_session(proj, _tmp / f"{proj.name}", full_=full_,
                                        state_=ss.strip_paths(estado_ui, proj))
                S["_session_file"] = (sfile.name, sfile.read_bytes())
                st.success(f"{sfile.name} · {sfile.stat().st_size / 1e6:.1f} MB")
                shutil.rmtree(_tmp, ignore_errors=True)
            except Exception as e:
                st.error(t('Could not save: {v1}').format(v1=e))
        if S.get("_session_file"):
            st.download_button(t('Download {v1}').format(v1=S['_session_file'][0]), S["_session_file"][1],
                               file_name=S["_session_file"][0], mime="application/zip")

        _n_exp = sum(1 for v in ss.catalog(proj).values() if v["has"])
        if _n_exp:
            if st.button(t('Download results... ({v1} available)').format(v1=_n_exp),
                         width="stretch"):
                S["_open_downloads"] = True
                st.rerun()
            st.caption(t("You pick what to include with checkboxes and a single ZIP is built; nothing is written to the project folder."))
with _menu_datos:
    st.markdown(t("**Loaded in this project**"))
    st.write(t('Prepared receptors: **{v1}**').format(v1=len(S['receptors'])))
    st.write(t('Co-crystallized controls: **{v1}**').format(v1=len(S['controls'])))
    st.write(t("Ligands: **{n}**").format(n=len(S["ligands"])))

with _menu_ayuda:
    st.markdown(t("**PoliScreen manual**"))
    for _sec, _topics in ayuda.SECTIONS.items():
        with st.expander(t(_sec)):
            for _title, _body in _topics:
                st.markdown(f"**{t(_title)}**")
                st.markdown(t(_body))
                st.markdown("")
    _how_to_cite()

with _menu_cfg:
    st.radio(t("Language"), list(LANGUAGES), key="_lang_pick", horizontal=True,
             help=t("Only the interface translates; results, data and file names stay as generated."))
    st.markdown(t("**Appearance**"))
    st.caption(t("The interface is dark by default. Start with POLISCREEN_THEME=light for the "
                 "light one; the theme is fixed at start so the brand and the page always agree."))
    # Seeded through the key, not through a default value: these keys are reassigned every pass to
    # survive changing stage, and a widget cannot take both without Streamlit dropping the default.
    S.setdefault("cfg_split", 0.50)
    S.setdefault("cfg_height", 580)
    st.slider(t("Split between tools and viewer"), 0.3, 0.7, step=0.02, key="cfg_split",
              help=t("Left gives more space to the viewer; right, to the tools."))
    st.slider(t("Panel height (px)"), 380, 1200, step=20, key="cfg_height")

if S.pop("_open_downloads", False):
    _downloads_dialog(proj)


def _batch_chemotypes():
    """(hay_vina, hay_peptidos) present in the batch. Only decides WHICH settings to show; the real
    routing is done by the pipeline per ligand. Erring here at most shows one extra settings
    section, never sends a ligand to the wrong engine."""
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
    """Docking parameters. They live in the Run stage, which is where they are used.

    Only the settings of the engine or engines that will take part are shown: Vina for small
    molecules, ADCP for peptides, both in a mixed batch. A panel with settings that do nothing
    confuses about which engine is being used. The values of the hidden sections are returned at
    their default, so RunConfig always receives a complete configuration."""
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
            # 8 by default, which is fast enough to explore with. It is not converged: docking
            # one ligand against 8HTB with five seeds, one of the five landed 0.29 kcal/mol from
            # the other four, enough to move the top of the ranking. At 24 that outlier is gone
            # (spread 0.11) and 32 was no better. The help text says so; the choice is the user's,
            # because most of a session is exploring and only the final run has to be converged.
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
        # gnina's network was trained on protein-small molecule complexes; on peptides its score is extrapolation.
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


def _stage_receptors():
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
        # A modified residue is a heteroatom in the PDB but belongs to the chain, not a cofactor.
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
                    # Older projects wrote control_cadena{c}; nothing parses the name back, and the
                    # controls are found by globbing control_*, so those files keep loading.
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
        # These keys deliberately avoid the persistent prefixes above: that loop reassigns each
        # value to keep it across reruns, and Streamlit forbids assigning a value to a button.
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

def _drug_mode():
    """Approved drugs as a ligand source, filtered by property.

    Deliberately separate from the reaction builder, and not merely a different button: nothing
    here is designed or made, so there is no synthesizability to judge and no reaction under which
    to be feasible. Asking whether an approved drug is "synthesizable by this reaction" is a
    category error, and the two must not share that column.

    They do share the run. Compounds land in the same input folder as everything else, so a
    screening of five products from the builder beside five drugs from here is one experiment, not
    two.
    """
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
        done_set = {nm for nm, _, _ in made}
        # The whole folder, not only what was just built, so drugs sit beside anything the
        # reaction builder or an upload already put there. That is the point of the section.
        S["ligands"] = [str(p) for p in sorted(d.iterdir()) if p.is_file()]
        rows_ = [{"name": nm, "smiles": c.get("smiles"), "source": "chembl",
                  "product": c.get("name"), "iupac_name": None, "feasibility": None}
                 for nm, c in zip(names_, chosen) if nm in done_set]
        _merge_ligand_meta(rows_)
        _mark_done("use_drugs", signature)
        _notify(t('{v0} compounds built and ready for step 3.').format(v0=len(made)))
        st.success(t("{v0} drugs added. The screening now has {v2} compounds in total, from "
                     "every source you have used.").format(v0=len(made), v2=len(S["ligands"])))


def _merge_ligand_meta(rows):
    """Adds these rows to the project's ligand table, keeping what other sources already wrote.

    Every source writes the same file, and one overwriting another is how a mixed run loses half
    its provenance: the compounds still dock, and the table can no longer say where they came from.
    Matched on `name`, which is the 3D file stem and therefore what the rest of the pipeline joins
    on; a repeated name is the same compound rebuilt, so the new row wins.
    """
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


def _peptide_mode():
    """Peptide design: a path independent of reaction synthesis. Kept apart
    (its own state keys) so it does not mix with the products of the chemical builder."""
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

    # The termini are chosen before generating: they change the net charge and the built structure.
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
    # Built with the termini chemistry used in the screening, or cyclization would not be visible.
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
        done_set = {nm for nm, _, _ in made}
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


def _stage_ligands():
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
        _drug_mode()
    elif modo == "Generate peptides":
        _peptide_mode()
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
                # Without OPSIN nothing can be checked against the structure, and the count said so
                # honestly -- "0 of 44 verified" -- while names for the simple esters appeared in
                # the table anyway, because those are composed rather than verified. Two true
                # statements that read as a contradiction. A composed name is a guess: on a
                # benzofuroxan the N-oxide can sit on either nitrogen and both parse, so an
                # unchecked name is exactly the one that would reach a paper wrong. The button is
                # closed rather than left to produce names nobody can vouch for.
                if not nm.available():
                    st.button(t("Name (IUPAC, verified with OPSIN)"), disabled=True)
                    # The reasoning lives in Help and in docs/INSTALL.md. An explanation this
                    # long belongs where someone goes looking for it, not in the way of the work.
                    st.caption(t("IUPAC naming is off in this build. Products are identified by "
                                 "their SMILES, which is what gets docked. See Help."))
                elif st.button(t("Name (IUPAC, verified with OPSIN)")):
                    with st.spinner(t("Naming and verifying by round-trip...")):
                        named = AdmelabBridge().name_esters(
                            [p["smiles"] for p in prods], [p.get("alcohol_smiles") or "" for p in prods],
                            acid_smiles=nuc, alcohol_names=[p.get("product") for p in prods], use_web=True)
                    by = {n["smiles"]: n for n in named}
                    # On a benzofuroxan the N-oxide can be named on the wrong nitrogen and still parse: the name is checked against the structure.
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
                done_set = {nm for nm, _, _ in made}
                S["ligands"] = [str(p) for _, p, _ in made]
                meta = pd.DataFrame([{"name": nm, "smiles": p.get("smiles"), "source": p.get("source", "?"),
                                      "product": p.get("product"), "iupac_name": p.get("iupac_name"),
                                      "feasibility": p.get("feasibility")}
                                     for (nm, p) in zip(names_, chosen_ones) if nm in done_set])
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

def _stage_run():
    """Two things that can be launched from here, kept apart because one needs engines this
    package does not ship."""
    tab_screening, tab_tunnels = st.tabs([t("Screening"), t("Transport tunnels")])
    with tab_tunnels:
        _run_tunnels()
    with tab_screening:
        _run_screening()


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
    for label, ok, var, how in rows:
        if ok:
            st.success(f"{label} — {t('ready')}")
        else:
            st.info(f"**{label}** — {t('not installed')}. {how}")
    return all(ok for _l, ok, _v, _h in rows)


def _run_tunnels():
    """CAVER and CaverDock, driven from the box PoliScreen already placed on the active site.

    Nothing here is required to use the rest of PoliScreen, and nothing here is bundled: the two
    engines are found on the machine or reported absent. Results are read back in the Results tab,
    which needs neither of them.
    """
    st.subheader(t("Transport tunnels"))
    st.caption(t("A docking score says how well a compound sits in the site. This asks whether it "
                 "can reach it: CAVER finds the routes through the protein, CaverDock costs one."))

    has_caver = cv.caver_available()
    has_dock = cv.caverdock_available()
    if not (has_caver and has_dock):
        _engine_status()

    # Tunnels already found stay readable and drawable whatever is installed now. Computing them
    # needs CAVER; looking at what was computed does not, and an engine uninstalled later is no
    # reason to hide a result that is already on disk.
    existing = S.get("tun_drawn")
    if existing and cv.clusters(existing):
        _tunnel_table(existing, cv.clusters(existing))

    recs = [Path(p) for p in S.get("receptors", []) if Path(p).exists()]
    if not recs:
        st.info(t("Stage a receptor in step 1 first."))
        return

    out_root = proj / lay.TUNNELS
    rec = st.selectbox(t("Receptor"), recs, format_func=_rname, key="tun_rec")

    # --- 1. the structure ----------------------------------------------------------------------
    st.markdown(t("**1 · What CAVER should look at**"))
    st.caption(t("Hydrogens make every atom effectively larger and close the narrow routes, so "
                 "they are removed. A receptor prepared for docking is not the right input here."))

    # The original download is offered first because it is the one that reproduces a CaverWeb run:
    # preparing a receptor for docking also drops atoms CAVER wants.
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

    # --- 2. where to measure from --------------------------------------------------------------
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

    # CAVER measures outwards from this point. Placed outside the protein it finds the outside:
    # one enormous tunnel, a wide bottleneck, and a warning nobody reads.
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
    # The table itself was drawn at the top, so it is there whether or not the engines are.
    S["tun_drawn"] = str(caver_out)

    # --- 2. the transport ----------------------------------------------------------------------
    st.markdown(t("**2 · Push a compound through one**"))
    if not has_dock:
        return
    # The two engines want different files and get different files, which is easy to miss and
    # changes what the numbers mean. CAVER measured the route on the stripped structure; CaverDock
    # docks with Vina's scoring, which needs the hydrogens and charges of the docking receptor.
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

    # Entering and leaving are two separate dockings of the same route and the barrier is not the
    # same in both; lower and upper bound are two different quantities. All four are usually
    # wanted, so all four are checkboxes rather than a choice of one.
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

    # Asking for both bounds is one job, not two. CaverDock computes the lower bound on its way to
    # the upper one and writes both profiles into the same folder -- measured on this project, the
    # separate lower-bound run took 7 min 33 s per compound and produced nothing the upper-bound
    # run had not already produced. Running it as well also puts two rows in the table for one
    # calculation, one of them marked as having no upper bound.
    bounds = ["ub"] if want_ub else (["lb"] if want_lb else [])
    jobs = [(cmp_, tun, d, b) for cmp_ in ligands_pick for tun in tunnels_pick
            for d in directions for b in bounds]

    if not jobs:
        st.info(t("Pick at least one tunnel, one compound, one direction and one bound."))
        return
    if want_ub and want_lb:
        st.caption(t("Both bounds is one calculation: the upper bound produces the lower one on "
                     "its way."))
    # An upper bound is roughly five times a lower one, and every combination is a separate docking.
    minutes = sum(28 if b == "ub" else 8 for _c, _t, _d, b in jobs)
    st.caption(t("{n} calculations, roughly {m} minutes in total.").format(n=len(jobs), m=minutes))

    if not cv.reproducible(int(cpus)):
        # CaverDock says this once, in the middle of a log nobody reads.
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
                    # One failed combination is not a reason to lose the ones that worked.
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


def _run_screening():
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
        xbox = {}
        # Control assigned by geometry: the box, the reference fingerprint and the redocking validation depend on it.
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
                # Kept in state because the rerun below wipes anything drawn now, and a run that
                # found nothing is exactly when the reason matters.
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
                        # Flexibility omitted: fpocket derives it from B-factors, which PDBFixer zeroes.
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
            # The ligand must be able to rotate in the box, not merely fit: otherwise the search is restricted to the orientations that enter.
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
                    # Docking is the long phase: it gets a bar instead of one line per job.
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


def _run_log_panel():
    """Log of the last run, preserved across tabs."""
    state_ = S.get("_log_state", "completo")
    with st.status("Screening completed" if state_ == "completo" else "The run failed",
                   state="complete" if state_ == "completo" else "error", expanded=False):
        for n, d in S.get("_log_run", []):
            st.write(f"**{n}** · {d}")
        _download_table(pd.DataFrame(S["_log_run"], columns=["stage", "detail"]),
                         "registro_corrida", key="log_run")

def _pleiotropic_summary(rk, targets_):
    """Which compound binds well in SEVERAL targets at once. A per-target summary is already given by
    the blocks below; this one looks for broad spectrum, not potency in a single one.

    Effectiveness is compared across targets because it is normalized against each one's control: it
    is the same 0-1 scale even though the energies are not. It is sorted by the Minimum effectiveness
    across targets —the weakest link—, which is what really measures 'good in all'; an average would
    let through a compound excellent in one target and useless in another."""
    import numpy as np
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


def _stage_results():
    """Two questions about the same target, kept apart because they are answered separately.

    The screening asks how well a compound sits in the site. The tunnels ask whether it can reach
    it -- from a calculation run outside PoliScreen, which is why it is a tab and not a step.
    """
    tab_screening, tab_tunnels = st.tabs([t("Screening"), t("Transport tunnels")])
    with tab_tunnels:
        _results_tunnels()
    with tab_screening:
        _results_screening()


def _with_tunnel_geometry(table, folder):
    """Join CAVER's own numbers for each tunnel onto the transport rows.

    Ea says what entering costs. Whether that means anything depends on the tunnel: one with
    nothing to cross costs nothing to cross. The geometry was computed when the routes were found,
    and it lives in the CAVER output beside the runs rather than inside them, so it is looked up
    here instead of arriving with the profile.
    """
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
    """One row per receptor, compound, tunnel and direction: the richest calculation of each.

    Asking for both bounds used to run two jobs, and CaverDock computes the lower bound on its way
    to the upper one, so the pair appears as two rows of the same numbers -- one of them marked as
    having no upper bound. New runs do not produce that, but the folders already on disk do, and
    two rows for one calculation is not something to leave for the reader to reconcile.
    """
    if table.empty or "flags" not in table:
        return table
    keys = [c for c in ("receptor", "ligand", "tunnel", "direction") if c in table]
    if not keys:
        return table
    # A row with an upper bound outranks one without; a refusal outranks silence about it.
    rank = table["flags"].fillna("").map(
        lambda f: 0 if "failed" in f.split() else 1 if "lower_bound_only" in f else 2)
    ordered = table.assign(_rank=rank).sort_values("_rank", ascending=False)
    return ordered.drop_duplicates(subset=keys, keep="first").drop(columns="_rank")


def _readable_transport(table):
    """The same rows, named for reading, with what happened said in words.

    A run that did not finish used to appear as a row of blanks, which reads as a result of zero
    rather than as an absence.
    """
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
    """Which route each compound takes, and by how much it beats its next one.

    A table sorted by Ea answers "what is easiest anywhere in this dataset". The question usually
    being asked is narrower: *this* compound, which way does it go. The margin matters as much as
    the winner -- two routes within a few tenths of a kcal/mol is not a preference, it is a tie,
    and CaverDock's own repeatability is of that order.
    """
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
    # Every calculation may have failed, and an empty frame has no column to sort on.
    return pd.DataFrame(rows).sort_values("Ea") if rows else pd.DataFrame()


def _route_preferences(table):
    """The per-compound answer, above the full table."""
    preference = route_preference(table)
    if preference.empty or preference["routes"].max() < 2:
        return                       # with one route each there is no preference to report
    with st.expander(t("Which route each compound takes"), expanded=True):
        shown = preference.rename(columns={
            "ligand": t("Compound"), "tunnel": t("Tunnel"), "direction": t("Direction"),
            "margin": t("Beats the next by"), "routes": t("Routes tried")})
        st.dataframe(shown, width="stretch", hide_index=True)
        st.caption(t("Lowest Ea per compound. A margin under about 0.5 kcal/mol is a tie, not a "
                     "preference: that is the order of CaverDock's own repeatability."))


def _results_tunnels():
    """CAVER and CaverDock output, read into the same kind of table as everything else.

    Nothing here runs an engine. The folder is asked for rather than discovered because the
    calculation does not live in the project: it comes back from CaverWeb, or from a local
    caverdock-run, and the same reader handles both.
    """
    st.subheader(t("Transport tunnels"))

    if not tn.available():
        st.info(t("Tunnel reading needs **caver-translate**, a separate package with no "
                  "dependencies of its own. The rest of PoliScreen works without it."))
        st.code(tn.INSTALL_HINT, language="bash")
        return

    # This project's own runs, without being asked for: PoliScreen wrote them and knows where.
    # The box is for the other case -- a CaverWeb download, which lives wherever it was unzipped.
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
    except Exception as e:                       # a half-written run is common; do not lose the tab
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


def _results_screening():
    st.subheader(t("Results"))
    meta_p, inter_p, dock_p = proj / "run.json", lay.artifact(proj, lay.INTERACTIONS_CSV), lay.artifact(proj, lay.DOCKING_CSV)
    if not (meta_p.exists() and inter_p.exists()):
        st.info(t("No results in this folder yet. Run step 3."))
    else:
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
            # Catalytic and secondary are exclusive roles; the engine gives priority to the catalytic one.
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

        # Verifies the docking only: compounds are compared with the crystallographic ligand, not with this pose.
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
        # Handed to the summary panel, which used to re-read ranking.csv from disk. That file was
        # written with the weights of the run; this table is recomputed with the weights on screen.
        # Move a slider and the two disagreed about the same compound, with nothing saying why.
        S["_rk_live"] = rk.copy()
        rk["Ki"] = rk["pred_ki_M"].map(_fmt_ki) if "pred_ki_M" in rk.columns else None
        # An axis weighted without data does not score, and must not be declared in Methods.
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

        # ADCP and Vina energies come from different scoring functions and are not comparable.
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
        # The control first, not whatever sorts first alphabetically. This diagram is read to judge
        # every other one -- the contacts it reproduces are the reference the ranking is built on --
        # so opening on 1-butanol meant a deliberate step before the panel said anything useful.
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
                import io as _io
                _b = _io.BytesIO(); fig_int.savefig(_b, format="png", dpi=160, bbox_inches="tight")
                _download_image(_b.getvalue(), f"interaccion_{cmp_}_pose{mod}", key=f"int_{R}_{cmp_}_{mod}")
            except Exception:
                pass
            st.caption(t("Green = reproduces a control interaction (same residue and same bond). Gray = extra contact or the same residue with a different bond type."))

        st.markdown("---")
        _how_to_cite()

def _tunnel_table(caver_out, found):
    """The tunnels, with a box to pick which ones are drawn and the colour they are drawn in.

    Six routes at once is a knot. The colour circle is the same device the cavity table uses, so a
    row is matched to what is on screen without reading a legend.
    """
    geometry = {}
    summary = Path(caver_out) / "summary.txt"
    if summary.exists() and tn.available():
        from caver_translate.parse import parse_tunnels
        geometry = {g.tunnel: g for g in parse_tunnels(summary)}

    numbers = [_tunnel_number(c) for c in found]

    # A row of real checkboxes rather than st.data_editor. The editor keeps the user's edits as a
    # diff and replays it over the frame it is given, and that frame is rebuilt from the state the
    # diff just produced: the two disagree for one run, which is felt as having to click twice.
    # st.checkbox has no diff -- its value IS session state -- so what it returns is what was just
    # done. It is also the only one of the two that AppTest can operate, so this table is testable
    # and the other was not.
    widths = [0.6, 1.6, 1.2, 1.2, 1.2, 1.2]
    head = st.columns(widths, vertical_alignment="bottom")
    for col, label in zip(head, ["", t("Tunnel"), t("Bottleneck (Å)"), t("Length (Å)"),
                                 t("Curvature"), t("Priority")]):
        col.caption(label)

    shown = []
    for i, n in enumerate(numbers):
        g = geometry.get(n)
        dot = vw.emoji_for_color(vw.TUNNEL_PALETTE[(n - 1) % len(vw.TUNNEL_PALETTE)])
        # Only the first is on to begin with; after that every row remembers what was done to it,
        # and unticking them all leaves the viewer empty, which is what unticking them all means.
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


def _tunnel_number(cluster) -> int:
    """The cluster's own number, which is its rank: tun_cl_003 and tun_cl_003_1 are both 3."""
    digits = re.search(r"tun_cl_0*(\d+)", Path(cluster).stem)
    return int(digits.group(1)) if digits else 0


def _tunnel_groups(only=None, root=None):
    """The tunnels found for this project, as sphere groups the viewer can draw.

    A CAVER tunnel cluster is a chain of spheres in a PDB, which is the same shape fpocket's alpha
    spheres arrive in, so nothing new has to be taught to the viewer. Opaque and saturated, unlike
    the cavities: a cavity is a volume to judge through, a tunnel is a route to follow, and a
    translucent tube through a ribbon is harder to trace, not easier.

    `only` is the set of tunnel numbers to draw. Six routes at once is a knot; the point of the
    picture is one of them. `root` defaults to the CAVER output this session last produced, and is
    passed in by the tests, which run outside the session this reads.
    """
    root = root if root is not None else S.get("tun_drawn")
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


def _viewer_height(reserve: int) -> int:
    """Height of the 3D viewer so the panel shows selectors, viewer and footer without scroll. The
    height of each stage's chrome is subtracted from the panel height (cfg_alto); the minimum keeps
    it from ending up tiny if the panel is very short."""
    return max(190, int(S.get("cfg_height", 580)) - reserve)


def _scene_height() -> int:
    """Height for a panel whose 3D scene is the point and whose other parts can be scrolled to.

    Subtracting the chrome works when everything has to fit at once. The transport panel also
    carries a plot, four numbers and a button, and subtracting all of that left the molecules in a
    letterbox. The container scrolls, so the scene takes most of the panel and the rest is below
    it -- which is also the order it gets read in.
    """
    return max(380, int(int(S.get("cfg_height", 580)) * 0.86))


def _viewer_panel(etapa: str):
    """Output of the active stage. Drawn after the tools panel, so it can read what the latter has
    just left in session_state (for example the search box)."""
    if etapa == "Receptors":
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

    elif etapa == "Ligands":
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
            # Reopening a folder restores the ligand files but not what was built in this session,
            # so the structures are redrawn from the table on disk: seeing what is about to be
            # docked is the point of this panel.
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

    elif etapa == "Run":
        boxes_ = S.get("_boxes") or {}
        cav_map = S.get("_cavities") or {}
        if boxes_:
            c1, c2, c3, c4 = st.columns([2, 1, 1, 1])
            rsel = c1.selectbox(t("Receptor"), list(boxes_), format_func=_rname,
                                key="vis_box_rec", label_visibility="collapsed")
            groups_by_receptor = cav_map.get(rsel)
            S.setdefault("vis_show_cav", bool(groups_by_receptor))
            ver_cav = c2.checkbox(t("Cavities"), key="vis_show_cav")
            # A CAVER tunnel is a chain of spheres, which is the shape the viewer already draws
            # cavities in. Both on by default: a route is read against the pockets it connects.
            # Absent and empty both mean none: before the table has been drawn nothing has been
            # chosen, and after every row is unticked nothing was chosen. Neither draws anything.
            routes = _tunnel_groups(set(S.get("tun_shown") or []))
            # Disabled only when there are none to draw at all. Basing it on the current selection
            # instead greys the control out the moment the last row is unticked, which reads as a
            # fault rather than as the empty view that was asked for.
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

    else:
        views = ["Summary", "3D complex"]
        if tn.available() and tn.runs_in(proj / lay.TUNNELS):
            views.append("Transport")
        view_ = st.radio(t("View"), views, horizontal=True,
                         format_func=t, key="vis_res_view", label_visibility="collapsed")
        if view_ == "3D complex":
            _complex_viewer()
        elif view_ == "Transport":
            _transport_viewer()
        else:
            _visual_summary()


def _transport_viewer():
    """One compound through one tunnel: the poses that matter, and the profile they came from.

    A trajectory is one pose per disc -- eighty of them here. Three are worth showing, and which
    three is a question about the energy, not about spacing: where it starts, where it is hardest,
    and where it ends. The plot beside them says where the reported numbers came from, so Ea and
    dE_BS do not have to be taken on trust.
    """
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
    # How many context poses fit is a property of the tunnel, so the range comes from its length
    # rather than from a number picked once and applied to every route.
    most = tn.most_extra(profile)
    if most:
        S.setdefault("vis_tun_poses", 3 + tn.suggested_extra(profile))
        # Clamped before the slider is drawn, not passed to it: this key is reassigned every pass
        # to survive changing stage, and a widget cannot take both. Another route's answer can sit
        # outside this one's range, and a slider handed a value past its maximum raises.
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

    # Scene and profile together: the poses are the route and the plot says why those poses. Split
    # across tabs, matching a molecule to its point on the curve costs a click each time.
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
        # Cofactors kept in the receptor are drawn: one sitting in the route is the reason a
        # barrier is where it is, and an unexplained wall reads as a fault in the calculation.
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


def _complex_viewer():
    """Receptor-ligand complex in 3D, stepping through compound and pose.

    The 2D diagram of the tools panel says WHICH contacts there are; this says WHERE. The same
    interaction table is read and the same keys are used, so that the pose chosen here and the one
    in the diagram are the same entity and not two numberings to be matched by hand.
    """
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


def _visual_summary():
    """Results summary meant to be read at a glance: the exhaustive detail stays in the tools panel;
    here goes what you show someone in ten seconds."""
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


_STAGE_FN = {"Receptors": _stage_receptors, "Ligands": _stage_ligands,
             "Run": _stage_run, "Results": _stage_results}

_HEIGHT = int(S.get("cfg_height", 580))
_rep = float(S.get("cfg_split", 0.50))
_izq, _der = st.columns([_rep, 1.0 - _rep], gap="medium")
with _izq:
    st.markdown(t("**Tools · {stage}**").format(stage=t(S["stage"])))
    with st.container(height=_HEIGHT, border=True, key="left_panel"):
        _STAGE_FN[S["stage"]]()
with _der:
    st.markdown(t("**Viewer**"))
    with st.container(height=_HEIGHT, border=True, key="right_panel"):
        _viewer_panel(S["stage"])

_nav = st.container(key="stage_bar")
_cols = _nav.columns(len(STAGES))
for _i, _e in enumerate(STAGES):
    if _cols[_i].button(t(_e), key=f"nav_{_e}", width="stretch",
                        type=("primary" if _e == S["stage"] else "secondary")):
        S["stage"] = _e
        st.rerun()

# Opened last, so it also catches what the active stage just did: a run finishing sets the notice
# while this script is already past the panels.
if S.get("_notice"):
    _headline, _detail = S.pop("_notice")
    _notice_dialog(_headline, _detail)
