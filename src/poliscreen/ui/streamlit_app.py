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

# Re-exports for backward compatibility and test contracts
from poliscreen.ui.common import (
    CITATIONS,
    ACKNOWLEDGMENTS,
    _confirm_delete,
    _controls_of,
    _download_image,
    _download_table,
    _empty_state,
    _fmt_ki,
    _forget_all_receptors,
    _forget_receptor,
    _how_to_cite,
    _rname,
    _to_smiles,
)
from poliscreen.ui.components.admet import (
    _render_adme,
    _scatter_dock_inter,
    _shade,
)
from poliscreen.ui.components.transport import (
    STATUS_TEXT,
    TRANSPORT_LABELS,
    _one_row_per_route,
    _readable_transport,
    _route_preferences,
    _tunnel_groups,
    _tunnel_number,
    _tunnel_table,
    _with_tunnel_geometry,
    route_preference,
)
from poliscreen.ui.views.receptors import (
    render_receptors_tools,
    render_receptors_viewer,
)
from poliscreen.ui.views.ligands import (
    render_ligands_tools,
    render_ligands_viewer,
)
from poliscreen.ui.views.run import (
    _batch_chemotypes,
    _docking_params,
    render_run_tools,
    render_run_viewer,
)
from poliscreen.ui.views.results import (
    render_results_tools,
    render_results_viewer,
)


def _notify(headline: str, detail: str = "") -> None:
    """Report the end of a long operation. Shown as a modal, the same mechanism as the deletion
    prompt: a toast clips anything longer than a line, and custom HTML does not survive
    Streamlit's sanitizing, so neither was ever seen."""
    st.session_state["_notice"] = (headline, detail)


@st.dialog("PoliScreen")
def _notice_dialog(headline: str, detail: str):
    st.markdown(f"### {headline}")
    if detail:
        st.caption(detail)
    if st.button(t("Close"), type="primary", width="stretch"):
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
            st.session_state[f"dl_{k}"] = bool(mark and v["reason"])
        else:
            st.session_state.setdefault(f"dl_{k}", bool(v["reason"]))

    rec_items = {k: v for k, v in has_.items() if v["reason"]}
    otros_items = {k: v for k, v in has_.items() if not v["reason"]}
    if rec_items:
        st.markdown(t("**Recommended**"))
        for k, v in rec_items.items():
            st.checkbox(f"{v['desc']} · {_human_size(v['bytes'])}", key=f"dl_{k}")
            st.caption(v["reason"])
    if otros_items:
        st.markdown(t("**Optional**"))
        for k, v in otros_items.items():
            st.checkbox(f"{v['desc']} · {_human_size(v['bytes'])}", key=f"dl_{k}")
            if v["regenerable"]:
                st.caption(t("Regenerated by running again; takes space in the package."))

    chosen_items = [k for k in has_ if st.session_state.get(f"dl_{k}")]
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
            st.session_state["_zip"] = (f"{Path(folder).name}_PoliScreen.zip", datos, included_items)
        except Exception as e:
            st.error(t('Could not build the package: {v1}').format(v1=e))
    if st.session_state.get("_zip"):
        name_, datos, included_items = st.session_state["_zip"]
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


def _stage_receptors():
    return render_receptors_tools(proj)


def _stage_ligands():
    return render_ligands_tools(proj)


def _stage_run():
    return render_run_tools(proj)


def _stage_results():
    return render_results_tools(proj)


def _viewer_panel(stage: str):
    if stage == "Receptors":
        render_receptors_viewer(proj)
    elif stage == "Ligands":
        render_ligands_viewer(proj)
    elif stage == "Run":
        render_run_viewer(proj)
    elif stage == "Results":
        render_results_viewer(proj)


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
