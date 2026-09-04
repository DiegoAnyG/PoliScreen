"""Shared UI helper functions, formatters, and download utilities for PoliScreen."""
from __future__ import annotations

import io
from pathlib import Path

import pandas as pd
import streamlit as st

from poliscreen import __version__
from poliscreen.core import layout as lay
from poliscreen.core import pipeline as pl
from poliscreen.core import receptor as rc
from poliscreen.core import screening as sc
from poliscreen.core import viewer as vw
from poliscreen.ui.i18n import t

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

_SUMMARY_STEPS = ("3d", "docking", "plip", "rescoring", "validation")


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


def _to_smiles(path):
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    suf = Path(path).suffix.lower()
    m = (Chem.MolFromMol2File(str(path)) if suf == ".mol2" else
         next(iter(Chem.SDMolSupplier(str(path))), None) if suf == ".sdf" else
         Chem.MolFromMolFile(str(path)) if suf == ".mol" else None)
    return Chem.MolToSmiles(m) if m is not None else None


def _how_to_cite(expanded_: bool = False):
    with st.expander(t("How to cite"), expanded=expanded_):
        st.markdown(CITATIONS)
    with st.expander(t("Acknowledgments")):
        st.markdown(t(ACKNOWLEDGMENTS))


def _download_table(df, name_: str, key: str):
    """Offers downloading a DataFrame as CSV or XLSX below the table."""
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
    """Offers downloading an image (PNG) below the viewer."""
    if not png:
        return
    c = st.columns([3, 1])
    c[0].caption(t("Download as:"))
    c[1].download_button(t("PNG"), png, file_name=f"{name_}.png", mime="image/png",
                         key=f"dl_png_{key}", width="stretch")


def _already_done(key_: str, firma) -> bool:
    """True if the action already ran with THESE parameters."""
    return st.session_state.get("_signature_" + key_) == firma


def _mark_done(key_: str, firma):
    st.session_state["_signature_" + key_] = firma


def _load_control_map(rec_dir: Path) -> dict:
    cmap = dict(st.session_state.get("_control_map") or {})
    if rec_dir.is_dir():
        mf = rec_dir / "control_map.json"
        if mf.exists():
            try:
                import json
                for k, v in json.loads(mf.read_text(encoding="utf-8")).items():
                    cmap.setdefault(sc.normalize_key(k), v)
            except Exception:
                pass
    return cmap


def _save_control_map(rec_dir: Path, cmap: dict) -> None:
    if rec_dir.is_dir():
        mf = rec_dir / "control_map.json"
        try:
            import json
            mf.write_text(json.dumps(cmap, indent=2), encoding="utf-8")
        except Exception:
            pass


def _controls_of(rec: Path, receptors: list, controls: list) -> list:
    """Controls belonging to a receptor, prioritizing manual/saved mapping then geometry."""
    cmap = _load_control_map(Path(rec).parent)
    assign = pl._assign_controls([Path(c) for c in controls],
                                 [Path(r) for r in receptors], cmap)
    target_keys = {rec.stem, sc.normalize_key(rec.stem)}
    return [c for c in controls
            if assign.get(sc.normalize_key(Path(c).stem)) in target_keys]


def _forget_receptor(path_str: str) -> None:
    """Deletes a prepared receptor together with its controls and the choices made for it."""
    S = st.session_state
    rec = Path(path_str)
    doomed = [path_str] + _controls_of(rec, S["receptors"], S["controls"])
    for f in doomed:
        Path(f).unlink(missing_ok=True)
    S["receptors"] = [p for p in S["receptors"] if p != path_str]
    S["controls"] = [p for p in S["controls"] if p not in doomed]
    S.setdefault("_forget_prep", []).append(sc.normalize_key(_rname(rec)))
    if S.get("_control_map"):
        for doomed_f in doomed:
            k = sc.normalize_key(Path(doomed_f).stem)
            S["_control_map"].pop(k, None)
        _save_control_map(rec.parent, S["_control_map"])
    if S.get("last_prepared") == path_str:
        S.pop("last_prepared", None)
        S.pop("last_original", None)


def _forget_all_receptors() -> None:
    """Empties the project's receptor folder: prepared structures and controls alike."""
    S = st.session_state
    first_p = Path(S["receptors"][0]) if S["receptors"] else (Path(S["controls"][0]) if S["controls"] else None)
    if first_p and first_p.parent.exists():
        (first_p.parent / "control_map.json").unlink(missing_ok=True)
    for f in list(S["receptors"]) + list(S["controls"]):
        Path(f).unlink(missing_ok=True)
    S["_forget_prep"] = [sc.normalize_key(_rname(p)) for p in S["receptors"]]
    S["receptors"], S["controls"] = [], []
    S.pop("_control_map", None)
    S.pop("last_prepared", None)
    S.pop("last_original", None)


def _notify(headline: str, detail: str = "") -> None:
    """Report the end of a long operation. Shown as a modal dialog."""
    st.session_state["_notice"] = (headline, detail)


def _run_summary(log) -> str:
    """One line per phase, with what each one processed."""
    seen, out = set(), []
    for name_, detail in reversed(list(log or [])):
        if name_ in _SUMMARY_STEPS and name_ not in seen and detail:
            seen.add(name_)
            out.append(f"{name_}: {detail}")
    return " · ".join(reversed(out))


def _empty_state(message: str):
    """Empty state of the viewer: the logo as a watermark and a line with what is missing."""
    msg = t(message)
    logo = vw.logo_path()
    if logo and logo.suffix.lower() != ".svg":
        c = st.columns([1, 2, 1])[1]
        c.image(str(logo), width="stretch")
        c.markdown(f"<div style='text-align:center;opacity:.75;font-size:.92rem'>{msg}</div>",
                   unsafe_allow_html=True)
        return
    brand = logo.read_text(errors="ignore") if logo else vw.logo_svg()
    st.markdown(
        f"<div style='display:flex;flex-direction:column;align-items:center;justify-content:center;"
        f"padding:3rem 1rem;opacity:.75;text-align:center'>{brand}"
        f"<div style='margin-top:.9rem;font-size:.92rem;opacity:.8'>{msg}</div></div>",
        unsafe_allow_html=True)


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


def _viewer_height(reserve: int) -> int:
    """Height of the 3D viewer so the panel shows selectors, viewer and footer without scroll."""
    return max(190, int(st.session_state.get("cfg_height", 580)) - reserve)


def _scene_height() -> int:
    """Height for a panel whose 3D scene takes most of the panel."""
    return max(380, int(int(st.session_state.get("cfg_height", 580)) * 0.86))
