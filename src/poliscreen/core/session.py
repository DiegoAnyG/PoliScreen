"""Work sessions and result export.

A `.poliscreen` session is a ZIP with a manifest. It lets you close the app and come back later to
reopen the analysis, change the weighting and re-examine the results **without repeating the
docking** and without having to remember the working folder path.

Two sizes:
  - light (default): configuration, tables and receptors. Enough to reopen, re-score and view
    diagrams. On the order of megabytes.
  - full: adds poses and complexes, so the 3D structures can also be re-examined and PLIP re-run.
    It can take hundreds of megabytes.
"""
from __future__ import annotations

import io
import json
import os
import re
import zipfile
from datetime import datetime, timezone
from pathlib import Path

from . import layout as lay
from typing import Optional, Sequence

FORMAT = 1
EXT = ".poliscreen"

_UNC_WSL = re.compile(r"^/{2,}(?:wsl\.localhost|wsl\$)/[^/]+", re.I)
_UNIDAD = re.compile(r"^([A-Za-z]):(?=/|$)")


def default_root() -> Path:
    """Where projects are created when the user types nothing.

    POLISCREEN_PROJECTS lets the container point this at its mounted volume: inside Docker the
    home directory is not persistent, so results written there disappear with the container.
    """
    env = os.environ.get("POLISCREEN_PROJECTS")
    return Path(env) if env else Path.home() / "poliscreen_proyectos"


def normalize_path(text_: str, base: Optional[Path] = None) -> tuple:
    """Turns what the user types into a usable POSIX path.

    Returns (path, notice). The notice is non-empty when the text is reinterpreted: it must be seen,
    because writing to another folder silently ruins the analysis and leaves the results where no one
    looks.
    """
    raw_text = (text_ or "").strip().strip('"').strip("'")
    if not raw_text:
        return default_root() / "demo", ""

    s = raw_text.replace("\\", "/")
    s = _UNC_WSL.sub("", s) or "/"
    m = _UNIDAD.match(s)
    if m:
        s = f"/mnt/{m.group(1).lower()}/{s[m.end():].lstrip('/')}".rstrip("/")
    s = re.sub(r"/{2,}", "/", s) or "/"

    p = Path(s).expanduser()
    if not p.is_absolute():
        first_ = p.parts[0] if p.parts else ""
        p = Path("/") / p if (Path("/") / first_).is_dir() else (base or Path.home()) / p

    if str(p) == raw_text:
        return p, ""
    if "\\" in raw_text:
        notice = (f"A Windows path was detected. PoliScreen runs inside Linux (WSL), "
                 f"so `{p}` will be used.")
    else:
        notice = f"Path adjusted to `{p}`."
    return p, notice

BASE_FILES = ("run.json", "ranking.csv", lay.SUMMARY_CSV, lay.INTERACTIONS_CSV,
                 lay.DOCKING_CSV, lay.VALIDATION_CSV, lay.ANALOGUES_CSV,
                 "ligands_meta.csv")
BASE_FOLDERS = (lay.RECEPTORS, lay.INPUT_LIGANDS)
HEAVY_FOLDERS = ("poses", lay.COMPLEXES, "xml_plip", "xtal")

EXPORTS = {
    "results_csv":    ("Full ranking with all metrics", "file", "ranking.csv"),
    "summary_csv":       ("Compact summary per compound", "file", lay.SUMMARY_CSV),
    "interactions_csv": ("Interaction matrix per pose (PLIP)", "file", lay.INTERACTIONS_CSV),
    "docking_csv":       ("Energies of all poses", "file", lay.DOCKING_CSV),
    "ligands_csv":      ("Ligand table: name, SMILES, IUPAC and provenance", "file", "ligands_meta.csv"),
    "validation_csv":    ("Redocking validation of the control", "file", lay.VALIDATION_CSV),
    lay.RECEPTORS:        ("Prepared receptors and co-crystallized controls", "folder", lay.RECEPTORS),
    "ligands_zip":      ("3D structures of the ligands (SDF)", "folder", lay.INPUT_LIGANDS),
    "complexes_zip":     ("Receptor-ligand complexes (PDB)", "folder", lay.COMPLEXES),
    "poses_zip":         ("Docking poses per model", "folder", "poses"),
    "methods":           ("Methods section: parameters and versions", "generado", "PoliScreen_Methods.md"),
}

RECOMMENDED = {
    "methods":           "Not in the folder: written on export",
    "results_csv":    "Main table: score, Ki, efficiency and confidence",
    "interactions_csv": "Contact fingerprint that underpins the score",
    "ligands_csv":      "Provenance of each ligand: SMILES and IUPAC",
    "validation_csv":    "Redocking RMSD: the first thing a reviewer checks",
    lay.RECEPTORS:        "Exact inputs to repeat the run",
    "ligands_zip":      "Exact inputs to repeat the run",
}
REGENERABLE = ("summary_csv", "docking_csv", "complexes_zip", "poses_zip")


def _peso(proj: Path, tipo: str, origin: str) -> tuple:
    """(bytes, n_files) of a catalog element."""
    if tipo == "generado":
        return 0, 1
    p = lay.artifact(proj, origin)
    if tipo == "file":
        return (p.stat().st_size, 1) if p.is_file() else (0, 0)
    if not p.is_dir():
        return 0, 0
    fs = [f for f in p.rglob("*") if f.is_file()]
    return sum(f.stat().st_size for f in fs), len(fs)


def catalog(proj) -> dict:
    """Everything exportable plus what the interface needs to decide: whether it exists, how much it
    weighs and whether it is worth downloading. {key: {desc, hay, bytes, n, motivo, regenerable}}."""
    proj = Path(proj)
    out = {}
    for key_, (desc, tipo, origin) in EXPORTS.items():
        b, n = _peso(proj, tipo, origin)
        has_ = (proj / "run.json").exists() if tipo == "generado" else n > 0
        out[key_] = {"desc": desc, "has": has_, "bytes": b, "n": n,
                      "reason": RECOMMENDED.get(key_, ""), "regenerable": key_ in REGENERABLE}
    return out


def package_bytes(proj, keys_: Sequence[str], methods_text: Optional[str] = None) -> tuple:
    """In-memory ZIP with the chosen elements. Returns (bytes, list of what was included).

    It writes nothing to the project folder: the download is the user's decision, not the app's.
    """
    proj = Path(proj)
    buf, included_items = io.BytesIO(), []
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for key_ in keys_:
            if key_ not in EXPORTS:
                continue
            _desc, tipo, origin = EXPORTS[key_]
            if tipo == "generado":
                if methods_text:
                    zf.writestr(origin, methods_text)
                    included_items.append(origin)
            elif tipo == "file":
                if _add_file(zf, lay.artifact(proj, origin), f"tablas/{origin}"):
                    included_items.append(f"tablas/{origin}")
            else:
                n = _add_dir(zf, lay.artifact(proj, origin), origin)
                if n:
                    included_items.append(f"{origin}/ ({n} files)")
        if included_items:
            stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
            zf.writestr("CONTENTS.txt",
                        f"PoliScreen · proyecto '{proj.name}' · exportado {stamp}\n\n"
                        + "\n".join(f"- {i}" for i in included_items) + "\n")
    return buf.getvalue(), included_items


def _add_file(zf: zipfile.ZipFile, src: Path, arc: str) -> bool:
    if src.exists() and src.is_file():
        zf.write(src, arc)
        return True
    return False


def _add_dir(zf: zipfile.ZipFile, src: Path, arc: str) -> int:
    n = 0
    if not src.is_dir():
        return 0
    for p in sorted(src.rglob("*")):
        if p.is_file():
            zf.write(p, f"{arc}/{p.relative_to(src).as_posix()}")
            n += 1
    return n


UI_STATE = "ui_state.json"
LEGACY_UI_STATE = "estado_ui.json"

# A session is meant to be shared, so it must not carry the folder layout of the machine that
# produced it. Paths inside the project are stored relative to this marker and expanded again on
# restore, which also makes a session portable between machines; anything outside the project is
# reduced to its file name.
PROJECT_MARK = "{project}"


def _map_paths(obj, fn):
    if isinstance(obj, str):
        return fn(obj)
    if isinstance(obj, list):
        return [_map_paths(v, fn) for v in obj]
    if isinstance(obj, dict):
        return {_map_paths(k, fn): _map_paths(v, fn) for k, v in obj.items()}
    return obj


def strip_paths(obj, proj):
    """Local paths out of a state about to be written into a session."""
    root = str(Path(proj)).rstrip("/")

    def fn(s):
        if s.startswith(root + "/"):
            return PROJECT_MARK + s[len(root):]
        if s.startswith("/") and len(s.split("/")) > 2 and Path(s).suffix:
            return Path(s).name          # outside the project: the name is all that is needed
        return s
    return _map_paths(obj, fn)


def restore_paths(obj, proj):
    """Inverse of strip_paths, against the folder the session is restored into."""
    root = str(Path(proj)).rstrip("/")
    return _map_paths(obj, lambda s: root + s[len(PROJECT_MARK):]
                      if s.startswith(PROJECT_MARK) else s)

LEGACY_WIDGET_KEYS = {
    "etapa": "stage", "modo_ligandos": "ligand_mode", "ultimo_preparado": "last_prepared",
    "ultimo_original": "last_original", "cfg_alto": "cfg_height", "cfg_reparto": "cfg_split",
    "nuc_smiles": "core_smiles", "rx_reaccion": "rx_reaction", "pep_ciclo": "pep_cyclic",
    "pep_clases": "pep_classes", "pep_cons": "pep_consecutive", "pep_entrada": "pep_input",
    "pep_excl": "pep_exclude", "pep_nac": "pep_n_acetyl", "pep_pre": "pep_prefix",
    "pep_suf": "pep_suffix", "pep_sinrep": "pep_no_repeat", "pep_txt": "pep_text",
    "pep_usag": "pep_use_g", "pep_usaq": "pep_use_q", "vis_ver_cav": "vis_show_cav",
    "vis_ver_rec": "vis_show_rec", "vis_res_sitio": "vis_res_site",
    "vis_res_vista": "vis_res_view", "vis_ejes_box": "vis_axes_box",
    "vis_ejes_rec": "vis_axes_rec", "vis_cx_sup": "vis_cx_surface",
}


def _upgrade_widget_keys(widgets: dict) -> dict:
    """Translates the widget keys of a session written before the rename."""
    if not isinstance(widgets, dict):
        return widgets
    return {LEGACY_WIDGET_KEYS.get(k, k): v for k, v in widgets.items()}


def read_state(archivo) -> dict:
    """UI state saved in the session (built products, core, reagents).
    Returns {} if the session does not carry it."""
    for member in (UI_STATE, LEGACY_UI_STATE):
        try:
            with zipfile.ZipFile(Path(archivo)) as zf:
                state = json.loads(zf.read(member).decode("utf-8"))
        except Exception:
            continue
        if isinstance(state, dict) and "widgets" in state:
            state["widgets"] = _upgrade_widget_keys(state["widgets"])
        return state
    return {}


def save_session(proj, dest, full_: bool = False, notes: str = "", state_: Optional[dict] = None) -> Path:
    """Packs the project folder into a .poliscreen file.

    completa=True includes poses and complexes (heavy); otherwise only configuration, tables and inputs.
    estado: snapshot of the interface (products, core, reagents) so that on restore the series does
    not have to be rebuilt by hand.
    """
    proj = Path(proj)
    dest = Path(dest)
    if dest.suffix != EXT:
        dest = dest.with_suffix(EXT)
    dest.parent.mkdir(parents=True, exist_ok=True)

    included: list = []
    with zipfile.ZipFile(dest, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for name_ in BASE_FILES:
            if _add_file(zf, lay.artifact(proj, name_), name_):
                included.append(name_)
        for folder in BASE_FOLDERS:
            n = _add_dir(zf, lay.artifact(proj, folder), folder)
            if n:
                included.append(f"{folder}/ ({n} files)")
        if full_:
            for folder in HEAVY_FOLDERS:
                n = _add_dir(zf, lay.artifact(proj, folder), folder)
                if n:
                    included.append(f"{folder}/ ({n} files)")

        if state_:
            zf.writestr(UI_STATE, json.dumps(state_, indent=2, ensure_ascii=False, default=str))
            included.append(UI_STATE)

        try:
            from .. import __version__ as ver
        except Exception:
            ver = "unknown"
        manifest = {
            "format": FORMAT,
            "poliscreen": ver,
            "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "project": proj.name,
            "full": bool(full_),
            "contents": included,
            "notes": notes,
        }
        zf.writestr("manifest.json", json.dumps(manifest, indent=2, ensure_ascii=False))
    return dest


def session_info(archivo) -> dict:
    """Reads the manifest without extracting anything: used to show what the session carries before
    opening it."""
    with zipfile.ZipFile(Path(archivo)) as zf:
        try:
            return json.loads(zf.read("manifest.json").decode("utf-8"))
        except KeyError:
            return {"format": None, "contents": zf.namelist()[:50]}


def load_session(archivo, dest_dir) -> Path:
    """Extracts a session into dest_dir and returns the project folder ready to use.

    Rejects paths that escape the destination (protection against tampered archives).
    """
    archivo = Path(archivo)
    dest = Path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archivo) as zf:
        for member in zf.infolist():
            dest_ = (dest / member.filename).resolve()
            if not str(dest_).startswith(str(dest.resolve())):
                raise ValueError(f"Ruta insegura en la sesion: {member.filename}")
        zf.extractall(dest)
    return dest
