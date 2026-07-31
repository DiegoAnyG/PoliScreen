"""Interfaz de PoliScreen. Envuelve el nucleo; no contiene ciencia propia.

Lanzar:  poliscreen ui      (o: streamlit run .../ui/streamlit_app.py)
"""
from __future__ import annotations

import json
import shutil
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from poliscreen import __version__
from poliscreen.core import adcp
from poliscreen.core import docking as dk
from poliscreen.core import ligands as lig
from poliscreen.core import pipeline as pl
from poliscreen.core import peptides as pp
from poliscreen.core import pockets as pk
from poliscreen.core import reactions as rx
from poliscreen.core import reagents as rg
from poliscreen.core import report as rp
from poliscreen.core import receptor as rc
from poliscreen.core import screening as sc
from poliscreen.core import session as ss
from poliscreen.core import validation as vl
from poliscreen.core import viewer as vw
from poliscreen.core.design import AdmelabBridge
from poliscreen.ui import ayuda


def _shade(df, col, value="tuyo", color="rgba(255,205,60,0.20)"):
    """Resalta las filas cuya columna `col` vale `value` (para marcar lo que aporto el usuario)."""
    if col not in df.columns:
        return df
    return df.style.apply(lambda r: [f"background-color: {color}" if str(r.get(col)) == value else ""
                                     for _ in r], axis=1)


def _scatter_dock_inter(sub):
    """Dispersión docking vs. calidad de interacción. Muestra el compromiso: arriba-derecha = bueno en ambos."""
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
    ax.set_xlabel("Docking (kcal/mol; more negative = better)")
    ax.set_ylabel("Interaction quality (0-1 vs. control)")
    ax.invert_xaxis()
    ax.grid(alpha=0.3)
    ax.set_title("Docking vs. quality · red = control · ideal: top-right")
    fig.tight_layout()
    return fig


def _render_adme(admet, items, keyp):
    """items: [(etiqueta, smiles)]. Muestra tabla resumen de todos + detalle por compuesto."""
    filas = []
    for lb, smi in items:
        r = admet.get(rg.inchikey(smi)) or {}
        filas.append({"compound": lb, "MW": r.get("MW"), "LogP": r.get("LogP"), "QED": r.get("QED"),
                      "LD50 (mg/kg)": r.get("LD50_mg_per_kg"), "GHS": r.get("GHS_category"),
                      "AMES": r.get("AMES"), "hERG": r.get("hERG"), "DILI": r.get("DILI")})
    st.markdown("**ADMET summary of all compounds**")
    st.dataframe(pd.DataFrame(filas), width="stretch", height=min(320, 60 + 34 * len(filas)))
    st.caption("AMES/hERG/DILI = toxicity probability (lower is better). LD50 in mg/kg (higher is better). Predicted on the WHOLE molecule (core + reagent), not the reagent alone.")
    labels = dict(items)
    sel = st.selectbox("View detail of", list(labels), key=f"adme_det_{keyp}")
    row = admet.get(rg.inchikey(labels[sel]))
    if not row:
        return
    ca, cb = st.columns([1, 1])
    ca.pyplot(rp.radar_fig(row, title=sel))
    cb.metric("Oral LD50 (mg/kg)", rp._f(row.get("LD50_mg_per_kg"), 0))
    cb.metric("GHS category", str(row.get("GHS_category") or "-"))
    cb.metric("QED", rp._f(row.get("QED")))
    cb.caption("Green = favorable · amber = intermediate · red = unfavorable.")
    st.warning(rp.AVISO_LD50)
    _col = {"bueno": "background-color:rgba(46,158,126,0.22)",
            "medio": "background-color:rgba(226,168,44,0.22)",
            "malo": "background-color:rgba(214,70,70,0.22)", "info": ""}
    for titulo, fs in rp.sections(row):
        st.markdown(f"**{titulo}**")
        dd = pd.DataFrame(fs, columns=["Propiedad", "Valor", "v"])
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
    """Ki en unidades legibles. Es informativa: no interviene en el ranking."""
    if ki is None or pd.isna(ki):
        return ""
    if ki < 1e-9:
        return f"{ki * 1e12:.1f} pM"
    if ki < 1e-6:
        return f"{ki * 1e9:.1f} nM"
    if ki < 1e-3:
        return f"{ki * 1e6:.2f} uM"
    return f"{ki * 1e3:.2f} mM"

CITAS = f"""**PoliScreen v{__version__}** — Anaya Guerrero DC.

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

AGRADECIMIENTOS = """PoliScreen builds on these open-source projects:

- **Streamlit** (Apache-2.0) — web interface
- **3Dmol.js / py3Dmol** (BSD-3-Clause) — 3D molecular viewer
- **pandas** (BSD-3-Clause) and **NumPy** (BSD-3-Clause) — data handling
- **Matplotlib** (Matplotlib/PSF license) — interaction diagrams
- **OpenPyXL** (MIT) — XLSX export
- **OPSIN** (MIT) — IUPAC name verification

The scientific tools the engine runs (Vina, ADCP, gnina, PLIP, RDKit, Open Babel,
fpocket, OpenMM) have their full citation in «How to cite»."""


def _como_citar(expandido: bool = False):
    with st.expander("How to cite", expanded=expandido):
        st.markdown(CITAS)
    with st.expander("Acknowledgments"):
        st.markdown(AGRADECIMIENTOS)


def _descargar_tabla(df, nombre: str, key: str):
    """Ofrece descargar un DataFrame como CSV o XLSX bajo la tabla."""
    import io
    c = st.columns([2, 1, 1])
    c[0].caption("Download as:")
    c[1].download_button("CSV", df.to_csv(index=False).encode("utf-8"),
                         file_name=f"{nombre}.csv", mime="text/csv", key=f"dl_csv_{key}",
                         use_container_width=True)
    buf = io.BytesIO()
    try:
        with pd.ExcelWriter(buf, engine="openpyxl") as w:
            df.to_excel(w, index=False, sheet_name="PoliScreen")
        c[2].download_button("XLSX", buf.getvalue(), file_name=f"{nombre}.xlsx", key=f"dl_xlsx_{key}",
                             mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                             use_container_width=True)
    except Exception:
        pass


def _descargar_imagen(png: bytes, nombre: str, key: str):
    """Ofrece descargar una imagen (PNG) bajo el visor. jpeg no aporta aquí: son diagramas con
    fondo plano, PNG conserva nitidez sin artefactos."""
    if not png:
        return
    c = st.columns([3, 1])
    c[0].caption("Download as:")
    c[1].download_button("PNG", png, file_name=f"{nombre}.png", mime="image/png",
                         key=f"dl_png_{key}", use_container_width=True)


def _ya_hecho(clave: str, firma) -> bool:
    """True si la acción ya se ejecuto con ESTOS parámetros. Sirve para deshabilitar el boton
    hasta que algo cambie: evita repetir un cálculo largo por error y las pulsaciones dobles."""
    return S.get("_firma_" + clave) == firma


def _marcar_hecho(clave: str, firma):
    S["_firma_" + clave] = firma


def _vacio(mensaje: str):
    """Estado vacio del visualizador: el logotipo como marca de agua y una línea con lo que falta.
    Si no hay logotipo instalado en assets/ se usa un glifo monocromo de reserva."""
    logo = vw.logo_path()
    if logo and logo.suffix.lower() != ".svg":
        c = st.columns([1, 2, 1])[1]
        c.image(str(logo), use_container_width=True)
        c.markdown(f"<div style='text-align:center;opacity:.75;font-size:.92rem'>{mensaje}</div>",
                   unsafe_allow_html=True)
        return
    marca = logo.read_text(errors="ignore") if logo else vw.logo_svg()
    st.markdown(
        f"<div style='display:flex;flex-direction:column;align-items:center;justify-content:center;"
        f"padding:3rem 1rem;opacity:.75;text-align:center'>{marca}"
        f"<div style='margin-top:.9rem;font-size:.92rem;opacity:.8'>{mensaje}</div></div>",
        unsafe_allow_html=True)


@st.dialog("Confirm deletion")
def _confirmar_borrado(carpeta):
    """Ventana modal antes de una acción destructiva: borrar resultados no tiene vuelta atras."""
    st.warning("This deletes the poses, complexes, PLIP XML and all result tables in this folder. **This action cannot be undone.**")
    st.caption(f"Carpeta: `{carpeta}`")
    st.caption("Prepared receptors, controls and input ligands are kept.")
    st.info("To keep this analysis, cancel and use File -> Save session first.")
    c1, c2 = st.columns(2)
    if c1.button("Yes, delete", type="primary", use_container_width=True):
        pl.clean(carpeta)
        st.session_state["_aviso"] = "Results deleted. Receptors, controls and ligands are kept."
        st.rerun()
    if c2.button("Cancel", use_container_width=True):
        st.rerun()


def _tam(b: int) -> str:
    return f"{b / 1e6:.1f} MB" if b >= 1e6 else (f"{b / 1e3:.0f} kB" if b else "—")


@st.dialog("Download results", width="large")
def _dialogo_descargas(carpeta):
    """Selector de exportaciones. Arma un único ZIP en memoria con lo que el usuario marque.

    La carpeta del proyecto ya contiene casi todos estos archivos; el paquete existe para llevarse
    el análisis a otra maquina, adjuntarlo a un articulo o archivarlo sin los intermedios pesados.
    Por eso cada elemento dice si vale la pena bajarlo y cuanto ocupa, en vez de ofrecerlos a ciegas.
    """
    cat = ss.catalogo(carpeta)
    hay = {k: v for k, v in cat.items() if v["hay"]}
    if not hay:
        st.info("Nothing to export in this folder yet.")
        return
    st.caption(f"Project folder: `{carpeta}`")
    st.caption("Everything here already lives in that folder, except the Methods section, written on export. The package is for moving the analysis elsewhere or attaching it to a manuscript without the heavy intermediates.")

    c1, c2 = st.columns(2)
    marcar = c1.button("Select the recommended", use_container_width=True)
    limpiar = c2.button("Clear all", use_container_width=True)
    for k, v in hay.items(): 
        if marcar or limpiar:
            S[f"dl_{k}"] = bool(marcar and v["motivo"])
        else:
            S.setdefault(f"dl_{k}", bool(v["motivo"]))

    rec = {k: v for k, v in hay.items() if v["motivo"]}
    otros = {k: v for k, v in hay.items() if not v["motivo"]}
    if rec:
        st.markdown("**Recommended**")
        for k, v in rec.items():
            st.checkbox(f"{v['desc']} · {_tam(v['bytes'])}", key=f"dl_{k}")
            st.caption(v["motivo"])
    if otros:
        st.markdown("**Optional**")
        for k, v in otros.items():
            st.checkbox(f"{v['desc']} · {_tam(v['bytes'])}", key=f"dl_{k}")
            if v["regenerable"]:
                st.caption("Regenerated by running again; takes space in the package.")

    elegidas = [k for k in hay if S.get(f"dl_{k}")]
    total = sum(hay[k]["bytes"] for k in elegidas)
    st.divider()
    if not elegidas:
        st.caption("You have not selected anything.")
        return
    st.caption(f"{len(elegidas)} item(s) · {_tam(total)} uncompressed")
    if st.button("Prepare package", type="primary", use_container_width=True):
        try:
            meta = carpeta / "run.json"
            mtxt = rp.methods_text(json.loads(meta.read_text())) if meta.exists() else None
            datos, incluidos = ss.paquete(carpeta, elegidas, methods_text=mtxt)
            S["_zip"] = (f"{Path(carpeta).name}_PoliScreen.zip", datos, incluidos)
        except Exception as e:
            st.error(f"Could not build the package: {e}")
    if S.get("_zip"):
        nombre, datos, incluidos = S["_zip"]
        st.download_button(f"Download {nombre} ({_tam(len(datos))})", datos,
                           file_name=nombre, mime="application/zip",
                           type="primary", use_container_width=True)
        st.caption("Contains: " + ", ".join(incluidos))


st.set_page_config(page_title="PoliScreen", layout="wide", initial_sidebar_state="collapsed")

st.markdown(
    """
    <style>
      /* La cabecera de Streamlit es una barra FIJA que cubre todo el ancho superior. Aunque se
         vea vacia, interceptaba los clics de la barra de menu. Se hace transparente a los
         eventos y solo los botones reales (el de tres puntos) vuelven a recibirlos; ademas el
         contenido se baja lo justo para no solaparse con ella. */
      [data-testid="stAppDeployButton"] { display: none; }
      header[data-testid="stHeader"] { background: transparent; pointer-events: none; }
      header[data-testid="stHeader"] button { pointer-events: auto; }
      .block-container { padding-top: 2.4rem; padding-bottom: 0.3rem;
                         padding-left: 1.6rem; padding-right: 1.6rem; }
      /* Adornos de Streamlit que no aportan: el boton de pantalla completa sobre cada imagen
         y el ancla que aparece junto a los titulos y solo enlaza a la propia pagina. */
      [data-testid="StyledFullScreenButton"] { display: none; }
      [data-testid="stHeaderActionElements"] { display: none; }
      h1 a.anchor-link, h2 a.anchor-link, h3 a.anchor-link,
      h4 a.anchor-link, h5 a.anchor-link { display: none; }
      /* Cabecera: linea tenue de separacion */
      .st-key-barra_menu { border-bottom: 1px solid rgba(128,128,128,.28); padding-bottom: .35rem; }
      /* Barra inferior de etapas: solo texto, plano y centrado */
      .st-key-barra_etapas { border-top: 1px solid rgba(128,128,128,.28); padding-top: .35rem; }
      .st-key-barra_etapas button { border: none !important; background: transparent !important;
                                    font-weight: 600; letter-spacing: .01em; }
      .st-key-barra_etapas button:hover { background: rgba(128,128,128,.14) !important; }
      /* Etapa activa: color de acento y subrayado. Un fondo gris tenue no se distinguia. */
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

ETAPAS = ["Receptors", "Ligands", "Run", "Results"]
S.setdefault("etapa", ETAPAS[0])

# Streamlit prohíbe modificar el estado de un widget Después de instanciarlo. Para cambiar la
# carpeta al restaurar una sesión se deja aquí la peticion y se aplica al principio de la pasada
# siguiente, antes de que exista el campo de texto.
if S.get("_proj_pendiente"):
    S["proj_dir"] = S.pop("_proj_pendiente")
    S["_proj_cargado"] = None       # fuerza el auto-descubrimiento de la carpeta nueva
# Valores de widgets restaurados de una sesión: se aplican antes de crear los widgets.
for _k, _v in (S.pop("_widgets_pendientes", None) or {}).items():
    S[_k] = _v

# Streamlit descarta del estado las claves de los widgets que no se dibujan en una pasada. Como
# solo se ejecuta la etapa activa, al cambiar de etapa se perderian el modo de ligandos y los
# parámetros de generación. Reasignarlas al principio de cada pasada las conserva.
_PREFIJOS_PERSISTENTES = ("pep_", "modo_", "cat_", "sec_", "rec_", "box_", "sitios_", "rx_",
                          "cx_", "cy_", "cz_", "sx_", "sy_", "sz_", "src_", "vis_", "cfg_")
for _k in [k for k in S.keys() if isinstance(k, str) and k.startswith(_PREFIJOS_PERSISTENTES)]:
    S[_k] = S[_k]

# --- Barra de menu superior. Los popover se comportan como menus desplegables clasicos; su
# contenido se ejecuta en cada pasada, así que las variables que definen (proj, seed...) quedan
# disponibles para el resto de la aplicación.
_barra = st.container(key="barra_menu")
_logo_f = vw.logo_path()
_wm = vw.wordmark_path()
# Marca a la izquierda y menus a la derecha. El icono, el nombre y la versión van en un único bloque
# flexbox con las imagenes incrustadas: es la única forma de mantenerlos pegados y centrados a la
# misma altura. Repartirlos en columnas los separaba y parecía que uno estaba más arriba que otro.
_izq, _hueco, _m1, _m2, _m3, _m4 = _barra.columns([3.6, 1.4, 1.1, 1.1, 1.6, 1.1],
                                                  vertical_alignment="center")


def _img_inline(ruta, alto, clase="", estilo=""):
    import base64
    mime = "svg+xml" if Path(ruta).suffix.lower() == ".svg" else Path(ruta).suffix.lower().lstrip(".")
    b64 = base64.b64encode(Path(ruta).read_bytes()).decode()
    return f"<img class='{clase}' src='data:image/{mime};base64,{b64}' style='height:{alto}px;{estilo}'>"


# Icono y logotipo son oscuros sobre fondo claro; en tema oscuro se invierten para contrastar. Con el
# tema resuelto en el servidor el filtro va en linea sobre la imagen, no por clase: Streamlit puede
# descartar el atributo class al sanear el HTML, pero conserva el style. Si el tema sigue al sistema
# la deteccion queda al navegador por prefers-color-scheme, unico caso donde hace falta la clase.
_INV = "filter:invert(1) hue-rotate(180deg) brightness(1.1);"
try:
    _tema = st.context.theme.type       # 'dark' | 'light' | None (sistema/indeterminado)
except Exception:
    _tema = None
_est_marca = _INV if _tema == "dark" else ""
_cls_img = "marca-sistema" if _tema is None else ""
_col_ver = "#e6e6e6" if _tema == "dark" else "#666"
st.markdown(f"""<style>
  @media (prefers-color-scheme: dark) {{ .marca-sistema {{ {_INV} }} }}
</style>""", unsafe_allow_html=True)

_marca = []
if _logo_f and _logo_f.suffix.lower() != ".svg":
    _marca.append(_img_inline(_logo_f, 72, _cls_img, _est_marca))
_marca.append(_img_inline(_wm, 56, _cls_img, _est_marca) if _wm
              else "<span style='font-size:2.1rem;font-weight:700'>PoliScreen</span>")
_marca.append(f"<span style='color:{_col_ver};font-size:.85rem;align-self:flex-end;"
              f"padding-bottom:.5rem'>v{__version__}</span>")
# margin-top negativo para subirlo y no tocar el filete inferior de la barra; gap generoso.
_izq.markdown("<div style='display:flex;align-items:center;gap:.7rem;margin:-.3rem 0 .4rem'>"
              + "".join(_marca) + "</div>", unsafe_allow_html=True)
_menu_archivo = _m1.popover("File", use_container_width=True)
_menu_datos = _m2.popover("Data", use_container_width=True)
_menu_cfg = _m3.popover("Settings", use_container_width=True)
_menu_ayuda = _m4.popover("Help", use_container_width=True)

with _menu_archivo:
    st.markdown("**Project**")
    S.setdefault("proj_dir", str(Path.home() / "poliscreen_proyectos" / "demo"))
    _escrito = st.text_input("Project folder", key="proj_dir",
                             help="Path inside Linux (WSL). If you paste a Windows path (`\\\\wsl.localhost\\...` or `C:\\...`) it is translated automatically.")
    proj, _aviso_ruta = ss.normalizar_ruta(_escrito)
    if _aviso_ruta:
        # Se corrige el propio campo en la pasada siguiente para que lo guardado en la sesión y lo
        # que el usuario ve coincidan con la carpeta que se usa de verdad.
        st.warning(_aviso_ruta)
        S["_proj_pendiente"] = str(proj)
    try:
        proj.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        st.error(f"Cannot create that folder: {e}")
        st.stop()
    st.caption(f"Results in `{proj}`")
    # Auto-descubrir lo YA preparado en disco: sobrevive a reinicios de Streamlit y a cambiar de carpeta.
    # Sin esto, la lista de receptores/controles/ligandos (session_state) se pierde al reiniciar y el
    # docking corre SIN el control aunque el archivo siga en la carpeta.
    if S.get("_proj_cargado") != str(proj):
        S["_proj_cargado"] = str(proj)
        # Se parte de cero y se rellena con lo que haya en la carpeta. Antes solo se reasignaba
        # cuando el subdirectorio existía, de modo que al pasar a un proyecto nuevo persistian los
        # receptores, controles y ligandos del anterior: una tanda podia arrastrar el control de
        # otra diana sin que nada lo indicara.
        rec_dir, lig_dir = proj / "receptores", proj / "ligandos_entrada"
        S["receptors"] = (sorted(str(p) for p in rec_dir.glob("*_listo.pdb"))
                          if rec_dir.exists() else [])
        # Un control puede ser un PDB cuando es una cadena peptídica extraida del cristal; se exige
        # el prefijo para no confundirlo con el receptor, que vive en la misma carpeta.
        S["controls"] = (sorted(str(p) for p in rec_dir.iterdir()
                                if p.suffix.lower() in (".sdf", ".mol2", ".mol")
                                or (p.suffix.lower() == ".pdb" and p.name.startswith("control_")))
                         if rec_dir.exists() else [])
        S["ligands"] = (sorted(str(p) for p in lig_dir.iterdir()
                               if p.suffix.lower() in (".sdf", ".mol2", ".mol", ".smi", ".pdbqt"))
                        if lig_dir.exists() else [])
    if S["receptors"] or S["controls"] or S["ligands"]:
        st.caption(f"Recovered from disk: {len(S['receptors'])} receptor(s), "
                   f"{len(S['controls'])} control(s), {len(S['ligands'])} ligand(s).")

    # --- Sesión de trabajo: abrir/guardar sin tener que escribir rutas ---
    with st.expander("Session and export"):
        sub = st.file_uploader("Open session (.poliscreen)", type=["poliscreen", "zip"],
                               help="Restores a previous analysis: tables, receptors and ligands. You can change the weighting without repeating the docking.")
        if sub is not None and st.button("Restore this session"):
            destino = Path.home() / "poliscreen_proyectos" / Path(sub.name).stem
            tmp_s = destino.parent / f"_{Path(sub.name).stem}.poliscreen"
            tmp_s.parent.mkdir(parents=True, exist_ok=True)
            tmp_s.write_bytes(sub.getvalue())
            try:
                ss.load_session(tmp_s, destino)
                info = ss.session_info(tmp_s)
                est = ss.leer_estado(tmp_s)          # productos, nucleo, reactivos, cavidades, registro
                for clave in ("products", "nuc_smiles", "reagents", "reag_info", "pockets",
                              "pep_seqs", "pep_aviso", "_log_run", "_log_estado"):
                    if est.get(clave) is not None:
                        S[clave] = est[clave]
                # Los valores de los widgets se aplican en la pasada siguiente, cuando aún no existen.
                S["_widgets_pendientes"] = est.get("widgets") or {}
                # No se puede tocar S["proj_dir"] aquí (el widget ya existe): se deja pendiente.
                S["_proj_pendiente"] = str(destino)
                st.success(f"Session '{info.get('proyecto', '?')}' restored in {destino}")
                st.rerun()
            except Exception as e:
                st.error(f"Could not restore the session: {e}")
            finally:
                tmp_s.unlink(missing_ok=True)

        completa = st.checkbox("Include poses and complexes (heavy session)", value=False,
                               help="Unchecked, the session is a few MB and enough to reopen and re-score. Checked, it also lets you re-examine 3D structures.")
        if st.button("Save session"):
            try:
                # La sesión replica el último estado del usuario. Los pesos que tenga en Resultados se
                # escriben en run.json antes de empaquetar, de modo que al restaurar los sliders
                # arranquen con ellos (leen de ahí) y no haya que reponderar.
                _rj = proj / "run.json"
                if S.get("_pesos_ui") and _rj.exists():
                    try:
                        _d = json.loads(_rj.read_text())
                        _d["weights"] = S["_pesos_ui"]
                        _rj.write_text(json.dumps(_d, indent=2))
                    except Exception:
                        pass
                # Se incluyen las cavidades, el registro de la corrida y la selección de la interfaz
                # (residuos catalíticos y secundarios entre otros): recalcular fpocket o reponderar
                # tras restaurar cuesta tiempo, y la caja quedaría sin su origen.
                claves_ui = ("products", "nuc_smiles", "reagents", "reag_info", "pockets",
                             "pep_seqs", "pep_aviso", "_log_run", "_log_estado")
                estado_ui = {k: S.get(k) for k in claves_ui if S.get(k) is not None}
                estado_ui["widgets"] = {k: S[k] for k in S.keys()
                                        if isinstance(k, str) and k.startswith(_PREFIJOS_PERSISTENTES)
                                        and isinstance(S[k], (str, int, float, bool, list))}
                # Se arma en una carpeta temporal y se conserva solo en memoria: la sesión es una
                # descarga, no un archivo más que ensucie la carpeta de resultados.
                _tmp = Path(tempfile.mkdtemp())
                sfile = ss.save_session(proj, _tmp / f"{proj.name}", completa=completa,
                                        estado=estado_ui)
                S["_sesion"] = (sfile.name, sfile.read_bytes())
                st.success(f"{sfile.name} · {sfile.stat().st_size / 1e6:.1f} MB")
                shutil.rmtree(_tmp, ignore_errors=True)
            except Exception as e:
                st.error(f"Could not save: {e}")
        if S.get("_sesion"):
            st.download_button(f"Download {S['_sesion'][0]}", S["_sesion"][1],
                               file_name=S["_sesion"][0], mime="application/zip")

        _n_exp = sum(1 for v in ss.catalogo(proj).values() if v["hay"])
        if _n_exp:
            # El modal se abre desde el flujo principal, no desde dentro del popover: un dialogo
            # invocado dentro de un contenedor desplegable se dibuja dentro de el y queda recortado.
            if st.button(f"Download results... ({_n_exp} available)",
                         use_container_width=True):
                S["_abrir_descargas"] = True
                st.rerun()
            st.caption("You pick what to include with checkboxes and a single ZIP is built; nothing is written to the project folder.")
with _menu_datos:
    st.markdown("**Loaded in this project**")
    st.write(f"Prepared receptors: **{len(S['receptors'])}**")
    st.write(f"Co-crystallized controls: **{len(S['controls'])}**")
    st.write(f"Ligands: **{len(S['ligands'])}**")
    st.caption("Auto-detected from the project folder; they survive restarting the app.")

with _menu_ayuda:
    st.markdown("**PoliScreen manual**")
    st.caption("Each section expands with the full detail.")
    for _sec, _temas in ayuda.SECCIONES.items():
        with st.expander(_sec):
            for _titulo, _cuerpo in _temas:
                st.markdown(f"**{_titulo}**")
                st.markdown(_cuerpo)
                st.markdown("")
    _como_citar()

with _menu_cfg:
    st.markdown("**Appearance**")
    st.caption("The light/dark theme is set in the top-right menu (⋮) -> Settings -> Theme. Leave it on \"Use system setting\" to follow the system.")
    st.slider("Split between tools and viewer", 0.3, 0.7, 0.46, 0.02,
              key="cfg_reparto",
              help="Left gives more space to the viewer; right, to the tools.")
    st.slider("Panel height (px)", 380, 900, 520, 20, key="cfg_alto")
    st.caption("Docking parameters are in step 3 (Run), next to the launch button.")

if S.pop("_abrir_descargas", False):
    _dialogo_descargas(proj)


def _quimiotipos_tanda():
    """(hay_vina, hay_peptidos) presentes en la tanda. Solo decide QUE ajustes mostrar; el enrutado
    real lo hace el pipeline por ligando. Errar aquí como mucho muestra una sección de ajustes de
    más, nunca envia un ligando al motor equivocado."""
    peps = S.get("modo_ligandos") == "Generate peptides" or bool(S.get("pep_seqs"))
    vina = False
    for c in S.get("controls", []):
        cl = str(c).lower()
        if cl.endswith(".pdb") and pp.secuencia_de_estructura(c):
            peps = True
        elif cl.endswith((".sdf", ".mol2", ".mol")):
            vina = True
    if S.get("products") or (S.get("modo_ligandos") in ("Build by reaction", "Upload ready ligands")
                             and S.get("ligands")):
        vina = True
    if not peps and not vina:
        vina = True                      # sin señal clara, por defecto moléculas pequeñas
    return vina, peps


def _params_docking():
    """Parámetros de acoplamiento. Viven en la etapa Ejecutar, que es donde se usan.

    Solo se muestran los ajustes del motor o motores que van a intervenir: Vina para moléculas
    pequeñas, ADCP para péptidos, ambos en una tanda mixta. Un panel con ajustes que no hacen nada
    confunde sobre que motor se esta usando. Los valores de las secciones ocultas se devuelven en
    su valor por defecto, de modo que RunConfig siempre recibe una configuración completa."""
    hay_vina, hay_peps = _quimiotipos_tanda()
    adcp_ok = adcp.available()
    usa_adcp = hay_peps and adcp_ok
    usa_vina = hay_vina or (hay_peps and not adcp_ok)      # sin ADCP, los péptidos caen a Vina

    exhaust, energy_range, ph, cpu, workers = 8, 3.0, 7.4, 1, 0
    adcp_pasos, adcp_reps = 250_000, 20

    with st.expander("Advanced docking settings"):
        st.caption("The defaults are fine for a first exploration. Raise them for a definitive screen.")
        if usa_vina and usa_adcp:
            st.info("**Mixed** run: small molecules dock with **Vina** and peptides with **ADCP**. Each engine's settings appear separately below.")
        elif usa_adcp:
            st.info("**Peptide** screening with **ADCP**. Only its settings are shown; Vina's do not apply.")
        elif hay_peps and not adcp_ok:
            st.warning("There are peptides and ADCP is not installed: they will dock with **Vina**, whose sampling does not cover that flexibility. Install it with scripts/get_adcp.sh.")

        seed = st.number_input("Seed", value=42, step=1,
                               help="Fixes the randomness: the same seed gives the same result in both engines.")
        n_poses = st.slider("Poses per ligand", 1, 20, 5,
                            help="Below 3 the confidence metric loses resolution.")

        if usa_vina:
            st.markdown("**Vina** — small molecules")
            exhaust = st.slider("Exhaustiveness", 8, 64, 8, 8,
                                help="Higher = finer and slower search.")
            energy_range = st.slider("Energy range (kcal/mol)", 1.0, 8.0, 3.0, 0.5,
                                     help="Energy window relative to the best pose for reporting alternative modes.")
            ph = st.slider("Protonation pH", 5.0, 9.0, 7.4, 0.1,
                           help="pH at which OpenBabel protonates before docking (physiological ≈ 7.4).")
            cpu = st.number_input("Threads per docking", 1, 16, 1,
                                  help="1 keeps the result reproducible. Raise it only if you do not mind.")
            workers = st.number_input("Dockings in parallel (0 = automatic)", 0, 32, 0)
            if cpu > 1:
                st.warning("With more than one thread per docking, Vina stops being deterministic.")

        if usa_adcp:
            st.markdown("**ADCP** — peptides")
            st.caption("It uses the machine's cores automatically and is reproducible with the seed; Vina's thread settings do not affect it.")
            adcp_pasos = st.select_slider(
                "Steps per replica", [50_000, 100_000, 250_000, 500_000, 1_000_000, 2_500_000],
                value=250_000, format_func=lambda v: f"{v // 1000} k",
                help="Length of each search. Raise it if the control does not recover its pose or if the energy keeps improving as you increase it.")
            adcp_reps = st.slider("Independent replicas", 4, 100, 20, 2,
                                  help="Parallel searches from different starting points. More replicas lower the chance of getting stuck in a local minimum.")
            st.caption(f"Approximate cost: {adcp_pasos * adcp_reps / 5e6:.0f}× that of an octapeptide "
                       "with the default values (~35 s on six threads).")

        st.markdown("**Second opinion (neural network)**")
        hay_gnina = dk.gnina_available()
        # La red de gnina se entreno con complejos proteína-molécula pequeña. Con péptidos sigue
        # ejecutandose y devuelve numeros, que es justo el peligro: son extrapolacion fuera de su
        # dominio de entrenamiento. Se avisa en vez de bloquear.
        rescnn = st.checkbox("Re-score the poses with gnina (CNN, GPU)", value=False,
                             disabled=not hay_gnina,
                             help="It does not re-dock: it keeps the poses and evaluates them with a neural network trained on crystallographic complexes. Adds independent evidence to the confidence metric.")
        if not hay_gnina:
            st.caption("gnina is not installed. It is optional: without it, confidence uses the other evidence.")
        elif rescnn and hay_peps:
            st.warning("You are screening peptides. gnina's network was trained on small-molecule complexes, so here it scores outside its domain: low values do not necessarily mean the pose is bad. Use it to compare, not as a criterion, and declare it in Methods.")
        elif rescnn:
            st.caption("The best pose of each compound is re-scored (~2 s per compound).")
    return dict(seed=int(seed), exhaustiveness=int(exhaust), n_poses=int(n_poses),
                energy_range=float(energy_range), ph=float(ph), cpu=int(cpu), workers=int(workers),
                rescoring_cnn=bool(rescnn), adcp_pasos=int(adcp_pasos),
                adcp_replicas=int(adcp_reps))


# Cada etapa es una función: solo se ejecuta la activa, no las cuatro como ocurria con st.tabs.
# ---------------------------------------------------------------- 1. receptor
def _etapa_receptores():
    st.subheader("Prepare a receptor")
    st.caption("Type a PDB identifier or upload your own file. Waters are removed, hydrogens added, and the original residue numbering is kept.")
    c1, c2 = st.columns([1, 2])
    pdb_id = c1.text_input("PDB identifier", placeholder="4D44")
    up = c2.file_uploader("...or upload a .pdb file", type=["pdb"])

    src = None
    if up is not None:
        src = proj / "receptores" / up.name
        src.parent.mkdir(parents=True, exist_ok=True)
        src.write_bytes(up.getvalue())
    elif pdb_id.strip() and st.button("Download from the PDB"):
        try:
            src = rc.fetch_pdb(pdb_id, proj / "receptores")
            S["src_pdb"] = str(src)
        except rc.ReceptorError as e:
            st.error(str(e))
    if src is None and S.get("src_pdb"):
        src = Path(S["src_pdb"])

    if src and src.exists():
        st.success(f"Structure loaded: {src.name}")
        info = rc.inspect(src)
        st.write(f"**{info.n_atoms}** atoms · chains **{', '.join(info.chains)}** · **{info.n_waters}** waters")
        if info.het:
            st.dataframe(pd.DataFrame([{"group": h.resname, "chain": h.chain, "number": h.resseq,
                                        "atoms": h.n_atoms, "key": h.key} for h in info.het]),
                         width="stretch", height=220)
        # Claves explicitas ligadas al archivo: sin ellas la selección de cadenas, cofactores y
        # control se pierde al cambiar de etapa, porque Streamlit descarta los widgets no dibujados.
        kb = sc.normalize_key(src.stem)
        c1, c2, c3 = st.columns(3)
        chains = c1.multiselect("Chains to keep", info.chains, default=info.chains[:1],
                                key=f"rec_chains_{kb}")
        # Un residuo modificado (fosfotirosina, selenometionina...) el PDB lo declara como
        # heteroatomo, pero no es un cofactor: pertenece a la cadena. Ofrecerlo entre los cofactores
        # invitaba a conservarlo por esa vía, que lo duplicaba sobre el que ya esta en la proteína.
        mods = rc.modified_residues(src)
        claves_mod = {(m.chain, m.resseq) for m in mods}
        keys = [h.key for h in info.het if (h.chain, str(h.resseq).strip()) not in claves_mod]
        keep = c2.multiselect("Keep (cofactors)", keys, key=f"rec_keep_{kb}",
                              help="A site cofactor, e.g. NADP.")
        # Un solo selector de control que admite las dos formas en que un ligando de referencia
        # puede estar en el cristal: como heterogrupo (una molécula pequeña) o como cadena (un
        # péptido). Las cadenas se marcan con un prefijo para distinguirlas de una clave de
        # heterogrupo, y se ofrecen solo las que no se conservan ya como receptor.
        _cad_libres = [c for c in info.chains if c not in chains]
        _opc_ctrl = keys + [f"chain:{c}" for c in _cad_libres]

        def _fmt_ctrl(o):
            if o.startswith("chain:"):
                c = o.split(":", 1)[1]
                n = sum(1 for _l in src.read_text(errors="ignore").splitlines()
                        if _l.startswith("ATOM") and _l[21] == c and _l[12:16].strip() == "CA")
                return f"Chain {c} · {n} residues (peptide)"
            return f"{o} (hetero group)"

        _sel_ctrl = c3.multiselect("Extract as control", _opc_ctrl, key=f"rec_extract_{kb}",
                                   format_func=_fmt_ctrl,
                                   help="The co-crystallized ligand that defines the reference fingerprint. It can be a hetero group or a peptide chain; both appear here.")
        extract = [o for o in _sel_ctrl if not o.startswith("chain:")]
        cad_ctrl = [o.split(":", 1)[1] for o in _sel_ctrl if o.startswith("chain:")]
        smiles = st.text_input("SMILES of the extracted ligand (optional)", key=f"rec_smiles_{kb}",
                               help="Fixes bond orders, which the PDB does not store.")
        keep_mod = []
        if mods:
            st.markdown("**Modified residues of the chain**")
            st.caption("Detected in the structure. Checked, they are kept with their modification; unchecked, they are replaced by the amino acid they derive from and **the modification is lost** — which is often the function, as in a phosphorylated activation loop.")
            keep_mod = st.multiselect("Keep with its modification",
                                      [m.key for m in mods],
                                      default=[m.key for m in mods],
                                      key=f"rec_mod_{kb}",
                                      format_func=lambda k: next(
                                          (m.label for m in mods if m.key == k), k))
        for _c in cad_ctrl:
            _n = sum(1 for _l in src.read_text(errors="ignore").splitlines()
                     if _l.startswith("ATOM") and _l[21] == _c and _l[12:16].strip() == "CA")
            if _n > pp.MAX_LARGO:
                st.warning(f"Chain {_c} has {_n} residues: too long to treat it "
                           "as a reference ligand.")
        firma_prep = (str(src), tuple(chains), tuple(keep), tuple(extract), tuple(cad_ctrl),
                      tuple(keep_mod), smiles)
        prep_hecho = _ya_hecho("prep_" + kb, firma_prep)
        if prep_hecho:
            st.caption("Receptor already prepared with this selection. Change something to prepare it again.")
        if st.button("Prepare receptor", type="primary", disabled=prep_hecho):
            with st.spinner("Preparing..."):
                dest = proj / "receptores" / f"{src.stem}_listo.pdb"
                rc.prepare(src, dest, keep_chains=chains or None, keep_het=keep, ph=7.4,
                           keep_modified=keep_mod, on_aviso=st.warning)
                if str(dest) not in S["receptors"]:
                    S["receptors"].append(str(dest))
                S["ultimo_preparado"] = str(dest)
                S["ultimo_original"] = str(src)
                for k in extract:
                    het = info.find(k)
                    p = rc.extract_ligand(src, het, proj / "receptores" / f"control_{het.resname}.sdf",
                                          smiles=smiles or None)
                    if str(p) not in S["controls"]:
                        S["controls"].append(str(p))
                for c in cad_ctrl:
                    p = rc.extract_chain(src, c, proj / "receptores" / f"control_cadena{c}.pdb",
                                         on_aviso=st.warning)
                    if str(p) not in S["controls"]:
                        S["controls"].append(str(p))
                    _seq = pp.secuencia_de_estructura(p)
                    if _seq:
                        st.info(f"Chain {c} extracted as control: `{_seq[0]}` "
                                f"({len(_seq[0])} residues)"
                                + (", it will be docked with ADCP." if adcp.available()
                                   and adcp.MIN_RESIDUOS <= len(_seq[0]) <= adcp.MAX_RESIDUOS
                                   else "."))
            _marcar_hecho("prep_" + kb, firma_prep)
            st.success(f"Done: {dest.name}")
            st.rerun()

    # Comprobación visual: que se fue y que se quedo
    if S.get("ultimo_preparado") and Path(S["ultimo_preparado"]).exists():
        st.markdown("---")
        st.subheader("Preparation check")
        antes = vw.resumen_estructura(S["ultimo_original"])
        despues = vw.resumen_estructura(S["ultimo_preparado"])
        # Todo como texto: la columna mezcla recuentos y listas de cadenas, y Arrow no puede
        # serializar una columna con enteros y cadenas a la vez.
        comp = pd.DataFrame([
            {"": "Atoms", "before": str(antes["atomos"]), "after": str(despues["atomos"])},
            {"": "Hydrogens", "before": str(antes["hidrogenos"]), "after": str(despues["hidrogenos"])},
            {"": "Waters", "before": str(antes["aguas"]), "after": str(despues["aguas"])},
            {"": "Chains", "before": ", ".join(antes["cadenas"]), "after": ", ".join(despues["cadenas"])},
            {"": "Hetero groups", "before": ", ".join(sorted(antes["heterogrupos"])) or "-",
             "after": ", ".join(sorted(despues["heterogrupos"])) or "-"},
        ])
        st.dataframe(comp, width="stretch", hide_index=True)
        if despues["aguas"] == 0 and despues["hidrogenos"] > 0:
            st.success("No waters and hydrogens added.")
        else:
            st.warning("Check: there should be 0 waters and hydrogens present.")

        st.caption("The structure is shown in the right panel; there you can change the view and style.")

    if S["receptors"]:
        st.write("**Receptores preparados:**", ", ".join(Path(p).name for p in S["receptors"]))
    if S["controls"]:
        st.write("**Controles:**", ", ".join(Path(p).name for p in S["controls"]))

# ---------------------------------------------------------------- 2b. péptidos
def _modo_peptidos():
    """Diseño de péptidos: vía independiente de la síntesis por reacción. Se mantiene aparte
    (S['péptidos']) para que no se mezcle con los productos del constructor químico."""
    st.caption("Peptides undergo no chemical reactions: they are built directly from the sequence. Between 1 and 20 residues.")
    entrada = st.radio("How to obtain the sequences", ["Generate library", "Write sequences"],
                       horizontal=True, key="pep_entrada")

    secuencias, aviso, problemas = [], "", []
    if entrada == "Write sequences":
        txt = st.text_area("One sequence per line, in one-letter code",
                           placeholder="KWKLFKKI\nGIGKFLHSAK\nRRWWRF", height=130, key="pep_txt")
        crudas = [s.strip().upper() for s in txt.splitlines() if s.strip()]
        malas = []
        for s in crudas:
            fuera = set(s) - set(pp.AMINOACIDOS)
            if fuera:
                malas.append(f"{s} (invalid symbols: {', '.join(sorted(fuera))})")
            elif not (pp.MIN_LONGITUD <= len(s) <= pp.MAX_LONGITUD):
                malas.append(f"{s} (length {len(s)}; the maximum is {pp.MAX_LONGITUD})")
            else:
                secuencias.append(s)
        if malas:
            st.warning("These lines are ignored: " + " · ".join(malas[:6]))
    else:
        c1, c2, c3 = st.columns(3)
        largo = c1.number_input("Residues per peptide", pp.MIN_LONGITUD, pp.MAX_LONGITUD, 7, key="pep_len")
        cuantos = c2.number_input("How many peptides", 1, 2000, 50, key="pep_n")
        semilla = c3.number_input("Seed", value=42, step=1, key="pep_seed",
                                  help="Same seed and same rules = same library.")
        with st.expander("Composition: which amino acids it may use", expanded=True):
            clases = st.multiselect("Allowed classes (empty = all 20)", list(pp.CLASES),
                                    format_func=lambda k: pp.CLASES[k], key="pep_clases")
            excl = st.multiselect("Exclude specific residues", sorted(pp.AMINOACIDOS),
                                  format_func=lambda a: f"{a} · {pp.AMINOACIDOS[a][0]}", key="pep_excl")
            alf = pp.alfabeto(incluir=clases, excluir_residuos=excl)
            st.caption(f"Resulting alphabet ({len(alf)}): {', '.join(alf) if alf else 'empty'}")
        with st.expander("Sequence rules"):
            r1, r2 = st.columns(2)
            sin_rep = r1.checkbox("No repeated residues", key="pep_sinrep")
            maxcons = r1.number_input("Max identical in a row (0 = no limit)", 0, 10, 0, key="pep_cons")
            maxres = r2.number_input("Max times per residue (0 = no limit)", 0, 20, 0, key="pep_maxres")
            pre = r2.text_input("Starts with", key="pep_pre", placeholder="p. ej. KK").upper()
            suf = r1.text_input("Ends with", key="pep_suf", placeholder="p. ej. GG").upper()
        with st.expander("Physicochemical filters"):
            st.caption("In antimicrobial peptides, positive net charge and moderate hydrophobicity are the traits most associated with activity.")
            f1, f2 = st.columns(2)
            usar_q = f1.checkbox("Filter by net charge", key="pep_usaq")
            q_rng = f1.slider("Net charge at pH 7.4", -10.0, 10.0, (2.0, 9.0), 0.5,
                              key="pep_q", disabled=not usar_q)
            usar_g = f2.checkbox("Filter by hydropathy (GRAVY)", key="pep_usag")
            g_rng = f2.slider("GRAVY", -4.5, 4.5, (-1.0, 1.0), 0.1, key="pep_g", disabled=not usar_g)

        reglas = pp.Reglas(longitud=int(largo), alfabeto=alf, sin_repetir=sin_rep,
                           max_consecutivos=int(maxcons), max_por_residuo=int(maxres),
                           prefijo=pre, sufijo=suf,
                           carga_min=q_rng[0] if usar_q else None,
                           carga_max=q_rng[1] if usar_q else None,
                           gravy_min=g_rng[0] if usar_g else None,
                           gravy_max=g_rng[1] if usar_g else None)
        problemas = reglas.validar()
        for p in problemas:
            st.error(p)

    # Los extremos se eligen ANTES de generar: cambian la estructura que se construye y, por tanto,
    # la carga neta y las estructuras que se muestran.
    st.markdown("---")
    st.markdown("**Terminus chemistry**")
    e1, e2, e3 = st.columns(3)
    n_ac = e1.checkbox("Acetylate N-terminus", key="pep_nac",
                       help="Protects against aminopeptidases.")
    c_am = e2.checkbox("Amidate C-terminus", key="pep_cam",
                       help="Removes the terminal negative charge: +1 net charge, which usually increases antimicrobial activity.")
    ciclo = e3.checkbox("Cyclize head-to-tail", key="pep_ciclo",
                        help="Rigidifies the peptide and greatly reduces degrees of freedom, which also makes docking more reliable.")

    if entrada == "Generate library" and not problemas:
        st.caption(f"Available combinatorial space: ~{reglas.espacio():.0f} sequences.")
        # La química de los extremos NO entra aquí: no cambia las secuencias, que es lo único que
        # produce este boton. Quien depende de ella es la construcción de los compuestos, y es ese
        # boton el que se rehabilita al cambiarla.
        firma = (int(largo), int(cuantos), int(semilla), tuple(alf), sin_rep, int(maxcons),
                 int(maxres), pre, suf, usar_q, q_rng, usar_g, g_rng)
        hecho = S.get("_pep_firma") == firma and S.get("pep_seqs")
        # Este boton NO se bloquea: generar es barato y determinista con la semilla, así que volver a
        # pulsarlo no tiene coste ni riesgo. El bloqueo se reserva para construir las estructuras 3D,
        # que es la operación cara. Aquí solo se informa de si la biblioteca ya esta al dia.
        if st.button("Generate library", type="primary"):
            with st.spinner("Generating sequences..."):
                secuencias, aviso = pp.generate(reglas, int(cuantos), seed=int(semilla))
            S["pep_seqs"] = secuencias
            S["pep_aviso"] = aviso
            S["_pep_firma"] = firma
        if hecho:
            st.caption(f"Library generated with these parameters ({len(S['pep_seqs'])} sequences).")
        secuencias = secuencias or S.get("pep_seqs", [])
        aviso = aviso or S.get("pep_aviso", "")

    if aviso:
        st.warning(aviso)
    if not secuencias:
        return

    # La propia secuencia hace de nombre: es única (se deduplica al generar), valida como nombre de
    # archivo y mucho más informativa que un correlativo cuando se leen los resultados.
    filas = [pp.propiedades(s, c_amida=c_am, n_acetil=n_ac, ciclico=ciclo) for s in secuencias]
    df = pd.DataFrame(filas)[["nombre", "secuencia", "longitud", "carga_neta", "gravy",
                              "momento_hidrofobico", "fraccion_hidrofobica", "indice_boman"]]
    st.dataframe(df, width="stretch", hide_index=True, height=260)
    _descargar_tabla(df, "peptidos", key="pep_tabla")
    st.caption("`momento_hidrofobico` measures amphipathicity (hydrophobic vs. polar face when folded into a helix); `indice_boman` estimates the tendency to bind other proteins: above 2.5 kcal/mol is considered promiscuous.")

    nivel, msg = pp.viabilidad_docking(int(df["longitud"].max()), n_peptidos=len(df),
                                       hay_adcp=adcp.available())
    (st.success if nivel == "bueno" else st.warning if nivel == "medio" else st.error)(
        f"**Docking of {len(df)} peptides of {int(df['longitud'].max())} residues:** {msg}")

    S["_pep_preview"] = [(f["nombre"], f["secuencia"]) for f in filas[:24]]
    # El visualizador construye las estructuras y necesita la misma química de extremos que se
    # usara en el cribado; si no, dibuja el péptido lineal sin proteger y el ciclo no se ve.
    S["_pep_quimica"] = (bool(n_ac), bool(c_am), bool(ciclo))
    _firma_pep = (tuple(secuencias), n_ac, c_am, ciclo)
    _pep_listo = _ya_hecho("usar_peptidos", _firma_pep)
    if st.button("Use these peptides in the screening", type="primary", disabled=_pep_listo,
                 help="Change the sequences or terminus chemistry to build them again."
                      if _pep_listo else None):
        with st.spinner(f"Building the 3D structure of {len(secuencias)} peptides..."):
            smiles, nombres, fallos = [], [], 0
            for f in filas:
                smi = pp.to_smiles(f["secuencia"], n_acetil=n_ac, c_amida=c_am, ciclico=ciclo)
                if smi:
                    smiles.append(smi); nombres.append(f["nombre"])
                else:
                    fallos += 1
            made = lig.materialize(smiles, proj / "ligandos_entrada", names=nombres)
        hechos = {nm for nm, _, _ in made}
        S["ligands"] = [str(p) for _, p, _ in made]
        meta = pd.DataFrame([{"name": f["nombre"], "smiles": smi, "fuente": "péptido",
                              "producto": f["secuencia"], "iupac_name": None,
                              "viabilidad": f"{f['longitud']} residuos · carga {f['carga_neta']}"}
                             for f, smi in zip(filas, smiles) if f["nombre"] in hechos])
        (proj / "ligands_meta.csv").write_text(meta.to_csv(index=False))
        _marcar_hecho("usar_peptidos", _firma_pep)
        # Se nombran los que se pierden en la construcción 3D. Antes solo se contaban los fallos de
        # SMILES, así que un péptido que no lograba encajarse en tres dimensiones desaparecía del
        # cribado sin dejar rastro y el recuento final no cuadraba con la biblioteca.
        perdidos = [n for n in nombres if n not in hechos]
        st.toast(f"{len(made)} peptides built and ready for step 3.", icon="🧬")
        st.success(f"{len(made)} peptides ready for step 3."
                   + (f" {fallos} could not be built." if fallos else ""))
        if perdidos:
            st.warning("Could not generate the 3D structure of: " + ", ".join(perdidos)
                       + ". They are long, flexible chains; try cyclizing them to rigidify them.")
    if _pep_listo:
        st.caption(f"{len(S['ligands'])} peptides built with these parameters.")


# ---------------------------------------------------------------- 2. ligandos
def _etapa_ligandos():
    st.subheader("What do you want to dock?")
    modo = st.radio("Source of the compounds",
                    ["Build by reaction", "Generate peptides", "Upload ready ligands"],
                    horizontal=True, key="modo_ligandos")

    if modo == "Generate peptides":
        _modo_peptidos()
    elif modo == "Build by reaction":
        S["lead"] = None
        izq, der = st.columns(2)
        with izq:
            st.markdown("#### Core (your starting molecule)")
            nuc = st.text_input("Core SMILES", value=S.get("nuc_smiles", ""),
                                placeholder="O=C(O)c1ccc2[n+]([O-])onc2c1")
            fnuc = st.file_uploader("...or core file", type=["sdf", "mol2", "mol"])
            if fnuc is not None:
                d = proj / "nucleo"
                d.mkdir(parents=True, exist_ok=True)
                fp = d / fnuc.name
                fp.write_bytes(fnuc.getvalue())
                nuc = _to_smiles(fp) or nuc
            S["nuc_smiles"] = nuc
            rxkey = st.selectbox("Reaction", list(rx.REACTIONS), key="rx_reaccion",
                                 format_func=lambda k: rx.get(k).nombre)
            reaction = rx.get(rxkey)
            st.caption(reaction.descripcion)
            if nuc:
                aplica = any(r.key == rxkey for r in rx.applicable(nuc))
                if not aplica:
                    st.warning(f"The core has no {reaction.lead_grupo}; this reaction does not apply.")
                else:
                    sitios = rx.lead_sites(nuc, reaction)
                    st.success(f"The core can undergo {reaction.nombre}: {len(sitios)} reactive site(s).")
                    idx = 0
                    if len(sitios) > 1:
                        idx = st.selectbox("Growth point", range(len(sitios)),
                                           format_func=lambda i: f"atoms {sitios[i]['atomos']}")
                    hl = sitios[idx]["atomos"] if sitios else []
                    # Se dibuja en el panel derecho, que es el visualizador.
                    S["_nucleo_png"] = vw.molecule_png_indexed(nuc, highlight=hl, size=420)
        with der:
            if reaction.kind == "coupling":
                st.markdown("#### Reagents that couple")
                usar_int = st.checkbox("Internal library", value=True,
                                       help=f"{len(rg.load_internal(reaction)) if nuc else 0} curated reagents.")
                ups = st.file_uploader("Your reagents (csv/xlsx with columns name and smiles · sdf · mol2 · smi)",
                                       type=["csv", "xlsx", "xls", "sdf", "mol2", "mol", "smi"],
                                       accept_multiple_files=True)
                with st.expander("What columns must my Excel/CSV have?"):
                    st.markdown(
                        "Two columns: a **name** and a **SMILES**. Accepted headers are:\n- Name: `name`, `nombre`, `compound`, `compuesto`, `Alcohol origen`, `Nombre clave`\n- SMILES: `smiles`, `smile`, `SMILES alcohol`\n\nDeduplicated by structure (InChIKey) and filtered to those bearing the reaction group (for esterification, an alcohol/phenol OH; acids and amines are discarded).")
                    st.dataframe(pd.DataFrame({"name": ["Bencílico", "Mentol", "Ciclohexanol"],
                                               "smiles": ["OCc1ccccc1", "CC(C)C1CCC(C)CC1O", "OC1CCCCC1"]}),
                                 width="stretch", hide_index=True)
                usar_pc = st.checkbox("Supplement with PubChem (experimental, needs internet)", value=False)
                pc_max = st.number_input("PubChem maximum", 5, 100, 25) if usar_pc else 25
                upaths = []
                if ups:
                    d = proj / "reactivos"; d.mkdir(parents=True, exist_ok=True)
                    for u in ups:
                        (d / u.name).write_bytes(u.getvalue()); upaths.append(str(d / u.name))
                if st.button("Gather reagents", type="primary"):
                    with st.spinner("Gathering and deduplicating..."):
                        reags, info = rg.build(reaction, use_internal=usar_int, user_paths=upaths,
                                               use_pubchem=usar_pc, pubchem_max=int(pc_max))
                    S["reagents"] = [(r.name, r.smiles, r.source, r.inchikey) for r in reags]
                    S["reag_info"] = info
                    if info.get("aviso_pubchem"):
                        st.warning(info["aviso_pubchem"])
                if S.get("reag_info"):
                    info = S["reag_info"]
                    st.write(f"**{info['total']} reagents** — " +
                             " · ".join(f"{k}: {v}" for k, v in info["por_fuente"].items()))
                    dfa = pd.DataFrame([{"nombre": n, "SMILES": s, "fuente": src} for n, s, src, ik in S["reagents"]])
                    st.dataframe(_shade(dfa, "fuente"), width="stretch", height=240)
                    st.caption("Highlighted = reagents you provided.")
            else:   # decoración: sin reactivo externo, usa sustituyentes internos
                st.markdown("#### Substituents")
                st.caption("Decoration uses small internal groups (F, Cl, CN, OMe...); you upload no reagents.")
                c1, c2 = st.columns(2)
                S["n_analogs"] = c1.number_input("How many analogues", 1, 200, S.get("n_analogs", 20))
                S["n_sub"] = c2.multiselect("Number of substitutions", [1, 2, 3], default=S.get("n_sub", [1]))
                S["use_ml"] = st.checkbox("Predict ADMET with AI (slower the first time)", value=S.get("use_ml", True))
                b = AdmelabBridge()
                if not b.available():
                    st.error("Cannot find the design engine (admelab).")
                elif nuc and st.button("Generate analogues"):
                    with st.spinner("Generating and predicting properties..."):
                        d = b.design(nuc, use_ml=bool(S["use_ml"]),
                                     n_substitutions=S.get("n_sub", [1]) or [1], max_rows=int(S["n_analogs"]))
                    S["products"] = [dict(producto=(r.get("name") or f"analogo{i + 1:03d}"), smiles=r["SMILES"],
                                          fuente="interno", sintetizable=True, viabilidad="decoración")
                                     for i, r in enumerate(d.rows) if r.get("SMILES")]
                    if d.n_generated < int(S["n_analogs"]):
                        st.warning(f"{d.n_generated} analogues generated: with {S.get('n_sub', [1])} "
                                   "substitution(s) the chemical space runs out there. Try 2.")

        # ---- productos (coupling: esterifica; decoración: ya generados arriba) ----
        if reaction.kind == "coupling" and nuc and S.get("reagents"):
            st.markdown("---")
            b = AdmelabBridge()
            if not b.available():
                st.error("Cannot find the reaction engine (admelab).")
            else:
                firma_p = (nuc, rxkey, tuple(sorted(ik for _n, _s, _src, ik in S["reagents"])))
                hecho_p = _ya_hecho("productos", firma_p) and S.get("products")
                if hecho_p:
                    st.caption(f"{len(S['products'])} products already built with this core and "
                               "these reagents. Change one to rebuild.")
                if st.button("Build products", type="primary", disabled=bool(hecho_p)):
                    alcs = [{"name": n, "smiles": s} for n, s, src, ik in S["reagents"]]
                    with st.spinner("Building the series..."):
                        prods = b.esterify(nuc, alcs, policy="preferred")
                    src_by_ik = {ik: src for n, s, src, ik in S["reagents"]}
                    for p in prods:
                        p["fuente"] = src_by_ik.get(rg.inchikey(p.get("alcohol_smiles", "") or ""), "?")
                        p["producto"] = p.get("alcohol")
                    _marcar_hecho("productos", firma_p)
                    S["products"] = prods

        prods = S.get("products")
        if prods:
            st.markdown("---")
            n_ok = sum(1 for p in prods if p.get("sintetizable"))
            st.info(f"{n_ok} of {len(prods)} products are synthesizable by this reaction.")
            if not any(p.get("alcohol_smiles") for p in prods):
                # El nombrado compone el nombre a partir de las dos piezas de la reacción (radical y
                # acilo) y lo verifica con OPSIN. En decoración no hay dos piezas: habria que nombrar
                # una molécula arbitraria desde su estructura, que es un problema distinto y sin
                # solución fiable offline.
                st.caption("The IUPAC name is only generated for coupling reactions, where it is composed from the two fragments and verified with OPSIN. In decoration, products are identified by their SMILES.")
            if any(p.get("alcohol_smiles") for p in prods):
                if st.button("Name (IUPAC, verified with OPSIN)"):
                    with st.spinner("Naming and verifying by round-trip..."):
                        named = AdmelabBridge().name_esters(
                            [p["smiles"] for p in prods], [p.get("alcohol_smiles") or "" for p in prods],
                            acid_smiles=nuc, alcohol_names=[p.get("producto") for p in prods], use_web=True)
                    by = {n["smiles"]: n for n in named}
                    for p in prods:
                        inf = by.get(p["smiles"], {})
                        p["iupac_name"] = inf.get("iupac_name"); p["iupac_verif"] = inf.get("verified")
                    S["products"] = prods
                    nver = sum(1 for p in prods if p.get("iupac_verif"))
                    st.success(f"{nver} of {len(prods)} with a verified IUPAC name. The rest keep their label "
                               "(the alcohol name); they are niche and OPSIN does not always cover them offline.")
            dfp = pd.DataFrame(prods)
            cols = [c for c in ("producto", "iupac_name", "fuente", "oh_type", "viabilidad", "sintetizable", "smiles")
                    if c in dfp.columns]
            st.dataframe(_shade(dfp[cols], "fuente"), width="stretch", height=320)
            st.caption("Highlighted = products with YOUR reagents. `sintetizable`=False are infeasible by this reaction.")

            st.caption("The 2D structures of the products are shown in the right panel so you can check the bond and stereochemistry.")

            with st.expander("ADMET report (predicts ~40 endpoints with AI for ALL at once)"):
                if st.button("Predict ADMET"):
                    with st.spinner("Predicting with ADMET-AI for all (the model downloads the first time)..."):
                        pr = AdmelabBridge().predict([p["smiles"] for p in prods], use_ml=True)
                    S["admet"] = {rg.inchikey(r.get("SMILES")): r for r in pr.rows}
                if S.get("admet"):
                    _render_adme(S["admet"], [(p.get("producto") or f"prod{i}", p["smiles"])
                                              for i, p in enumerate(prods)], keyp="lig")

            c_sel1, c_sel2 = st.columns(2)
            solo_ok = c_sel1.checkbox("Dock only the synthesizable ones", value=True)
            incluir_nuc = c_sel2.checkbox("Add the bare core (reference)", value=True,
                                          help="Docks the unesterified core as a baseline: reveals how much activity the scaffold contributes on its own, apart from the tail.")
            firma_uso = (tuple(p.get("smiles") for p in prods), solo_ok, bool(incluir_nuc and nuc))
            usado = _ya_hecho("usar_productos", firma_uso)
            if usado:
                st.caption(f"These products are already loaded for the screening ({len(S['ligands'])} "
                           "compounds). Change the selection to regenerate them.")
            if st.button("Use these products in the screening", type="primary", disabled=usado):
                elegidos = [p for p in prods if (p.get("sintetizable") or not solo_ok)]
                # El nucleo desnudo entra como un candidato más (mismo pipeline), etiquetado aparte.
                if incluir_nuc and nuc:
                    elegidos = [dict(producto="nucleo_libre", smiles=nuc, fuente="núcleo",
                                     iupac_name=None, viabilidad="referencia (sin esterificar)",
                                     sintetizable=True)] + elegidos
                nombres = [lig.safe_name(p.get("producto") or f"prod{i}") for i, p in enumerate(elegidos)]
                # A la carpeta de ENTRADA, no a 'ligands': esa la borra la limpieza de cada corrida.
                with st.spinner(f"Generating 3D of {len(elegidos)} compounds..."):
                    made = lig.materialize([p["smiles"] for p in elegidos], proj / "ligandos_entrada", names=nombres)
                hechos = {nm for nm, _, _ in made}
                S["ligands"] = [str(p) for _, p, _ in made]
                meta = pd.DataFrame([{"name": nm, "smiles": p.get("smiles"), "fuente": p.get("fuente", "?"),
                                      "producto": p.get("producto"), "iupac_name": p.get("iupac_name"),
                                      "viabilidad": p.get("viabilidad")}
                                     for (nm, p) in zip(nombres, elegidos) if nm in hechos])
                (proj / "ligands_meta.csv").write_text(meta.to_csv(index=False))
                _marcar_hecho("usar_productos", firma_uso)
                extra = " (includes the bare core as reference)" if (incluir_nuc and nuc) else ""
                st.toast(f"{len(made)} compounds built and ready for step 3.", icon="⚗️")
                st.success(f"{len(made)} compounds ready for step 3{extra}.")

    else:  # Subir ligandos listos
        S["lead"] = None
        ups = st.file_uploader("Upload ligands", type=["mol2", "sdf", "mol", "smi"], accept_multiple_files=True)
        if ups:
            d = proj / "ligandos_entrada"
            d.mkdir(parents=True, exist_ok=True)
            for u in ups:
                (d / u.name).write_bytes(u.getvalue())
            S["ligands"] = [str(p) for p in sorted(d.iterdir())]
            # Se leen los SMILES de los archivos subidos y se escribe ligands_meta.csv. Sin esto no
            # hay estructura de la que derivar ADME, eficiencia de ligando, SAscore ni PAINS: esas
            # columnas salian vacias para todo lo que no fuese el control.
            smap = sc.build_smiles_map(str(d))
            filas = []
            for p in S["ligands"]:
                nombre = Path(p).stem
                smi = smap.get(sc.normalize_key(nombre))
                filas.append({"name": nombre, "smiles": smi, "fuente": "subido",
                              "producto": nombre, "iupac_name": None, "viabilidad": None})
            sin_smiles = sum(1 for f in filas if not f["smiles"])
            (proj / "ligands_meta.csv").write_text(pd.DataFrame(filas).to_csv(index=False))
            if sin_smiles:
                st.warning(f"{sin_smiles} of {len(filas)} ligands gave no readable structure; "
                           "there will be no ADMET or descriptors for them.")
        if S["ligands"]:
            st.write(f"**{len(S['ligands'])} ligands:** " + ", ".join(Path(p).name for p in S["ligands"][:8]))
            ml = proj / "ligands_meta.csv"
            if ml.exists():
                mdf = pd.read_csv(ml)
                con = mdf["smiles"].notna().sum() if "smiles" in mdf.columns else 0
                st.caption(f"Structure read from {con} of {len(mdf)}: allows computing ADMET, "
                           "ligand efficiency, SAscore and PAINS alerts.")
                items = [(r["name"], r["smiles"]) for _, r in mdf.iterrows() if pd.notna(r.get("smiles"))]
                if items:
                    with st.expander("ADMET report of the uploaded ligands"):
                        if st.button("Predict ADMET", key="adme_subidos"):
                            with st.spinner("Predicting with ADMET-AI..."):
                                pr = AdmelabBridge().predict([s for _, s in items], use_ml=True)
                            S["admet"] = {**(S.get("admet") or {}),
                                          **{rg.inchikey(r.get("SMILES")): r for r in pr.rows}}
                        if S.get("admet"):
                            _render_adme(S["admet"], items, keyp="sub")

    cup = st.file_uploader("Controls (co-crystallized ligand)", type=["mol2", "sdf", "mol"],
                           accept_multiple_files=True,
                           help="If you already extracted it in step 1, no need to upload anything.")
    if cup:
        d = proj / "receptores"
        d.mkdir(parents=True, exist_ok=True)
        nuevos = []
        for u in cup:
            p = d / u.name
            try:
                p.write_bytes(u.getvalue())
            except Exception as e:
                st.error(f"Could not save {u.name}: {e}")
                continue
            if str(p) not in S["controls"]:
                S["controls"].append(str(p))
                nuevos.append(u.name)
        # Sin este aviso, subir un control parecía no hacer nada: el archivo se guardaba pero no
        # había ninguna senal de que se hubiera registrado.
        if nuevos:
            st.success(f"Control(s) loaded: {', '.join(nuevos)}. "
                       f"Total controls: {len(S['controls'])}.")
        st.caption("Loaded controls and those extracted in step 1 are docked alongside the ligands and define the reference fingerprint. With several receptors, assign each to its receptor in step 3.")

# ---------------------------------------------------------------- 3. ejecutar
def _etapa_ejecutar():
    st.subheader("Run the screening")
    recs = [Path(p) for p in S["receptors"]]
    ctrls = [Path(p) for p in S["controls"]]
    ligs = [Path(p) for p in S["ligands"]]
    st.write(f"Receptores: **{len(recs)}** · Controles: **{len(ctrls)}** · "
             + (f"Lead: `{S.get('lead')}`" if S.get("lead") else f"Ligandos: **{len(ligs)}**"))
    if recs and not ctrls:
        st.warning("No control loaded. The control is docked alongside the ligands and defines the reference; without it there is no baseline or validation. Extract the co-crystallized one in step 1 (or upload it below). If you already extracted it, check it is in the project's `receptores/` folder.")

    boxes = {}
    site_boxes = {}   # docking hibrido: receptor -> [(etiqueta, Box)] con bolsillos adicionales
    if recs:
        st.markdown("**Search box** — where to search inside the protein.")
        st.caption("Most reliable centered on the co-crystallized ligand: it marks the real site. The geometric center or a cofactor point elsewhere.")
        xbox = {}
        # La asignación control->receptor se resuelve sola por geometría: cada control cae en el
        # espacio de SU receptor. Solo se pide intervencion cuando algun control no se puede ubicar
        # (subido suelto, sin coincidir con ninguna estructura cargada). De la asignación dependen la
        # caja centrada en el control, la huella de referencia y la validación por redocking.
        manual = {}
        if len(recs) > 1 and ctrls:
            auto = pl._assign_controls(ctrls, recs, {})
            sin_ubicar = [c for c in ctrls if sc.normalize_key(c.stem) not in auto]
            if sin_ubicar:
                with st.expander("Assign the missing controls to their receptor", expanded=True):
                    st.caption("These controls could not be placed by geometry; indicate which receptor each belongs to.")
                    _rec_labels = ["(none)"] + [r.stem for r in recs]
                    for c in sin_ubicar:
                        sel = st.selectbox(f"Control \"{c.stem}\"", _rec_labels,
                                           key=f"ctrlrec_{sc.normalize_key(c.stem)}")
                        if sel != "(none)":
                            manual[sc.normalize_key(c.stem)] = sel
            else:
                _pares = ", ".join(f"{c.stem} → {auto[sc.normalize_key(c.stem)]}" for c in ctrls)
                st.caption(f"Controls assigned automatically by geometry: {_pares}.")
        asignacion = pl._assign_controls(ctrls, recs, manual)
        S["_control_map"] = manual
        S.setdefault("pockets", {})
        for _i_r, r in enumerate(recs):
            # Separador y titulo destacado por receptor: con varios apilados era difícil saber donde
            # empezaban los ajustes de cada uno. Se evita envolver el bloque en un contenedor para no
            # reindentar todo el cuerpo del bucle, que es fuente de errores.
            st.divider()
            _th, _tv = st.columns([3, 1], vertical_alignment="center")
            _th.markdown(f"### ▸ {r.name}")
            # Enfocar el visor en el receptor que se esta editando, en vez de quedarse en otro.
            if _tv.button("Show in the viewer", key=f"verrec_{r.name}", use_container_width=True):
                S["vis_box_rec"] = str(r)
                st.rerun()
            grupos = dk.hetero_groups(r)
            ctrl = next((c for c in ctrls if asignacion.get(sc.normalize_key(c.stem)) == r.stem), None)
            # detección de pockets (fpocket), bajo demanda
            b1, b2 = st.columns([1, 3])
            ya_pk = bool(S["pockets"].get(str(r)))
            if b1.button("Detect pockets", key=f"pk_{r.name}", type="primary",
                         disabled=not pk.fpocket_available() or ya_pk,
                         help="Cavities already detected for this receptor." if ya_pk else None):
                S["vis_box_rec"] = str(r)              # el visor sigue al receptor recien detectado
                with st.spinner("Searching cavities with fpocket..."):
                    S["pockets"][str(r)] = pk.detect(r)
                # Se muestran en el visor sin que haya que activarlo: es el motivo de detectarlas.
                S["vis_ver_cav"] = True
                st.rerun()
            pkts = S["pockets"].get(str(r), [])
            if not pkts and not pk.fpocket_available():
                b2.caption("fpocket not installed: `conda install -n cribado -c conda-forge fpocket`.")
            # fuente de la caja. Al cambiarla, el visor pasa a este receptor: se esta editando, y
            # ver el cambio sobre otra estructura confundiría. Se hace por callback porque corre
            # antes de instanciar el selector del visor, que es cuando su valor aún se puede fijar.
            pk_opts = {p["label"]: p for p in pkts}
            opts = ([f"Center on the control ({ctrl.name})"] if ctrl else []) \
                + list(pk_opts.keys()) + ["Automatic"] + list(grupos.keys())
            pick = st.selectbox("Box source", opts, key=f"box_{r.name}",
                                on_change=lambda rr=str(r): S.__setitem__("vis_box_rec", rr))

            if pkts:
                cats_r = {x.lower() for x in S.get(f"cat_{r.stem}", [])}
                elegida = pk_opts.get(pick, {}).get("n")
                # Todos los sitios que entraran al docking: el principal y los del hibrido. El
                # multiselect del hibrido se dibuja más abajo, así que se lee su valor guardado.
                hib_labels = set(S.get(f"sitios_{r.name}") or [])
                usados = {elegida} | {p["n"] for p in pkts if p["label"] in hib_labels}
                usados.discard(None)
                filas, cav = [], []
                for i, p in enumerate(pkts[:8]):
                    en_uso = p["n"] in usados
                    # Amarillo solo el sitio principal; los del hibrido mantienen su color para
                    # poder distinguirlos entre si, pero se dibujan igual de opacos por estar en uso.
                    color = (vw.COLOR_ELEGIDA if p["n"] == elegida
                             else vw.PALETA_CAVIDADES[i % len(vw.PALETA_CAVIDADES)])
                    cav.append({"alpha": p.get("alpha_xyz"), "color": color, "elegida": en_uso})
                    pr = p.get("props", {})
                    resid = p.get("residues") or []
                    hay_cat = bool(cats_r and {x.lower() for x in resid} & cats_r)
                    fila = {
                        "Color": vw.emoji_de_color(color), "Pocket": p["n"],
                        "Used": ("principal" if p["n"] == elegida else "híbrido") if en_uso else "",
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
                        # Flexibility se omite: fpocket la deriva de los B-factors, que la preparación
                        # con PDBFixer deja a 0; siempre saldría vacia y confundiría.
                        "Residues": ", ".join(resid[:14]) + ("…" if len(resid) > 14 else ""),
                    }
                    # La columna catalitica solo aparece si el usuario ya designo residuos ancla;
                    # sin esa referencia no hay forma de decidirlo y un "?" no informa de nada.
                    if cats_r:
                        fila["Catalytic"] = "sí" if hay_cat else "no"
                    filas.append(fila)
                # Por receptor: con más de uno, un único dict global lo sobrescribia en cada
                # iteración y el visor acababa mostrando (o borrando) las cavidades del último.
                S.setdefault("_cavidades", {})[str(r)] = cav
                dfp = pd.DataFrame(filas)
                st.dataframe(
                    dfp, width="stretch", hide_index=True, height=230,
                    column_config={"Color": st.column_config.TextColumn(
                        "Color", width="small",
                        help="Color the cavity is drawn with in the viewer")})
                _descargar_tabla(pd.DataFrame([{"Pocket": p["n"], **p.get("props", {}),
                                                "Residues": ", ".join(p.get("residues") or [])}
                                               for p in pkts]),
                                 f"cavidades_{r.stem}", key=f"cav_{r.name}")
                st.caption("All cavities are drawn at once in the right panel. The one **used for docking** is highlighted and more opaque. `Cavity` is its real extent; `Box` is the search region, with a 14 Å minimum because below that a ligand would not fit (marked `*` when that minimum was applied).")
                with st.expander("All properties fpocket computes"):
                    st.dataframe(pd.DataFrame([{"pocket": p["n"], **p.get("props", {})} for p in pkts]),
                                 width="stretch", hide_index=True)
            else:
                S.get("_cavidades", {}).pop(str(r), None)

            if ctrl and pick.startswith("Center on the control"):
                base = dk.box_from_file(ctrl)
            elif pick in pk_opts:
                base = pk.pocket_box(pk_opts[pick])
            elif pick in grupos:
                base = grupos[pick]
            else:
                base = dk.auto_box(r)
            # Se siembran los campos con la caja base al cambiar de fuente O si aún no existen.
            # Sin la segunda condición, al volver a esta pestana los number_input arrancaban en 0
            # y la caja saltaba fuera de la proteína (el widget, sin valor guardado, usa su mínimo).
            if S.get(f"src_{r.name}") != pick or f"cx_{r.name}" not in S:
                S[f"src_{r.name}"] = pick
                for k, v in (("cx", base.cx), ("cy", base.cy), ("cz", base.cz),
                             ("sx", base.sx), ("sy", base.sy), ("sz", base.sz)):
                    S[f"{k}_{r.name}"] = float(v)
            gc, gs = st.columns(2)
            with gc:
                st.markdown("**Center** — where the box sits (Å)")
                st.caption("Moves the box through space. The axes are shown in the viewer on the right.")
                cc = st.columns(3)
                cx = cc[0].number_input("← X →", step=1.0, key=f"cx_{r.name}", format="%.1f",
                                        help="Left / right (red axis).")
                cy = cc[1].number_input("↓ Y ↑", step=1.0, key=f"cy_{r.name}", format="%.1f",
                                        help="Down / up (green axis).")
                cz = cc[2].number_input("⊙ Z ⊗", step=1.0, key=f"cz_{r.name}", format="%.1f",
                                        help="Into / out of the screen (blue axis).")
            with gs:
                st.markdown("**Size** — how much the box spans (Å)")
                st.caption("Grows or shrinks each side. If the ligand does not fit, Vina fails.")
                cs = st.columns(3)
                sx = cs[0].number_input("width X", min_value=6.0, step=1.0, key=f"sx_{r.name}", format="%.1f")
                sy = cs[1].number_input("height Y", min_value=6.0, step=1.0, key=f"sy_{r.name}", format="%.1f")
                sz = cs[2].number_input("depth Z", min_value=6.0, step=1.0, key=f"sz_{r.name}", format="%.1f")
            boxes[str(r)] = dk.Box(cx, cy, cz, sx, sy, sz)
            # El panel derecho se dibuja después que el izquierdo, así que puede leer esto.
            S.setdefault("_boxes", {})[str(r)] = boxes[str(r)].as_dict()
            st.caption("The box is drawn over the receptor in the right panel.")
            # Un ligando tiene que poder GIRAR dentro de la caja, no solo caber quieto. Con una
            # caja del tamaño justo la busqueda queda limitada a las orientaciones que entran, y el
            # resultado parece un mal acoplamiento cuando es una restricción geométrica involuntaria.
            # Se mide contra los ligandos (que se acoplan en TODOS los receptores) y el control de
            # ESTE receptor; los controles de otras dianas no entran en esta caja, así que un control
            # grande de otra diana no debe disparar el aviso aquí.
            _lig_este = list(S["ligands"]) + ([str(ctrl)] if ctrl else [])
            _minimo = dk.caja_minima(_lig_este)
            if _minimo and min(sx, sy, sz) < _minimo:
                st.warning(
                    f"The largest ligand is **{_minimo - 4:.0f} Å** on its major axis and the box "
                    f"is {min(sx, sy, sz):.0f} Å on its shortest side. It fits, but cannot "
                    f"reorient: the search is restricted to the orientations that fit. "
                    f"Raise all three sides to at least **{_minimo:.0f} Å**.")

            # Docking hibrido: acoplar el mismo ligando en varios bolsillos del mismo receptor.
            if pkts:
                # El bolsillo que ya es la caja principal no se ofrece: elegirlo creaba un sitio
                # duplicado que luego había que colapsar.
                disponibles_hib = [p["label"] for p in pkts if p["label"] != pick]
                extra = st.multiselect(
                    "Also dock in other pockets (hybrid docking)",
                    disponibles_hib, key=f"sitios_{r.name}",
                    help="Each chosen pocket is docked separately and gets its own ranking. Reveals whether a compound prefers the catalytic site or slips into an allosteric one.")
                if extra:
                    lst = [("principal", boxes[str(r)])]
                    for lab in extra:
                        pdd = next((p for p in pkts if p["label"] == lab), None)
                        if pdd:
                            lst.append((f"Pk{pdd['n']}", pk.pocket_box(pdd)))
                    site_boxes[str(r)] = lst
                    st.caption(f"Hybrid docking: {len(lst)} sites (main + {len(extra)} pocket(s)).")
                S[f"_hib_sel_{r.name}"] = set(extra)

    params = _params_docking()
    c1, c2 = st.columns([2, 1])
    reuse = c1.checkbox("Reuse previous calculations from this folder", value=False,
                        help="Off, each run recomputes everything. Enable only if nothing has changed: reusing poses made with another box gives false results.")
    if c2.button("Delete this folder's results"):
        _confirmar_borrado(proj)
    st.caption(f"Everything is saved in `{proj}` — poses, complexes, PLIP XML and the CSV tables.")

    if not recs:
        st.info("Prepare at least one receptor in step 1.")
    elif not (S.get("lead") or ligs):
        st.info("Choose compounds in step 2.")
    else:
        # El boton se bloquea mientras la configuración no cambie: un cribado dura minutos y
        # volver a lanzarlo por olvido o por una doble pulsacion desperdicia el trabajo.
        firma = (tuple(sorted(str(x) for x in recs)), tuple(sorted(str(x) for x in ligs)),
                 tuple(sorted(str(x) for x in ctrls)),
                 tuple(sorted((k, tuple(v.as_dict().values())) for k, v in boxes.items())),
                 tuple(sorted((k, len(v)) for k, v in site_boxes.items())),
                 tuple(sorted(params.items())), reuse, str(proj))
        hecho = _ya_hecho("run", firma)
        if st.button("Run", type="primary", disabled=hecho,
                     help="Already run with this configuration. Change something to launch again."
                          if hecho else None):
            cfg = pl.RunConfig(receptors=recs, out_dir=proj, lead=S.get("lead") or None, ligands=ligs,
                               controls=ctrls, boxes=boxes, site_boxes=site_boxes,
                               control_map=S.get("_control_map") or {},
                               n_analogs=int(S.get("n_analogs", 20)),
                               n_substitutions=S.get("n_sub", [1]) or [1], use_ml=bool(S.get("use_ml", True)),
                               reuse=reuse, **params)
            # El registro se guarda además de mostrarse: contiene los avisos que explican una corrida
            # (cobertura por sitio, fallos por ligando, validación del control) y se perdía al
            # cambiar de pestana, justo cuando hacen falta para interpretar los resultados.
            S["_log_run"] = []

            def _paso(n, d):
                S["_log_run"].append((n, d))
                st.write(f"**{n}** · {d}")

            with st.status("Running...", expanded=True) as status:
                try:
                    pl.run(cfg, on_step=_paso)
                    _marcar_hecho("run", firma)
                    status.update(label="Screening completed", state="complete")
                    S["_log_estado"] = "completo"
                    st.toast("Screening completed. Go to the results tab.", icon="✅")
                    st.success("Done. Go to the results tab.")
                except Exception as e:
                    status.update(label="Failed", state="error")
                    S["_log_estado"] = "error"
                    S["_log_run"].append(("error", str(e)))
                    st.toast("The screening failed. Check the message.", icon="⚠️")
                    st.error(str(e))
        elif S.get("_log_run"):
            _registro_corrida()
        if hecho:
            st.caption("Screening completed with this configuration. Change a parameter to enable the button again.")


def _registro_corrida():
    """Registro de la última corrida, conservado entre pestanas."""
    estado = S.get("_log_estado", "completo")
    with st.status("Screening completed" if estado == "completo" else "The run failed",
                   state="complete" if estado == "completo" else "error", expanded=False):
        for n, d in S.get("_log_run", []):
            st.write(f"**{n}** · {d}")
        _descargar_tabla(pd.DataFrame(S["_log_run"], columns=["etapa", "detalle"]),
                         "registro_corrida", key="log_run")

def _resumen_pleiotropico(rk, dianas):
    """Que compuesto une bien en VARIAS dianas a la vez. Un resumen por diana ya lo dan los bloques
    de abajo; este busca amplio espectro, no potencia en una sola.

    La efectividad se compara entre dianas porque esta normalizada contra el control de cada una: es
    la misma escala 0-1 aunque las energías no lo sean. Se ordena por la efectividad Mínima entre
    dianas —el eslabon más debil—, que es lo que de verdad mide 'bueno en todas'; una media dejaría
    pasar un compuesto excelente en una diana e inútil en otra."""
    import numpy as np
    sub = rk[rk.get("is_control", 0) != 1].copy() if "is_control" in rk.columns else rk.copy()
    sub = sub[pd.notna(sub.get("efectividad_pct"))]
    if sub.empty:
        return
    sub["_diana"] = sub["receptor"].map(sc.base_of)
    # Un compuesto puede tener varios bolsillos por diana (hibrido): se queda su mejor efectividad.
    mejor = (sub.groupby(["compound", "_diana"])["efectividad_pct"].max()
             .reset_index())
    piv = mejor.pivot(index="compound", columns="_diana", values="efectividad_pct")
    presente_en_todas = piv.dropna()
    st.markdown("### Pleiotropic summary — activity across several targets")
    if presente_en_todas.empty:
        st.caption("No compound docked in all targets; there is no broad-spectrum comparison.")
        return
    presente_en_todas = presente_en_todas.assign(
        **{"minimum": presente_en_todas.min(axis=1).round(1),
           "mean": presente_en_todas.mean(axis=1).round(1)})
    tabla = presente_en_todas.sort_values("minimum", ascending=False).reset_index()
    tabla.columns = ["compound"] + [f"{c} (%)" if c in dianas else c for c in tabla.columns[1:]]
    st.caption("Effectiveness (%) of each compound in each target, ordered by the **minimum** across targets: broad-spectrum ones on top. Only those docked in all.")
    st.dataframe(tabla.round(1), width="stretch", hide_index=True,
                 height=min(340, 60 + 34 * len(tabla)))
    _descargar_tabla(tabla, "pleiotropico", key="pleio")
    mejor_amplio = tabla.iloc[0]
    st.success(f"Best broad-spectrum: **{mejor_amplio['compound']}** "
               f"(minimum {mejor_amplio['minimum']:.0f} % across {len(dianas)} targets).")
    st.divider()


# ---------------------------------------------------------------- 4. resultados
def _etapa_resultados():
    st.subheader("Results")
    meta_p, inter_p, dock_p = proj / "run.json", proj / "interacciones.csv", proj / "resultados_docking.csv"
    if not (meta_p.exists() and inter_p.exists()):
        st.info("No results in this folder yet. Run step 3.")
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

        st.markdown("**Catalytic / anchor residues**. The score also rewards the quality of the pocket's other interactions.")
        st.caption("Auto-suggested from the directional interactions of the crystallographic ligand. Edit them if you know your target's real catalytic site.")
        cat, sec = {}, {}
        cols = st.columns(max(1, len(ref_info)))
        for i, R in enumerate(sorted(ref_info)):
            opciones = sorted(set(ref_info[R].get("residues", [])) | set(pocket_res_map.get(R, [])),
                              key=lambda r: (sc.resnum(r), r))
            sugeridos = ref_info[R].get("autocat", [])
            prev = ([x for x in meta.get("catalytic", {}).get(R, []) if x in opciones]
                    or [x for x in sugeridos if x in opciones])
            cat[R] = cols[i].multiselect(f"{R}  ·  ref: {ref_info[R].get('src', '?')}", opciones,
                                         default=prev, key=f"cat_{R}")
            # Un residuo ya designado catalitico no puede ofrecerse también como secundario: son
            # roles excluyentes y el motor da prioridad al catalitico, así que elegirlo dos veces
            # solo confunde.
            libres = [x for x in opciones if x not in cat[R]]
            # Si un residuo ya elegido como secundario se promueve a catalitico, se retira del
            # estado del otro selector para que no quede un valor fuera de las opciones.
            k_sec = f"sec_{R}"
            if k_sec in S:
                S[k_sec] = [x for x in S[k_sec] if x in libres]
            prev_s = [x for x in meta.get("secondary", {}).get(R, []) if x in libres]
            sec[R] = cols[i].multiselect(f"{R} · secondary (bonus, not required)", libres,
                                         default=prev_s, key=f"sec_{R}")
        st.caption("**Gate** (catalytic): mandatory, missing them is penalized. **Secondary**: known pocket anchors that add more than an ordinary contact (×w_sec) but are not required.")

        # Validación: si el control no recupera su postura, el resto no es fiable
        val_p = proj / "validacion_redocking.csv"
        if val_p.exists():
            val = pd.read_csv(val_p)
            msg = vl.resumen(val)
            (st.success if not msg.startswith("WARNING") else st.error)(msg)
            with st.expander("Redocking validation detail"):
                st.dataframe(val, width="stretch", hide_index=True)
                st.caption("RMSD against the co-crystallized ligand. Valid below 2 Å.")

        # Validación de INTERACCIONES: cuanto reproduce el control DOCKEADO la huella del Cristalográfico.
        # Es solo verificación del docking; los compuestos se comparan con el cristalográfico, no con esta pose.
        crystal = meta.get("crystal_feats", {})
        if crystal:
            filas = []
            for R in sorted(crystal):
                cks = {ck for ck, rc in cassign.items() if rc == R} or ckeys
                sub = inter[(inter["receptor"] == R) & (inter["ckey"].isin(cks))]
                if sub.empty:
                    continue
                s = sub["name"].map(lambda n: dscore.get(sc.pose_key(n), float("nan")))
                best = sub.loc[s.idxmin()] if s.notna().any() else sub.iloc[0]
                pose_feats = [c for c in icols if best[c] > 0]
                rec = sc.fp_recovery(crystal[R], pose_feats)
                filas.append({"receptor": R, "recovery": rec["recovery"], "Tanimoto": rec["tanimoto"],
                              "reproduced": f"{rec['shared']}/{rec['ref_n']}", "extra (non-crystal)": rec["extra"]})
            if filas:
                st.markdown("**Interaction validation** — docked control vs. crystallographic ligand.")
                st.dataframe(pd.DataFrame(filas), width="stretch", hide_index=True)
                st.caption("`recovery` = fraction of the crystallographic interactions reproduced by the control's docked pose; `Tanimoto` also includes the extra contacts docking adds. ")

        st.markdown("**Weighting**")
        mw = meta.get("weights", {})
        metric_afin = st.radio(
            "Affinity-axis metric", ["dock", "le"], horizontal=True,
            index=1 if str(mw.get("dock_metric", "dock")).lower() == "le" else 0,
            format_func=lambda m: "Raw score (kcal/mol)" if m == "dock" else "Ligand efficiency (LE)",
            help="Vina's raw score favors large molecules (size bias). LE = -ΔG/heavy atoms corrects it. Recommended if your library varies widely in size; both columns are reported.")
        c1, c2, c3, c4 = st.columns(4)
        w_dock = c1.slider("Docking weight", 0.0, 1.0, float(mw.get("dock", 0.5)), 0.05)
        w_inter = c2.slider("Interactions weight", 0.0, 1.0, float(mw.get("inter", 0.5)), 0.05)
        w_adme = c3.slider("ADME weight", 0.0, 1.0, float(mw.get("adme", 0.0)), 0.05,
                           help="Physicochemical (drug-likeness) quality of the compound. Guards against rewarding only large/greasy molecules.")
        w_tox = c4.slider("Toxicity weight", 0.0, 1.0, float(mw.get("tox", 0.0)), 0.05,
                          help="Requires ADMET predicted (Ligands tab); otherwise this axis is ignored.")
        c5, c6, c7 = st.columns(3)
        w_cat = c5.slider("Catalytic-residue weight", 1.0, 6.0, float(mw.get("w_cat", 3.0)), 0.5,
                          help="How much an interaction with a catalytic (gate) residue is worth vs. an ordinary pocket one.")
        w_sec = c6.slider("Secondary-residue weight", 1.0, 3.0, float(mw.get("w_sec", 1.5)), 0.25,
                          help="How much an interaction with a SECONDARY anchor is worth vs. an ordinary pocket contact (×1).")
        cat_gate = c7.slider("Catalytic strictness", 0.0, 1.0, float(mw.get("cat_gate", 0.5)), 0.05,
                             help="0 = missing a catalytic residue is not penalized; 1 = missing all nullifies the score.")
        _axw = {"docking": w_dock, "interaction": w_inter, "ADME": w_adme, "tox": w_tox}
        _tot = sum(_axw.values())
        if _tot > 0:
            st.caption("Real contribution of each axis: "
                       + " · ".join(f"{k} {v / _tot * 100:.0f}%" for k, v in _axw.items() if v > 0))
        else:
            st.warning("All axis weights are 0: there will be no score. Raise at least one.")
        with st.expander("Weights by interaction type (advanced)"):
            st.caption("Merit value per type (0-1). Default: salt bridge > H-bond > π > halogen > hydrophobic. Literature-guided; adjust to your judgment.")
            st.caption("`water` (water-mediated bridges) only matters if you keep water molecules when preparing the receptor. In the usual flow they are removed, so this weight has no effect.")
            tw = {}; tcols = st.columns(4)
            for j, (tk, tv) in enumerate(sc.TYPE_WEIGHTS.items()):
                tw[tk] = tcols[j % 4].number_input(
                    tk, 0.0, 1.0, float((mw.get("type_weights") or {}).get(tk, tv)), 0.05, key=f"tw_{tk}")
        w = dict(pl.DEFAULT_WEIGHTS)
        w.update(dock=w_dock, inter=w_inter, adme=w_adme, tox=w_tox, dock_metric=metric_afin,
                 w_cat=w_cat, w_sec=w_sec, cat_gate=cat_gate, type_weights=tw)
        # Pesos que el usuario tiene ahora mismo en la interfaz. Al guardar la sesión se escriben en
        # run.json para que, al restaurar, los sliders arranquen con ellos y no haya que reponderar.
        S["_pesos_ui"] = w

        # SMILES para que se calculen ADME, SAscore, PAINS y Lipinski en la tabla
        smap = sc.build_smiles_map(str(proj / "ligandos_entrada"))
        for _d in {str(Path(c).parent) for c in meta.get("controls", [])}:  # SMILES de los controles
            for _k, _v in sc.build_smiles_map(_d).items():
                smap.setdefault(_k, _v)
        ml0 = proj / "ligands_meta.csv"
        if ml0.exists():
            _m0 = pd.read_csv(ml0)
            for _n, _s in zip(_m0.get("name", []), _m0.get("smiles", [])):
                k0 = sc.normalize_key(_n)
                if k0 and pd.notna(_s) and k0 not in smap:
                    smap[k0] = str(_s)
        # Confianza: estabilidad de pose precalculada (run.json) + fiabilidad de la diana (validación)
        pose_stab = meta.get("pose_stability", {})
        reliable_map = {}
        _vp = proj / "validacion_redocking.csv"
        if _vp.exists():
            _v = pd.read_csv(_vp)
            for _, _r in _v.iterrows():
                _rm = pd.to_numeric(pd.Series([_r.get("rmsd_min_A")]), errors="coerce").iloc[0]
                if pd.notna(_r.get("diana")):
                    reliable_map[str(_r["diana"])] = bool(pd.notna(_rm) and _rm < 2.0)
        rk = sc.compute_ranking(inter, dc, ckeys, cassign, ref_info, icols, dscore, cat, w,
                                smiles_map=smap, pocket_res_map=pocket_res_map, sec_map=sec,
                                pose_stability=pose_stab, reliable_map=reliable_map)
        rk["Ki"] = rk["pred_ki_M"].map(_fmt_ki) if "pred_ki_M" in rk.columns else None
        # Honestidad de ejes: uno con peso pero SIN datos no puntua; avisar para no declararlo en Methods.
        faltan = []
        if w_adme > 0 and pd.to_numeric(rk.get("adme"), errors="coerce").isna().all():
            faltan.append("ADME")
        if w_tox > 0 and pd.to_numeric(rk.get("ld50_mgkg"), errors="coerce").isna().all():
            faltan.append("toxicidad")
        if faltan:
            st.warning(f"You weight **{' and '.join(faltan)}** but there is no data for that axis in this run: it is ignored "
                       "in the score. Predict ADMET first, or lower its weight to 0 so Methods does not declare it.")

        with st.expander("Export Methods (for the paper)"):
            st.caption("Parameters, box, weights, reference and exact software versions. Reproducibility ready to paste into the Methods section.")
            metodos = rp.methods_text(meta, weights=w, catalytic=cat, secondary=sec)
            st.download_button("Download Methods.md", metodos, file_name="PoliScreen_Methods.md",
                               mime="text/markdown")
            st.code(metodos, language="markdown")

        # Motor por compuesto. Cuando hay péptidos, a unos ligandos los acopla ADCP y el control
        # sigue en Vina; sus energías proceden de funciones de puntuación distintas y no son
        # comparables entre si. La columna tiene que estar a la vista y la mezcla anunciada ANTES de
        # que nadie ordene la tabla por la energía y lea una diferencia que solo es de escala.
        motores = {}
        _dock_p = proj / "resultados_docking.csv"
        if _dock_p.exists():
            try:
                _dd = pd.read_csv(_dock_p)
                if "motor" in _dd.columns:
                    motores = {sc.normalize_key(c): m for c, m in
                               zip(_dd["compound_name"], _dd["motor"]) if pd.notna(m)}
            except Exception:
                motores = {}
        if motores:
            rk["motor"] = rk["compound"].map(lambda c: motores.get(sc.normalize_key(c), ""))
            _distintos = sorted({m for m in rk["motor"] if m})
            if len(_distintos) > 1:
                st.warning(
                    f"This table mixes **{' and '.join(_distintos)}**. Their energies come from "
                    "different scoring functions and are not comparable: `best_dock`, "
                    "`pKi` and `LE` only make sense within each engine. To compare across "
                    "them use `efectividad_pct`, computed from the contacts and independent "
                    "of the engine.")

        # Columnas fijas: las esenciales siempre visibles y el resto en el desplegable de abajo.
        # El selector manual sobraba: la tabla ya se puede ordenar y ampliar desde su propia cabecera.
        elegidas = [c for c in ("compound", "IUPAC", "motor", "pose", "best_dock", "pKi", "LE",
                                "best_inter", "cat_coverage", "efectividad_pct", "percentil",
                                "confidence", "consenso", "key_interaction", "sa_score", "pains",
                                "tipo")
                    if c in rk.columns]

        # procedencia (sombreado) y nombre IUPAC, leidos del metadata de ligandos
        meta_lig = proj / "ligands_meta.csv"
        tuyos, iupac_map, nombre_real = set(), {}, {}
        if meta_lig.exists():
            m = pd.read_csv(meta_lig)
            tuyos = {sc.normalize_key(n) for n, f in zip(m.get("name", []), m.get("fuente", []))
                     if str(f) == "tuyo"}
            if "iupac_name" in m.columns:
                iupac_map = {sc.normalize_key(n): v for n, v in zip(m.get("name", []), m.get("iupac_name", []))
                             if pd.notna(v) and str(v).strip()}
            # El nombre de las poses se pasa a minusculas para emparejar archivos, lo que estropea
            # los identificadores sensibles a mayusculas (una secuencia peptídica se escribe SIEMPRE
            # en mayusculas). Se recupera el nombre tal como se creo.
            nombre_real = {sc.normalize_key(n): n for n in m.get("name", []) if pd.notna(n)}
        if nombre_real:
            rk["compound"] = rk["compound"].map(lambda c: nombre_real.get(sc.normalize_key(c), c))
        if iupac_map:
            rk["IUPAC"] = rk["compound"].map(lambda c: iupac_map.get(sc.normalize_key(c), ""))
        if tuyos:
            rk["fuente"] = rk["compound"].map(lambda c: "tuyo" if sc.normalize_key(c) in tuyos else "")
            st.caption("Highlighted rows = compounds made with reagents you provided.")

        if any("~" in str(x) for x in rk["receptor"].unique()):
            st.info("**Hybrid docking**: each block is a different pocket of the same receptor. Compare a compound's effectiveness across sites to see where it prefers to bind.")
        # Con más de un receptor (dianas distintas, no bolsillos del mismo), un resumen pleiotropico
        # ANTES de los bloques por receptor: que compuesto es bueno en TODOS a la vez, que es lo que
        # buscaria quien quiere un ligando de amplio espectro y no uno específico de una diana.
        _dianas = sorted({sc.base_of(x) for x in rk["receptor"].unique()})
        if len(_dianas) > 1:
            _resumen_pleiotropico(rk, _dianas)
        for R in sorted(rk["receptor"].unique()):
            sub = rk[rk["receptor"] == R].copy()
            _et = f"{R.split('~')[0]} · sitio **{R.split('~', 1)[1]}**" if "~" in str(R) else f"**{R}**"
            _refsrc = meta.get("site_reference", {}).get(R) or ref_info.get(R, {}).get("src", "?")
            st.markdown(f"{_et} · interaction reference: `{_refsrc}`")
            noc = sub[sub["is_control"] != 1]
            if not noc.empty:
                m1, m2, m3, m4 = st.columns(4)
                try:
                    bd = noc.loc[pd.to_numeric(noc["best_dock"], errors="coerce").idxmin()]
                    m1.metric("Best docking", str(bd["compound"])[:18], f"{bd['best_dock']:.2f} kcal/mol",
                              delta_color="inverse")
                except Exception:
                    pass
                try:
                    bi = noc.loc[pd.to_numeric(noc["best_inter"], errors="coerce").idxmax()]
                    m2.metric("Best interaction quality", str(bi["compound"])[:18], f"{bi['best_inter']:.2f}")
                except Exception:
                    pass
                try:
                    be = noc.loc[pd.to_numeric(noc["efectividad_pct"], errors="coerce").idxmax()]
                    m3.metric("Best effectiveness", str(be["compound"])[:18], f"{be['efectividad_pct']:.0f} %")
                except Exception:
                    pass
                try:
                    bc = noc.loc[pd.to_numeric(noc["confidence"], errors="coerce").idxmax()]
                    m4.metric("Highest confidence", str(bc["compound"])[:18], f"{bc['confidence']:.2f}")
                except Exception:
                    pass
            vista = sub[elegidas]
            st.dataframe(_shade(vista.assign(fuente=sub.get("fuente", "")), "fuente") if tuyos else vista,
                         width="stretch", height=min(400, 60 + 34 * len(sub)))
            _descargar_tabla(vista, f"ranking_{R}", key=f"rk_{R}")
            g1, g2 = st.columns(2)
            ch = sub.dropna(subset=["efectividad_pct"]).set_index("compound")["efectividad_pct"]
            if not ch.empty:
                g1.bar_chart(ch, height=260)
            fig = _scatter_dock_inter(sub)
            if fig:
                g2.pyplot(fig)
        st.download_button("Download ranking (CSV)", rk.to_csv(index=False).encode(), "ranking.csv")

        # reporte ADMET sobre TODOS los ligandos del ranking (compuestos + nucleo + control), a elección
        # del usuario. smap ya incluye compuestos, nucleo y controles.
        items_all = [(c, smap[sc.normalize_key(c)]) for c in rk["compound"].unique()
                     if sc.normalize_key(c) in smap and pd.notna(smap.get(sc.normalize_key(c)))]
        if items_all:
            with st.expander("ADMET report (compounds + core + control, those you choose)"):
                nombres_all = [c for c, _ in items_all]
                elegidos_adme = st.multiselect("Which ligands to predict ADMET for?", nombres_all,
                                               default=nombres_all, key="adme_sel_res")
                items = [(c, s) for c, s in items_all if c in elegidos_adme]
                if st.button("Predict ADMET", key="pred_res") and items:
                    with st.spinner(f"Predicting with ADMET-AI for {len(items)} ligand(s)..."):
                        pr = AdmelabBridge().predict([s for _, s in items], use_ml=True)
                    S["admet"] = {**(S.get("admet") or {}), **{rg.inchikey(r.get("SMILES")): r for r in pr.rows}}
                if S.get("admet") and items:
                    _render_adme(S["admet"], items, keyp="res")

        st.markdown("---")
        st.markdown("**Interaction diagram** of a specific pose.")
        d1, d2, d3 = st.columns(3)
        R = d1.selectbox("Receptor", sorted(inter["receptor"].unique()))
        sr = inter[inter["receptor"] == R]
        cmp_ = d2.selectbox("Compound", sorted(sr["compound"].unique()))
        scmp = sr[sr["compound"] == cmp_]
        mods = sorted({sc.model_of(n) for n in scmp["name"]})
        mod = d3.selectbox("Pose", mods)
        row = scmp[scmp["name"].apply(lambda n: sc.model_of(n) == mod)]
        if not row.empty:
            referencia = ref_info.get(R, {}).get("feats", [])
            fig_int = sc.draw_2d(row.iloc[0], f"{R} · {cmp_} · pose {mod}", reference=referencia)
            # Tamaño nativo (moderado): estirarlo al ancho del panel lo hacía enorme; la figura ya se
            # dibuja pequeña para caber en pantalla sin desplazarse.
            st.pyplot(fig_int, use_container_width=False)
            try:
                import io as _io
                _b = _io.BytesIO(); fig_int.savefig(_b, format="png", dpi=160, bbox_inches="tight")
                _descargar_imagen(_b.getvalue(), f"interaccion_{cmp_}_pose{mod}", key=f"int_{R}_{cmp_}_{mod}")
            except Exception:
                pass
            st.caption("Green = reproduces a control interaction (same residue and same bond). Gray = extra contact or the same residue with a different bond type.")

        st.markdown("---")
        _como_citar()

# ---------------------------------------------------------------- visualizador (panel derecho)
def _alto_visor(reserva: int) -> int:
    """Altura del visor 3D para que el panel muestre selectores, visor y pie sin scroll. Se resta a la
    altura del panel (cfg_alto) el alto del cromo de cada etapa; el mínimo evita que quede diminuto si
    el panel es muy bajo."""
    return max(190, int(S.get("cfg_alto", 520)) - reserva)


def _visualizador(etapa: str):
    """Salida de la etapa activa. Se dibuja después del panel de herramientas, así que puede leer
    lo que este acaba de dejar en session_state (por ejemplo la caja de busqueda)."""
    if etapa == "Receptors":
        preparados = [p for p in S.get("receptors", []) if Path(p).exists()]
        if preparados:
            # Con varios receptores preparados, un selector para ver cada uno. Por defecto, el
            # último preparado; tras restaurar una sesión no existe, y se cae al último de la lista.
            _ult = S.get("ultimo_preparado")
            if len(preparados) > 1:
                _idx = preparados.index(_ult) if _ult in preparados else len(preparados) - 1
                rsel = st.selectbox("Receptor", preparados, index=_idx,
                                    format_func=lambda p: Path(p).name, key="vis_rec_sel")
            else:
                rsel = preparados[0]

            # El original de CADA receptor, no solo el del último preparado: es el archivo de partida
            # que queda junto al preparado (X.pdb frente a X_listo.pdb). El control es el ASIGNADO a
            # este receptor, no el último cargado, para que se vea el par correcto.
            _stem = Path(rsel).stem
            _orig = next((str(o) for o in [Path(rsel).with_name(_stem[:-6] + ".pdb")
                                           if _stem.endswith("_listo") else None,
                                           Path(S["ultimo_original"]) if rsel == _ult
                                           and S.get("ultimo_original") else None]
                          if o and o.exists()), None)
            _ctrls = [Path(c) for c in S["controls"]]
            _asig = pl._assign_controls(_ctrls, [Path(p) for p in preparados], S.get("_control_map") or {})
            _ctrl = next((str(c) for c in _ctrls
                          if _asig.get(sc.normalize_key(c.stem)) == _stem), None)

            c1, c2 = st.columns([3, 1])
            opciones = ["Prepared"] + (["Original"] if _orig else []) \
                + (["With its control"] if _ctrl else [])
            cual = c1.radio("View", opciones, horizontal=True, key="vis_ver_rec",
                            label_visibility="collapsed")
            ejes = c2.checkbox("XYZ axes", value=True, key="vis_ejes_rec")
            try:
                receptor = _orig if (cual == "Original" and _orig) else rsel
                ligando = _ctrl if cual == "With its control" else None
                _h = _alto_visor(250)      # selector + radio + posible pie del control
                components.html(vw.view_html(receptor=receptor, ligando=ligando,
                                             mostrar_aguas=False, ejes=ejes, alto=_h), height=_h + 12)
                if cual == "With its control" and _ctrl:
                    st.caption(f"Control of this receptor: `{Path(_ctrl).stem}`.")
            except Exception as e:
                st.error(f"Could not draw the structure: {e}")
        else:
            _vacio("Prepare a receptor and it will appear here in 3D.")

    elif etapa == "Ligands":
        prods = S.get("products")
        nuc_png = S.get("_nucleo_png")
        peps = S.get("_pep_preview")
        if peps and S.get("modo_ligandos") == "Generate peptides":
            st.markdown("**Generated sequences**")
            # La secuencia se colorea por clase de residuo: de un vistazo se ve si el péptido es
            # anfipatico (bloques hidrofóbicos y cationicos alternados) o uniforme.
            leyenda = [("#3d7ea6", "hydrophobic"), ("#b5453c", "charge +"),
                       ("#3f7d4e", "charge -"), ("#7a6ba8", "polar"), ("#8a8a8a", "G/P")]
            def _color(a):
                cls = pp.AMINOACIDOS[a][1]
                if "cargado_pos" in cls: return "#b5453c"
                if "cargado_neg" in cls: return "#3f7d4e"
                if "hidrofobico" in cls: return "#3d7ea6"
                if "especial" in cls:    return "#8a8a8a"
                return "#7a6ba8"
            html = []
            for nom, seq in peps:
                letras = "".join(
                    f"<span style='display:inline-block;width:1.35em;text-align:center;"
                    f"background:{_color(a)}22;color:{_color(a)};border-radius:3px;margin:1px;"
                    f"font-weight:600'>{a}</span>" for a in seq)
                html.append(f"<div style='margin:.35rem 0'><code style='opacity:.6'>{nom}</code> "
                            f"<span style='font-family:monospace;font-size:1.05rem'>{letras}</span></div>")
            st.markdown("".join(html), unsafe_allow_html=True)
            st.markdown(" ".join(f"<span style='color:{c};font-size:.8rem'>■ {n}</span>"
                                 for c, n in leyenda), unsafe_allow_html=True)
            # Rejilla con TODAS las estructuras: verificar una sola no dice nada de la biblioteca.
            # Se construyen con la MISMA química de extremos que se usara en el cribado; dibujarlas
            # siempre lineales haría invisible el efecto de ciclar, acetilar o amidar.
            _nac, _cam, _cic = S.get("_pep_quimica", (False, False, False))
            _firma_grid = (tuple(s for _, s in peps), _nac, _cam, _cic)
            if S.get("_pep_grid_firma") != _firma_grid:
                smis, etiquetas = [], []
                for nom, seq in peps:
                    s = pp.to_smiles(seq, n_acetil=_nac, c_amida=_cam, ciclico=_cic)
                    if s:
                        smis.append(s); etiquetas.append(nom)
                S["_pep_grid"] = vw.grid_png(smis, legends=etiquetas, cols=3, sub=250) if smis else None
                S["_pep_grid_n"] = len(smis)
                S["_pep_grid_firma"] = _firma_grid
            if S.get("_pep_grid"):
                quim = ("head-to-tail cyclized" if _cic else
                        ", ".join(filter(None, ["N-acetylated" if _nac else "",
                                                "C-amidated" if _cam else ""])) or "free termini")
                st.image(S["_pep_grid"],
                         caption=f"Structure of {S.get('_pep_grid_n', 0)} peptides · {quim}.")
        elif prods:
            png = vw.grid_png([p.get("smiles") for p in prods],
                              legends=[str(p.get("producto") or "") for p in prods])
            if png:
                st.image(png, caption=f"{len(prods)} products built. "
                                      "Check the ester bond and stereochemistry.")
        elif nuc_png:
            st.image(nuc_png, caption="Core with atom indices; in color, the reactive site.")
        elif S["ligands"]:
            st.success(f"{len(S['ligands'])} ligands ready to dock.")
            st.caption(", ".join(Path(p).stem for p in S["ligands"][:20]))
        else:
            _vacio("Build or upload ligands and you will see their structures here.")

    elif etapa == "Run":
        cajas = S.get("_boxes") or {}
        cav_map = S.get("_cavidades") or {}
        if cajas:
            c1, c2, c3 = st.columns([2, 1, 1])
            rsel = c1.selectbox("Receptor", list(cajas), format_func=lambda p: Path(p).name,
                                key="vis_box_rec", label_visibility="collapsed")
            # Cada receptor tiene sus propias cavidades; se dibujan las del que esta seleccionado.
            grupos_r = cav_map.get(rsel)
            ver_cav = c2.checkbox("Cavities", value=bool(grupos_r), key="vis_ver_cav")
            ejes = c3.checkbox("XYZ axes", value=True, key="vis_ejes_box")
            grupos = grupos_r if (ver_cav and grupos_r) else None
            try:
                _h = _alto_visor(210)      # fila de selector/casillas + pie de la caja
                components.html(vw.view_html(receptor=rsel, caja=cajas[rsel], cavidades=grupos,
                                             mostrar_aguas=False, ejes=ejes, alto=_h), height=_h + 12)
                b = cajas[rsel]
                st.caption(f"Caja (malva): centro ({b['cx']}, {b['cy']}, {b['cz']}) · "
                           f"{b['sx']} × {b['sy']} × {b['sz']} Å"
                           + (f" · {len(grupos)} cavities; the one used is highlighted." if grupos else ""))
            except Exception as e:
                st.error(f"No pude dibujar: {e}")
        else:
            _vacio("Define the search box and it will be drawn here over the receptor.")

    else:  # Resultados
        vista = st.radio("View", ["Summary", "3D complex"], horizontal=True,
                         key="vis_res_vista", label_visibility="collapsed")
        if vista == "3D complex":
            _visor_complejo()
        else:
            _resumen_visual()


def _visor_complejo():
    """Complejo receptor-ligando en 3D, recorriendo compuesto y pose.

    El diagrama 2D del panel de herramientas dice QUE contactos hay; esto dice DONDE. Se lee la
    misma tabla de interacciones y se usan las mismas claves, de modo que la pose elegida aquí y la
    del diagrama son la misma entidad y no dos numeraciones que haya que casar a mano.
    """
    # La MISMA carpeta normalizada que usa el panel de herramientas. Leer aquí el texto crudo del
    # campo (S['proj_dir']) hacia que, si la ruta se había normalizado, el visor buscara en otra
    # carpeta y mostrara el estado vacio mientras las herramientas si encontraban los resultados.
    proj_p = proj
    inter_p = proj_p / "interacciones.csv"
    if not inter_p.exists():
        _vacio("Run a screening and you can browse the 3D complexes here.")
        return
    t = pd.read_csv(inter_p)
    if "name" not in t.columns or t.empty:
        _vacio("The interactions table is empty.")
        return
    t = t.assign(_rec=t["name"].map(sc.receptor_from_name),
                 _cmp=t["name"].map(sc.compound_from_pose_name),
                 _mod=t["name"].map(sc.model_of))

    # Los controles se acoplan contra todos los receptores, pero cada uno pertenece a una sola
    # diana. Se ocultan de las demas para no ofrecer combinaciones control-receptor que no tienen
    # sentido; la asignación se lee de run.json, escrita en la corrida.
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
            return True                       # no es control: visible en todos
        return _asig.get(ck) == sc.base_of(rec)   # control: solo en su receptor

    c1, c2, c3 = st.columns([1.3, 1.7, 0.9])
    R = c1.selectbox("Receptor", sorted(t["_rec"].unique()), key="vis_cx_rec")
    sr = t[t["_rec"] == R]
    _cmps = sorted(c for c in sr["_cmp"].unique() if _visible(R, c))
    C = c2.selectbox("Compound", _cmps or sorted(sr["_cmp"].unique()), key="vis_cx_cmp")
    scmp = sr[sr["_cmp"] == C]
    M = c3.selectbox("Pose", sorted(scmp["_mod"].unique()), key="vis_cx_pose")
    o1, o2 = st.columns(2)
    sup = o1.checkbox("Show the surface", value=False, key="vis_cx_sup",
                      help="Translucent molecular surface of the receptor. With the ribbon alone you cannot tell whether the ligand is inside the cavity or resting outside.")
    het = o2.checkbox("Cofactors and hetero", value=True, key="vis_cx_het")
    fila = scmp[scmp["_mod"] == M]
    if fila.empty:
        _vacio("There is no pose with that combination.")
        return

    # El nombre de la tabla es el del complejo fusionado; la pose suelta es el mismo sin el prefijo.
    # Se prefiere la pose suelta porque así el ligando se dibuja aparte y en otro color; el complejo
    # fusionado sirve de reserva cuando la carpeta de poses se ha limpiado.
    nombre = str(fila.iloc[0]["name"])
    pose_f = proj_p / "poses" / f"{nombre.removeprefix('Complejo_')}.pdb"
    rec_f = next((p for p in (proj_p / "receptores").glob(f"{R}.*")
                  if p.suffix.lower() in (".pdb", ".pdbqt")), None)

    # Pie de una sola linea abajo; asi selectores, casillas, visor y resumen caben sin scroll.
    _h = _alto_visor(300)
    try:
        if pose_f.exists() and rec_f is not None:
            html = vw.view_html(receptor=rec_f, ligando=pose_f, mostrar_aguas=False,
                                mostrar_hetero=het, superficie=sup, alto=_h)
        else:
            fus = proj_p / "Complejos_Fusionados" / f"{nombre}.pdb"
            if not fus.exists():
                _vacio("Cannot find this pose's files in the project folder.")
                return
            html = vw.view_html(receptor=fus, mostrar_aguas=False, mostrar_hetero=het,
                                superficie=sup, alto=_h)
        components.html(html, height=_h + 12)
    except Exception as e:
        st.error(f"Could not draw the complex: {e}")
        return

    # Contactos de ESTA pose, para que lo que se ve en 3D y lo que dice la tabla coincidan. Todo en un
    # pie de una linea: la descripcion larga y la lista completa gastaban cuatro renglones y forzaban el
    # scroll. El detalle exhaustivo de contactos esta en el diagrama 2D del panel de herramientas.
    r0 = fila.iloc[0]
    feats = [c for c in fila.columns
             if "_" in str(c) and str(c).rsplit("_", 1)[-1] in sc.TYPE_STYLE
             and pd.notna(r0[c]) and r0[c] > 0]
    _cab = f"{C} · pose {M} · {R}"
    if feats:
        _cab += f" — **{len(feats)} contactos**"
    st.caption(_cab)


def _resumen_visual():
    """Resumen de resultados pensado para leerse de un vistazo: el detalle exhaustivo queda en el
    panel de herramientas; aquí va lo que se ensena a alguien en diez segundos."""
    rk_p = proj / "ranking.csv"          # misma carpeta normalizada que el panel de herramientas
    if not rk_p.exists():
        _vacio("Run a screening and the results summary will appear here.")
        return
    rk = pd.read_csv(rk_p)
    # Con docking hibrido hay un bloque por sitio: sin selector solo se veria el primero.
    sitios = sorted(rk["receptor"].unique()) if "receptor" in rk.columns else []
    if len(sitios) > 1:
        sel = st.selectbox("Site", sitios, key="vis_res_sitio",
                           format_func=lambda s: s.split("~", 1)[1] if "~" in str(s) else s)
        rk = rk[rk["receptor"] == sel]
        st.caption(f"Summary of site **{sel}**. Switch site to compare where each compound binds best.")
    # Sin la comprobación explícita, cuando falta la columna rk.get() devuelve un escalar y la
    # comparación produce un booleano suelto: pandas lo interpreta como nombre de columna y falla.
    noc = (rk[rk["is_control"] != 1] if "is_control" in rk.columns else rk).copy()
    if noc.empty or "efectividad_pct" not in noc.columns:
        st.info("No compounds to summarize yet.")
        return

    ef = pd.to_numeric(noc["efectividad_pct"], errors="coerce")
    conf = pd.to_numeric(noc.get("confidence"), errors="coerce") if "confidence" in noc else None
    superan = int((ef >= 105).sum())
    fiables = int((conf >= 0.5).sum()) if conf is not None else 0

    # Titular: el mejor compuesto por efectividad
    mejor = noc.loc[ef.idxmax()]
    enc, img = st.columns([2, 1])
    enc.markdown(f"### {str(mejor['compound'])[:38]}")
    enc.caption("Highest-effectiveness compound")
    # Estructura del ganador: pequeña, solo para reconocerlo de un vistazo.
    smi_top = None
    ml = proj / "ligands_meta.csv"
    if ml.exists():
        _m = pd.read_csv(ml)
        if {"name", "smiles"} <= set(_m.columns):
            _k = sc.normalize_key(mejor["compound"])
            _hit = _m[_m["name"].map(lambda n: sc.normalize_key(n) == _k)]
            if not _hit.empty:
                smi_top = str(_hit.iloc[0]["smiles"])
    if smi_top:
        png = vw.molecule_png(smi_top, size=170)
        if png:
            img.image(png)
    k1, k2, k3 = st.columns(3)
    k1.metric("Effectiveness", f"{ef.max():.0f} %")
    if "best_dock" in noc.columns:
        k2.metric("Affinity", f"{float(mejor['best_dock']):.1f}", "kcal/mol", delta_color="off")
    if conf is not None and pd.notna(mejor.get("confidence")):
        k3.metric("Confidence", f"{float(mejor['confidence']):.2f}")

    st.divider()
    c1, c2, c3 = st.columns(3)
    c1.metric("Compounds", len(noc))
    c2.metric("Beat the control", superan, help="Effectiveness ≥ 105 % vs. the crystallographic ligand.")
    if conf is not None:
        c3.metric("Confidence ≥ 0.5", fiables, help="Concordant evidence: reliable result.")

    # Podio: los cinco mejores, con barra proporcional para leerlos de un vistazo
    st.markdown("**Top five**")
    top = noc.assign(_ef=ef).nlargest(5, "_ef")
    tope = max(float(top["_ef"].max()), 1e-9)
    for i, (_, fila) in enumerate(top.iterrows(), 1):
        nombre = str(fila["compound"])[:34]
        val = float(fila["_ef"])
        cf = float(fila["confidence"]) if (conf is not None and pd.notna(fila.get("confidence"))) else None
        etiqueta = f"**{i}. {nombre}** — {val:.0f} %"
        if cf is not None:
            etiqueta += f" · confianza {cf:.2f}"
        # La efectividad puede ser negativa (un compuesto que reproduce menos contactos que la
        # penalización de referencia), y st.progress solo admite [0, 1]: un valor fuera de rango
        # lanzaba una excepción que tiraba el panel entero del visualizador.
        st.progress(max(0.0, min(val / tope, 1.0)), text=etiqueta)

    if conf is not None and fiables == 0:
        st.warning("No compound reaches confidence 0.5. With few poses the metric loses resolution: raise \"Poses per ligand\" in step 3 and run again.")

    # Descarga del resumen que se esta viendo, no de todo el ranking.
    cols_res = [c for c in ("compound", "best_dock", "pKi", "LE", "best_inter", "efectividad_pct",
                            "percentil", "confidence", "cnn_score", "consenso") if c in noc.columns]
    _descargar_tabla(noc[cols_res].sort_values("efectividad_pct", ascending=False),
                     "resumen_" + str(noc["receptor"].iloc[0]).replace("~", "_"), key="resumen_vis")


# ---------------------------------------------------------------- composicion de la pantalla
_ETAPA_FN = {"Receptors": _etapa_receptores, "Ligands": _etapa_ligandos,
             "Run": _etapa_ejecutar, "Results": _etapa_resultados}

# Altura fija en ambos paneles: cada uno hace su propio scroll y la pagina en conjunto no se
# desplaza, de modo que la cabecera y la barra de etapas quedan siempre a la vista.
_ALTO = int(S.get("cfg_alto", 520))
if S.get("_aviso"):
    st.success(S.pop("_aviso"))
_rep = float(S.get("cfg_reparto", 0.46))
_izq, _der = st.columns([_rep, 1.0 - _rep], gap="medium")
with _izq:
    st.markdown(f"**Tools · {S['etapa']}**")
    with st.container(height=_ALTO, border=True, key="panel_izq"):
        _ETAPA_FN[S["etapa"]]()
with _der:
    st.markdown("**Viewer**")
    with st.container(height=_ALTO, border=True, key="panel_der"):
        _visualizador(S["etapa"])

# Barra de etapas: botones planos, solo texto, repartidos por todo el ancho. Se dibuja al final,
# pero su valor ya esta en session_state al principio de la pasada, así que los paneles de arriba
# muestran la etapa correcta.
_nav = st.container(key="barra_etapas")
_cols = _nav.columns(len(ETAPAS))
for _i, _e in enumerate(ETAPAS):
    if _cols[_i].button(_e, key=f"nav_{_e}", use_container_width=True,
                        type=("primary" if _e == S["etapa"] else "secondary")):
        S["etapa"] = _e
        st.rerun()
