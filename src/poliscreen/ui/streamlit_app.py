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
    ax.set_xlabel("Docking (kcal/mol; más negativo = mejor)")
    ax.set_ylabel("Calidad de interacción (0-1 vs. control)")
    ax.invert_xaxis()
    ax.grid(alpha=0.3)
    ax.set_title("Docking vs. calidad · rojo = control · ideal: arriba-derecha")
    fig.tight_layout()
    return fig


def _render_adme(admet, items, keyp):
    """items: [(etiqueta, smiles)]. Muestra tabla resumen de todos + detalle por compuesto."""
    filas = []
    for lb, smi in items:
        r = admet.get(rg.inchikey(smi)) or {}
        filas.append({"compuesto": lb, "MW": r.get("MW"), "LogP": r.get("LogP"), "QED": r.get("QED"),
                      "LD50 (mg/kg)": r.get("LD50_mg_per_kg"), "GHS": r.get("GHS_category"),
                      "AMES": r.get("AMES"), "hERG": r.get("hERG"), "DILI": r.get("DILI")})
    st.markdown("**Resumen ADMET de todos los compuestos**")
    st.dataframe(pd.DataFrame(filas), width="stretch", height=min(320, 60 + 34 * len(filas)))
    st.caption("AMES/hERG/DILI = probabilidad de toxicidad (menor mejor). LD50 en mg/kg (mayor mejor). "
               "Se predice sobre la molécula COMPLETA (núcleo + reactivo), no el reactivo solo.")
    labels = dict(items)
    sel = st.selectbox("Ver detalle de", list(labels), key=f"adme_det_{keyp}")
    row = admet.get(rg.inchikey(labels[sel]))
    if not row:
        return
    ca, cb = st.columns([1, 1])
    ca.pyplot(rp.radar_fig(row, title=sel))
    cb.metric("LD50 oral (mg/kg)", rp._f(row.get("LD50_mg_per_kg"), 0))
    cb.metric("Categoría GHS", str(row.get("GHS_category") or "-"))
    cb.metric("QED", rp._f(row.get("QED")))
    cb.caption("Verde = favorable · ámbar = intermedio · rojo = desfavorable.")
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

**Cita también las herramientas que PoliScreen ejecuta.** No hacerlo es una omisión frecuente en
revisión por pares:

- **AutoDock Vina 1.2** — Eberhardt J, Santos-Martins D, Tillack AF, Forli S. *AutoDock Vina 1.2.0:
  New Docking Methods, Expanded Force Field, and Python Bindings.* J Chem Inf Model. 2021;61(8):3891–3898.
- **AutoDock Vina (original)** — Trott O, Olson AJ. J Comput Chem. 2010;31(2):455–461.
- **AutoDock CrankPep (ADCP)** (si acoplaste péptidos) — Zhang Y, Sanner MF. *Docking Flexible
  Cyclic Peptides with AutoDock CrankPep.* J Chem Theory Comput. 2019;15(10):5161–5168.
- **AGFR / AutoGridFR** (preparación de la diana para ADCP) — Zhang Y, Forli S, Omelchenko A,
  Sanner MF. *AutoGridFR.* J Comput Chem. 2019;40(32):2882–2891.
- **GNINA** (si re-puntuaste con la red neuronal) — McNutt AT, et al. *GNINA 1.0: molecular docking
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
- **ADMET-AI** (si usaste la predicción ADMET) — Swanson K, et al. Bioinformatics. 2024;40(7):btae416.

Las versiones exactas de cada herramienta se exportan con **Archivo → Métodos**, para que la
sección de Métodos del artículo sea reproducible."""

AGRADECIMIENTOS = """PoliScreen se apoya en estos proyectos de código abierto:

- **Streamlit** (Apache-2.0) — interfaz web
- **3Dmol.js / py3Dmol** (BSD-3-Clause) — visor molecular tridimensional
- **pandas** (BSD-3-Clause) y **NumPy** (BSD-3-Clause) — manejo de datos
- **Matplotlib** (licencia Matplotlib/PSF) — diagramas de interacción
- **OpenPyXL** (MIT) — exportación a XLSX
- **OPSIN** (MIT) — verificación de nombres IUPAC

Las herramientas científicas que ejecuta el motor (Vina, ADCP, gnina, PLIP, RDKit, Open Babel,
fpocket, OpenMM) tienen su cita completa en «Cómo citar»."""


def _como_citar(expandido: bool = False):
    with st.expander("Cómo citar", expanded=expandido):
        st.markdown(CITAS)
    with st.expander("Agradecimientos"):
        st.markdown(AGRADECIMIENTOS)


def _descargar_tabla(df, nombre: str, key: str):
    """Ofrece descargar un DataFrame como CSV o XLSX bajo la tabla."""
    import io
    c = st.columns([2, 1, 1])
    c[0].caption("Descargar como:")
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
    c[0].caption("Descargar como:")
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


@st.dialog("Confirmar borrado")
def _confirmar_borrado(carpeta):
    """Ventana modal antes de una acción destructiva: borrar resultados no tiene vuelta atras."""
    st.warning("Se borrarán las poses, los complejos, los XML de PLIP y todas las tablas de "
               "resultados de esta carpeta. **Esta acción no se puede deshacer.**")
    st.caption(f"Carpeta: `{carpeta}`")
    st.caption("Se conservan los receptores preparados, los controles y los ligandos de entrada.")
    st.info("Si quieres conservar este análisis, cancela y usa antes Archivo → Guardar sesión.")
    c1, c2 = st.columns(2)
    if c1.button("Sí, borrar", type="primary", use_container_width=True):
        pl.clean(carpeta)
        st.session_state["_aviso"] = "Resultados borrados. Se conservan receptores, controles y ligandos."
        st.rerun()
    if c2.button("Cancelar", use_container_width=True):
        st.rerun()


def _tam(b: int) -> str:
    return f"{b / 1e6:.1f} MB" if b >= 1e6 else (f"{b / 1e3:.0f} kB" if b else "—")


@st.dialog("Descargar resultados", width="large")
def _dialogo_descargas(carpeta):
    """Selector de exportaciones. Arma un único ZIP en memoria con lo que el usuario marque.

    La carpeta del proyecto ya contiene casi todos estos archivos; el paquete existe para llevarse
    el análisis a otra maquina, adjuntarlo a un articulo o archivarlo sin los intermedios pesados.
    Por eso cada elemento dice si vale la pena bajarlo y cuanto ocupa, en vez de ofrecerlos a ciegas.
    """
    cat = ss.catalogo(carpeta)
    hay = {k: v for k, v in cat.items() if v["hay"]}
    if not hay:
        st.info("Todavía no hay nada que exportar en esta carpeta.")
        return
    st.caption(f"Carpeta del proyecto: `{carpeta}`")
    st.caption("Todo esto ya está en esa carpeta, salvo la sección de Métodos, que se redacta al "
               "exportar. El paquete sirve para llevarte el análisis a otro sitio o adjuntarlo a "
               "un manuscrito sin arrastrar los archivos intermedios.")

    c1, c2 = st.columns(2)
    marcar = c1.button("Marcar los recomendados", use_container_width=True)
    limpiar = c2.button("Desmarcar todo", use_container_width=True)
    for k, v in hay.items(): 
        if marcar or limpiar:
            S[f"dl_{k}"] = bool(marcar and v["motivo"])
        else:
            S.setdefault(f"dl_{k}", bool(v["motivo"]))

    rec = {k: v for k, v in hay.items() if v["motivo"]}
    otros = {k: v for k, v in hay.items() if not v["motivo"]}
    if rec:
        st.markdown("**Recomendado**")
        for k, v in rec.items():
            st.checkbox(f"{v['desc']} · {_tam(v['bytes'])}", key=f"dl_{k}")
            st.caption(v["motivo"])
    if otros:
        st.markdown("**Opcional**")
        for k, v in otros.items():
            st.checkbox(f"{v['desc']} · {_tam(v['bytes'])}", key=f"dl_{k}")
            if v["regenerable"]:
                st.caption("Se regenera volviendo a ejecutar; ocupa espacio en el paquete.")

    elegidas = [k for k in hay if S.get(f"dl_{k}")]
    total = sum(hay[k]["bytes"] for k in elegidas)
    st.divider()
    if not elegidas:
        st.caption("No has marcado nada.")
        return
    st.caption(f"{len(elegidas)} elemento(s) · {_tam(total)} sin comprimir")
    if st.button("Preparar paquete", type="primary", use_container_width=True):
        try:
            meta = carpeta / "run.json"
            mtxt = rp.methods_text(json.loads(meta.read_text())) if meta.exists() else None
            datos, incluidos = ss.paquete(carpeta, elegidas, methods_text=mtxt)
            S["_zip"] = (f"{Path(carpeta).name}_PoliScreen.zip", datos, incluidos)
        except Exception as e:
            st.error(f"No pude armar el paquete: {e}")
    if S.get("_zip"):
        nombre, datos, incluidos = S["_zip"]
        st.download_button(f"Descargar {nombre} ({_tam(len(datos))})", datos,
                           file_name=nombre, mime="application/zip",
                           type="primary", use_container_width=True)
        st.caption("Contiene: " + ", ".join(incluidos))


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

ETAPAS = ["Receptores", "Ligandos", "Ejecutar", "Resultados"]
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
_menu_archivo = _m1.popover("Archivo", use_container_width=True)
_menu_datos = _m2.popover("Datos", use_container_width=True)
_menu_cfg = _m3.popover("Configuración", use_container_width=True)
_menu_ayuda = _m4.popover("Ayuda", use_container_width=True)

with _menu_archivo:
    st.markdown("**Proyecto**")
    S.setdefault("proj_dir", str(Path.home() / "poliscreen_proyectos" / "demo"))
    _escrito = st.text_input("Carpeta del proyecto", key="proj_dir",
                             help="Ruta dentro de Linux (WSL). Si pegas una ruta de Windows "
                                  "(`\\\\wsl.localhost\\...` o `C:\\...`) se traduce sola.")
    proj, _aviso_ruta = ss.normalizar_ruta(_escrito)
    if _aviso_ruta:
        # Se corrige el propio campo en la pasada siguiente para que lo guardado en la sesión y lo
        # que el usuario ve coincidan con la carpeta que se usa de verdad.
        st.warning(_aviso_ruta)
        S["_proj_pendiente"] = str(proj)
    try:
        proj.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        st.error(f"No puedo crear esa carpeta: {e}")
        st.stop()
    st.caption(f"Resultados en `{proj}`")
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
        st.caption(f"Recuperado del disco: {len(S['receptors'])} receptor(es), "
                   f"{len(S['controls'])} control(es), {len(S['ligands'])} ligando(s).")

    # --- Sesión de trabajo: abrir/guardar sin tener que escribir rutas ---
    with st.expander("Sesión y exportación"):
        sub = st.file_uploader("Abrir sesión (.poliscreen)", type=["poliscreen", "zip"],
                               help="Restaura un análisis anterior: tablas, receptores y ligandos. "
                                    "Puedes cambiar la ponderación sin repetir el docking.")
        if sub is not None and st.button("Restaurar esta sesión"):
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
                st.success(f"Sesión '{info.get('proyecto', '?')}' restaurada en {destino}")
                st.rerun()
            except Exception as e:
                st.error(f"No pude restaurar la sesión: {e}")
            finally:
                tmp_s.unlink(missing_ok=True)

        completa = st.checkbox("Incluir poses y complejos (sesión pesada)", value=False,
                               help="Sin marcar, la sesión pesa unos pocos MB y basta para reabrir "
                                    "y repuntuar. Marcada, permite además reexaminar estructuras 3D.")
        if st.button("Guardar sesión"):
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
                st.error(f"No pude guardar: {e}")
        if S.get("_sesion"):
            st.download_button(f"Descargar {S['_sesion'][0]}", S["_sesion"][1],
                               file_name=S["_sesion"][0], mime="application/zip")

        _n_exp = sum(1 for v in ss.catalogo(proj).values() if v["hay"])
        if _n_exp:
            # El modal se abre desde el flujo principal, no desde dentro del popover: un dialogo
            # invocado dentro de un contenedor desplegable se dibuja dentro de el y queda recortado.
            if st.button(f"Descargar resultados… ({_n_exp} disponibles)",
                         use_container_width=True):
                S["_abrir_descargas"] = True
                st.rerun()
            st.caption("Eliges con casillas qué incluir y se arma un solo ZIP; nada se escribe "
                       "en la carpeta del proyecto.")
with _menu_datos:
    st.markdown("**Cargado en este proyecto**")
    st.write(f"Receptores preparados: **{len(S['receptors'])}**")
    st.write(f"Controles co-cristalizados: **{len(S['controls'])}**")
    st.write(f"Ligandos: **{len(S['ligands'])}**")
    st.caption("Se detectan solos desde la carpeta del proyecto; sobreviven a reiniciar la aplicación.")

with _menu_ayuda:
    st.markdown("**Manual de PoliScreen**")
    st.caption("Cada sección se despliega con el detalle completo.")
    for _sec, _temas in ayuda.SECCIONES.items():
        with st.expander(_sec):
            for _titulo, _cuerpo in _temas:
                st.markdown(f"**{_titulo}**")
                st.markdown(_cuerpo)
                st.markdown("")
    _como_citar()

with _menu_cfg:
    st.markdown("**Apariencia**")
    st.caption("El tema claro/oscuro se cambia en el menú de la esquina superior derecha (⋮) → "
               "Settings → Theme. Puedes dejarlo en «Use system setting» para que siga al sistema.")
    st.slider("Reparto entre herramientas y visualizador", 0.3, 0.7, 0.46, 0.02,
              key="cfg_reparto",
              help="Hacia la izquierda, más espacio para el visualizador; hacia la derecha, para las "
                   "herramientas.")
    st.slider("Alto de los paneles (px)", 380, 900, 520, 20, key="cfg_alto")
    st.caption("Los parámetros de docking están en el paso 3 (Ejecutar), junto al botón de lanzar.")

if S.pop("_abrir_descargas", False):
    _dialogo_descargas(proj)


def _quimiotipos_tanda():
    """(hay_vina, hay_peptidos) presentes en la tanda. Solo decide QUE ajustes mostrar; el enrutado
    real lo hace el pipeline por ligando. Errar aquí como mucho muestra una sección de ajustes de
    más, nunca envia un ligando al motor equivocado."""
    peps = S.get("modo_ligandos") == "Generar péptidos" or bool(S.get("pep_seqs"))
    vina = False
    for c in S.get("controls", []):
        cl = str(c).lower()
        if cl.endswith(".pdb") and pp.secuencia_de_estructura(c):
            peps = True
        elif cl.endswith((".sdf", ".mol2", ".mol")):
            vina = True
    if S.get("products") or (S.get("modo_ligandos") in ("Construir por reacción", "Subir ligandos listos")
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

    with st.expander("Ajustes avanzados de docking"):
        st.caption("Los valores por defecto sirven para una primera exploración. Súbelos para "
                   "un cribado definitivo.")
        if usa_vina and usa_adcp:
            st.info("Tanda **mixta**: las moléculas pequeñas se acoplan con **Vina** y los péptidos "
                    "con **ADCP**. Abajo aparecen los ajustes de cada motor por separado.")
        elif usa_adcp:
            st.info("Cribado de **péptidos** con **ADCP**. Se muestran solo sus ajustes; los de Vina "
                    "no intervienen.")
        elif hay_peps and not adcp_ok:
            st.warning("Hay péptidos y ADCP no está instalado: se acoplarán con **Vina**, cuyo "
                       "muestreo no cubre esa flexibilidad. Instálalo con scripts/get_adcp.sh.")

        seed = st.number_input("Semilla", value=42, step=1,
                               help="Fija el azar: la misma semilla da el mismo resultado en ambos motores.")
        n_poses = st.slider("Poses por ligando", 1, 20, 5,
                            help="Con menos de 3 la métrica de confianza pierde resolución.")

        if usa_vina:
            st.markdown("**Vina** — moléculas pequeñas")
            exhaust = st.slider("Exhaustividad", 8, 64, 8, 8,
                                help="Más alto = búsqueda más fina y más lenta.")
            energy_range = st.slider("Rango de energía (kcal/mol)", 1.0, 8.0, 3.0, 0.5,
                                     help="Ventana de energía respecto a la mejor pose para reportar modos alternativos.")
            ph = st.slider("pH de protonación", 5.0, 9.0, 7.4, 0.1,
                           help="pH al que OpenBabel protona antes de acoplar (fisiológico ≈ 7.4).")
            cpu = st.number_input("Hilos por acoplamiento", 1, 16, 1,
                                  help="1 mantiene el resultado reproducible. Súbelo solo si no te importa.")
            workers = st.number_input("Acoplamientos en paralelo (0 = automático)", 0, 32, 0)
            if cpu > 1:
                st.warning("Con más de un hilo por acoplamiento, Vina deja de ser determinista.")

        if usa_adcp:
            st.markdown("**ADCP** — péptidos")
            st.caption("Usa automáticamente los núcleos de la máquina y es reproducible con la "
                       "semilla; los ajustes de hilos de Vina no le afectan.")
            adcp_pasos = st.select_slider(
                "Pasos por réplica", [50_000, 100_000, 250_000, 500_000, 1_000_000, 2_500_000],
                value=250_000, format_func=lambda v: f"{v // 1000} k",
                help="Longitud de cada búsqueda. Súbelo si el control no recupera su postura o si "
                     "la energía sigue mejorando al aumentarlo.")
            adcp_reps = st.slider("Réplicas independientes", 4, 100, 20, 2,
                                  help="Búsquedas paralelas desde puntos distintos. Más réplicas "
                                       "reducen la probabilidad de quedarse en un mínimo local.")
            st.caption(f"Coste aproximado: {adcp_pasos * adcp_reps / 5e6:.0f}× el de un octapéptido "
                       "con los valores por defecto (~35 s con seis hilos).")

        st.markdown("**Segunda opinión (red neuronal)**")
        hay_gnina = dk.gnina_available()
        # La red de gnina se entreno con complejos proteína-molécula pequeña. Con péptidos sigue
        # ejecutandose y devuelve numeros, que es justo el peligro: son extrapolacion fuera de su
        # dominio de entrenamiento. Se avisa en vez de bloquear.
        rescnn = st.checkbox("Re-puntuar las poses con gnina (CNN, GPU)", value=False,
                             disabled=not hay_gnina,
                             help="No vuelve a acoplar: mantiene las poses y las evalúa con una red "
                                  "neuronal entrenada sobre complejos cristalográficos. Añade una "
                                  "evidencia independiente a la métrica de confianza.")
        if not hay_gnina:
            st.caption("gnina no está instalado. Es opcional: sin él, la confianza usa las otras evidencias.")
        elif rescnn and hay_peps:
            st.warning("Estás cribando péptidos. La red de gnina se entrenó con complejos de "
                       "molécula pequeña, así que aquí puntúa fuera de su dominio: sus valores "
                       "bajos no significan necesariamente que la pose sea mala. Úsala para "
                       "comparar, no como criterio, y decláralo en Métodos.")
        elif rescnn:
            st.caption("Se re-puntúa la mejor pose de cada compuesto (~2 s por compuesto).")
    return dict(seed=int(seed), exhaustiveness=int(exhaust), n_poses=int(n_poses),
                energy_range=float(energy_range), ph=float(ph), cpu=int(cpu), workers=int(workers),
                rescoring_cnn=bool(rescnn), adcp_pasos=int(adcp_pasos),
                adcp_replicas=int(adcp_reps))


# Cada etapa es una función: solo se ejecuta la activa, no las cuatro como ocurria con st.tabs.
# ---------------------------------------------------------------- 1. receptor
def _etapa_receptores():
    st.subheader("Preparar un receptor")
    st.caption("Escribe un identificador del PDB o sube tu propio archivo. Se quitan las aguas, "
               "se añaden hidrógenos y se conserva la numeración original de los residuos.")
    c1, c2 = st.columns([1, 2])
    pdb_id = c1.text_input("Identificador del PDB", placeholder="4D44")
    up = c2.file_uploader("...o sube un archivo .pdb", type=["pdb"])

    src = None
    if up is not None:
        src = proj / "receptores" / up.name
        src.parent.mkdir(parents=True, exist_ok=True)
        src.write_bytes(up.getvalue())
    elif pdb_id.strip() and st.button("Descargar del PDB"):
        try:
            src = rc.fetch_pdb(pdb_id, proj / "receptores")
            S["src_pdb"] = str(src)
        except rc.ReceptorError as e:
            st.error(str(e))
    if src is None and S.get("src_pdb"):
        src = Path(S["src_pdb"])

    if src and src.exists():
        st.success(f"Estructura cargada: {src.name}")
        info = rc.inspect(src)
        st.write(f"**{info.n_atoms}** átomos · cadenas **{', '.join(info.chains)}** · **{info.n_waters}** aguas")
        if info.het:
            st.dataframe(pd.DataFrame([{"grupo": h.resname, "cadena": h.chain, "número": h.resseq,
                                        "átomos": h.n_atoms, "clave": h.key} for h in info.het]),
                         width="stretch", height=220)
        # Claves explicitas ligadas al archivo: sin ellas la selección de cadenas, cofactores y
        # control se pierde al cambiar de etapa, porque Streamlit descarta los widgets no dibujados.
        kb = sc.normalize_key(src.stem)
        c1, c2, c3 = st.columns(3)
        chains = c1.multiselect("Cadenas a conservar", info.chains, default=info.chains[:1],
                                key=f"rec_chains_{kb}")
        # Un residuo modificado (fosfotirosina, selenometionina...) el PDB lo declara como
        # heteroatomo, pero no es un cofactor: pertenece a la cadena. Ofrecerlo entre los cofactores
        # invitaba a conservarlo por esa vía, que lo duplicaba sobre el que ya esta en la proteína.
        mods = rc.modified_residues(src)
        claves_mod = {(m.chain, m.resseq) for m in mods}
        keys = [h.key for h in info.het if (h.chain, str(h.resseq).strip()) not in claves_mod]
        keep = c2.multiselect("Conservar (cofactores)", keys, key=f"rec_keep_{kb}",
                              help="Un cofactor del sitio, por ejemplo NADP.")
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
                return f"Cadena {c} · {n} residuos (péptido)"
            return f"{o} (heterogrupo)"

        _sel_ctrl = c3.multiselect("Extraer como control", _opc_ctrl, key=f"rec_extract_{kb}",
                                   format_func=_fmt_ctrl,
                                   help="El ligando co-cristalizado que define la huella de "
                                        "referencia. Puede ser un heterogrupo o una cadena "
                                        "peptídica; ambos aparecen aquí.")
        extract = [o for o in _sel_ctrl if not o.startswith("chain:")]
        cad_ctrl = [o.split(":", 1)[1] for o in _sel_ctrl if o.startswith("chain:")]
        smiles = st.text_input("SMILES del ligando extraído (opcional)", key=f"rec_smiles_{kb}",
                               help="Corrige los órdenes de enlace, que el PDB no guarda.")
        keep_mod = []
        if mods:
            st.markdown("**Residuos modificados de la cadena**")
            st.caption("Detectados en la estructura. Marcados, se conservan con su modificación; "
                       "desmarcados, se sustituyen por el aminoácido del que derivan y **se pierde "
                       "la modificación** — que a menudo es la función, como en un bucle de "
                       "activación fosforilado.")
            keep_mod = st.multiselect("Conservar con su modificación",
                                      [m.key for m in mods],
                                      default=[m.key for m in mods],
                                      key=f"rec_mod_{kb}",
                                      format_func=lambda k: next(
                                          (m.label for m in mods if m.key == k), k))
        for _c in cad_ctrl:
            _n = sum(1 for _l in src.read_text(errors="ignore").splitlines()
                     if _l.startswith("ATOM") and _l[21] == _c and _l[12:16].strip() == "CA")
            if _n > pp.MAX_LARGO:
                st.warning(f"La cadena {_c} tiene {_n} residuos: demasiado larga para tratarla "
                           "como ligando de referencia.")
        firma_prep = (str(src), tuple(chains), tuple(keep), tuple(extract), tuple(cad_ctrl),
                      tuple(keep_mod), smiles)
        prep_hecho = _ya_hecho("prep_" + kb, firma_prep)
        if prep_hecho:
            st.caption("Receptor ya preparado con esta selección. Cambia algo para volver a prepararlo.")
        if st.button("Preparar receptor", type="primary", disabled=prep_hecho):
            with st.spinner("Preparando..."):
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
                        st.info(f"Cadena {c} extraída como control: `{_seq[0]}` "
                                f"({len(_seq[0])} residuos)"
                                + (", se acoplará con ADCP." if adcp.available()
                                   and adcp.MIN_RESIDUOS <= len(_seq[0]) <= adcp.MAX_RESIDUOS
                                   else "."))
            _marcar_hecho("prep_" + kb, firma_prep)
            st.success(f"Listo: {dest.name}")
            st.rerun()

    # Comprobación visual: que se fue y que se quedo
    if S.get("ultimo_preparado") and Path(S["ultimo_preparado"]).exists():
        st.markdown("---")
        st.subheader("Comprobación de la preparación")
        antes = vw.resumen_estructura(S["ultimo_original"])
        despues = vw.resumen_estructura(S["ultimo_preparado"])
        # Todo como texto: la columna mezcla recuentos y listas de cadenas, y Arrow no puede
        # serializar una columna con enteros y cadenas a la vez.
        comp = pd.DataFrame([
            {"": "Átomos", "antes": str(antes["atomos"]), "después": str(despues["atomos"])},
            {"": "Hidrógenos", "antes": str(antes["hidrogenos"]), "después": str(despues["hidrogenos"])},
            {"": "Aguas", "antes": str(antes["aguas"]), "después": str(despues["aguas"])},
            {"": "Cadenas", "antes": ", ".join(antes["cadenas"]), "después": ", ".join(despues["cadenas"])},
            {"": "Heterogrupos", "antes": ", ".join(sorted(antes["heterogrupos"])) or "-",
             "después": ", ".join(sorted(despues["heterogrupos"])) or "-"},
        ])
        st.dataframe(comp, width="stretch", hide_index=True)
        if despues["aguas"] == 0 and despues["hidrogenos"] > 0:
            st.success("Sin aguas y con hidrógenos añadidos.")
        else:
            st.warning("Revisa: deberían quedar 0 aguas y haber hidrógenos.")

        st.caption("La estructura se muestra en el panel derecho; ahí puedes cambiar la vista y el estilo.")

    if S["receptors"]:
        st.write("**Receptores preparados:**", ", ".join(Path(p).name for p in S["receptors"]))
    if S["controls"]:
        st.write("**Controles:**", ", ".join(Path(p).name for p in S["controls"]))

# ---------------------------------------------------------------- 2b. péptidos
def _modo_peptidos():
    """Diseño de péptidos: vía independiente de la síntesis por reacción. Se mantiene aparte
    (S['péptidos']) para que no se mezcle con los productos del constructor químico."""
    st.caption("Los péptidos no se someten a reacciones químicas: se construyen directamente a "
               "partir de la secuencia. Entre 1 y 20 residuos.")
    entrada = st.radio("Cómo obtener las secuencias", ["Generar biblioteca", "Escribir secuencias"],
                       horizontal=True, key="pep_entrada")

    secuencias, aviso, problemas = [], "", []
    if entrada == "Escribir secuencias":
        txt = st.text_area("Una secuencia por línea, en código de una letra",
                           placeholder="KWKLFKKI\nGIGKFLHSAK\nRRWWRF", height=130, key="pep_txt")
        crudas = [s.strip().upper() for s in txt.splitlines() if s.strip()]
        malas = []
        for s in crudas:
            fuera = set(s) - set(pp.AMINOACIDOS)
            if fuera:
                malas.append(f"{s} (símbolos no válidos: {', '.join(sorted(fuera))})")
            elif not (pp.MIN_LONGITUD <= len(s) <= pp.MAX_LONGITUD):
                malas.append(f"{s} (longitud {len(s)}; el máximo es {pp.MAX_LONGITUD})")
            else:
                secuencias.append(s)
        if malas:
            st.warning("Se ignoran estas líneas: " + " · ".join(malas[:6]))
    else:
        c1, c2, c3 = st.columns(3)
        largo = c1.number_input("Residuos por péptido", pp.MIN_LONGITUD, pp.MAX_LONGITUD, 7, key="pep_len")
        cuantos = c2.number_input("Cuántos péptidos", 1, 2000, 50, key="pep_n")
        semilla = c3.number_input("Semilla", value=42, step=1, key="pep_seed",
                                  help="Misma semilla y mismas reglas = misma biblioteca.")
        with st.expander("Composición: qué aminoácidos puede usar", expanded=True):
            clases = st.multiselect("Clases permitidas (vacío = los 20)", list(pp.CLASES),
                                    format_func=lambda k: pp.CLASES[k], key="pep_clases")
            excl = st.multiselect("Excluir residuos concretos", sorted(pp.AMINOACIDOS),
                                  format_func=lambda a: f"{a} · {pp.AMINOACIDOS[a][0]}", key="pep_excl")
            alf = pp.alfabeto(incluir=clases, excluir_residuos=excl)
            st.caption(f"Alfabeto resultante ({len(alf)}): {', '.join(alf) if alf else 'vacío'}")
        with st.expander("Reglas de secuencia"):
            r1, r2 = st.columns(2)
            sin_rep = r1.checkbox("Sin residuos repetidos", key="pep_sinrep")
            maxcons = r1.number_input("Máx. iguales seguidos (0 = sin límite)", 0, 10, 0, key="pep_cons")
            maxres = r2.number_input("Máx. veces por residuo (0 = sin límite)", 0, 20, 0, key="pep_maxres")
            pre = r2.text_input("Empieza por", key="pep_pre", placeholder="p. ej. KK").upper()
            suf = r1.text_input("Termina en", key="pep_suf", placeholder="p. ej. GG").upper()
        with st.expander("Filtros fisicoquímicos"):
            st.caption("En péptidos antimicrobianos, la carga neta positiva y una hidrofobicidad "
                       "moderada son los rasgos más asociados a la actividad.")
            f1, f2 = st.columns(2)
            usar_q = f1.checkbox("Filtrar por carga neta", key="pep_usaq")
            q_rng = f1.slider("Carga neta a pH 7.4", -10.0, 10.0, (2.0, 9.0), 0.5,
                              key="pep_q", disabled=not usar_q)
            usar_g = f2.checkbox("Filtrar por hidropatía (GRAVY)", key="pep_usag")
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
    st.markdown("**Química de los extremos**")
    e1, e2, e3 = st.columns(3)
    n_ac = e1.checkbox("Acetilar extremo N", key="pep_nac",
                       help="Protege frente a aminopeptidasas.")
    c_am = e2.checkbox("Amidar extremo C", key="pep_cam",
                       help="Elimina la carga negativa terminal: +1 en la carga neta, lo que "
                            "suele aumentar la actividad antimicrobiana.")
    ciclo = e3.checkbox("Ciclar cabeza-cola", key="pep_ciclo",
                        help="Rigidiza el péptido y reduce mucho los grados de libertad, "
                             "lo que además hace más fiable el acoplamiento.")

    if entrada == "Generar biblioteca" and not problemas:
        st.caption(f"Espacio combinatorio disponible: ~{reglas.espacio():.0f} secuencias.")
        # La química de los extremos NO entra aquí: no cambia las secuencias, que es lo único que
        # produce este boton. Quien depende de ella es la construcción de los compuestos, y es ese
        # boton el que se rehabilita al cambiarla.
        firma = (int(largo), int(cuantos), int(semilla), tuple(alf), sin_rep, int(maxcons),
                 int(maxres), pre, suf, usar_q, q_rng, usar_g, g_rng)
        hecho = S.get("_pep_firma") == firma and S.get("pep_seqs")
        # Este boton NO se bloquea: generar es barato y determinista con la semilla, así que volver a
        # pulsarlo no tiene coste ni riesgo. El bloqueo se reserva para construir las estructuras 3D,
        # que es la operación cara. Aquí solo se informa de si la biblioteca ya esta al dia.
        if st.button("Generar biblioteca", type="primary"):
            with st.spinner("Generando secuencias..."):
                secuencias, aviso = pp.generate(reglas, int(cuantos), seed=int(semilla))
            S["pep_seqs"] = secuencias
            S["pep_aviso"] = aviso
            S["_pep_firma"] = firma
        if hecho:
            st.caption(f"Biblioteca generada con estos parámetros ({len(S['pep_seqs'])} secuencias).")
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
    st.caption("`momento_hidrofobico` mide la anfipaticidad (cara hidrofóbica frente a cara polar "
               "al plegarse en hélice); `indice_boman` estima la tendencia a unirse a otras "
               "proteínas: por encima de 2.5 kcal/mol se considera promiscuo.")

    nivel, msg = pp.viabilidad_docking(int(df["longitud"].max()), n_peptidos=len(df),
                                       hay_adcp=adcp.available())
    (st.success if nivel == "bueno" else st.warning if nivel == "medio" else st.error)(
        f"**Acoplamiento de {len(df)} péptidos de {int(df['longitud'].max())} residuos:** {msg}")

    S["_pep_preview"] = [(f["nombre"], f["secuencia"]) for f in filas[:24]]
    # El visualizador construye las estructuras y necesita la misma química de extremos que se
    # usara en el cribado; si no, dibuja el péptido lineal sin proteger y el ciclo no se ve.
    S["_pep_quimica"] = (bool(n_ac), bool(c_am), bool(ciclo))
    _firma_pep = (tuple(secuencias), n_ac, c_am, ciclo)
    _pep_listo = _ya_hecho("usar_peptidos", _firma_pep)
    if st.button("Usar estos péptidos en el cribado", type="primary", disabled=_pep_listo,
                 help="Cambia las secuencias o la química de los extremos para volver a construirlos."
                      if _pep_listo else None):
        with st.spinner(f"Construyendo la estructura 3D de {len(secuencias)} péptidos..."):
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
        st.toast(f"{len(made)} péptidos construidos y listos para el paso 3.", icon="🧬")
        st.success(f"{len(made)} péptidos listos para el paso 3."
                   + (f" {fallos} no se pudieron construir." if fallos else ""))
        if perdidos:
            st.warning("No se pudo generar la estructura 3D de: " + ", ".join(perdidos)
                       + ". Son cadenas largas y flexibles; prueba a ciclarlas para rigidizarlas.")
    if _pep_listo:
        st.caption(f"{len(S['ligands'])} péptidos construidos con estos parámetros.")


# ---------------------------------------------------------------- 2. ligandos
def _etapa_ligandos():
    st.subheader("¿Qué quieres acoplar?")
    modo = st.radio("Origen de los compuestos",
                    ["Construir por reacción", "Generar péptidos", "Subir ligandos listos"],
                    horizontal=True, key="modo_ligandos")

    if modo == "Generar péptidos":
        _modo_peptidos()
    elif modo == "Construir por reacción":
        S["lead"] = None
        izq, der = st.columns(2)
        with izq:
            st.markdown("#### Núcleo (tu molécula de partida)")
            nuc = st.text_input("SMILES del núcleo", value=S.get("nuc_smiles", ""),
                                placeholder="O=C(O)c1ccc2[n+]([O-])onc2c1")
            fnuc = st.file_uploader("...o archivo del núcleo", type=["sdf", "mol2", "mol"])
            if fnuc is not None:
                d = proj / "nucleo"
                d.mkdir(parents=True, exist_ok=True)
                fp = d / fnuc.name
                fp.write_bytes(fnuc.getvalue())
                nuc = _to_smiles(fp) or nuc
            S["nuc_smiles"] = nuc
            rxkey = st.selectbox("Reacción", list(rx.REACTIONS), key="rx_reaccion",
                                 format_func=lambda k: rx.get(k).nombre)
            reaction = rx.get(rxkey)
            st.caption(reaction.descripcion)
            if nuc:
                aplica = any(r.key == rxkey for r in rx.applicable(nuc))
                if not aplica:
                    st.warning(f"El núcleo no tiene {reaction.lead_grupo}; esta reacción no aplica.")
                else:
                    sitios = rx.lead_sites(nuc, reaction)
                    st.success(f"El núcleo puede sufrir {reaction.nombre}: {len(sitios)} sitio(s) reactivo(s).")
                    idx = 0
                    if len(sitios) > 1:
                        idx = st.selectbox("Punto de crecimiento", range(len(sitios)),
                                           format_func=lambda i: f"átomos {sitios[i]['atomos']}")
                    hl = sitios[idx]["atomos"] if sitios else []
                    # Se dibuja en el panel derecho, que es el visualizador.
                    S["_nucleo_png"] = vw.molecule_png_indexed(nuc, highlight=hl, size=420)
        with der:
            if reaction.kind == "coupling":
                st.markdown("#### Reactivos que se unen")
                usar_int = st.checkbox("Biblioteca interna", value=True,
                                       help=f"{len(rg.load_internal(reaction)) if nuc else 0} reactivos curados.")
                ups = st.file_uploader("Tus reactivos (csv/xlsx con columnas name y smiles · sdf · mol2 · smi)",
                                       type=["csv", "xlsx", "xls", "sdf", "mol2", "mol", "smi"],
                                       accept_multiple_files=True)
                with st.expander("¿Qué columnas debe tener mi Excel/CSV?"):
                    st.markdown(
                        "Dos columnas: una de **nombre** y una de **SMILES**. Las cabeceras aceptadas son:\n"
                        "- Nombre: `name`, `nombre`, `compound`, `compuesto`, `Alcohol origen`, `Nombre clave`\n"
                        "- SMILES: `smiles`, `smile`, `SMILES alcohol`\n\n"
                        "Se deduplica por estructura (InChIKey) y se filtra a los que tienen el grupo de la "
                        "reacción (para esterificación, un OH de alcohol/fenol; los ácidos y aminas se descartan).")
                    st.dataframe(pd.DataFrame({"name": ["Bencílico", "Mentol", "Ciclohexanol"],
                                               "smiles": ["OCc1ccccc1", "CC(C)C1CCC(C)CC1O", "OC1CCCCC1"]}),
                                 width="stretch", hide_index=True)
                usar_pc = st.checkbox("Complementar con PubChem (experimental, requiere internet)", value=False)
                pc_max = st.number_input("Máximo de PubChem", 5, 100, 25) if usar_pc else 25
                upaths = []
                if ups:
                    d = proj / "reactivos"; d.mkdir(parents=True, exist_ok=True)
                    for u in ups:
                        (d / u.name).write_bytes(u.getvalue()); upaths.append(str(d / u.name))
                if st.button("Reunir reactivos", type="primary"):
                    with st.spinner("Reuniendo y deduplicando..."):
                        reags, info = rg.build(reaction, use_internal=usar_int, user_paths=upaths,
                                               use_pubchem=usar_pc, pubchem_max=int(pc_max))
                    S["reagents"] = [(r.name, r.smiles, r.source, r.inchikey) for r in reags]
                    S["reag_info"] = info
                    if info.get("aviso_pubchem"):
                        st.warning(info["aviso_pubchem"])
                if S.get("reag_info"):
                    info = S["reag_info"]
                    st.write(f"**{info['total']} reactivos** — " +
                             " · ".join(f"{k}: {v}" for k, v in info["por_fuente"].items()))
                    dfa = pd.DataFrame([{"nombre": n, "SMILES": s, "fuente": src} for n, s, src, ik in S["reagents"]])
                    st.dataframe(_shade(dfa, "fuente"), width="stretch", height=240)
                    st.caption("Resaltado = reactivos que aportaste tú.")
            else:   # decoración: sin reactivo externo, usa sustituyentes internos
                st.markdown("#### Sustituyentes")
                st.caption("La decoración usa grupos pequeños internos (F, Cl, CN, OMe...); no subes reactivos.")
                c1, c2 = st.columns(2)
                S["n_analogs"] = c1.number_input("Cuántos análogos", 1, 200, S.get("n_analogs", 20))
                S["n_sub"] = c2.multiselect("Nº de sustituciones", [1, 2, 3], default=S.get("n_sub", [1]))
                S["use_ml"] = st.checkbox("Predecir ADMET con IA (más lento la 1ª vez)", value=S.get("use_ml", True))
                b = AdmelabBridge()
                if not b.available():
                    st.error("No encuentro el motor de diseño (admelab).")
                elif nuc and st.button("Generar análogos"):
                    with st.spinner("Generando y prediciendo propiedades..."):
                        d = b.design(nuc, use_ml=bool(S["use_ml"]),
                                     n_substitutions=S.get("n_sub", [1]) or [1], max_rows=int(S["n_analogs"]))
                    S["products"] = [dict(producto=(r.get("name") or f"analogo{i + 1:03d}"), smiles=r["SMILES"],
                                          fuente="interno", sintetizable=True, viabilidad="decoración")
                                     for i, r in enumerate(d.rows) if r.get("SMILES")]
                    if d.n_generated < int(S["n_analogs"]):
                        st.warning(f"Se generaron {d.n_generated} análogos: con {S.get('n_sub', [1])} "
                                   "sustitución(es) el espacio químico se agota ahí. Prueba con 2.")

        # ---- productos (coupling: esterifica; decoración: ya generados arriba) ----
        if reaction.kind == "coupling" and nuc and S.get("reagents"):
            st.markdown("---")
            b = AdmelabBridge()
            if not b.available():
                st.error("No encuentro el motor de reacción (admelab).")
            else:
                firma_p = (nuc, rxkey, tuple(sorted(ik for _n, _s, _src, ik in S["reagents"])))
                hecho_p = _ya_hecho("productos", firma_p) and S.get("products")
                if hecho_p:
                    st.caption(f"{len(S['products'])} productos ya construidos con este núcleo y "
                               "estos reactivos. Cambia alguno para volver a construir.")
                if st.button("Construir productos", type="primary", disabled=bool(hecho_p)):
                    alcs = [{"name": n, "smiles": s} for n, s, src, ik in S["reagents"]]
                    with st.spinner("Construyendo la serie..."):
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
            st.info(f"{n_ok} de {len(prods)} productos son sintetizables por esta reacción.")
            if not any(p.get("alcohol_smiles") for p in prods):
                # El nombrado compone el nombre a partir de las dos piezas de la reacción (radical y
                # acilo) y lo verifica con OPSIN. En decoración no hay dos piezas: habria que nombrar
                # una molécula arbitraria desde su estructura, que es un problema distinto y sin
                # solución fiable offline.
                st.caption("El nombre IUPAC solo se genera en reacciones de acoplamiento, donde se "
                           "compone a partir de los dos fragmentos y se verifica con OPSIN. En "
                           "decoración los productos se identifican por su SMILES.")
            if any(p.get("alcohol_smiles") for p in prods):
                if st.button("Nombrar (IUPAC, verificado con OPSIN)"):
                    with st.spinner("Nombrando y verificando por round-trip..."):
                        named = AdmelabBridge().name_esters(
                            [p["smiles"] for p in prods], [p.get("alcohol_smiles") or "" for p in prods],
                            acid_smiles=nuc, alcohol_names=[p.get("producto") for p in prods], use_web=True)
                    by = {n["smiles"]: n for n in named}
                    for p in prods:
                        inf = by.get(p["smiles"], {})
                        p["iupac_name"] = inf.get("iupac_name"); p["iupac_verif"] = inf.get("verified")
                    S["products"] = prods
                    nver = sum(1 for p in prods if p.get("iupac_verif"))
                    st.success(f"{nver} de {len(prods)} con nombre IUPAC verificado. El resto conserva su etiqueta "
                               "(nombre del alcohol); son de nicho y OPSIN no siempre los cubre offline.")
            dfp = pd.DataFrame(prods)
            cols = [c for c in ("producto", "iupac_name", "fuente", "oh_type", "viabilidad", "sintetizable", "smiles")
                    if c in dfp.columns]
            st.dataframe(_shade(dfp[cols], "fuente"), width="stretch", height=320)
            st.caption("Resaltado = productos con TUS reactivos. `sintetizable`=False son inviables por esta reacción.")

            st.caption("Las estructuras 2D de los productos se muestran en el panel derecho para "
                       "que verifiques el enlace y la estereoquímica.")

            with st.expander("Reporte ADMET (predice ~40 endpoints con IA para TODOS a la vez)"):
                if st.button("Predecir ADMET"):
                    with st.spinner("Prediciendo con ADMET-AI para todos (la 1ª vez descarga el modelo)..."):
                        pr = AdmelabBridge().predict([p["smiles"] for p in prods], use_ml=True)
                    S["admet"] = {rg.inchikey(r.get("SMILES")): r for r in pr.rows}
                if S.get("admet"):
                    _render_adme(S["admet"], [(p.get("producto") or f"prod{i}", p["smiles"])
                                              for i, p in enumerate(prods)], keyp="lig")

            c_sel1, c_sel2 = st.columns(2)
            solo_ok = c_sel1.checkbox("Acoplar solo los sintetizables", value=True)
            incluir_nuc = c_sel2.checkbox("Añadir el núcleo solo (referencia)", value=True,
                                          help="Acopla el núcleo sin esterificar como línea base: revela cuánta "
                                               "actividad aporta el scaffold por sí mismo, aparte de la cola.")
            firma_uso = (tuple(p.get("smiles") for p in prods), solo_ok, bool(incluir_nuc and nuc))
            usado = _ya_hecho("usar_productos", firma_uso)
            if usado:
                st.caption(f"Estos productos ya están cargados para el cribado ({len(S['ligands'])} "
                           "compuestos). Cambia la selección para volver a generarlos.")
            if st.button("Usar estos productos en el cribado", type="primary", disabled=usado):
                elegidos = [p for p in prods if (p.get("sintetizable") or not solo_ok)]
                # El nucleo desnudo entra como un candidato más (mismo pipeline), etiquetado aparte.
                if incluir_nuc and nuc:
                    elegidos = [dict(producto="nucleo_libre", smiles=nuc, fuente="núcleo",
                                     iupac_name=None, viabilidad="referencia (sin esterificar)",
                                     sintetizable=True)] + elegidos
                nombres = [lig.safe_name(p.get("producto") or f"prod{i}") for i, p in enumerate(elegidos)]
                # A la carpeta de ENTRADA, no a 'ligands': esa la borra la limpieza de cada corrida.
                with st.spinner(f"Generando 3D de {len(elegidos)} compuestos..."):
                    made = lig.materialize([p["smiles"] for p in elegidos], proj / "ligandos_entrada", names=nombres)
                hechos = {nm for nm, _, _ in made}
                S["ligands"] = [str(p) for _, p, _ in made]
                meta = pd.DataFrame([{"name": nm, "smiles": p.get("smiles"), "fuente": p.get("fuente", "?"),
                                      "producto": p.get("producto"), "iupac_name": p.get("iupac_name"),
                                      "viabilidad": p.get("viabilidad")}
                                     for (nm, p) in zip(nombres, elegidos) if nm in hechos])
                (proj / "ligands_meta.csv").write_text(meta.to_csv(index=False))
                _marcar_hecho("usar_productos", firma_uso)
                extra = " (incluye el núcleo solo como referencia)" if (incluir_nuc and nuc) else ""
                st.toast(f"{len(made)} compuestos construidos y listos para el paso 3.", icon="⚗️")
                st.success(f"{len(made)} compuestos listos para el paso 3{extra}.")

    else:  # Subir ligandos listos
        S["lead"] = None
        ups = st.file_uploader("Sube ligandos", type=["mol2", "sdf", "mol", "smi"], accept_multiple_files=True)
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
                st.warning(f"{sin_smiles} de {len(filas)} ligandos no dieron una estructura legible; "
                           "para ellos no habrá ADMET ni descriptores.")
        if S["ligands"]:
            st.write(f"**{len(S['ligands'])} ligandos:** " + ", ".join(Path(p).name for p in S["ligands"][:8]))
            ml = proj / "ligands_meta.csv"
            if ml.exists():
                mdf = pd.read_csv(ml)
                con = mdf["smiles"].notna().sum() if "smiles" in mdf.columns else 0
                st.caption(f"Estructura leída de {con} de {len(mdf)}: permite calcular ADMET, "
                           "eficiencia de ligando, SAscore y alertas PAINS.")
                items = [(r["name"], r["smiles"]) for _, r in mdf.iterrows() if pd.notna(r.get("smiles"))]
                if items:
                    with st.expander("Reporte ADMET de los ligandos subidos"):
                        if st.button("Predecir ADMET", key="adme_subidos"):
                            with st.spinner("Prediciendo con ADMET-AI..."):
                                pr = AdmelabBridge().predict([s for _, s in items], use_ml=True)
                            S["admet"] = {**(S.get("admet") or {}),
                                          **{rg.inchikey(r.get("SMILES")): r for r in pr.rows}}
                        if S.get("admet"):
                            _render_adme(S["admet"], items, keyp="sub")

    cup = st.file_uploader("Controles (ligando co-cristalizado)", type=["mol2", "sdf", "mol"],
                           accept_multiple_files=True,
                           help="Si ya lo extrajiste en el paso 1, no hace falta subir nada.")
    if cup:
        d = proj / "receptores"
        d.mkdir(parents=True, exist_ok=True)
        nuevos = []
        for u in cup:
            p = d / u.name
            try:
                p.write_bytes(u.getvalue())
            except Exception as e:
                st.error(f"No pude guardar {u.name}: {e}")
                continue
            if str(p) not in S["controls"]:
                S["controls"].append(str(p))
                nuevos.append(u.name)
        # Sin este aviso, subir un control parecía no hacer nada: el archivo se guardaba pero no
        # había ninguna senal de que se hubiera registrado.
        if nuevos:
            st.success(f"Control(es) cargado(s): {', '.join(nuevos)}. "
                       f"Total de controles: {len(S['controls'])}.")
        st.caption("Los controles cargados y los extraídos en el paso 1 se acoplan junto a los "
                   "ligandos y definen la huella de referencia. Con varios receptores, asigna cada "
                   "uno a su receptor en el paso 3.")

# ---------------------------------------------------------------- 3. ejecutar
def _etapa_ejecutar():
    st.subheader("Ejecutar el cribado")
    recs = [Path(p) for p in S["receptors"]]
    ctrls = [Path(p) for p in S["controls"]]
    ligs = [Path(p) for p in S["ligands"]]
    st.write(f"Receptores: **{len(recs)}** · Controles: **{len(ctrls)}** · "
             + (f"Líder: `{S.get('lead')}`" if S.get("lead") else f"Ligandos: **{len(ligs)}**"))
    if recs and not ctrls:
        st.warning("No hay ningún control cargado. El control se dockea junto a los ligandos y define la "
                   "referencia; sin él no hay línea base ni validación. Extrae el co-cristalizado en el paso 1 "
                   "(o súbelo abajo). Si ya lo extrajiste, revisa que esté en la carpeta `receptores/` del proyecto.")

    boxes = {}
    site_boxes = {}   # docking hibrido: receptor -> [(etiqueta, Box)] con bolsillos adicionales
    if recs:
        st.markdown("**Caja de búsqueda** — dónde buscar dentro de la proteína.")
        st.caption("Lo más fiable es centrarla en el ligando co-cristalizado: marca el sitio real. "
                   "El centro geométrico o un cofactor apuntan a otro lugar.")
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
                with st.expander("Asignar los controles que faltan a su receptor", expanded=True):
                    st.caption("Estos controles no se pudieron ubicar por geometría; indica a qué "
                               "receptor pertenece cada uno.")
                    _rec_labels = ["(ninguno)"] + [r.stem for r in recs]
                    for c in sin_ubicar:
                        sel = st.selectbox(f"Control «{c.stem}»", _rec_labels,
                                           key=f"ctrlrec_{sc.normalize_key(c.stem)}")
                        if sel != "(ninguno)":
                            manual[sc.normalize_key(c.stem)] = sel
            else:
                _pares = ", ".join(f"{c.stem} → {auto[sc.normalize_key(c.stem)]}" for c in ctrls)
                st.caption(f"Controles asignados automáticamente por geometría: {_pares}.")
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
            if _tv.button("Ver en el visor", key=f"verrec_{r.name}", use_container_width=True):
                S["vis_box_rec"] = str(r)
                st.rerun()
            grupos = dk.hetero_groups(r)
            ctrl = next((c for c in ctrls if asignacion.get(sc.normalize_key(c.stem)) == r.stem), None)
            # detección de pockets (fpocket), bajo demanda
            b1, b2 = st.columns([1, 3])
            ya_pk = bool(S["pockets"].get(str(r)))
            if b1.button("Detectar pockets", key=f"pk_{r.name}", type="primary",
                         disabled=not pk.fpocket_available() or ya_pk,
                         help="Cavidades ya detectadas para este receptor." if ya_pk else None):
                S["vis_box_rec"] = str(r)              # el visor sigue al receptor recien detectado
                with st.spinner("Buscando cavidades con fpocket..."):
                    S["pockets"][str(r)] = pk.detect(r)
                # Se muestran en el visor sin que haya que activarlo: es el motivo de detectarlas.
                S["vis_ver_cav"] = True
                st.rerun()
            pkts = S["pockets"].get(str(r), [])
            if not pkts and not pk.fpocket_available():
                b2.caption("fpocket no instalado: `conda install -n cribado -c conda-forge fpocket`.")
            # fuente de la caja. Al cambiarla, el visor pasa a este receptor: se esta editando, y
            # ver el cambio sobre otra estructura confundiría. Se hace por callback porque corre
            # antes de instanciar el selector del visor, que es cuando su valor aún se puede fijar.
            pk_opts = {p["label"]: p for p in pkts}
            opts = ([f"Centrar en el control ({ctrl.name})"] if ctrl else []) \
                + list(pk_opts.keys()) + ["Automática"] + list(grupos.keys())
            pick = st.selectbox("Fuente de la caja", opts, key=f"box_{r.name}",
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
                        "Usada": ("principal" if p["n"] == elegida else "híbrido") if en_uso else "",
                        "Drogabilidad": p["druggability"], "Score": p["score"],
                        "Volumen (Å³)": round(p.get("volume") or 0),
                        "Cavidad (Å)": "%.0f×%.0f×%.0f" % (p.get("ex", 0), p.get("ey", 0), p.get("ez", 0)),
                        "Caja (Å)": ("%.0f×%.0f×%.0f" % (p["sx"], p["sy"], p["sz"])
                                     + (" *" if p.get("minimo_aplicado") else "")),
                        "Esferas-α": p.get("spheres"),
                        "Hidrofobicidad": pr.get("Hydrophobicity score"),
                        "Polaridad": pr.get("Polarity score"),
                        "Carga": pr.get("Charge score"),
                        "SASA apolar": pr.get("Apolar SASA"),
                        # Flexibility se omite: fpocket la deriva de los B-factors, que la preparación
                        # con PDBFixer deja a 0; siempre saldría vacia y confundiría.
                        "Residuos": ", ".join(resid[:14]) + ("…" if len(resid) > 14 else ""),
                    }
                    # La columna catalitica solo aparece si el usuario ya designo residuos ancla;
                    # sin esa referencia no hay forma de decidirlo y un "?" no informa de nada.
                    if cats_r:
                        fila["Catalítica"] = "sí" if hay_cat else "no"
                    filas.append(fila)
                # Por receptor: con más de uno, un único dict global lo sobrescribia en cada
                # iteración y el visor acababa mostrando (o borrando) las cavidades del último.
                S.setdefault("_cavidades", {})[str(r)] = cav
                dfp = pd.DataFrame(filas)
                st.dataframe(
                    dfp, width="stretch", hide_index=True, height=230,
                    column_config={"Color": st.column_config.TextColumn(
                        "Color", width="small",
                        help="Color con el que la cavidad se dibuja en el visor")})
                _descargar_tabla(pd.DataFrame([{"Pocket": p["n"], **p.get("props", {}),
                                                "Residuos": ", ".join(p.get("residues") or [])}
                                               for p in pkts]),
                                 f"cavidades_{r.stem}", key=f"cav_{r.name}")
                st.caption("Todas las cavidades se dibujan a la vez en el panel derecho. La **usada para "
                           "el docking** va resaltada y más opaca. `Cavidad` es su extensión real; `Caja` es "
                           "la de búsqueda, con un mínimo de 14 Å porque por debajo no cabría un ligando "
                           "(marcado con `*` cuando se ha aplicado ese mínimo).")
                with st.expander("Todas las propiedades que calcula fpocket"):
                    st.dataframe(pd.DataFrame([{"pocket": p["n"], **p.get("props", {})} for p in pkts]),
                                 width="stretch", hide_index=True)
            else:
                S.get("_cavidades", {}).pop(str(r), None)

            if ctrl and pick.startswith("Centrar en el control"):
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
                st.markdown("**Centro** — dónde se sitúa la caja (Å)")
                st.caption("Mueve la caja por el espacio. Los ejes se ven en el visor de la derecha.")
                cc = st.columns(3)
                cx = cc[0].number_input("← X →", step=1.0, key=f"cx_{r.name}", format="%.1f",
                                        help="Izquierda / derecha (eje rojo).")
                cy = cc[1].number_input("↓ Y ↑", step=1.0, key=f"cy_{r.name}", format="%.1f",
                                        help="Abajo / arriba (eje verde).")
                cz = cc[2].number_input("⊙ Z ⊗", step=1.0, key=f"cz_{r.name}", format="%.1f",
                                        help="Hacia dentro / hacia fuera de la pantalla (eje azul).")
            with gs:
                st.markdown("**Tamaño** — cuánto abarca la caja (Å)")
                st.caption("Agranda o encoge cada lado. Si el ligando no cabe, Vina falla.")
                cs = st.columns(3)
                sx = cs[0].number_input("ancho X", min_value=6.0, step=1.0, key=f"sx_{r.name}", format="%.1f")
                sy = cs[1].number_input("alto Y", min_value=6.0, step=1.0, key=f"sy_{r.name}", format="%.1f")
                sz = cs[2].number_input("fondo Z", min_value=6.0, step=1.0, key=f"sz_{r.name}", format="%.1f")
            boxes[str(r)] = dk.Box(cx, cy, cz, sx, sy, sz)
            # El panel derecho se dibuja después que el izquierdo, así que puede leer esto.
            S.setdefault("_boxes", {})[str(r)] = boxes[str(r)].as_dict()
            st.caption("La caja se dibuja sobre el receptor en el panel derecho.")
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
                    f"El ligando más grande mide **{_minimo - 4:.0f} Å** en su eje mayor y la caja "
                    f"tiene {min(sx, sy, sz):.0f} Å en su lado menor. Cabe, pero no puede "
                    f"reorientarse: la búsqueda queda restringida a las orientaciones que entran. "
                    f"Sube los tres lados a **{_minimo:.0f} Å** como mínimo.")

            # Docking hibrido: acoplar el mismo ligando en varios bolsillos del mismo receptor.
            if pkts:
                # El bolsillo que ya es la caja principal no se ofrece: elegirlo creaba un sitio
                # duplicado que luego había que colapsar.
                disponibles_hib = [p["label"] for p in pkts if p["label"] != pick]
                extra = st.multiselect(
                    "Acoplar también en otros bolsillos (docking híbrido)",
                    disponibles_hib, key=f"sitios_{r.name}",
                    help="Cada bolsillo elegido se acopla por separado y su ranking sale aparte. Revela si "
                         "un compuesto prefiere el sitio catalítico o se cuela en uno alostérico.")
                if extra:
                    lst = [("principal", boxes[str(r)])]
                    for lab in extra:
                        pdd = next((p for p in pkts if p["label"] == lab), None)
                        if pdd:
                            lst.append((f"Pk{pdd['n']}", pk.pocket_box(pdd)))
                    site_boxes[str(r)] = lst
                    st.caption(f"Docking híbrido: {len(lst)} sitios (principal + {len(extra)} bolsillo(s)).")
                S[f"_hib_sel_{r.name}"] = set(extra)

    params = _params_docking()
    c1, c2 = st.columns([2, 1])
    reuse = c1.checkbox("Reutilizar cálculos previos de esta carpeta", value=False,
                        help="Desactivado, cada ejecución recalcula todo. Actívalo solo si nada ha cambiado: "
                             "reutilizar poses hechas con otra caja da resultados falsos.")
    if c2.button("Borrar resultados de esta carpeta"):
        _confirmar_borrado(proj)
    st.caption(f"Todo se guarda en `{proj}` — poses, complejos, XML de PLIP y las tablas CSV.")

    if not recs:
        st.info("Prepara al menos un receptor en el paso 1.")
    elif not (S.get("lead") or ligs):
        st.info("Elige compuestos en el paso 2.")
    else:
        # El boton se bloquea mientras la configuración no cambie: un cribado dura minutos y
        # volver a lanzarlo por olvido o por una doble pulsacion desperdicia el trabajo.
        firma = (tuple(sorted(str(x) for x in recs)), tuple(sorted(str(x) for x in ligs)),
                 tuple(sorted(str(x) for x in ctrls)),
                 tuple(sorted((k, tuple(v.as_dict().values())) for k, v in boxes.items())),
                 tuple(sorted((k, len(v)) for k, v in site_boxes.items())),
                 tuple(sorted(params.items())), reuse, str(proj))
        hecho = _ya_hecho("run", firma)
        if st.button("Ejecutar", type="primary", disabled=hecho,
                     help="Ya ejecutado con esta configuración. Cambia algo para volver a lanzarlo."
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

            with st.status("Ejecutando...", expanded=True) as status:
                try:
                    pl.run(cfg, on_step=_paso)
                    _marcar_hecho("run", firma)
                    status.update(label="Cribado completado", state="complete")
                    S["_log_estado"] = "completo"
                    st.toast("Cribado completado. Ve a la pestaña de resultados.", icon="✅")
                    st.success("Listo. Ve a la pestaña de resultados.")
                except Exception as e:
                    status.update(label="Falló", state="error")
                    S["_log_estado"] = "error"
                    S["_log_run"].append(("error", str(e)))
                    st.toast("El cribado falló. Revisa el mensaje.", icon="⚠️")
                    st.error(str(e))
        elif S.get("_log_run"):
            _registro_corrida()
        if hecho:
            st.caption("Cribado completado con esta configuración. Modifica algún parámetro "
                       "para habilitar de nuevo el botón.")


def _registro_corrida():
    """Registro de la última corrida, conservado entre pestanas."""
    estado = S.get("_log_estado", "completo")
    with st.status("Cribado completado" if estado == "completo" else "La corrida falló",
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
    st.markdown("### Resumen pleiotrópico — actividad en varias dianas")
    if presente_en_todas.empty:
        st.caption("Ningún compuesto se acopló en todas las dianas; no hay comparación de amplio espectro.")
        return
    presente_en_todas = presente_en_todas.assign(
        **{"mínima": presente_en_todas.min(axis=1).round(1),
           "media": presente_en_todas.mean(axis=1).round(1)})
    tabla = presente_en_todas.sort_values("mínima", ascending=False).reset_index()
    tabla.columns = ["compuesto"] + [f"{c} (%)" if c in dianas else c for c in tabla.columns[1:]]
    st.caption("Efectividad (%) de cada compuesto en cada diana, ordenado por la **mínima** entre "
               "dianas: arriba, los de amplio espectro. Solo los que se acoplaron en todas.")
    st.dataframe(tabla.round(1), width="stretch", hide_index=True,
                 height=min(340, 60 + 34 * len(tabla)))
    _descargar_tabla(tabla, "pleiotropico", key="pleio")
    mejor_amplio = tabla.iloc[0]
    st.success(f"Mejor amplio espectro: **{mejor_amplio['compuesto']}** "
               f"(mínima {mejor_amplio['mínima']:.0f} % entre {len(dianas)} dianas).")
    st.divider()


# ---------------------------------------------------------------- 4. resultados
def _etapa_resultados():
    st.subheader("Resultados")
    meta_p, inter_p, dock_p = proj / "run.json", proj / "interacciones.csv", proj / "resultados_docking.csv"
    if not (meta_p.exists() and inter_p.exists()):
        st.info("Todavía no hay resultados en esta carpeta. Ejecuta el paso 3.")
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

        st.markdown("**Residuos catalíticos / ancla**. El puntaje "
                    "premia además la calidad de las demás interacciones del pocket.")
        st.caption("Sugeridos automáticamente desde las interacciones direccionales del ligando "
                   "cristalográfico. Edítalos si conoces el sitio catalítico real de tu diana.")
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
            sec[R] = cols[i].multiselect(f"{R} · secundarios (bonificación, no obligatorios)", libres,
                                         default=prev_s, key=f"sec_{R}")
        st.caption("**Gate** (catalíticos): obligatorios, se penaliza faltarlos. **Secundarios**: anclas "
                   "conocidas del bolsillo que suman más que un contacto cualquiera (×w_sec) pero no se exigen.")

        # Validación: si el control no recupera su postura, el resto no es fiable
        val_p = proj / "validacion_redocking.csv"
        if val_p.exists():
            val = pd.read_csv(val_p)
            msg = vl.resumen(val)
            (st.success if not msg.startswith("ATENCION") else st.error)(msg)
            with st.expander("Detalle de la validación por redocking"):
                st.dataframe(val, width="stretch", hide_index=True)
                st.caption("RMSD frente al ligando co-cristalizado. Válido por debajo de 2 Å.")

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
                filas.append({"receptor": R, "recuperación": rec["recovery"], "Tanimoto": rec["tanimoto"],
                              "reproducidas": f"{rec['shared']}/{rec['ref_n']}", "extra (no cristal)": rec["extra"]})
            if filas:
                st.markdown("**Validación de interacciones** — control dockeado vs. ligando cristalográfico.")
                st.dataframe(pd.DataFrame(filas), width="stretch", hide_index=True)
                st.caption("`recuperación` = fracción de las interacciones del cristalográfico que la pose dockeada "
                           "del control reproduce; `Tanimoto` incluye además los contactos extra que añade el docking. ")

        st.markdown("**Ponderación**")
        mw = meta.get("weights", {})
        metric_afin = st.radio(
            "Métrica del eje de afinidad", ["dock", "le"], horizontal=True,
            index=1 if str(mw.get("dock_metric", "dock")).lower() == "le" else 0,
            format_func=lambda m: "Score crudo (kcal/mol)" if m == "dock" else "Eficiencia de ligando (LE)",
            help="El score crudo de Vina premia moléculas grandes (sesgo de tamaño). LE = -ΔG/átomos pesados "
                 "corrige ese sesgo. Recomendado si tu biblioteca varía mucho de tamaño; reporta ambas columnas.")
        c1, c2, c3, c4 = st.columns(4)
        w_dock = c1.slider("Peso docking", 0.0, 1.0, float(mw.get("dock", 0.5)), 0.05)
        w_inter = c2.slider("Peso interacciones", 0.0, 1.0, float(mw.get("inter", 0.5)), 0.05)
        w_adme = c3.slider("Peso ADME", 0.0, 1.0, float(mw.get("adme", 0.0)), 0.05,
                           help="Calidad fisicoquímica (drug-likeness) del compuesto. Guarda contra premiar "
                                "solo moléculas grandes/grasas.")
        w_tox = c4.slider("Peso toxicidad", 0.0, 1.0, float(mw.get("tox", 0.0)), 0.05,
                          help="Requiere haber predicho ADMET (pestaña Ligandos); si no, este eje se ignora.")
        c5, c6, c7 = st.columns(3)
        w_cat = c5.slider("Peso de residuo catalítico", 1.0, 6.0, float(mw.get("w_cat", 3.0)), 0.5,
                          help="Cuánto vale una interacción con un residuo catalítico (gate) frente a uno normal del pocket.")
        w_sec = c6.slider("Peso de residuo secundario", 1.0, 3.0, float(mw.get("w_sec", 1.5)), 0.25,
                          help="Cuánto vale una interacción con un ancla SECUNDARIA frente a un contacto de pocket normal (×1).")
        cat_gate = c7.slider("Exigencia catalítica", 0.0, 1.0, float(mw.get("cat_gate", 0.5)), 0.05,
                             help="0 = no penaliza faltar a un catalítico; 1 = faltar a todos anula el puntaje.")
        _axw = {"docking": w_dock, "interacción": w_inter, "ADME": w_adme, "tox": w_tox}
        _tot = sum(_axw.values())
        if _tot > 0:
            st.caption("Contribución real de cada eje: "
                       + " · ".join(f"{k} {v / _tot * 100:.0f}%" for k, v in _axw.items() if v > 0))
        else:
            st.warning("Todos los pesos de eje están a 0: no habrá puntaje. Sube al menos uno.")
        with st.expander("Pesos por tipo de interacción (avanzado)"):
            st.caption("Valor de mérito por tipo (0-1). Por defecto: salino > H-bond > π > halógeno > hidrofóbica. "
                       "Guiado por literatura; ajústalo a tu criterio.")
            st.caption("`water` (puentes mediados por agua) solo interviene si conservas moléculas de "
                       "agua al preparar el receptor. En el flujo habitual se eliminan, así que este "
                       "peso no tiene efecto.")
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
            st.warning(f"Pesas **{' y '.join(faltan)}** pero no hay datos de ese eje en esta corrida: se ignora "
                       "en el puntaje. Predice ADMET antes, o baja su peso a 0 para que el Methods no lo declare.")

        with st.expander("Exportar Métodos (para el paper)"):
            st.caption("Parámetros, caja, pesos, referencia y versiones exactas del software. "
                       "Reproducibilidad lista para pegar en la sección de Métodos.")
            metodos = rp.methods_text(meta, weights=w, catalytic=cat, secondary=sec)
            st.download_button("Descargar Methods.md", metodos, file_name="PoliScreen_Methods.md",
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
                    f"Esta tabla mezcla **{' y '.join(_distintos)}**. Sus energías salen de "
                    "funciones de puntuación distintas y no son comparables entre sí: `best_dock`, "
                    "`pKi` y `LE` solo tienen sentido dentro de cada motor. Para comparar entre "
                    "ellos usa `efectividad_pct`, que se calcula sobre los contactos y no depende "
                    "del motor.")

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
            st.caption("Filas resaltadas = compuestos hechos con reactivos que aportaste tú.")

        if any("~" in str(x) for x in rk["receptor"].unique()):
            st.info("**Docking híbrido**: cada bloque es un bolsillo distinto del mismo receptor. "
                    "Compara la efectividad de un compuesto entre sitios para ver dónde prefiere unirse.")
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
            st.markdown(f"{_et} · referencia de interacciones: `{_refsrc}`")
            noc = sub[sub["is_control"] != 1]
            if not noc.empty:
                m1, m2, m3, m4 = st.columns(4)
                try:
                    bd = noc.loc[pd.to_numeric(noc["best_dock"], errors="coerce").idxmin()]
                    m1.metric("Mejor docking", str(bd["compound"])[:18], f"{bd['best_dock']:.2f} kcal/mol",
                              delta_color="inverse")
                except Exception:
                    pass
                try:
                    bi = noc.loc[pd.to_numeric(noc["best_inter"], errors="coerce").idxmax()]
                    m2.metric("Mejor calidad de interacción", str(bi["compound"])[:18], f"{bi['best_inter']:.2f}")
                except Exception:
                    pass
                try:
                    be = noc.loc[pd.to_numeric(noc["efectividad_pct"], errors="coerce").idxmax()]
                    m3.metric("Mejor efectividad", str(be["compound"])[:18], f"{be['efectividad_pct']:.0f} %")
                except Exception:
                    pass
                try:
                    bc = noc.loc[pd.to_numeric(noc["confidence"], errors="coerce").idxmax()]
                    m4.metric("Mayor confianza", str(bc["compound"])[:18], f"{bc['confidence']:.2f}")
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
        st.download_button("Descargar ranking (CSV)", rk.to_csv(index=False).encode(), "ranking.csv")

        # reporte ADMET sobre TODOS los ligandos del ranking (compuestos + nucleo + control), a elección
        # del usuario. smap ya incluye compuestos, nucleo y controles.
        items_all = [(c, smap[sc.normalize_key(c)]) for c in rk["compound"].unique()
                     if sc.normalize_key(c) in smap and pd.notna(smap.get(sc.normalize_key(c)))]
        if items_all:
            with st.expander("Reporte ADMET (compuestos + núcleo + control, los que elijas)"):
                nombres_all = [c for c, _ in items_all]
                elegidos_adme = st.multiselect("¿De qué ligandos predecir ADMET?", nombres_all,
                                               default=nombres_all, key="adme_sel_res")
                items = [(c, s) for c, s in items_all if c in elegidos_adme]
                if st.button("Predecir ADMET", key="pred_res") and items:
                    with st.spinner(f"Prediciendo con ADMET-AI para {len(items)} ligando(s)..."):
                        pr = AdmelabBridge().predict([s for _, s in items], use_ml=True)
                    S["admet"] = {**(S.get("admet") or {}), **{rg.inchikey(r.get("SMILES")): r for r in pr.rows}}
                if S.get("admet") and items:
                    _render_adme(S["admet"], items, keyp="res")

        st.markdown("---")
        st.markdown("**Diagrama de interacciones** de una pose concreta.")
        d1, d2, d3 = st.columns(3)
        R = d1.selectbox("Receptor", sorted(inter["receptor"].unique()))
        sr = inter[inter["receptor"] == R]
        cmp_ = d2.selectbox("Compuesto", sorted(sr["compound"].unique()))
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
            st.caption("Verde = reproduce una interacción del control (mismo residuo y mismo enlace). "
                       "Gris = contacto de más o el mismo residuo con otro tipo de enlace.")

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
    if etapa == "Receptores":
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
            opciones = ["Preparado"] + (["Original"] if _orig else []) \
                + (["Con su control"] if _ctrl else [])
            cual = c1.radio("Ver", opciones, horizontal=True, key="vis_ver_rec",
                            label_visibility="collapsed")
            ejes = c2.checkbox("Ejes XYZ", value=True, key="vis_ejes_rec")
            try:
                receptor = _orig if (cual == "Original" and _orig) else rsel
                ligando = _ctrl if cual == "Con su control" else None
                _h = _alto_visor(250)      # selector + radio + posible pie del control
                components.html(vw.view_html(receptor=receptor, ligando=ligando,
                                             mostrar_aguas=False, ejes=ejes, alto=_h), height=_h + 12)
                if cual == "Con su control" and _ctrl:
                    st.caption(f"Control de este receptor: `{Path(_ctrl).stem}`.")
            except Exception as e:
                st.error(f"No pude dibujar la estructura: {e}")
        else:
            _vacio("Prepara un receptor y aparecerá aquí en 3D.")

    elif etapa == "Ligandos":
        prods = S.get("products")
        nuc_png = S.get("_nucleo_png")
        peps = S.get("_pep_preview")
        if peps and S.get("modo_ligandos") == "Generar péptidos":
            st.markdown("**Secuencias generadas**")
            # La secuencia se colorea por clase de residuo: de un vistazo se ve si el péptido es
            # anfipatico (bloques hidrofóbicos y cationicos alternados) o uniforme.
            leyenda = [("#3d7ea6", "hidrofóbico"), ("#b5453c", "carga +"),
                       ("#3f7d4e", "carga −"), ("#7a6ba8", "polar"), ("#8a8a8a", "G/P")]
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
                quim = ("ciclado cabeza-cola" if _cic else
                        ", ".join(filter(None, ["N-acetilado" if _nac else "",
                                                "C-amidado" if _cam else ""])) or "extremos libres")
                st.image(S["_pep_grid"],
                         caption=f"Estructura de {S.get('_pep_grid_n', 0)} péptidos · {quim}.")
        elif prods:
            png = vw.grid_png([p.get("smiles") for p in prods],
                              legends=[str(p.get("producto") or "") for p in prods])
            if png:
                st.image(png, caption=f"{len(prods)} productos construidos. "
                                      "Verifica el enlace éster y la estereoquímica.")
        elif nuc_png:
            st.image(nuc_png, caption="Núcleo con índices de átomo; en color, el sitio reactivo.")
        elif S["ligands"]:
            st.success(f"{len(S['ligands'])} ligandos listos para acoplar.")
            st.caption(", ".join(Path(p).stem for p in S["ligands"][:20]))
        else:
            _vacio("Construye o sube ligandos y verás aquí sus estructuras.")

    elif etapa == "Ejecutar":
        cajas = S.get("_boxes") or {}
        cav_map = S.get("_cavidades") or {}
        if cajas:
            c1, c2, c3 = st.columns([2, 1, 1])
            rsel = c1.selectbox("Receptor", list(cajas), format_func=lambda p: Path(p).name,
                                key="vis_box_rec", label_visibility="collapsed")
            # Cada receptor tiene sus propias cavidades; se dibujan las del que esta seleccionado.
            grupos_r = cav_map.get(rsel)
            ver_cav = c2.checkbox("Cavidades", value=bool(grupos_r), key="vis_ver_cav")
            ejes = c3.checkbox("Ejes XYZ", value=True, key="vis_ejes_box")
            grupos = grupos_r if (ver_cav and grupos_r) else None
            try:
                _h = _alto_visor(210)      # fila de selector/casillas + pie de la caja
                components.html(vw.view_html(receptor=rsel, caja=cajas[rsel], cavidades=grupos,
                                             mostrar_aguas=False, ejes=ejes, alto=_h), height=_h + 12)
                b = cajas[rsel]
                st.caption(f"Caja (malva): centro ({b['cx']}, {b['cy']}, {b['cz']}) · "
                           f"{b['sx']} × {b['sy']} × {b['sz']} Å"
                           + (f" · {len(grupos)} cavidades; la usada va resaltada." if grupos else ""))
            except Exception as e:
                st.error(f"No pude dibujar: {e}")
        else:
            _vacio("Define la caja de búsqueda y se dibujará aquí sobre el receptor.")

    else:  # Resultados
        vista = st.radio("Vista", ["Resumen", "Complejo 3D"], horizontal=True,
                         key="vis_res_vista", label_visibility="collapsed")
        if vista == "Complejo 3D":
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
        _vacio("Ejecuta un cribado y aquí podrás recorrer los complejos en 3D.")
        return
    t = pd.read_csv(inter_p)
    if "name" not in t.columns or t.empty:
        _vacio("La tabla de interacciones está vacía.")
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
    C = c2.selectbox("Compuesto", _cmps or sorted(sr["_cmp"].unique()), key="vis_cx_cmp")
    scmp = sr[sr["_cmp"] == C]
    M = c3.selectbox("Pose", sorted(scmp["_mod"].unique()), key="vis_cx_pose")
    o1, o2 = st.columns(2)
    sup = o1.checkbox("Mostrar la superficie", value=False, key="vis_cx_sup",
                      help="Superficie molecular translúcida del receptor. Con la cinta sola no se "
                           "distingue si el ligando está dentro de la cavidad o apoyado por fuera.")
    het = o2.checkbox("Cofactores y hetero", value=True, key="vis_cx_het")
    fila = scmp[scmp["_mod"] == M]
    if fila.empty:
        _vacio("No hay ninguna pose con esa combinación.")
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
                _vacio("No encuentro los archivos de esta pose en la carpeta del proyecto.")
                return
            html = vw.view_html(receptor=fus, mostrar_aguas=False, mostrar_hetero=het,
                                superficie=sup, alto=_h)
        components.html(html, height=_h + 12)
    except Exception as e:
        st.error(f"No pude dibujar el complejo: {e}")
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
        _vacio("Ejecuta un cribado y aquí verás el resumen de resultados.")
        return
    rk = pd.read_csv(rk_p)
    # Con docking hibrido hay un bloque por sitio: sin selector solo se veria el primero.
    sitios = sorted(rk["receptor"].unique()) if "receptor" in rk.columns else []
    if len(sitios) > 1:
        sel = st.selectbox("Sitio", sitios, key="vis_res_sitio",
                           format_func=lambda s: s.split("~", 1)[1] if "~" in str(s) else s)
        rk = rk[rk["receptor"] == sel]
        st.caption(f"Resumen del sitio **{sel}**. Cambia de sitio para comparar dónde une mejor cada compuesto.")
    # Sin la comprobación explícita, cuando falta la columna rk.get() devuelve un escalar y la
    # comparación produce un booleano suelto: pandas lo interpreta como nombre de columna y falla.
    noc = (rk[rk["is_control"] != 1] if "is_control" in rk.columns else rk).copy()
    if noc.empty or "efectividad_pct" not in noc.columns:
        st.info("Sin compuestos que resumir todavía.")
        return

    ef = pd.to_numeric(noc["efectividad_pct"], errors="coerce")
    conf = pd.to_numeric(noc.get("confidence"), errors="coerce") if "confidence" in noc else None
    superan = int((ef >= 105).sum())
    fiables = int((conf >= 0.5).sum()) if conf is not None else 0

    # Titular: el mejor compuesto por efectividad
    mejor = noc.loc[ef.idxmax()]
    enc, img = st.columns([2, 1])
    enc.markdown(f"### {str(mejor['compound'])[:38]}")
    enc.caption("Compuesto con mayor efectividad")
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
    k1.metric("Efectividad", f"{ef.max():.0f} %")
    if "best_dock" in noc.columns:
        k2.metric("Afinidad", f"{float(mejor['best_dock']):.1f}", "kcal/mol", delta_color="off")
    if conf is not None and pd.notna(mejor.get("confidence")):
        k3.metric("Confianza", f"{float(mejor['confidence']):.2f}")

    st.divider()
    c1, c2, c3 = st.columns(3)
    c1.metric("Compuestos", len(noc))
    c2.metric("Superan el control", superan, help="Efectividad ≥ 105 % respecto al cristalográfico.")
    if conf is not None:
        c3.metric("Confianza ≥ 0.5", fiables, help="Evidencias concordantes: resultado fiable.")

    # Podio: los cinco mejores, con barra proporcional para leerlos de un vistazo
    st.markdown("**Los cinco mejores**")
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
        st.warning("Ningún compuesto alcanza confianza 0.5. Con pocas poses la métrica pierde "
                   "resolución: sube «Poses por ligando» en el paso 3 y vuelve a ejecutar.")

    # Descarga del resumen que se esta viendo, no de todo el ranking.
    cols_res = [c for c in ("compound", "best_dock", "pKi", "LE", "best_inter", "efectividad_pct",
                            "percentil", "confidence", "cnn_score", "consenso") if c in noc.columns]
    _descargar_tabla(noc[cols_res].sort_values("efectividad_pct", ascending=False),
                     "resumen_" + str(noc["receptor"].iloc[0]).replace("~", "_"), key="resumen_vis")


# ---------------------------------------------------------------- composicion de la pantalla
_ETAPA_FN = {"Receptores": _etapa_receptores, "Ligandos": _etapa_ligandos,
             "Ejecutar": _etapa_ejecutar, "Resultados": _etapa_resultados}

# Altura fija en ambos paneles: cada uno hace su propio scroll y la pagina en conjunto no se
# desplaza, de modo que la cabecera y la barra de etapas quedan siempre a la vista.
_ALTO = int(S.get("cfg_alto", 520))
if S.get("_aviso"):
    st.success(S.pop("_aviso"))
_rep = float(S.get("cfg_reparto", 0.46))
_izq, _der = st.columns([_rep, 1.0 - _rep], gap="medium")
with _izq:
    st.markdown(f"**Herramientas · {S['etapa']}**")
    with st.container(height=_ALTO, border=True, key="panel_izq"):
        _ETAPA_FN[S["etapa"]]()
with _der:
    st.markdown("**Visualizador**")
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
