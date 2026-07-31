"""Sesiones de trabajo y exportación de resultados.

Una sesión `.poliscreen` es un ZIP con manifiesto. Permite cerrar la aplicación y volver después
para reabrir el análisis, cambiar la ponderacion y reexaminar los resultados **sin repetir el
docking** y sin tener que recordar la ruta de la carpeta de trabajo.

Dos tamaños:
  - ligera (por defecto): configuración, tablas y receptores. Basta para reabrir, repuntuar y
    ver diagramas. Del orden de megabytes.
  - completa: anade poses y complejos, de modo que también se pueden reexaminar las estructuras
    3D y reejecutar PLIP. Puede ocupar cientos de megabytes.
"""
from __future__ import annotations

import io
import json
import re
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Sequence

FORMATO = 1
EXT = ".poliscreen"

# --- Rutas escritas a mano -------------------------------------------------------------------
# En Linux la barra invertida es un carácter válido de nombre, no un separador. Pegar una ruta del
# Explorador de Windows crea así una única carpeta "home\diego\proyecto\..." en vez del árbol
# esperado. Se traducen las tres formas de ruta de Windows: UNC de WSL, unidad con letra y sin barra.
_UNC_WSL = re.compile(r"^/{2,}(?:wsl\.localhost|wsl\$)/[^/]+", re.I)
_UNIDAD = re.compile(r"^([A-Za-z]):(?=/|$)")


def normalizar_ruta(texto: str, base: Optional[Path] = None) -> tuple:
    """Convierte lo que el usuario escribe en una ruta POSIX utilizable.

    Devuelve (ruta, aviso). El aviso no va vacío cuando el texto se reinterpreta: debe verse, porque
    escribir en otra carpeta arruina el análisis en silencio y deja los resultados donde nadie mira.
    """
    crudo = (texto or "").strip().strip('"').strip("'")
    if not crudo:
        return Path.home() / "poliscreen_proyectos" / "demo", ""

    s = crudo.replace("\\", "/")
    s = _UNC_WSL.sub("", s) or "/"                  # \\wsl.localhost\Distro\home\x  ->  /home/x
    m = _UNIDAD.match(s)
    if m:                                            # C:\Users\...  ->  /mnt/c/Users/...
        s = f"/mnt/{m.group(1).lower()}/{s[m.end():].lstrip('/')}".rstrip("/")
    s = re.sub(r"/{2,}", "/", s) or "/"

    p = Path(s).expanduser()
    if not p.is_absolute():
        # Sin barra inicial suele ser una ruta absoluta a la que se le quitó el prefijo UNC a mano
        # ("home/diego/..."): si su primer tramo es un directorio de la raíz, era absoluta.
        primero = p.parts[0] if p.parts else ""
        p = Path("/") / p if (Path("/") / primero).is_dir() else (base or Path.home()) / p

    if str(p) == crudo:
        return p, ""
    if "\\" in crudo:
        aviso = (f"Se detectó una ruta de Windows. PoliScreen trabaja dentro de Linux (WSL), "
                 f"así que se usará `{p}`.")
    else:
        aviso = f"Ruta ajustada a `{p}`."
    return p, aviso

# Tablas y configuración: el nucleo reproducible de una corrida.
ARCHIVOS_BASE = ("run.json", "ranking.csv", "resumen.csv", "interacciones.csv",
                 "resultados_docking.csv", "validacion_redocking.csv", "analogos.csv",
                 "ligands_meta.csv")
# Entradas del usuario: receptores preparados, controles y ligandos.
CARPETAS_BASE = ("receptores", "ligandos_entrada")
# Derivados pesados: solo en la sesión completa.
CARPETAS_PESADAS = ("poses", "Complejos_Fusionados", "xml_plip", "xtal")

# Catalogo de exportaciones sueltas. Cada entrada: (descripción, tipo, origen).
EXPORTS = {
    "resultados_csv":    ("Ranking completo con todas las métricas", "archivo", "ranking.csv"),
    "resumen_csv":       ("Resumen compacto por compuesto", "archivo", "resumen.csv"),
    "interacciones_csv": ("Matriz de interacciones por pose (PLIP)", "archivo", "interacciones.csv"),
    "docking_csv":       ("Energías de todas las poses", "archivo", "resultados_docking.csv"),
    "ligandos_csv":      ("Tabla de ligandos: nombre, SMILES, IUPAC y procedencia", "archivo", "ligands_meta.csv"),
    "validacion_csv":    ("Validación por redocking del control", "archivo", "validacion_redocking.csv"),
    "receptores":        ("Receptores preparados y controles co-cristalizados", "carpeta", "receptores"),
    "ligandos_zip":      ("Estructuras 3D de los ligandos (SDF)", "carpeta", "ligandos_entrada"),
    "complejos_zip":     ("Complejos receptor-ligando (PDB)", "carpeta", "Complejos_Fusionados"),
    "poses_zip":         ("Poses de docking por modelo", "carpeta", "poses"),
    "methods":           ("Sección de Métodos: parámetros y versiones", "generado", "PoliScreen_Methods.md"),
}

# Qué conviene descargar. La carpeta ya contiene casi todo, así que exportar solo importa en dos
# casos, que son los marcados: lo que no existe como archivo hasta exportarlo, y el subconjunto
# mínimo para reproducir y publicar la corrida sin arrastrar los intermedios pesados (regenerables).
RECOMENDADO = {
    "methods":           "No está en la carpeta: se redacta al exportar",
    "resultados_csv":    "Tabla principal: puntuación, Ki, eficiencia y confianza",
    "interacciones_csv": "Huella de contactos que sustenta la puntuación",
    "ligandos_csv":      "Procedencia de cada ligando: SMILES e IUPAC",
    "validacion_csv":    "RMSD del redocking: lo primero que revisa un evaluador",
    "receptores":        "Entradas exactas para repetir la corrida",
    "ligandos_zip":      "Entradas exactas para repetir la corrida",
}
# Lo que se puede rehacer a partir de lo anterior; se ofrece, pero pesa y no hace falta guardarlo.
REGENERABLE = ("resumen_csv", "docking_csv", "complejos_zip", "poses_zip")


def _peso(proj: Path, tipo: str, origen: str) -> tuple:
    """(bytes, n_archivos) de un elemento del catalogo."""
    if tipo == "generado":
        return 0, 1
    p = proj / origen
    if tipo == "archivo":
        return (p.stat().st_size, 1) if p.is_file() else (0, 0)
    if not p.is_dir():
        return 0, 0
    fs = [f for f in p.rglob("*") if f.is_file()]
    return sum(f.stat().st_size for f in fs), len(fs)


def catalogo(proj) -> dict:
    """Todo lo exportable con lo que la interfaz necesita para decidir: existe, cuanto pesa y si
    merece la pena bajarlo. {clave: {desc, hay, bytes, n, motivo, regenerable}}."""
    proj = Path(proj)
    out = {}
    for clave, (desc, tipo, origen) in EXPORTS.items():
        b, n = _peso(proj, tipo, origen)
        hay = (proj / "run.json").exists() if tipo == "generado" else n > 0
        out[clave] = {"desc": desc, "hay": hay, "bytes": b, "n": n,
                      "motivo": RECOMENDADO.get(clave, ""), "regenerable": clave in REGENERABLE}
    return out


def paquete(proj, claves: Sequence[str], methods_text: Optional[str] = None) -> tuple:
    """ZIP en memoria con los elementos elegidos. Devuelve (bytes, lista de lo incluido).

    No escribe nada en la carpeta del proyecto: la descarga la decide el usuario, no la aplicación.
    """
    proj = Path(proj)
    buf, incluidos = io.BytesIO(), []
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for clave in claves:
            if clave not in EXPORTS:
                continue
            _desc, tipo, origen = EXPORTS[clave]
            if tipo == "generado":
                if methods_text:
                    zf.writestr(origen, methods_text)
                    incluidos.append(origen)
            elif tipo == "archivo":
                if _add_file(zf, proj / origen, f"tablas/{origen}"):
                    incluidos.append(f"tablas/{origen}")
            else:
                n = _add_dir(zf, proj / origen, origen)
                if n:
                    incluidos.append(f"{origen}/ ({n} archivos)")
        if incluidos:
            sello = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
            zf.writestr("CONTENIDO.txt",
                        f"PoliScreen · proyecto '{proj.name}' · exportado {sello}\n\n"
                        + "\n".join(f"- {i}" for i in incluidos) + "\n")
    return buf.getvalue(), incluidos


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


ESTADO_UI = "estado_ui.json"


def leer_estado(archivo) -> dict:
    """Estado de la interfaz guardado en la sesión (productos construidos, nucleo, reactivos).
    Devuelve {} si la sesión no lo trae."""
    try:
        with zipfile.ZipFile(Path(archivo)) as zf:
            return json.loads(zf.read(ESTADO_UI).decode("utf-8"))
    except Exception:
        return {}


def save_session(proj, dest, completa: bool = False, notas: str = "", estado: Optional[dict] = None) -> Path:
    """Empaqueta la carpeta de proyecto en un archivo .poliscreen.

    completa=True incluye poses y complejos (pesado); si no, solo configuración, tablas y entradas.
    estado: instantanea de la interfaz (productos, nucleo, reactivos) para que al restaurar no haya
    que reconstruir la serie a mano.
    """
    proj = Path(proj)
    dest = Path(dest)
    if dest.suffix != EXT:
        dest = dest.with_suffix(EXT)
    dest.parent.mkdir(parents=True, exist_ok=True)

    incluido: list = []
    with zipfile.ZipFile(dest, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for nombre in ARCHIVOS_BASE:
            if _add_file(zf, proj / nombre, nombre):
                incluido.append(nombre)
        for carpeta in CARPETAS_BASE:
            n = _add_dir(zf, proj / carpeta, carpeta)
            if n:
                incluido.append(f"{carpeta}/ ({n} archivos)")
        if completa:
            for carpeta in CARPETAS_PESADAS:
                n = _add_dir(zf, proj / carpeta, carpeta)
                if n:
                    incluido.append(f"{carpeta}/ ({n} archivos)")

        if estado:
            zf.writestr(ESTADO_UI, json.dumps(estado, indent=2, ensure_ascii=False, default=str))
            incluido.append(ESTADO_UI)

        try:
            from .. import __version__ as ver
        except Exception:
            ver = "desconocida"
        manifiesto = {
            "formato": FORMATO,
            "poliscreen": ver,
            "creado": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "proyecto": proj.name,
            "completa": bool(completa),
            "contenido": incluido,
            "notas": notas,
        }
        zf.writestr("manifest.json", json.dumps(manifiesto, indent=2, ensure_ascii=False))
    return dest


def session_info(archivo) -> dict:
    """Lee el manifiesto sin extraer nada: sirve para mostrar que trae la sesión antes de abrirla."""
    with zipfile.ZipFile(Path(archivo)) as zf:
        try:
            return json.loads(zf.read("manifest.json").decode("utf-8"))
        except KeyError:
            return {"formato": None, "contenido": zf.namelist()[:50]}


def load_session(archivo, dest_dir) -> Path:
    """Extrae una sesión en dest_dir y devuelve la carpeta de proyecto lista para usar.

    Rechaza rutas que escapen del destino (protección frente a archivos manipulados).
    """
    archivo = Path(archivo)
    dest = Path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archivo) as zf:
        for miembro in zf.infolist():
            destino = (dest / miembro.filename).resolve()
            if not str(destino).startswith(str(dest.resolve())):
                raise ValueError(f"Ruta insegura en la sesion: {miembro.filename}")
        zf.extractall(dest)
    return dest
