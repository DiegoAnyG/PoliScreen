"""Visor 3D. Devuelve HTML autocontenido, así que sirve en cualquier interfaz web.

Ver la estructura es la única forma barata de comprobar que la preparación hizo lo que dice:
que se fueron las aguas, que el cofactor sigue ahi y que el ligando cae donde debe.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

WATERS = {"HOH", "WAT", "H2O", "DOD"}
COLOR_PROTEINA = "spectrum"

# py3Dmol carga 3Dmol desde un CDN, lo que ata el visor a la red: si la descarga tarda o se bloquea,
# el panel sale en blanco sin mensaje, de forma intermitente. Con una copia local se trabaja sin
# conexión y se fija la versión usada en un análisis publicado.
_JS_LOCAL = Path(__file__).resolve().parent.parent / "assets" / "3Dmol-min.js"
_JS_CACHE: Optional[str] = None


def _incrustar_js(html: str) -> str:
    """Antepone la librería incrustada para que el cargador remoto de py3Dmol no llegue a usarse.

    py3Dmol llama a loadScriptAsync() y guarda la promesa en $3Dmolpromise, protegida por
    `if (typeof $3Dmolpromise === 'undefined')`. Definir esa variable ya resuelta evita la descarga
    sin tener que reescribir el código generado, que cambia entre versiones.
    """
    global _JS_CACHE
    if not _JS_LOCAL.exists() or "$3Dmolpromise" not in html:
        return html
    if _JS_CACHE is None:
        _JS_CACHE = _JS_LOCAL.read_text(errors="ignore")
    return (f"<script>{_JS_CACHE}</script>\n"
            "<script>var $3Dmolpromise = Promise.resolve();</script>\n") + html


def js_local_disponible() -> bool:
    return _JS_LOCAL.exists()


def _leer(path) -> tuple:
    p = Path(path)
    fmt = {".pdb": "pdb", ".pdbqt": "pdb", ".sdf": "sdf", ".mol": "sdf", ".mol2": "mol2"}.get(p.suffix.lower())
    if fmt is None:
        raise ValueError(f"Formato no soportado para visualizar: {p.suffix}")
    return p.read_text(errors="ignore"), fmt


def _box_dims(caja) -> Optional[dict]:
    """Normaliza una caja (Box o dict con cx..sz) a {centro, dimensiones} para py3Dmol."""
    if caja is None:
        return None
    g = caja.as_dict() if hasattr(caja, "as_dict") else dict(caja)
    try:
        return {"center": {"x": float(g["cx"]), "y": float(g["cy"]), "z": float(g["cz"])},
                "dimensions": {"w": float(g["sx"]), "h": float(g["sy"]), "d": float(g["sz"])}}
    except (KeyError, TypeError, ValueError):
        return None


# Colores para distinguir cavidades a la vez. El primero se reserva a la elegida para acoplar.
COLOR_ELEGIDA = "#f4d35e"
EMOJI_ELEGIDA = "🟡"
# Cada color va emparejado con un circulo del mismo tono: la tabla puede mostrar el circulo y
# el usuario reconoce la cavidad en el visor sin leer un código hexadecimal.
PALETA_CAVIDADES = ["#4C9BE8", "#B06FD6", "#54B87A", "#E8823C", "#D65C7A", "#B4885E", "#E8E8E8"]
EMOJI_CAVIDADES = ["🔵", "🟣", "🟢", "🟠", "🔴", "🟤", "⚪"]


def emoji_de_color(color: str) -> str:
    """Circulo de color equivalente al usado en el visor 3D."""
    if color == COLOR_ELEGIDA:
        return EMOJI_ELEGIDA
    try:
        return EMOJI_CAVIDADES[PALETA_CAVIDADES.index(color)]
    except ValueError:
        return "⚫"


def _ejes_xyz(v, receptor, largo: float = 10.0):
    """Dibuja los ejes X, Y, Z en una esquina. Sin ellos no se sabe qué significa mover cx, cy o cz;
    con ellos la caja se ajusta mirando, no adivinando."""
    try:
        import numpy as np
        pts = np.array(_coords_pdb(receptor))
        if pts.size == 0:
            return
        o = pts.min(0) - 6.0
        for d, color, nombre in ((0, "#E05C5C", "X"), (1, "#5CE07A", "Y"), (2, "#5C8CE0", "Z")):
            fin = o.copy()
            fin[d] += largo
            v.addArrow({"start": {"x": float(o[0]), "y": float(o[1]), "z": float(o[2])},
                        "end": {"x": float(fin[0]), "y": float(fin[1]), "z": float(fin[2])},
                        "color": color, "radius": 0.28})
            v.addLabel(nombre, {"position": {"x": float(fin[0]), "y": float(fin[1]), "z": float(fin[2])},
                                "fontSize": 13, "fontColor": color, "backgroundOpacity": 0.0})
    except Exception:
        pass


def _coords_pdb(path) -> list:
    pts = []
    for l in Path(path).read_text(errors="ignore").splitlines():
        if l.startswith(("ATOM", "HETATM")):
            try:
                pts.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
            except ValueError:
                continue
    return pts


def view_html(receptor=None, ligando=None, alto: int = 480, ancho="100%",
              estilo_proteina: str = "cartoon", mostrar_aguas: bool = True,
              mostrar_hetero: bool = True, caja=None, pocket_spheres=None,
              cavidades=None, opacidad: float = 0.95, ejes: bool = False,
              superficie: bool = False) -> str:
    """HTML de un visor con receptor, ligando o ambos. Cualquiera puede omitirse.

    caja: si se pasa (Box o dict cx,cy,cz,sx,sy,sz), dibuja la caja de docking en malla.
    pocket_spheres: lista [(x,y,z,radio)] de esferas-alfa (fpocket); se renderizan translucidas
    como el VOLUMEN de la cavidad (isosuperficie tipo CaverWeb), no una caja rectangular.
    superficie: superficie molecular del receptor, translucida. Con la cinta sola no se distingue
    si un ligando queda dentro de la cavidad o apoyado por fuera, porque la cinta solo dibuja el
    esqueleto y deja pasar la vista entre las cadenas laterales. Es la comprobación visual
    inmediata para un péptido, que por su tamaño puede quedar mayormente expuesto.
    """
    import py3Dmol
    v = py3Dmol.view(width=ancho, height=alto)
    if receptor is not None:
        texto, fmt = _leer(receptor)
        v.addModel(texto, fmt)
        v.setStyle({"model": -1}, {estilo_proteina: {"color": COLOR_PROTEINA}})
        if mostrar_hetero:
            # Heteroatomos que no son agua: cofactores, iones, ligandos co-cristalizados
            v.addStyle({"model": -1, "hetflag": True},
                       {"stick": {"radius": 0.18, "colorscheme": "greenCarbon"}})
        if mostrar_aguas:
            for w in WATERS:
                v.addStyle({"model": -1, "resn": w}, {"sphere": {"radius": 0.28, "color": "lightblue"}})
    if ligando is not None:
        texto, fmt = _leer(ligando)
        v.addModel(texto, fmt)
        # Bolas y varillas siempre, también cuando el ligando es un péptido. Un péptido llega en
        # formato PDB con residuos de aminoacido, así que dibujado como proteína saldría en cinta y
        # sería indistinguible del receptor; lo que interesa de un ligando son sus átomos.
        v.setStyle({"model": -1}, {"stick": {"radius": 0.22, "colorscheme": "yellowCarbon"},
                                   "sphere": {"radius": 0.42, "colorscheme": "yellowCarbon"}})
    dims = _box_dims(caja)
    if dims:
        # Malla gruesa y opaca: con líneas finas y translucidas la caja se perdía sobre el cartoon.
        v.addBox({**dims, "color": "#FF3DDA", "opacity": 1.0, "wireframe": True, "linewidth": 6.0})
    # Cavidades: lista de {alpha, color, elegida}. Se dibujan todas a la vez, cada una de un color,
    # y la elegida para acoplar se resalta más opaca para distinguirla de un vistazo.
    grupos = list(cavidades or [])
    if pocket_spheres:
        grupos.append({"alpha": pocket_spheres, "color": COLOR_ELEGIDA, "elegida": True})
    for g in grupos:
        col = g.get("color", COLOR_ELEGIDA)
        op = opacidad if g.get("elegida") else max(0.55, opacidad * 0.70)
        for s in (g.get("alpha") or []):
            try:
                x, y, z, r = float(s[0]), float(s[1]), float(s[2]), float(s[3])
            except (IndexError, ValueError, TypeError):
                continue
            v.addSphere({"center": {"x": x, "y": y, "z": z}, "radius": r,
                         "color": col, "opacity": op})
    # El encuadre se hace SOBRE EL RECEPTOR (modelo 0), no sobre la escena completa: si se
    # incluyeran los ejes o las cavidades, el centro de rotacion se desplazaria fuera de la
    # proteína. Los ejes y la caja son objetos de dibujo; no intervienen en ningun cálculo.
    if receptor is not None:
        v.zoomTo({"model": 0})
    else:
        v.zoomTo()
    if ejes and receptor is not None:
        _ejes_xyz(v, receptor)
    # La superficie se anade al final y solo sobre el modelo del receptor: incluir el ligando la
    # englobaria y dejaría de verse si esta dentro o fuera, que es justo lo que se quiere mirar.
    if superficie and receptor is not None:
        import py3Dmol as _p3d
        v.addSurface(_p3d.VDW, {"opacity": 0.62, "color": "#b9c6d6"}, {"model": 0})
    return _incrustar_js(v._make_html())


def grid_png(smiles_list, legends=None, cols: int = 4, sub: int = 220) -> Optional[bytes]:
    """Rejilla 2D de varias moléculas (para verificar de un vistazo que los productos son correctos)."""
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import Draw
        RDLogger.DisableLog("rdApp.*")
        mols, legs = [], []
        for i, smi in enumerate(smiles_list):
            m = Chem.MolFromSmiles(str(smi)) if smi else None
            if m is None:
                continue
            mols.append(m)
            legs.append(str(legends[i]) if legends and i < len(legends) else "")
        if not mols:
            return None
        img = Draw.MolsToGridImage(mols, molsPerRow=max(1, cols), subImgSize=(sub, sub),
                                   legends=legs, useSVG=False)
        import io
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return buf.getvalue()
    except Exception:
        return None


ASSETS = Path(__file__).resolve().parent.parent / "assets"


def logo_path() -> Optional[Path]:
    """Logo de la aplicación, si existe en el paquete. Se busca por nombre para que baste con
    depositar el archivo en src/poliscreen/assets/ sin tocar código."""
    for nombre in ("logo.png", "logo.svg", "logo.webp", "logo.jpg", "logo.jpeg"):
        p = ASSETS / nombre
        if p.exists() and p.stat().st_size > 0:
            return p
    return None


def wordmark_path() -> Optional[Path]:
    """Logotipo tipográfico (el nombre escrito), si está en assets/.

    Se busca por nombre, como el icono, para dejarlo sin tocar código: incrustar una tipografía
    ataría el proyecto a un archivo de fuente y su licencia; una imagen la aporta quien la tiene.
    """
    for nombre in ("titulo.svg", "titulo.png", "titulo.webp", "titulo.jpg"):
        p = ASSETS / nombre
        if p.exists() and p.stat().st_size > 0:
            return p
    return None


def logo_svg(color: str = "#6b7280", size: int = 150) -> str:
    """Marca de agua para los visores vacios: un glifo molecular monocromo. Un solo color para que
    no compita con el contenido y funcione igual sobre fondo claro u oscuro."""
    anillo = " ".join(f"{60 + 26 * __import__('math').cos(__import__('math').radians(a)):.1f},"
                      f"{60 + 26 * __import__('math').sin(__import__('math').radians(a)):.1f}"
                      for a in range(0, 360, 60))
    return f"""
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 120 120" width="{size}" height="{size}"
     fill="none" stroke="{color}" stroke-width="2.4" opacity="0.55">
  <polygon points="{anillo}" stroke-linejoin="round"/>
  <line x1="86" y1="60" x2="106" y2="48"/>
  <line x1="34" y1="60" x2="14" y2="72"/>
  <line x1="60" y1="86" x2="60" y2="108"/>
  <circle cx="106" cy="48" r="5.5" fill="{color}" stroke="none"/>
  <circle cx="14"  cy="72" r="5.5" fill="{color}" stroke="none"/>
  <circle cx="60"  cy="108" r="5.5" fill="{color}" stroke="none"/>
  <circle cx="60"  cy="60" r="3" fill="{color}" stroke="none" opacity="0.6"/>
</svg>"""


def molecule_png(smiles: str, size: int = 300) -> Optional[bytes]:
    """Dibujo 2D de una molécula a partir de su SMILES."""
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        m = Chem.MolFromSmiles(smiles)
        return _png(m, size) if m is not None else None
    except Exception:
        return None


def _png(mol, size: int) -> bytes:
    from rdkit.Chem.Draw import rdMolDraw2D
    d = rdMolDraw2D.MolDraw2DCairo(size, size)
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    return d.GetDrawingText()


def molecule_png_indexed(smiles: str, highlight=None, size: int = 360) -> Optional[bytes]:
    """Dibujo 2D con los indices de átomo y, resaltado, el sitio reactivo. Base del selector de sitio."""
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem.Draw import rdMolDraw2D
        RDLogger.DisableLog("rdApp.*")
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return None
        d = rdMolDraw2D.MolDraw2DCairo(size, size)
        d.drawOptions().addAtomIndices = True
        hl = [int(a) for a in (highlight or []) if a is not None and int(a) < m.GetNumAtoms()]
        rdMolDraw2D.PrepareAndDrawMolecule(d, m, highlightAtoms=hl)
        d.FinishDrawing()
        return d.GetDrawingText()
    except Exception:
        return None


def resumen_estructura(path) -> dict:
    """Cuenta átomos, aguas, heterogrupos y cadenas. Sirve para el antes y después de preparar."""
    p = Path(path)
    cadenas, het, n_atomos, n_aguas, n_h = set(), {}, 0, 0, 0
    for l in p.read_text(errors="ignore").splitlines():
        if l.startswith(("ATOM", "HETATM")):
            n_atomos += 1
            if (l[76:78].strip() or "") == "H":
                n_h += 1
            rn = l[17:20].strip()
            if l.startswith("ATOM"):
                cadenas.add(l[21].strip() or "_")
            elif rn in WATERS:
                n_aguas += 1
            else:
                het[rn] = het.get(rn, 0) + 1
    return {"atomos": n_atomos, "hidrogenos": n_h, "aguas": n_aguas,
            "cadenas": sorted(cadenas), "heterogrupos": het}
