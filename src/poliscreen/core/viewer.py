"""3D viewer. Returns self-contained HTML, so it works in any web interface.

Seeing the structure is the only cheap way to check that preparation did what it says: that the
waters are gone, that the cofactor is still there and that the ligand lands where it should.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

WATERS = {"HOH", "WAT", "H2O", "DOD"}
PROTEIN_COLOR = "spectrum"

_JS_LOCAL = Path(__file__).resolve().parent.parent / "assets" / "3Dmol-min.js"
_JS_CACHE: Optional[str] = None


# This travels inside an iframe of its own, and a browser gives every document an 8 px body
# margin. py3Dmol sizes its div in pixels, so those 16 px land outside the frame: a scrollbar
# appears and the viewer sits inset on all four sides, which is what makes the structure look
# squashed into a strip. Nothing in this document wants a margin.
_RESET = "<style>html,body{margin:0;padding:0}</style>\n"


def _embed_js(html: str) -> str:
    """Prepends the embedded library so py3Dmol's remote loader is never used.

    py3Dmol calls loadScriptAsync() and stores the promise in $3Dmolpromise, guarded by
    `if (typeof $3Dmolpromise === 'undefined')`. Defining that variable already resolved avoids the
    download without having to rewrite the generated code, which changes between versions.
    """
    global _JS_CACHE
    if not _JS_LOCAL.exists() or "$3Dmolpromise" not in html:
        return _RESET + html
    if _JS_CACHE is None:
        _JS_CACHE = _JS_LOCAL.read_text(errors="ignore")
    return _RESET + (f"<script>{_JS_CACHE}</script>\n"
                     "<script>var $3Dmolpromise = Promise.resolve();</script>\n") + html


def local_js_available() -> bool:
    return _JS_LOCAL.exists()


def _read(path) -> tuple:
    p = Path(path)
    fmt = {".pdb": "pdb", ".pdbqt": "pdb", ".sdf": "sdf", ".mol": "sdf", ".mol2": "mol2"}.get(p.suffix.lower())
    if fmt is None:
        raise ValueError(f"Unsupported format to visualize: {p.suffix}")
    return p.read_text(errors="ignore"), fmt


def _box_dims(box_) -> Optional[dict]:
    """Normalizes a box (Box or dict with cx..sz) to {center, dimensions} for py3Dmol."""
    if box_ is None:
        return None
    g = box_.as_dict() if hasattr(box_, "as_dict") else dict(box_)
    try:
        return {"center": {"x": float(g["cx"]), "y": float(g["cy"]), "z": float(g["cz"])},
                "dimensions": {"w": float(g["sx"]), "h": float(g["sy"]), "d": float(g["sz"])}}
    except (KeyError, TypeError, ValueError):
        return None


CHOSEN_COLOR = "#f4d35e"
CHOSEN_EMOJI = "🟡"
CAVITY_PALETTE = ["#4C9BE8", "#B06FD6", "#54B87A", "#E8823C", "#D65C7A", "#B4885E", "#E8E8E8"]
CAVITY_EMOJI = ["🔵", "🟣", "🟢", "🟠", "🔴", "🟤", "⚪"]

# Tunnels are told apart from cavities by being flat and saturated rather than soft: they are a
# route to follow, not a volume to judge, and a muted translucent tube is hard to trace through a
# ribbon. The order matches CAVER's own, so a tunnel keeps its colour between here and PyMOL.
TUNNEL_PALETTE = ["#E03131", "#2F9E44", "#1971C2", "#F1C40F",
                  "#E8590C", "#9C36B5", "#8D6E4A", "#DEE2E6"]
TUNNEL_EMOJI = ["🔴", "🟢", "🔵", "🟡", "🟠", "🟣", "🟤", "⚪"]


def emoji_for_color(color: str) -> str:
    """Colored circle equivalent to the one used in the 3D viewer."""
    if color == CHOSEN_COLOR:
        return CHOSEN_EMOJI
    for palette, emoji in ((CAVITY_PALETTE, CAVITY_EMOJI), (TUNNEL_PALETTE, TUNNEL_EMOJI)):
        try:
            return emoji[palette.index(color)]
        except ValueError:
            continue
    return "⚫"


def _xyz_axes(v, receptor, largo: float = 10.0):
    """Draws the X, Y, Z axes in a corner. Without them you cannot tell what moving cx, cy or cz means;
    with them the box is adjusted by looking, not guessing."""
    try:
        import numpy as np
        pts = np.array(_coords_pdb(receptor))
        if pts.size == 0:
            return
        o = pts.min(0) - 6.0
        for d, color, name_ in ((0, "#E05C5C", "X"), (1, "#5CE07A", "Y"), (2, "#5C8CE0", "Z")):
            fin = o.copy()
            fin[d] += largo
            v.addArrow({"start": {"x": float(o[0]), "y": float(o[1]), "z": float(o[2])},
                        "end": {"x": float(fin[0]), "y": float(fin[1]), "z": float(fin[2])},
                        "color": color, "radius": 0.28})
            v.addLabel(name_, {"position": {"x": float(fin[0]), "y": float(fin[1]), "z": float(fin[2])},
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


def view_html(receptor=None, ligand_=None, height_: int = 480, width_="100%",
              protein_style: str = "cartoon", show_waters: bool = True,
              show_hetero: bool = True, box_=None, pocket_spheres=None,
              cavities=None, opacity: float = 0.95, axes_: bool = False,
              surface: bool = False) -> str:
    """HTML of a viewer with receptor, ligand or both. Any of them can be omitted.

    box_: Box or dict cx..sz, drawn as a wireframe.
    pocket_spheres: [(x,y,z,radius)] alpha spheres (fpocket), drawn as the cavity volume.
    surface: molecular surface. The ribbon alone does not show whether a ligand sits inside the
    cavity or rests outside it.
    """
    import py3Dmol
    v = py3Dmol.view(width=width_, height=height_)
    if receptor is not None:
        text_, fmt = _read(receptor)
        v.addModel(text_, fmt)
        v.setStyle({"model": -1}, {protein_style: {"color": PROTEIN_COLOR}})
        if show_hetero:
            v.addStyle({"model": -1, "hetflag": True},
                       {"stick": {"radius": 0.18, "colorscheme": "greenCarbon"}})
        if show_waters:
            for w in WATERS:
                v.addStyle({"model": -1, "resn": w}, {"sphere": {"radius": 0.28, "color": "lightblue"}})
    if ligand_ is not None:
        text_, fmt = _read(ligand_)
        v.addModel(text_, fmt)
        v.setStyle({"model": -1}, {"stick": {"radius": 0.22, "colorscheme": "yellowCarbon"},
                                   "sphere": {"radius": 0.42, "colorscheme": "yellowCarbon"}})
    dims = _box_dims(box_)
    if dims:
        # Edges only: a filled box, however faint, hides the alpha spheres of the cavities it is
        # meant to frame. Thick lines make it readable without occluding anything.
        v.addBox({**dims, "color": "#FF3DDA", "opacity": 1.0, "wireframe": True, "linewidth": 10.0})
    groups_ = list(cavities or [])
    if pocket_spheres:
        groups_.append({"alpha": pocket_spheres, "color": CHOSEN_COLOR, "chosen": True})
    for g in groups_:
        col = g.get("color", CHOSEN_COLOR)
        # A group may fix its own opacity: cavities are judged through, so they stay translucent,
        # and a tunnel is followed, which a translucent tube makes harder rather than easier.
        op = g.get("opacity")
        if op is None:
            op = opacity if g.get("chosen") else max(0.55, opacity * 0.70)
        for s in (g.get("alpha") or []):
            try:
                x, y, z, r = float(s[0]), float(s[1]), float(s[2]), float(s[3])
            except (IndexError, ValueError, TypeError):
                continue
            v.addSphere({"center": {"x": x, "y": y, "z": z}, "radius": r,
                         "color": col, "opacity": op})
    if receptor is not None:
        v.zoomTo({"model": 0})
    else:
        v.zoomTo()
    if axes_ and receptor is not None:
        _xyz_axes(v, receptor)
    if surface and receptor is not None:
        import py3Dmol as _p3d
        v.addSurface(_p3d.VDW, {"opacity": 0.62, "color": "#b9c6d6"}, {"model": 0})
    return _embed_js(v._make_html())


def grid_png(smiles_list, legends=None, cols: int = 4, sub: int = 220) -> Optional[bytes]:
    """2D grid of several molecules (to check at a glance that the products are correct)."""
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
    """Application logo, if present in the package. It is looked up by name so that dropping the file
    into src/poliscreen/assets/ is enough, without touching code."""
    for name_ in ("logo.png", "logo.svg", "logo.webp", "logo.jpg", "logo.jpeg"):
        p = ASSETS / name_
        if p.exists() and p.stat().st_size > 0:
            return p
    return None


def wordmark_path() -> Optional[Path]:
    """Typographic wordmark (the written name), if in assets/.

    It is looked up by name, like the icon, to leave it code-free: embedding a typeface would tie the
    project to a font file and its license; an image is provided by whoever has it.
    """
    for name_ in ("titulo.svg", "titulo.png", "titulo.webp", "titulo.jpg"):
        p = ASSETS / name_
        if p.exists() and p.stat().st_size > 0:
            return p
    return None


def logo_svg(color: str = "#6b7280", size: int = 150) -> str:
    """Watermark for empty viewers: a monochrome molecular glyph. A single color so it does not
    compete with the content and works the same on a light or dark background."""
    ring = " ".join(f"{60 + 26 * __import__('math').cos(__import__('math').radians(a)):.1f},"
                      f"{60 + 26 * __import__('math').sin(__import__('math').radians(a)):.1f}"
                      for a in range(0, 360, 60))
    return f"""
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 120 120" width="{size}" height="{size}"
     fill="none" stroke="{color}" stroke-width="2.4" opacity="0.55">
  <polygon points="{ring}" stroke-linejoin="round"/>
  <line x1="86" y1="60" x2="106" y2="48"/>
  <line x1="34" y1="60" x2="14" y2="72"/>
  <line x1="60" y1="86" x2="60" y2="108"/>
  <circle cx="106" cy="48" r="5.5" fill="{color}" stroke="none"/>
  <circle cx="14"  cy="72" r="5.5" fill="{color}" stroke="none"/>
  <circle cx="60"  cy="108" r="5.5" fill="{color}" stroke="none"/>
  <circle cx="60"  cy="60" r="3" fill="{color}" stroke="none" opacity="0.6"/>
</svg>"""


def molecule_png(smiles: str, size: int = 300) -> Optional[bytes]:
    """2D drawing of a molecule from its SMILES."""
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
    """2D drawing with atom indices and, highlighted, the reactive site. Basis of the site selector."""
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


def structure_summary(path) -> dict:
    """Counts atoms, waters, hetero groups and chains. Used for the before and after of preparation."""
    p = Path(path)
    chains_, het, n_atoms_, n_waters_, n_h = set(), {}, 0, 0, 0
    for l in p.read_text(errors="ignore").splitlines():
        if l.startswith(("ATOM", "HETATM")):
            n_atoms_ += 1
            if (l[76:78].strip() or "") == "H":
                n_h += 1
            rn = l[17:20].strip()
            if l.startswith("ATOM"):
                chains_.add(l[21].strip() or "_")
            elif rn in WATERS:
                n_waters_ += 1
            else:
                het[rn] = het.get(rn, 0) + 1
    return {"atomos": n_atoms_, "hidrogenos": n_h, "aguas": n_waters_,
            "chains": sorted(chains_), "heterogrupos": het}
