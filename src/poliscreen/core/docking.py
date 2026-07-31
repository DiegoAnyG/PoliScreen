"""Docking: preparación de archivos, caja de búsqueda y motores.

Reproducible: semilla fija y un hilo por acoplamiento (Vina no es determinista con multihilo). El
paralelismo va entre acoplamientos independientes, no dentro. Sin dianas ni rutas cableadas.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence

WATERS = {"HOH", "WAT", "H2O", "DOD"}


class DockingError(RuntimeError):
    pass


# ---------------------------------------------------------------- caja
@dataclass
class Box:
    cx: float
    cy: float
    cz: float
    sx: float
    sy: float
    sz: float

    @classmethod
    def around(cls, points, pad: float = 10.0, lo: float = 16.0, hi: float = 30.0) -> "Box":
        import numpy as np
        a = np.asarray(points, dtype=float)
        c = a.mean(0)
        size = float(min(max(float((a.max(0) - a.min(0)).max()) + pad, lo), hi))
        return cls(*[round(float(v), 2) for v in c], size, size, size)

    def args(self) -> list:
        return ["--center_x", str(self.cx), "--center_y", str(self.cy), "--center_z", str(self.cz),
                "--size_x", str(self.sx), "--size_y", str(self.sy), "--size_z", str(self.sz)]

    def as_dict(self) -> dict:
        return dict(cx=self.cx, cy=self.cy, cz=self.cz, sx=self.sx, sy=self.sy, sz=self.sz)


def _coords(pdb, het_only=False):
    pts = []
    for line in Path(pdb).read_text(errors="ignore").splitlines():
        tag = line[:6].strip()
        if tag in (("HETATM",) if het_only else ("ATOM", "HETATM")):
            if het_only and line[17:20].strip() in WATERS:
                continue
            try:
                pts.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                continue
    return pts


def hetero_groups(pdb) -> dict:
    """Ligandos co-cristalizados (HETATM no agua) -> etiqueta: Box centrada en ellos.

    Sirve para centrar la caja en el sitio correcto cuando una proteína tiene varios.
    """
    groups = {}
    for line in Path(pdb).read_text(errors="ignore").splitlines():
        if not line.startswith("HETATM"):
            continue
        rn = line[17:20].strip()
        if rn in WATERS:
            continue
        key = (rn, (line[21].strip() or "_"), line[22:26].strip())
        try:
            groups.setdefault(key, []).append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        except ValueError:
            continue
    return {f"{k[0]} {k[1]}:{k[2]} ({len(v)} at.)": Box.around(v) for k, v in groups.items()}


def coords_from_file(path) -> list:
    """Coordenadas de cualquier formato habitual de ligando o estructura."""
    p = Path(path)
    suf = p.suffix.lower()
    text = p.read_text(errors="ignore")
    pts = []
    if suf in (".pdb", ".pdbqt", ".ent"):
        for l in text.splitlines():
            if l.startswith(("ATOM", "HETATM")):
                try:
                    pts.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
                except ValueError:
                    continue
    elif suf == ".mol2":
        dentro = False
        for l in text.splitlines():
            if l.startswith("@<TRIPOS>ATOM"):
                dentro = True
                continue
            if l.startswith("@<TRIPOS>"):
                dentro = False
            if dentro:
                f = l.split()
                if len(f) >= 5:
                    try:
                        pts.append((float(f[2]), float(f[3]), float(f[4])))
                    except ValueError:
                        continue
    elif suf in (".sdf", ".mol"):
        lines = text.splitlines()
        if len(lines) > 3:
            try:
                n = int(lines[3][0:3])
            except ValueError:
                n = 0
            for l in lines[4:4 + n]:
                f = l.split()
                if len(f) >= 4:
                    try:
                        pts.append((float(f[0]), float(f[1]), float(f[2])))
                    except ValueError:
                        continue
    if not pts:
        raise DockingError(f"No pude leer coordenadas de {p.name}.")
    return pts


def residues_in_box(receptor_pdb, box, pad: float = 0.0) -> set:
    """Residuos (etiqueta 'Tyr157', como los nombra PLIP) con al menos un átomo dentro de la caja.
    Define el 'pocket' objetivo: el ranking premia contactos productivos con CUALQUIERA de estos,
    no solo con los que toca el control."""
    hx, hy, hz = box.sx / 2 + pad, box.sy / 2 + pad, box.sz / 2 + pad
    res = set()
    for l in Path(receptor_pdb).read_text(errors="ignore").splitlines():
        if not l.startswith(("ATOM", "HETATM")):
            continue
        try:
            x, y, z = float(l[30:38]), float(l[38:46]), float(l[46:54])
        except ValueError:
            continue
        if abs(x - box.cx) <= hx and abs(y - box.cy) <= hy and abs(z - box.cz) <= hz:
            rt, rn = l[17:20].strip(), l[22:26].strip()
            if rt and rn and rt not in WATERS:
                res.add(f"{rt.capitalize()}{rn}")
    return res


def box_from_file(path, pad: float = 10.0, lo: float = 16.0, hi: float = 30.0) -> Box:
    """Caja centrada en un ligando concreto: lo más fiable cuando se tiene el co-cristalizado."""
    return Box.around(coords_from_file(path), pad=pad, lo=lo, hi=hi)


def load_boxes_xlsx(path, receptors: Sequence) -> dict:
    """Carga cajas de un xlsx con columnas receptor,cx,cy,cz,sx,sy,sz (formato de MolModa).

    Empareja por nombre de archivo sin extensión; ignora filas incompletas.
    """
    import pandas as pd
    t = pd.read_excel(path)
    cols = {str(c).lower().strip(): c for c in t.columns}
    need = ("receptor", "cx", "cy", "cz", "sx", "sy", "sz")
    if any(k not in cols for k in need):
        raise DockingError(f"Al xlsx le faltan columnas. Se requieren: {', '.join(need)}.")
    porclave = {}
    for _, row in t.iterrows():
        try:
            vals = {k: float(row[cols[k]]) for k in need[1:]}
        except (TypeError, ValueError):
            continue
        if any(v != v for v in vals.values()):
            continue
        porclave[Path(str(row[cols["receptor"]])).stem.lower()] = Box(**vals)
    out = {}
    for r in receptors:
        b = porclave.get(Path(r).stem.lower())
        if b:
            out[str(r)] = b
    return out


def auto_box(pdb) -> Box:
    """Caja automática: centrada en el ligando co-cristalizado si lo hay, si no en la proteína."""
    het = _coords(pdb, het_only=True)
    if len(het) >= 3:
        return Box.around(het)
    pts = _coords(pdb)
    if not pts:
        raise DockingError(f"{Path(pdb).name}: sin coordenadas legibles.")
    return Box.around(pts, pad=0.0, lo=24.0, hi=24.0)


# ---------------------------------------------------------------- preparación
def to_pdbqt(src, dst, receptor: bool = False, ph: float = 7.4) -> bool:
    """Convierte a pdbqt con OpenBabel. Receptor: rigido; ligando: cargas Gasteiger."""
    src, dst = Path(src), Path(dst)
    if src.suffix.lower() == ".pdbqt":
        shutil.copy(src, dst)
        return dst.exists()
    dst.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["obabel", str(src), "-O", str(dst)]
    # -r (ligando) conserva solo el fragmento conectado mayor: un ligando extraído de un PDB puede
    # salir fragmentado o duplicado y dar un pdbqt multi-molécula que Vina rechaza.
    cmd += ["-xr", "-p", str(ph)] if receptor else ["-r", "-p", str(ph), "--partialcharge", "gasteiger"]
    if src.suffix.lower() == ".smi":
        cmd += ["--gen3d"]
    subprocess.run(cmd, capture_output=True, text=True)
    return dst.exists() and dst.stat().st_size > 0


def scores_from_pdbqt(path) -> list:
    """Energías de las poses, leidas de las líneas REMARK del propio pdbqt."""
    out = []
    for line in Path(path).read_text(errors="ignore").splitlines():
        if line.startswith("REMARK VINA RESULT"):
            try:
                out.append(float(line.split(":")[1].split()[0]))
            except Exception:
                pass
    return out


def split_models(pose_pdbqt, out_dir) -> list:
    """Separa el pdbqt multi-pose en un PDB por modelo."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = Path(pose_pdbqt).stem
    existing = sorted(out_dir.glob(f"{stem}-model*.pdb"))
    if existing:
        return existing
    subprocess.run(["obabel", str(pose_pdbqt), "-opdb", "-O", str(out_dir / f"{stem}-model.pdb"), "-m"],
                   capture_output=True)
    return sorted(out_dir.glob(f"{stem}-model*.pdb"))


# ---------------------------------------------------------------- motores
class VinaEngine:
    """AutoDock Vina en CPU. Determinista con semilla fija y un hilo."""

    name = "vina"

    def __init__(self, exe: Optional[str] = None, cpu: int = 1, seed: int = 42,
                 exhaustiveness: int = 24, n_poses: int = 10, energy_range: float = 3.0):
        self.exe = exe or os.environ.get("POLISCREEN_VINA") or shutil.which("vina")
        self.cpu, self.seed = cpu, seed
        self.exhaustiveness, self.n_poses, self.energy_range = exhaustiveness, n_poses, energy_range

    def available(self) -> bool:
        return bool(self.exe) and Path(self.exe).exists()

    def dock(self, receptor_pdbqt, ligand_pdbqt, box: Box, out_pdbqt) -> list:
        if not self.available():
            raise DockingError("No encuentro el ejecutable de vina. Define POLISCREEN_VINA o ponlo en el PATH.")
        Path(out_pdbqt).parent.mkdir(parents=True, exist_ok=True)
        cmd = ([self.exe, "--receptor", str(receptor_pdbqt), "--ligand", str(ligand_pdbqt)]
               + box.args()
               + ["--exhaustiveness", str(self.exhaustiveness), "--num_modes", str(self.n_poses),
                  "--energy_range", str(self.energy_range),
                  "--seed", str(self.seed), "--cpu", str(self.cpu), "--out", str(out_pdbqt)])
        r = subprocess.run(cmd, capture_output=True, text=True)
        if not Path(out_pdbqt).exists():
            raise DockingError(f"vina no produjo salida. Causa probable: caja fuera del sitio o receptor mal "
                               f"preparado. stderr: {(r.stderr or '')[-300:]}")
        return scores_from_pdbqt(out_pdbqt)


class GninaEngine(VinaEngine):
    """gnina: motor derivado de Vina cuya puntuación la da una red neuronal convolucional
    entrenada sobre complejos cristalograficos, con aceleracion por GPU.

    Su interes en PoliScreen no es la velocidad sino ser una función de puntuación INDEPENDIENTE
    de la de Vina: dos metodos distintos que coinciden son más creibles que uno solo.
    """

    name = "gnina"

    def __init__(self, exe: Optional[str] = None, cpu: int = 1, seed: int = 42,
                 exhaustiveness: int = 24, n_poses: int = 10, energy_range: float = 3.0,
                 usar_gpu: bool = True):
        super().__init__(exe=exe or gnina_exe(), cpu=cpu, seed=seed,
                         exhaustiveness=exhaustiveness, n_poses=n_poses, energy_range=energy_range)
        self.usar_gpu = usar_gpu

    def dock(self, receptor_pdbqt, ligand_pdbqt, box: Box, out_pdbqt) -> list:
        if not self.available():
            raise DockingError("No encuentro gnina. Define POLISCREEN_GNINA o ponlo en el PATH.")
        Path(out_pdbqt).parent.mkdir(parents=True, exist_ok=True)
        cmd = ([self.exe, "-r", str(receptor_pdbqt), "-l", str(ligand_pdbqt)]
               + box.args()
               + ["--exhaustiveness", str(self.exhaustiveness), "--num_modes", str(self.n_poses),
                  "--seed", str(self.seed), "--cpu", str(self.cpu), "-o", str(out_pdbqt)])
        if not self.usar_gpu:
            cmd.append("--no_gpu")
        r = subprocess.run(cmd, capture_output=True, text=True)
        if not Path(out_pdbqt).exists():
            raise DockingError(f"gnina no produjo salida. stderr: {(r.stderr or '')[-300:]}")
        return scores_from_pdbqt(out_pdbqt)


def gnina_exe() -> Optional[str]:
    """Ruta al ejecutable de gnina: variable de entorno, envoltorio local, PATH o binario suelto.

    El envoltorio va antes que el binario: el ejecutable oficial no es autocontenido (necesita cuDNN
    y librerías de CUDA en LD_LIBRARY_PATH) y aborta al arrancar sin ellas; el envoltorio las
    configura. Priorizar el binario lo daba por disponible y fallaba solo al usarlo.
    """
    env = os.environ.get("POLISCREEN_GNINA")
    if env and Path(env).exists():
        return env
    base = Path.home() / "poliscreen_tools"
    for nombre in ("gnina-run", "gnina"):
        local = base / nombre
        if local.exists() and os.access(local, os.X_OK):
            return str(local)
    return shutil.which("gnina")


def gnina_available() -> bool:
    return gnina_exe() is not None


_RE_PUNTUACION = re.compile(r"^\s*(Affinity|CNNscore|CNNaffinity|CNN_VS)\s*:?\s*([-\d.]+)", re.M)


def rescore_poses(receptor, poses: Sequence, usar_gpu: bool = True, timeout: int = 900) -> dict:
    """Re-puntúa con gnina poses ya generadas por Vina, sin volver a buscar.

    Segunda opinión barata y reproducible: el muestreo sigue siendo el de Vina; solo cambia la
    función que evalúa cada pose. Devuelve {nombre_pose: {cnn_score, cnn_affinity, affinity}}, con
    cnn_score como probabilidad de pose correcta (0-1) y cnn_affinity en unidades de pK.
    """
    exe = gnina_exe()
    out = {}
    if not exe or not poses:
        return out
    for p in poses:
        p = Path(p)
        if not p.exists():
            continue
        cmd = [exe, "-r", str(receptor), "-l", str(p), "--score_only"]
        if not usar_gpu:
            cmd.append("--no_gpu")
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            continue
        vals = {}
        for clave, valor in _RE_PUNTUACION.findall(r.stdout or ""):
            try:
                vals.setdefault(clave, float(valor))
            except ValueError:
                pass
        if vals:
            out[p.stem] = {"affinity": vals.get("Affinity"),
                           "cnn_score": vals.get("CNNscore"),
                           "cnn_affinity": vals.get("CNNaffinity")}
    return out


ENGINES = {"vina": VinaEngine, "gnina": GninaEngine}


# ---------------------------------------------------------------- lote
def extension_maxima(archivos: Sequence) -> float:
    """Eje mayor del mayor de los ligandos, en ángstrom. 0 si no se puede medir.

    Dimensiona la caja: el ligando no solo debe caber, debe poder reorientarse. Una caja apenas
    mayor limita la búsqueda a las orientaciones que entran, y el mal resultado parece del
    acoplamiento cuando es una restricción geométrica impuesta sin querer.
    """
    mayor = 0.0
    for f in archivos:
        try:
            pts = coords_from_file(f)
        except Exception:
            continue
        if not pts:
            continue
        ejes = [max(p[i] for p in pts) - min(p[i] for p in pts) for i in range(3)]
        mayor = max(mayor, max(ejes))
    return mayor


def caja_minima(archivos: Sequence, margen: float = 4.0) -> float:
    """Lado mínimo recomendable para que el mayor de esos ligandos pueda girar dentro de la caja."""
    ext = extension_maxima(archivos)
    return round(ext + margen, 1) if ext else 0.0


def memoria_disponible_gb() -> Optional[float]:
    """GB realmente disponibles, leidos de /proc/meminfo. None si no se puede saber."""
    try:
        for l in Path("/proc/meminfo").read_text().splitlines():
            if l.startswith("MemAvailable:"):
                return int(l.split()[1]) / 1048576.0
    except Exception:
        pass
    return None


TORSDOF_LIMITE = 15


def torsdof(pdbqt) -> int:
    """Grados de libertad torsionales que Vina declara para un ligando ya preparado."""
    try:
        for l in Path(pdbqt).read_text(errors="ignore").splitlines():
            if l.startswith("TORSDOF"):
                return int(l.split()[1])
    except Exception:
        pass
    return 0


def coste_memoria_gb(box: "Box", tors: int = 6) -> float:
    """Memoria estimada de un acoplamiento, en GB.

    Dos sumandos: los mapas de afinidad, que crecen con el volumen de la caja (rejilla de 0,375 Å
    por tipo de átomo), y el árbol de conformaciones, que domina en cuanto el ligando deja de ser
    pequeño y se modela cuadrático en los grados de libertad, porque los estados vivos durante la
    búsqueda crecen mucho más deprisa que las torsiones.

    Cota generosa a propósito: pasarse cuesta tiempo; quedarse corto hace que el sistema mate el
    proceso a media tanda y se pierda todo el trabajo.
    """
    puntos = ((box.sx / 0.375) + 1) * ((box.sy / 0.375) + 1) * ((box.sz / 0.375) + 1)
    mapas = puntos * 22 * 4 / 1e9 * 3.0          # ~22 tipos de átomo, float de 4 bytes
    conformaciones = 0.05 * max(1, tors) ** 2 / 36.0
    return max(0.30, mapas + conformaciones)


def paralelismo_seguro(cajas: Sequence["Box"], ligandos_pdbqt: Sequence,
                       reserva_gb: float = 2.0) -> int:
    """Cuántos acoplamientos lanzar a la vez sin agotar la memoria.

    Repartir solo por núcleos es lo que hace que el sistema mate el proceso: muchos hilos con poca
    RAM es corriente (WSL toma por defecto la mitad del equipo, y 16 núcleos conviven con 7 GB). Se
    usa el límite más restrictivo entre núcleos y memoria disponible, con margen para el intérprete.
    """
    por_nucleos = max(1, (os.cpu_count() or 2) - 2)
    libre = memoria_disponible_gb()
    if not libre:
        return por_nucleos
    cajas = list(cajas) or [Box(0, 0, 0, 24, 24, 24)]
    mayor = max(cajas, key=lambda b: b.sx * b.sy * b.sz)
    tors = max([torsdof(l) for l in ligandos_pdbqt if l] or [6])
    coste = coste_memoria_gb(mayor, tors)
    por_memoria = max(1, int((libre - reserva_gb) / coste))
    return max(1, min(por_nucleos, por_memoria))


def pose_name(receptor_stem: str, ligand_stem: str) -> str:
    return f"docking_{receptor_stem}_compounds_a_{ligand_stem.lower()}"


def dock_batch(receptors: Sequence, ligands: Sequence, boxes: dict, work_dir,
               engine=None, workers: int = 0, on_progress=None, ph: float = 7.4, targets=None) -> list:
    """Acopla cada ligando contra cada sitio. Reanudable: salta lo que ya existe en disco.

    Un sitio es (ruta_receptor, id_sitio, Box). Por defecto, uno por receptor derivado de `boxes`.
    Con `targets` se acopla el mismo ligando en varios bolsillos del mismo receptor (híbrido); el
    id_sitio los distingue y aparece como 'receptor' en los resultados, separando el ranking por
    sitio. ph protona receptor y ligando al pasar a pdbqt. Devuelve una fila por pose.
    """
    engine = engine or VinaEngine()
    work = Path(work_dir)
    prep, poses = work / "prep", work / "poses"
    for d in (prep, poses):
        d.mkdir(parents=True, exist_ok=True)

    if targets is None:
        targets = [(r, Path(r).stem, boxes[str(r)]) for r in receptors if str(r) in boxes]

    # Un pdbqt por archivo de receptor, aunque varios sitios lo compartan.
    rec_pdbqt = {}
    for rpath, _sid, _box in targets:
        if str(rpath) in rec_pdbqt:
            continue
        dst = prep / f"{Path(rpath).stem}.pdbqt"
        rec_pdbqt[str(rpath)] = dst if (dst.exists() or to_pdbqt(rpath, dst, receptor=True, ph=ph)) else None

    # Ligandos a pdbqt ANTES de paralelizar: con varios sitios el mismo ligando entra en varias
    # tareas, y dos hilos escribiéndolo a la vez lo dejarían a medias, perdiéndolo en ese sitio.
    lig_pdbqt = {}
    for l in ligands:
        dst = prep / f"{Path(l).stem}.pdbqt"
        lig_pdbqt[str(l)] = dst if (dst.exists() or to_pdbqt(l, dst, ph=ph)) else None

    # Ligandos demasiado flexibles para Vina, apartados ANTES de lanzar nada: además de dar una pose
    # sin sentido, pueden agotar la memoria y tumbar la tanda entera, perdiendo el trabajo ya hecho.
    errors_previos = []
    for k, v in list(lig_pdbqt.items()):
        if v is None:
            continue
        t = torsdof(v)
        if t > TORSDOF_LIMITE:
            lig_pdbqt[k] = None
            errors_previos.append((Path(k).stem,
                                   f"{t} grados de libertad torsionales, por encima del límite "
                                   f"practicable de {TORSDOF_LIMITE} para Vina; usa ADCP"))

    tasks = []
    for rpath, sid, box in targets:
        if rec_pdbqt[str(rpath)] is None:
            continue
        for l in ligands:
            if lig_pdbqt.get(str(l)) is None:
                continue          # sin preparar o apartado por flexibilidad; ya esta reportado
            base = pose_name(sid, Path(l).stem)
            out = poses / f"{base}.pdbqt"
            if out.exists() and list(poses.glob(f"{base}-model*.pdb")):
                continue
            tasks.append((rpath, l, box, base, out))

    def run(t):
        rpath, l, box, base, out = t
        lp = lig_pdbqt.get(str(l))
        if lp is None:
            return (base, "no se pudo preparar el ligando")
        try:
            engine.dock(rec_pdbqt[str(rpath)], lp, box, out)
            split_models(out, poses)
            return (base, None)
        except Exception as e:
            return (base, str(e)[:200])

    n = workers if workers > 0 else paralelismo_seguro(
        [b for _r, _s, b in targets], [v for v in lig_pdbqt.values() if v])
    errors = list(errors_previos)
    if tasks:
        with ThreadPoolExecutor(max_workers=n) as ex:
            for i, (base, err) in enumerate(ex.map(run, tasks), 1):
                if err:
                    errors.append((base, err))
                if on_progress:
                    on_progress(i, len(tasks), base, err)

    rows = []
    for p in sorted(poses.glob("*.pdbqt")):
        m = re.search(r"compounds_a_(.+)$", p.stem)
        if not m:
            continue
        rec = re.search(r"docking_(.+?)_compounds_a_", p.stem)
        for i, s in enumerate(scores_from_pdbqt(p), 1):
            rows.append({"receptor": rec.group(1) if rec else "receptor",
                         "pose_name": f"{p.stem}-model{i}",
                         "compound_name": m.group(1),
                         "docking_score": s})
    return rows, errors
