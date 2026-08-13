"""Docking: file preparation, search box and engines.

Reproducible: fixed seed and one thread per docking (Vina is not deterministic with multithreading).
Parallelism goes between independent dockings, not inside. No hard-wired targets or paths.
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
    """Co-crystallized ligands (non-water HETATM) -> label: Box centered on them.

    Used to center the box on the correct site when a protein has several.
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
    """Coordinates from any common ligand or structure format."""
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
        inside = False
        for l in text.splitlines():
            if l.startswith("@<TRIPOS>ATOM"):
                inside = True
                continue
            if l.startswith("@<TRIPOS>"):
                inside = False
            if inside:
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
    """Residues (label 'Tyr157', as PLIP names them) with at least one atom inside the box.
    Defines the objective 'pocket': the ranking rewards productive contacts with ANY of these,
    not only with the ones the control touches."""
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
    """Box centered on a specific ligand: the most reliable when the co-crystallized one is available."""
    return Box.around(coords_from_file(path), pad=pad, lo=lo, hi=hi)


def load_boxes_xlsx(path, receptors: Sequence) -> dict:
    """Loads boxes from an xlsx with columns receptor,cx,cy,cz,sx,sy,sz (MolModa format).

    Matches by file name without extension; ignores incomplete rows.
    """
    import pandas as pd
    t = pd.read_excel(path)
    cols = {str(c).lower().strip(): c for c in t.columns}
    need = ("receptor", "cx", "cy", "cz", "sx", "sy", "sz")
    if any(k not in cols for k in need):
        raise DockingError(f"The xlsx is missing columns. Required: {', '.join(need)}.")
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
    """Automatic box: centered on the co-crystallized ligand if there is one, otherwise on the protein."""
    het = _coords(pdb, het_only=True)
    if len(het) >= 3:
        return Box.around(het)
    pts = _coords(pdb)
    if not pts:
        raise DockingError(f"{Path(pdb).name}: no readable coordinates.")
    return Box.around(pts, pad=0.0, lo=24.0, hi=24.0)


# Bumped whenever to_pdbqt changes how the file is made. The conversion was cached by existence
# alone, so a pdbqt written before a fix was reused for ever: on a project folder carried over from
# an earlier version the receptor reaching Vina was still the randomly re-protonated one, and the
# ranking moved with it while the code on disk was correct. Nothing said so, because from the
# outside a cached file and a fresh one are the same file.
#   1  the original conversion, which let obabel re-protonate the receptor at random
#   2  hydrogens kept as prepared -- see to_pdbqt
PDBQT_RECIPE = 2
_STAMP = "REMARK POLISCREEN PDBQT RECIPE "


def pdbqt_is_current(path) -> bool:
    """True when this pdbqt was written by the conversion in use now.

    An unstamped file predates the stamp, and therefore comes from the conversion that randomised
    the hydrogens, so it is rebuilt. Rebuilding costs about a second for a receptor; being wrong
    costs a ranking nobody can reproduce.
    """
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return False
    # The stamp is appended, so it is at the end; the whole file is already in memory from the read
    # and a receptor is a few thousand lines, so there is nothing to gain by stopping early.
    for line in reversed(path.read_text(errors="ignore").splitlines()):
        if line.startswith(_STAMP):
            return line[len(_STAMP):].strip() == str(PDBQT_RECIPE)
    return False


def _write_stamp(path) -> None:
    """A REMARK, which every pdbqt reader ignores and Vina has always skipped."""
    with Path(path).open("a", encoding="utf-8") as fh:
        fh.write(f"{_STAMP}{PDBQT_RECIPE}\n")


def has_hydrogens(path) -> bool:
    """True if the structure already carries explicit hydrogens."""
    for line in Path(path).read_text(errors="ignore").splitlines():
        if line[:6] in ("ATOM  ", "HETATM"):
            element = line[76:78].strip() or line[12:16].strip()
            if element[:1] == "H":
                return True
    return False


def to_pdbqt(src, dst, receptor: bool = False, ph: float = 7.4) -> bool:
    """Converts to pdbqt with OpenBabel. Receptor: rigid; ligand: Gasteiger charges."""
    src, dst = Path(src), Path(dst)
    if src.suffix.lower() == ".pdbqt":
        shutil.copy(src, dst)
        if dst.exists():
            _write_stamp(dst)
        return dst.exists()
    dst.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["obabel", str(src), "-O", str(dst)]
    # Keeps the largest fragment: a ligand extracted from a PDB can come out fragmented and Vina rejects it.
    if receptor:
        cmd += ["-xr"]
        # obabel -p STRIPS every hydrogen and re-adds them at RANDOM positions, so the same prepared
        # receptor produced a different pdbqt on every run: about 25 polar hydrogens moved, some by
        # 1.8 A. Those are exactly the atoms Vina types HD and PLIP reads its hydrogen bonds off, so
        # the poses and the interaction fingerprint moved with them. It is the second half of the
        # fault seeded out of receptor.prepare(), and the reason two machines holding a
        # byte-identical prepared receptor still disagreed on the ranking. What comes out of
        # prepare() is already protonated at this pH by PDBFixer, with a better pKa model than
        # obabel's and a fixed seed, so redoing it here could only lose that. Only a receptor that
        # arrives with no hydrogens at all still needs obabel to put them in.
        if not has_hydrogens(src):
            cmd += ["-p", str(ph)]
    else:
        cmd += ["-r", "-p", str(ph), "--partialcharge", "gasteiger"]
    if src.suffix.lower() == ".smi":
        cmd += ["--gen3d"]
    subprocess.run(cmd, capture_output=True, text=True)
    made = dst.exists() and dst.stat().st_size > 0
    if made:
        _write_stamp(dst)
    return made


def scores_from_pdbqt(path) -> list:
    """Pose energies, read from the REMARK lines of the pdbqt itself."""
    out = []
    for line in Path(path).read_text(errors="ignore").splitlines():
        if line.startswith("REMARK VINA RESULT"):
            try:
                out.append(float(line.split(":")[1].split()[0]))
            except Exception:
                pass
    return out


def split_models(pose_pdbqt, out_dir) -> list:
    """Splits the multi-pose pdbqt into one PDB per model."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = Path(pose_pdbqt).stem
    existing = sorted(out_dir.glob(f"{stem}-model*.pdb"))
    if existing:
        return existing
    subprocess.run(["obabel", str(pose_pdbqt), "-opdb", "-O", str(out_dir / f"{stem}-model.pdb"), "-m"],
                   capture_output=True)
    return sorted(out_dir.glob(f"{stem}-model*.pdb"))


class VinaEngine:
    """AutoDock Vina on CPU. Deterministic with a fixed seed and one thread."""

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
            raise DockingError(f"vina produced no output. Probable cause: box off the site or poorly "
                               f"prepared receptor. stderr: {(r.stderr or '')[-300:]}")
        return scores_from_pdbqt(out_pdbqt)


class GninaEngine(VinaEngine):
    """gnina: a Vina-derived engine whose scoring is given by a convolutional neural network trained on
    crystallographic complexes, with GPU acceleration.

    Its interest in PoliScreen is not speed but being a scoring function INDEPENDENT of Vina's: two
    different methods that agree are more credible than one alone.
    """

    name = "gnina"

    def __init__(self, exe: Optional[str] = None, cpu: int = 1, seed: int = 42,
                 exhaustiveness: int = 24, n_poses: int = 10, energy_range: float = 3.0,
                 use_gpu: bool = True):
        super().__init__(exe=exe or gnina_exe(), cpu=cpu, seed=seed,
                         exhaustiveness=exhaustiveness, n_poses=n_poses, energy_range=energy_range)
        self.use_gpu = use_gpu

    def dock(self, receptor_pdbqt, ligand_pdbqt, box: Box, out_pdbqt) -> list:
        if not self.available():
            raise DockingError("No encuentro gnina. Define POLISCREEN_GNINA o ponlo en el PATH.")
        Path(out_pdbqt).parent.mkdir(parents=True, exist_ok=True)
        cmd = ([self.exe, "-r", str(receptor_pdbqt), "-l", str(ligand_pdbqt)]
               + box.args()
               + ["--exhaustiveness", str(self.exhaustiveness), "--num_modes", str(self.n_poses),
                  "--seed", str(self.seed), "--cpu", str(self.cpu), "-o", str(out_pdbqt)])
        if not self.use_gpu:
            cmd.append("--no_gpu")
        r = subprocess.run(cmd, capture_output=True, text=True)
        if not Path(out_pdbqt).exists():
            raise DockingError(f"gnina produced no output. stderr: {(r.stderr or '')[-300:]}")
        return scores_from_pdbqt(out_pdbqt)


def gnina_exe() -> Optional[str]:
    """Path to the gnina executable: environment variable, local wrapper, PATH or standalone binary.

    The wrapper goes before the binary: the official executable is not self-contained (it needs cuDNN
    and CUDA libraries in LD_LIBRARY_PATH) and aborts at startup without them; the wrapper sets them
    up. Prioritizing the binary reported it as available and failed only when used.
    """
    env = os.environ.get("POLISCREEN_GNINA")
    if env and Path(env).exists():
        return env
    base = Path.home() / "poliscreen_tools"
    for name_ in ("gnina-run", "gnina"):
        local = base / name_
        if local.exists() and os.access(local, os.X_OK):
            return str(local)
    return shutil.which("gnina")


def gnina_available() -> bool:
    return gnina_exe() is not None


_RE_SCORE = re.compile(r"^\s*(Affinity|CNNscore|CNNaffinity|CNN_VS)\s*:?\s*([-\d.]+)", re.M)


def rescore_poses(receptor, poses: Sequence, use_gpu: bool = True, timeout: int = 900) -> dict:
    """Re-scores with gnina poses already generated by Vina, without searching again.

    A cheap, reproducible second opinion: the sampling is still Vina's; only the function that
    evaluates each pose changes. Returns {pose_name: {cnn_score, cnn_affinity, affinity}}, with
    cnn_score as the probability of a correct pose (0-1) and cnn_affinity in pK units.
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
        if not use_gpu:
            cmd.append("--no_gpu")
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            continue
        vals = {}
        for key_, value_ in _RE_SCORE.findall(r.stdout or ""):
            try:
                vals.setdefault(key_, float(value_))
            except ValueError:
                pass
        if vals:
            out[p.stem] = {"affinity": vals.get("Affinity"),
                           "cnn_score": vals.get("CNNscore"),
                           "cnn_affinity": vals.get("CNNaffinity")}
    return out


ENGINES = {"vina": VinaEngine, "gnina": GninaEngine}


def max_extent(archivos: Sequence) -> float:
    """Longest axis of the largest ligand, in angstrom. 0 if it cannot be measured.

    Sizes the box: the ligand must not only fit, it must be able to reorient. A box barely larger
    limits the search to the orientations that fit, and the poor result looks like docking when it is
    an unintended geometric restriction.
    """
    largest = 0.0
    for f in archivos:
        try:
            pts = coords_from_file(f)
        except Exception:
            continue
        if not pts:
            continue
        axes_ = [max(p[i] for p in pts) - min(p[i] for p in pts) for i in range(3)]
        largest = max(largest, max(axes_))
    return largest


def min_box(archivos: Sequence, margin: float = 4.0) -> float:
    """Minimum recommended side so the largest of those ligands can rotate inside the box."""
    ext = max_extent(archivos)
    return round(ext + margin, 1) if ext else 0.0


def available_memory_gb() -> Optional[float]:
    """GB actually available, read from /proc/meminfo. None if it cannot be known."""
    try:
        for l in Path("/proc/meminfo").read_text().splitlines():
            if l.startswith("MemAvailable:"):
                return int(l.split()[1]) / 1048576.0
    except Exception:
        pass
    return None


TORSDOF_LIMIT = 15


def torsdof(pdbqt) -> int:
    """Torsional degrees of freedom Vina declares for an already prepared ligand."""
    try:
        for l in Path(pdbqt).read_text(errors="ignore").splitlines():
            if l.startswith("TORSDOF"):
                return int(l.split()[1])
    except Exception:
        pass
    return 0


def memory_cost_gb(box: "Box", tors: int = 6) -> float:
    """Estimated memory of a docking, in GB.

    Two terms: the affinity maps, which grow with the box volume (0.375 A grid per atom type), and the
    conformation tree, which dominates as soon as the ligand stops being small and is modeled quadratic
    in the degrees of freedom, because the states alive during the search grow much faster than the
    torsions.

    A generous bound on purpose: overshooting costs time; undershooting makes the system kill the
    process mid-run and all the work is lost.
    """
    points_ = ((box.sx / 0.375) + 1) * ((box.sy / 0.375) + 1) * ((box.sz / 0.375) + 1)
    maps = points_ * 22 * 4 / 1e9 * 3.0
    conformaciones = 0.05 * max(1, tors) ** 2 / 36.0
    return max(0.30, maps + conformaciones)


def safe_parallelism(boxes_: Sequence["Box"], ligand_pdbqts: Sequence,
                       reserve_gb: float = 2.0) -> int:
    """How many dockings to launch at once without exhausting memory.

    Splitting only by cores is what makes the system kill the process: many threads with little RAM is
    common (WSL takes half the machine by default, and 16 cores coexist with 7 GB). The more
    restrictive limit between cores and available memory is used, with margin for the interpreter.

    Half the cores, not all of them: Vina saturates whatever it is given, and on a laptop a
    sustained all-core run throttles, which ends up slower than fewer jobs at full clock. Raise it
    with `workers` on a machine with cooling to spare.
    """
    by_cores = max(1, (os.cpu_count() or 2) // 2)
    libre = available_memory_gb()
    if not libre:
        return by_cores
    boxes_ = list(boxes_) or [Box(0, 0, 0, 24, 24, 24)]
    largest = max(boxes_, key=lambda b: b.sx * b.sy * b.sz)
    tors = max([torsdof(l) for l in ligand_pdbqts if l] or [6])
    coste = memory_cost_gb(largest, tors)
    by_memory = max(1, int((libre - reserve_gb) / coste))
    return max(1, min(by_cores, by_memory))


def pose_name(receptor_stem: str, ligand_stem: str) -> str:
    return f"docking_{receptor_stem}_compounds_a_{ligand_stem.lower()}"


def dock_batch(receptors: Sequence, ligands: Sequence, boxes: dict, work_dir,
               engine=None, workers: int = 0, on_progress=None, ph: float = 7.4, targets=None) -> list:
    """Docks each ligand against each site. Resumable: skips what already exists on disk.

    A site is (receptor_path, site_id, Box). By default, one per receptor derived from `boxes`. With
    `targets` the same ligand is docked in several pockets of the same receptor (hybrid); the site_id
    tells them apart and appears as 'receptor' in the results, separating the ranking per site. ph
    protonates receptor and ligand when converting to pdbqt. Returns one row per pose.
    """
    engine = engine or VinaEngine()
    work = Path(work_dir)
    prep, poses = work / "prep", work / "poses"
    for d in (prep, poses):
        d.mkdir(parents=True, exist_ok=True)

    if targets is None:
        targets = [(r, Path(r).stem, boxes[str(r)]) for r in receptors if str(r) in boxes]

    # Inputs rebuilt on this run. Their poses were docked against a different file and cannot be
    # reused: invalidating the pdbqt without invalidating what was docked with it would move the
    # staleness one level down and leave it just as invisible.
    rebuilt = set()

    rec_pdbqt = {}
    for rpath, _sid, _box in targets:
        if str(rpath) in rec_pdbqt:
            continue
        dst = prep / f"{Path(rpath).stem}.pdbqt"
        # pdbqt_is_current, not exists: a file from an older conversion has to be rebuilt, or a
        # project carried over from a previous version keeps docking against the receptor that
        # version produced, silently, for ever.
        if pdbqt_is_current(dst):
            rec_pdbqt[str(rpath)] = dst
        else:
            rec_pdbqt[str(rpath)] = dst if to_pdbqt(rpath, dst, receptor=True, ph=ph) else None
            rebuilt.add(str(rpath))

    lig_pdbqt = {}
    for l in ligands:
        dst = prep / f"{Path(l).stem}.pdbqt"
        if pdbqt_is_current(dst):
            lig_pdbqt[str(l)] = dst
        else:
            lig_pdbqt[str(l)] = dst if to_pdbqt(l, dst, ph=ph) else None
            rebuilt.add(str(l))

    # Ligands too flexible for Vina are set aside: their pose would be meaningless and they can exhaust memory.
    previous_errors = []
    for k, v in list(lig_pdbqt.items()):
        if v is None:
            continue
        t = torsdof(v)
        if t > TORSDOF_LIMIT:
            lig_pdbqt[k] = None
            previous_errors.append((Path(k).stem,
                                   f"{t} torsional degrees of freedom, above the practicable "
                                   f"limit of {TORSDOF_LIMIT} for Vina; use ADCP"))

    tasks = []
    for rpath, sid, box in targets:
        if rec_pdbqt[str(rpath)] is None:
            continue
        for l in ligands:
            if lig_pdbqt.get(str(l)) is None:
                continue
            base = pose_name(sid, Path(l).stem)
            out = poses / f"{base}.pdbqt"
            stale = str(rpath) in rebuilt or str(l) in rebuilt
            if not stale and out.exists() and list(poses.glob(f"{base}-model*.pdb")):
                continue
            if stale:
                # The split models are regenerated from the pose file, and split_models returns
                # early when they are already there, so they have to go with it.
                for old in [out] + list(poses.glob(f"{base}-model*.pdb")):
                    old.unlink(missing_ok=True)
            tasks.append((rpath, l, box, base, out))

    def run(t):
        rpath, l, box, base, out = t
        lp = lig_pdbqt.get(str(l))
        if lp is None:
            return (base, "the ligand could not be prepared")
        try:
            engine.dock(rec_pdbqt[str(rpath)], lp, box, out)
            split_models(out, poses)
            return (base, None)
        except Exception as e:
            return (base, str(e)[:200])

    n = workers if workers > 0 else safe_parallelism(
        [b for _r, _s, b in targets], [v for v in lig_pdbqt.values() if v])
    errors = list(previous_errors)
    if tasks:
        if on_progress:
            # done=0 announces the plan: how many jobs and how many at a time, so what the engine
            # actually does is visible instead of guessed at.
            on_progress(0, len(tasks), f"workers={n}", None)
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
