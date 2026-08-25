"""Running CAVER and CaverDock, when the machine has them.

Neither ships with PoliScreen, for different reasons, and they are found rather than assumed:

- **CAVER** is GPL-3 and cross-platform, so it *could* be bundled. It is not, for the reason OPSIN
  is not: it is a Java program, and a JVM is about 120 MB per installer and 500 MB unpacked.
  ``POLISCREEN_CAVER`` pointing at ``caver.jar`` (or the folder holding it), plus a Java on
  ``PATH``, turns it on.
- **CaverDock** is distributed as a Linux Apptainer image under an academic licence, so it cannot
  be redistributed at all -- the same position as ADCP. ``POLISCREEN_CAVERDOCK`` points at the
  ``.sif``.

Reading the results needs none of this; see ``core/tunnels.py``. This module is only for computing
them, and it carries two workarounds that each decide whether the answer is right:

**cd-analysis throws the tunnel away.** ``pycaverdock/bin/analysis.py`` extends the tunnel with
``tunnel = tunnel.extend(...)``, and ``extend()`` returns *only the new discs*. The tunnel is
discarded and the ligand is docked through the two-angstrom straight extrapolation past the mouth:
10 discs of constant radius instead of 68, finishing in 29 seconds instead of minutes, with nothing
in the log to say so. ``--skip-tunnel-extension`` uses the whole discretised tunnel, which is also
what CaverWeb docks through. cd-screening writes ``tunnel + tunnel.extend(...)`` and is unaffected.

**More than two MPI processes discard the seed.** CaverDock says so once, mid-run:
``Seed does not guarantee deterministic results when more than 2 MPI proccesses are used!`` -- and
it refuses to run a tunnel with one. Two is therefore both the minimum and the only reproducible
setting, and is the default. The same shape as Vina, which is deterministic on one thread and not
on several.

The command-line equivalent of this module is ``caverdock-run``, in the caver-translate repository.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional

# What CaverWeb runs, so a local run can be compared with one. Every one of these is overridable;
# they are defaults, not decisions made for the user.
DEFAULTS = {
    "probe_radius": 0.9,
    "shell_depth": 4.0,
    "shell_radius": 3.0,
    "clustering_threshold": 3.5,
    "max_distance": 3.0,
    "desired_radius": 5.0,
    # Without this CAVER writes only data/clusters_timeless and no data/clusters, which is the
    # directory everything downstream names. It costs nothing on a single structure and it is what
    # CaverWeb sets, so a local run and a server run leave the same folders behind.
    "save_dynamics_visualization": "yes",
}

# Version 3.0.2 reproduces CaverWeb to every digit on 8HTB; 3.0.3 BETA returns eight tunnels where
# the server returns six, because it swapped the ordinary Voronoi diagram for an additively
# weighted one. Nothing here pins a version -- whatever is installed is used -- but a run that has
# to line up with a CaverWeb job should be 3.0.2.
CAVER_HEAP = "4g"

# Below this CaverDock refuses to run a tunnel at all; above it the seed stops meaning anything.
MPI_REPRODUCIBLE = 2


class CaverError(RuntimeError):
    pass


def _from_env(name: str) -> Optional[Path]:
    value = os.environ.get(name)
    if not value:
        return None
    p = Path(value).expanduser()
    return p if p.exists() else None


def caver_jar() -> Optional[Path]:
    """caver.jar, from the environment, a java on PATH's neighbourhood, or the usual unpack sites.

    POLISCREEN_CAVER may name the jar or the folder holding it: nobody remembers which, and both
    are unambiguous.
    """
    env = _from_env("POLISCREEN_CAVER")
    if env:
        if env.is_file():
            return env
        jar = env / "caver.jar"
        if jar.exists():
            return jar
        found = next(iter(sorted(env.rglob("caver.jar"))), None)
        if found:
            return found

    # A CAVER download unpacks to caver_<version>/caver/caver.jar. Look wherever a person would
    # plausibly have put it, rather than at one path that only exists on one machine.
    for base in (Path.home(), Path.home() / "poliscreen_tools", Path.cwd()):
        try:
            if not base.exists():
                continue
        except OSError:
            continue
        for pattern in ("caver*/caver/caver.jar", "*/caver*/caver/caver.jar", "caver*/caver.jar"):
            found = next(iter(sorted(base.glob(pattern))), None)
            if found:
                return found
    return None


def caverdock_image() -> Optional[Path]:
    """The CaverDock .sif, from the environment or from where an image would plausibly sit."""
    env = _from_env("POLISCREEN_CAVERDOCK")
    if env:
        if env.is_file():
            return env
        found = next(iter(sorted(env.glob("caverdock*.sif"))), None)
        if found:
            return found

    for base in (Path.home(), Path.home() / "poliscreen_tools", Path.cwd()):
        try:
            if not base.exists():
                continue
        except OSError:
            continue
        for pattern in ("caverdock*.sif", "*/caverdock*.sif"):
            found = next(iter(sorted(base.glob(pattern))), None)
            if found:
                return found
    return None


def java_exe() -> Optional[str]:
    return shutil.which("java")


def apptainer_exe() -> Optional[str]:
    return shutil.which("apptainer") or shutil.which("singularity")


def caver_available() -> bool:
    return caver_jar() is not None and java_exe() is not None


def caverdock_available() -> bool:
    return caverdock_image() is not None and apptainer_exe() is not None


def reproducible(cpus: int) -> bool:
    """Whether a run with this many MPI processes can be repeated from its seed."""
    return cpus <= MPI_REPRODUCIBLE


WATER_NAMES = ("HOH", "WAT", "DOD", "H2O")


def hetero_groups(pdb) -> list:
    """The non-water hetero residues in a structure, with how many atoms each has.

    They decide the answer: a cofactor sitting in a channel closes it, and removing one that is
    genuinely there opens a route the protein does not have. Which is why this is asked rather
    than assumed.
    """
    counts = {}
    for line in Path(pdb).read_text(errors="ignore").splitlines():
        if line.startswith("HETATM"):
            name = line[17:20].strip()
            if name and name not in WATER_NAMES:
                counts[name] = counts.get(name, 0) + 1
    return sorted(counts.items())


def prepare_for_caver(pdb, out_path, keep_hetero=(), keep_waters: bool = False) -> Path:
    """The same structure, as CAVER needs to see it.

    Docking and tunnel search want opposite things from a receptor, and feeding one the other's
    file is the quietest way to get a wrong answer:

    **Hydrogens are removed.** CAVER measures a tunnel as the space left between van der Waals
    spheres, and a receptor prepared for docking carries thousands of added hydrogens that make
    every atom effectively larger. Measured on 8HTB: the docking-ready file (4449 atoms, 2237 of
    them hydrogens, GDP and Ca retained) yields **three** tunnels; the bare heavy-atom protein
    yields **six**, and the three that disappear are the narrow ones -- bottleneck radii 0.94, 1.00
    and 1.26 A. Nothing warns you; the run simply reports fewer routes.

    **Heteroatoms are chosen, not inherited.** Waters go unless asked for, and each cofactor or
    ligand is kept only if named. Keeping GDP is a statement that the channel is blocked in the
    physiological state; removing it is a statement that it is not. Both are legitimate; guessing
    is not.
    """
    keep = {str(h).strip().upper() for h in keep_hetero}
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    kept = []
    for line in Path(pdb).read_text(errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            # Column 77-78 is the element; the atom name is the fallback for files that omit it.
            element = line[76:78].strip().upper() or line[12:16].strip().lstrip("0123456789")[:1]
            if element == "H" or element == "D":
                continue
            if line.startswith("HETATM"):
                name = line[17:20].strip().upper()
                if name in WATER_NAMES:
                    if not keep_waters:
                        continue
                elif name not in keep:
                    continue
            kept.append(line)
        elif line.startswith(("TER", "END")):
            continue                      # rewritten below, so a filtered file stays well formed
    out.write_text("\n".join(kept + ["TER", "END", ""]), encoding="utf-8")
    return out


def centroid(points) -> tuple:
    """The middle of a set of atoms, which is how a ligand or a set of residues becomes a point."""
    pts = [tuple(float(v) for v in p[:3]) for p in points]
    if not pts:
        raise CaverError("No atoms to take a starting point from.")
    n = len(pts)
    return tuple(round(sum(p[i] for p in pts) / n, 3) for i in range(3))


def atoms_of(path, residues=()) -> list:
    """(x, y, z) of every atom, or only of the named residues.

    Residues are given as CAVER and PoliScreen both write them, name plus number: ASP199, Leu209.
    """
    wanted = {str(r).strip().upper() for r in residues}
    out = []
    for line in Path(path).read_text(errors="ignore").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if wanted:
            label = f"{line[17:20].strip()}{line[22:26].strip()}".upper()
            if label not in wanted:
                continue
        try:
            out.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        except ValueError:
            continue
    return out


def ligand_atoms(path) -> list:
    """(x, y, z) of a small molecule, in any of the formats a control arrives as.

    Each format is read by its own structure rather than by looking for lines that happen to start
    with three numbers. A bond in an SDF is written ``  2  1  1  0`` -- four numeric fields, which
    a permissive reader takes for a coordinate. On a 33-atom control that turned 33 atoms into 69
    "atoms" and moved the centre from (-14.97, -13.73, 18.64), which is inside the site, to
    (0.10, 2.04, 9.57), which is outside the protein. CAVER started there and reported the outside
    world as one enormous tunnel.
    """
    p = Path(path)
    suffix = p.suffix.lower()
    if suffix in (".pdb", ".pdbqt", ".ent"):
        return atoms_of(p)

    lines = p.read_text(errors="ignore").splitlines()
    if suffix in (".sdf", ".mol"):
        # Line 4 is the counts line: atoms in the first three columns, bonds in the next three.
        if len(lines) < 4:
            return []
        try:
            count = int(lines[3][:3])
        except ValueError:
            return []
        out = []
        for line in lines[4:4 + count]:
            try:
                out.append((float(line[0:10]), float(line[10:20]), float(line[20:30])))
            except ValueError:
                continue
        return out

    if suffix == ".mol2":
        out, inside = [], False
        for line in lines:
            if line.startswith("@<TRIPOS>"):
                inside = line.strip() == "@<TRIPOS>ATOM"
                continue
            if inside:
                fields = line.split()
                if len(fields) >= 5:
                    try:
                        out.append((float(fields[2]), float(fields[3]), float(fields[4])))
                    except ValueError:
                        continue
        return out

    raise CaverError(f"Cannot read coordinates from a {suffix or 'nameless'} file: {p.name}")


def inside_structure(point, structure, margin: float = 6.0) -> bool:
    """Whether a starting point has protein around it, rather than open space.

    CAVER measures outwards from the point, so one placed outside the protein finds the outside:
    a single enormous "tunnel" with a wide bottleneck and no meaning. It says so in a warning
    nobody reads, and reports the result either way.
    """
    x, y, z = start_point(point)
    near = 0
    for ax, ay, az in atoms_of(structure):
        if abs(ax - x) < margin and abs(ay - y) < margin and abs(az - z) < margin:
            near += 1
            if near >= 12:            # a dozen atoms within a few angstroms is a pocket, not space
                return True
    return False


def tunnel_spheres(cluster_pdb) -> list:
    """(x, y, z, radius) of one tunnel, for drawing it.

    A CAVER tunnel is written as a chain of spheres in a PDB, with the radius where the occupancy
    would be. Which is the same shape fpocket's alpha spheres arrive in, so the viewer already
    knows how to draw it.
    """
    out = []
    for line in Path(cluster_pdb).read_text(errors="ignore").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        radius = 0.0
        for start, end in ((54, 60), (60, 66)):
            try:
                value = float(line[start:end])
            except ValueError:
                continue
            if value > 0:
                radius = value
                break
        if radius > 0:
            out.append((x, y, z, radius))
    return out


def start_point(source) -> tuple:
    """A starting point from whatever the caller has: a Box, a point, or atoms to average.

    Tunnels are measured from one point outwards, and where that point sits changes the answer.
    The search box centre is a reasonable default because it was already placed on the site, but
    it is a cube's middle -- a ligand or a set of catalytic residues says it more precisely.
    """
    if source is None:
        raise CaverError("No starting point given.")
    if hasattr(source, "cx"):
        return (source.cx, source.cy, source.cz)
    seq = list(source)
    if len(seq) == 3 and all(isinstance(v, (int, float)) for v in seq):
        return tuple(float(v) for v in seq)
    return centroid(seq)


def write_config(out_dir, start, **overrides) -> Path:
    """The CAVER config, with wherever the tunnels are to be measured from.

    `start` is a Box, a point, or the atoms to take the middle of -- see start_point.
    """
    x, y, z = start_point(start)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    values = dict(DEFAULTS)
    values.update({k: v for k, v in overrides.items() if v is not None})

    lines = [
        "# Written by PoliScreen.",
        f"starting_point_coordinates {x} {y} {z}",
        "",
    ]
    lines += [f"{k} {v}" for k, v in values.items()]
    lines += [
        "visualization_tunnel_sampling_step 0.5",
        "compute_tunnel_residues yes",
        # Seeded like every other stage that touches coordinates: it fixes which structure is
        # sampled, never what it scores.
        "seed 1",
        "generate_trajectory yes",
        "",
    ]
    path = out / "config.txt"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def find_tunnels(structures, start, out_dir, config=None, heap: str = CAVER_HEAP,
                 timeout: int = 3600) -> Path:
    """Run CAVER on one structure, or on several, and return the folder holding its output.

    CAVER reads a *folder* of PDB files. One file is a static structure. Several are the snapshots
    of a trajectory, and CAVER clusters the tunnels across all of them -- which is what a dynamic
    tunnel is, and the only reason to pass more than one. Either is accepted here, so the caller
    does not have to know that a folder is what CAVER wants.

    The structure should have come through prepare_for_caver: a receptor prepared for docking has
    hydrogens on it, and they close the narrow tunnels without saying so.
    """
    jar = caver_jar()
    if jar is None:
        raise CaverError("CAVER is not installed. It is a Java program and is not shipped: point "
                         "POLISCREEN_CAVER at caver.jar (or the folder holding it).")
    if java_exe() is None:
        raise CaverError("CAVER needs a Java runtime and there is no java on PATH.")

    paths = [Path(structures)] if isinstance(structures, (str, Path)) else [Path(p) for p in structures]
    missing = [p for p in paths if not p.exists()]
    if missing:
        raise CaverError(f"No such structure: {missing[0]}")
    if not paths:
        raise CaverError("No structures to search.")

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    folder = out / "structures"
    folder.mkdir(exist_ok=True)
    # Numbered, because CAVER orders the snapshots by file name and a trajectory read out of order
    # clusters tunnels across time steps that never followed one another.
    width = len(str(len(paths)))
    for i, p in enumerate(paths):
        name = f"{p.stem}.pdb" if len(paths) == 1 else f"{i:0{width}d}_{p.stem}.pdb"
        shutil.copy2(p, folder / name)
    structures = folder

    conf = Path(config) if config else write_config(out, start)
    home = jar.parent
    result = out / "out"
    cmd = [java_exe(), f"-Xmx{heap}", "-cp", str(home / "lib"), "-jar", str(jar),
           "-home", str(home), "-pdb", str(structures), "-conf", str(conf), "-out", str(result)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if not clusters(result):
        tail = (r.stderr or r.stdout or "")[-400:]
        raise CaverError(f"CAVER produced no tunnel clusters: {tail}")
    return result


def clusters(caver_out) -> list:
    """The tunnel cluster files CAVER wrote, in priority order -- which is their numbering.

    Two directories can hold them and which one appears depends on a config flag rather than on
    anything about the run: without save_dynamics_visualization CAVER writes only
    clusters_timeless, named tun_cl_003_1.pdb instead of tun_cl_003.pdb. Both are read, because a
    config that came from somewhere else is not ours to assume anything about.
    """
    data = Path(caver_out) / "data"
    for name in ("clusters", "clusters_timeless"):
        folder = data / name
        if folder.is_dir():
            found = sorted(folder.glob("tun_cl_*.pdb"))
            if found:
                return found
    return []


def transport_command(image, workdir, receptor, ligand, tunnel, direction: str, bound: str,
                      cpus: int, seed: int, exhaustiveness: int = 1) -> list:
    """The cd-analysis invocation, kept separate so the two flags that matter can be tested.

    --skip-tunnel-extension is not a preference. Without it the tunnel is replaced by its own
    two-angstrom extension and the ligand is docked through that instead; see the module docstring.
    """
    return [
        apptainer_exe() or "apptainer", "exec", "--bind", f"{workdir}:/work", str(image),
        "cd-analysis",
        "--receptor", f"/work/{Path(receptor).name}",
        "--ligand", f"/work/{Path(ligand).name}",
        "--tunnel", f"/work/{Path(tunnel).name}",
        "--direction", "IN" if direction == "in" else "OUT",
        "--trajectory-type", "UPPERBOUND" if bound == "ub" else "LOWERBOUND",
        "--workdir", "/work",
        "--name", "analysis",
        "--threads", str(cpus),
        "--seed", str(seed),
        "--exhaustiveness", str(exhaustiveness),
        "--skip-tunnel-extension",
    ]


def run_name(receptor, ligand, tunnel, direction: str, bound: str) -> str:
    """cd-screening's own naming, which caver-translate already reads.

    Nothing inside a CaverDock output records what was calculated, so the folder name has to.
    Several runs written into one folder therefore come out as a single table, with the
    combinations not yet run counted as missing.
    """
    long = "lowerbound" if bound == "lb" else "upperbound"
    return (f"r{Path(receptor).stem}-l{Path(ligand).stem}-t{Path(tunnel).stem}"
            f"-d{direction}-{long}")


def transport(receptor, ligand, tunnel, out_dir, direction: str = "out", bound: str = "lb",
              cpus: int = MPI_REPRODUCIBLE, seed: int = 42, exhaustiveness: int = 1,
              timeout: int = 7200) -> Path:
    """Push one ligand down one tunnel and return the folder holding the result.

    The three inputs are copied in beside the output, so a finished run is self-contained: what
    went in sits next to what came out, which is what makes it reproducible a year later.
    """
    image = caverdock_image()
    if image is None:
        raise CaverError("CaverDock is not installed. It is a Linux image under an academic "
                         "licence and is not shipped: point POLISCREEN_CAVERDOCK at the .sif.")
    if apptainer_exe() is None:
        raise CaverError("CaverDock runs from an Apptainer image and there is no apptainer or "
                         "singularity on PATH.")
    if cpus < MPI_REPRODUCIBLE:
        raise CaverError(f"CaverDock needs at least {MPI_REPRODUCIBLE} MPI processes for a tunnel.")

    run = Path(out_dir) / run_name(receptor, ligand, tunnel, direction, bound)
    run.mkdir(parents=True, exist_ok=True)
    for src in (receptor, ligand, tunnel):
        shutil.copy2(src, run / Path(src).name)

    cmd = transport_command(image, run, receptor, ligand, tunnel, direction, bound, cpus, seed,
                            exhaustiveness)
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if not next(iter(run.glob("analysis-lb.pdbqt")), None):
        tail = (r.stderr or r.stdout or "")[-500:]
        raise CaverError(f"CaverDock produced no trajectory: {tail}")
    return run
