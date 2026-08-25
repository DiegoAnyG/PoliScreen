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


def write_config(out_dir, box, **overrides) -> Path:
    """The CAVER config, with the search box centre as the starting point.

    That centre is the one thing PoliScreen already knows and CAVER most needs: tunnels are
    measured from the active site outwards, and the box was put there for the docking.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    values = dict(DEFAULTS)
    values.update({k: v for k, v in overrides.items() if v is not None})

    lines = [
        "# Written by PoliScreen. The starting point is the centre of the search box.",
        f"starting_point_coordinates {box.cx} {box.cy} {box.cz}",
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


def find_tunnels(structure, box, out_dir, config=None, heap: str = CAVER_HEAP,
                 timeout: int = 3600) -> Path:
    """Run CAVER on one structure and return the folder holding its output.

    CAVER reads a *folder* of PDB files: one is a static structure, many are the snapshots of a
    trajectory and it clusters the tunnels across them. A single structure is copied into a folder
    of its own so the caller does not have to know that.
    """
    jar = caver_jar()
    if jar is None:
        raise CaverError("CAVER is not installed. It is a Java program and is not shipped: point "
                         "POLISCREEN_CAVER at caver.jar (or the folder holding it).")
    if java_exe() is None:
        raise CaverError("CAVER needs a Java runtime and there is no java on PATH.")

    structure = Path(structure)
    if not structure.exists():
        raise CaverError(f"No such structure: {structure}")

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    structures = out / "structures"
    structures.mkdir(exist_ok=True)
    shutil.copy2(structure, structures / f"{structure.stem}.pdb")

    conf = Path(config) if config else write_config(out, box)
    home = jar.parent
    result = out / "out"
    cmd = [java_exe(), f"-Xmx{heap}", "-cp", str(home / "lib"), "-jar", str(jar),
           "-home", str(home), "-pdb", str(structures), "-conf", str(conf), "-out", str(result)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if not (result / "data" / "clusters").is_dir():
        tail = (r.stderr or r.stdout or "")[-400:]
        raise CaverError(f"CAVER produced no tunnel clusters: {tail}")
    return result


def clusters(caver_out) -> list:
    """The tunnel cluster files CAVER wrote, in priority order -- which is their numbering."""
    folder = Path(caver_out) / "data" / "clusters"
    return sorted(folder.glob("tun_cl_*.pdb")) if folder.is_dir() else []


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
