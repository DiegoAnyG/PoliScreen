"""AutoDock CrankPep (ADCP): peptide-specific docking.

Vina treats the ligand as an independent torsion tree and its sampling degrades fast with
flexibility: on saFtsZ, a 5-residue peptide (23 rotatable bonds) no longer finishes in two minutes
and a 10-residue one does not converge. ADCP models the peptide with rotamers and samples
conformation and position at once, the only reasonable way to dock something so flexible.

It does not replace Vina for small molecules, and below MIN_RESIDUOS its sampler does not apply
(see that constant). Cost on 8HTB, octapeptide, six threads: 35 s (250,000 steps x 10 replicas) and
208 s (1,000,000 x 20), with the energy improving monotonically (-13.2 -> -17.4 -> -19.8 kcal/mol),
a sign that the sampling converges.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Sequence

SUITE_DIR = Path.home() / "poliscreen_tools" / "adfrsuite"
MAX_RESIDUES = 20
# Sampler lower bound: the crankshaft needs backbone between its pivots. Below 5 residues Vina is practicable.
MIN_RESIDUES = 5


def _suite_root() -> Optional[Path]:
    """Folder of the installed ADFR suite.

    Several installations are allowed and the most recent one with the executable is chosen: trying
    different versions is common (not all work), and looking for a fixed name would force exporting a
    variable every session.
    """
    env = os.environ.get("POLISCREEN_ADCP")
    if env and Path(env).exists():
        p = Path(env)
        return p if p.is_dir() else p.parent.parent

    cand = []
    base = SUITE_DIR.parent
    if base.exists():
        for d in base.glob("adfrsuite*"):
            cand.append(d)
            cand.extend(d.glob("ADFRsuite-*"))
    validos = [d for d in cand if (d / "bin" / "adcp").exists()]
    if validos:
        return max(validos, key=lambda d: d.stat().st_mtime)

    exe = shutil.which("adcp")
    return Path(exe).parent.parent if exe else None


def bin_dir() -> Optional[Path]:
    r = _suite_root()
    b = (r / "bin") if r else None
    return b if b and (b / "adcp").exists() else None


def available() -> bool:
    return bin_dir() is not None


def _env() -> dict:
    """PATH and LD_LIBRARY_PATH that ADCP needs.

    The suite binaries require libgomp (OpenMP), which is not always on the system but is in the
    conda environment; its lib is added so as not to depend on a privileged installation.
    """
    env = dict(os.environ)
    b = bin_dir()
    if b:
        env["PATH"] = str(b) + os.pathsep + env.get("PATH", "")
    libs = [p for p in (Path(os.environ.get("CONDA_PREFIX", "")) / "lib",) if p.exists()]
    if libs:
        env["LD_LIBRARY_PATH"] = os.pathsep.join(
            [str(p) for p in libs] + ([env["LD_LIBRARY_PATH"]] if env.get("LD_LIBRARY_PATH") else []))
    return env


def _rama_file() -> Optional[Path]:
    """Ramachandran probability table that the sampling engine expects in its working directory.
    It ships with the suite, but the binary does not look for it in its own installation."""
    r = _suite_root()
    if not r:
        return None
    cand = r / "CCSBpckgs" / "ADCP" / "ramaprob.data"
    if cand.exists():
        return cand
    return next(iter(r.rglob("ramaprob.data")), None)


class AdcpError(RuntimeError):
    pass


def prepare_target(receptor_pdb, box, out_dir, name_: str = "target", timeout: int = 900) -> Optional[Path]:
    """Prepares the target with AGFR and returns the .trg file that ADCP consumes.

    ADCP does not accept a PDB directly: it needs the affinity maps precomputed in a box, which is
    what AGFR generates. The same search box as the rest of PoliScreen is reused, so both engines
    explore exactly the same region.
    """
    b = bin_dir()
    if b is None:
        raise AdcpError("ADCP is not installed. Run scripts/get_adcp.sh.")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    env = _env()

    recq = out_dir / f"{name_}.pdbqt"
    if not recq.exists():
        r = subprocess.run([str(b / "prepare_receptor"), "-r", str(receptor_pdb), "-o", str(recq)],
                           capture_output=True, text=True, env=env, timeout=timeout)
        if not recq.exists():
            raise AdcpError(f"prepare_receptor failed: {(r.stderr or r.stdout or '')[-300:]}")

    dest_ = out_dir / name_
    cmd = [str(b / "agfr"), "-r", str(recq), "-o", str(dest_),
           "-b", f"user {box.cx} {box.cy} {box.cz} {box.sx} {box.sy} {box.sz}"]
    r = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=timeout)
    trg = next(iter(sorted(out_dir.glob(f"{name_}*.trg"))), None)
    if trg is None:
        raise AdcpError(f"AGFR did not generate the target: {(r.stderr or r.stdout or '')[-300:]}")
    return trg


_RE_ENERGY = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s", re.M)
_RE_RANK = re.compile(r"_ranked_(\d+)\.pdb$")


def dock_peptide(sequence_: str, target_trg, out_dir, name_: Optional[str] = None,
                 n_steps: int = 250000, n_replicas: int = 10, seed_: int = 42,
                 n_cores: int = 1, helix: bool = False, cyclic: bool = False,
                 cystine: bool = False, timeout: int = 3600) -> dict:
    """Docks a peptide from its SEQUENCE. Returns {'energy', 'poses', 'output'}.

    ADCP builds the conformation during docking, so no prior 3D structure is needed. In the sequence,
    an UPPERCASE letter means that residue starts from a helical conformation and lowercase from coil;
    useful when the peptide is known to be helical.

    ciclico: head-to-tail closure through the backbone. The sampler itself applies it, restricting the
    conformation to the cycle instead of docking the linear peptide; passing it matters because a
    cycle has far fewer degrees of freedom and a different binding mode from its open analogue.
    cistina: closure by a disulfide bridge between two cysteines.
    """
    b = bin_dir()
    if b is None:
        raise AdcpError("ADCP is not installed.")
    seq = "".join(ch for ch in sequence_ if ch.isalpha())
    if not MIN_RESIDUES <= len(seq) <= MAX_RESIDUES:
        raise AdcpError(f"ADCP docks from {MIN_RESIDUES} to {MAX_RESIDUES} residues; "
                        f"{len(seq)} were requested. Below {MIN_RESIDUES} use Vina.")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    name_ = name_ or seq

    trg = Path(target_trg)
    local_trg = out_dir / trg.name
    if not local_trg.exists():
        shutil.copy(trg, local_trg)

    rama = _rama_file()
    if rama and not (out_dir / rama.name).exists():
        shutil.copy(rama, out_dir / rama.name)

    cmd = [str(b / "adcp"), "-t", trg.name,
           "-s", seq.upper() if helix else seq.lower(),
           "-N", str(n_replicas), "-n", str(n_steps),
           "-o", name_, "-c", str(n_cores), "-S", str(seed_)]
    if cyclic:
        cmd.append("-cyc")
    if cystine:
        cmd.append("-cys")
    r = subprocess.run(cmd, capture_output=True, text=True, env=_env(),
                       cwd=str(out_dir), timeout=timeout)
    output = (r.stdout or "") + (r.stderr or "")

    poses = sorted((p for p in out_dir.glob(f"{name_}_ranked_*.pdb") if _tiene_atomos(p)),
                   key=lambda p: int(_RE_RANK.search(p.name).group(1)))
    energias = [float(m.group(2)) for m in _RE_ENERGY.finditer(output)]
    energias.sort()
    return {"sequence": seq, "energy": energias[0] if energias else None,
            "energies": energias[:len(poses)] or None,
            "poses": [str(p) for p in poses], "output": output[-1500:],
            "ok": bool(poses)}


def _tiene_atomos(pdb: Path) -> bool:
    try:
        return any(l.startswith(("ATOM", "HETATM")) for l in pdb.read_text(errors="ignore").splitlines())
    except Exception:
        return False


def diagnostic(receptor_pdb, box, timeout: int = 900) -> tuple:
    """(works, message). Checks end to end that ADCP produces usable poses.

    It is tested with a pentapeptide, not a shorter one: below MIN_RESIDUOS the sampler always aborts,
    so a shorter test would diagnose a correct installation as broken.
    """
    import tempfile
    if not available():
        return False, "ADCP is not installed. Run scripts/get_adcp.sh."
    tmp = Path(tempfile.mkdtemp())
    try:
        trg = prepare_target(receptor_pdb, box, tmp, name_="diag")
    except Exception as e:
        return False, f"AGFR could not prepare the target: {str(e)[:160]}"
    try:
        r = dock_peptide("KWKLF", trg, tmp, name_="diag_pep", n_steps=50000, n_replicas=2,
                         timeout=timeout)
    except Exception as e:
        return False, f"ADCP could not dock: {str(e)[:160]}"
    if r["ok"]:
        return True, (f"ADCP operational: {len(r['poses'])} poses, best energy "
                      f"{r['energy']:.1f} kcal/mol.")
    if "too short chains" in r["output"]:
        return False, (f"The sampler rejected the test sequence as too short; ADCP needs at "
                       f"least {MIN_RESIDUES} residues.")
    return False, "ADCP produced no poses with coordinates."


def dock_batch(sequences: Sequence[str], target_trg, out_dir, on_progress=None, **kw) -> list:
    """Docks several sequences against the same target. Sequential: ADCP already parallelizes internally."""
    rows_ = []
    for i, seq in enumerate(sequences, 1):
        try:
            rows_.append(dock_peptide(seq, target_trg, out_dir, **kw))
        except Exception as e:
            rows_.append({"sequence": seq, "energy": None, "poses": [],
                          "output": str(e)[:300], "ok": False})
        if on_progress:
            on_progress(i, len(sequences), seq)
    return rows_


MAX_CLOSURE = 1.8


def _close_ring(pdb) -> Optional[float]:
    """C(last)-N(first) distance; if it is a bond length, adds a CONECT.

    ADCP does not write the closing bond, so the viewer —which infers bonds by distance— draws the
    peptide open even when the ring is formed; the CONECT shows it as the cycle it is. It is only
    added if the distance justifies it: forcing it between distant atoms would draw a long stick and
    mislead about the real geometry.
    """
    pdb = Path(pdb)
    n_serial = c_serial = None
    n_xyz = c_xyz = None
    n_res = []
    try:
        lines_ = pdb.read_text(errors="ignore").splitlines()
    except Exception:
        return None
    for l in lines_:
        if l.startswith("ATOM"):
            n_res.append(int(l[22:26]))
    if not n_res:
        return None
    a, b = min(n_res), max(n_res)
    for l in lines_:
        if not l.startswith("ATOM"):
            continue
        name_, resi = l[12:16].strip(), int(l[22:26])
        try:
            serial = int(l[6:11])
        except ValueError:
            continue
        xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        if resi == a and name_ == "N":
            n_serial, n_xyz = serial, xyz
        elif resi == b and name_ == "C":
            c_serial, c_xyz = serial, xyz
    if not (n_xyz and c_xyz):
        return None
    d = sum((x - y) ** 2 for x, y in zip(n_xyz, c_xyz)) ** 0.5
    if d <= MAX_CLOSURE and n_serial and c_serial:
        body_ = [l for l in lines_ if not l.strip().startswith(("END", "CONECT"))]
        body_.append(f"CONECT{c_serial:>5}{n_serial:>5}")
        body_.append(f"CONECT{n_serial:>5}{c_serial:>5}")
        body_.append("END")
        pdb.write_text("\n".join(body_) + "\n")
    return d


def dock_sites(targets, peptides: dict, out_dir, receptor_by_site: dict,
                n_poses: Optional[int] = None, on_progress=None, **kw) -> tuple:
    """Docks peptides with ADCP; returns (rows, errors) in the shape of docking.dock_batch.

    It names the poses with the rest of the application's convention so that complex fusion, PLIP and
    scoring treat them like any engine's. ADCP's energy is not comparable with Vina's, so each row
    declares its engine.

    targets: [(receptor_path, site_id, box)]. peptidos: {name: (sequence, cyclic)}.
    """
    out_dir = Path(out_dir)
    poses_dir = out_dir / "poses"
    poses_dir.mkdir(parents=True, exist_ok=True)
    work_dir_ = out_dir / "adcp"
    work_dir_.mkdir(parents=True, exist_ok=True)

    rows_, errors_ = [], []
    total = len(targets) * max(1, len(peptides))
    done_ = 0
    for rpath, sid, box in targets:
        try:
            trg = prepare_target(receptor_by_site.get(sid, rpath), box, work_dir_ / sid, name_=sid)
        except Exception as e:
            errors_.append((sid, f"AGFR could not prepare the target: {str(e)[:160]}"))
            continue
        for name_, (seq, cyclic) in peptides.items():
            done_ += 1
            base = f"docking_{sid}_compounds_a_{name_}"
            try:
                r = dock_peptide(seq, trg, work_dir_ / sid / name_, name_=name_,
                                 cyclic=cyclic, **kw)
            except Exception as e:
                errors_.append((base, str(e)[:200]))
                if on_progress:
                    on_progress(done_, total, base, str(e)[:200])
                continue
            if not r["ok"]:
                errors_.append((base, "ADCP produced no poses"))
                if on_progress:
                    on_progress(done_, total, base, "no poses")
                continue

            # -N are search replicas, not poses: trimming them would worsen the result.
            output_poses = r["poses"][:n_poses] if n_poses else r["poses"]
            energias = (r.get("energies") or [r["energy"]] * len(r["poses"]))[:len(output_poses)]
            notice = None
            abiertas = 0
            for i, (p, e) in enumerate(zip(output_poses, energias), 1):
                dest_ = poses_dir / f"{base}-model{i}.pdb"
                dest_.write_bytes(Path(p).read_bytes())
                if cyclic:
                    # The ring is closed only if the geometry justifies it: ADCP does not always close small cycles.
                    d = _close_ring(dest_)
                    if d is None or d > MAX_CLOSURE:
                        abiertas += 1
                    if i == 1 and d is not None and d > MAX_CLOSURE:
                        notice = (f"cyclization was requested and the pose comes out with the termini {d:.2f} A apart "
                                 f"(an amide bond is 1.33 A): treat the geometry with caution")
                rows_.append({"receptor": sid, "pose_name": f"{base}-model{i}",
                              "compound_name": name_, "docking_score": e, "engine": "adcp"})
            if cyclic and abiertas == len(output_poses):
                notice = (f"ADCP did not close the ring in any of the {abiertas} poses; cycles of "
                         f"few residues stay strained and their termini do not reach each other. "
                         f"The structure used for descriptors is cyclic; the docked pose is not.")
            if notice:
                errors_.append((base, notice))
            if on_progress:
                on_progress(done_, total, base, None)
    return rows_, errors_
