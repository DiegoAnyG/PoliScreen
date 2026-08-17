"""What actually went into the docking, and what produced it.

Two machines running the same commit reported different rankings, and every explanation offered
for it -- a stale cache, a platform, a version -- was a hypothesis nobody could check, because
"the results differ" says nothing about *where* they start to differ. This turns that into a diff.

Each stage is hashed in the order the pipeline uses it, so the first line that differs is the
stage that caused it:

    prepared receptor   differs -> preparation: openmm, pdbfixer, or the source structure
    docking input       differs -> the conversion to pdbqt, or openbabel
    pose                differs -> the docking engine itself, on identical input

That last case is worth naming, because it is not a bug and no seed can remove it: two builds of
the same version of a docking engine, compiled by different compilers, can sum floating point in a
different order and land on a slightly different pose. If the inputs hash the same and the poses do
not, that is the finding -- and it belongs in the Methods section rather than in a bug tracker.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Iterable, Optional

from . import layout as lay

# In pipeline order. The label is what a reader compares, so it says which stage it belongs to
# rather than which folder it happens to live in.
STAGES = (
    ("source-receptor", (lay.RECEPTORS,), ("*.pdb",)),
    ("ligand-3d", (lay.INPUT_LIGANDS,), ("*.sdf", "*.mol", "*.mol2")),
    ("docking-input", ("prep",), ("*.pdbqt",)),
    ("pose", ("poses",), ("*.pdbqt",)),
)


def sha256(path) -> str:
    """Content hash, with line endings normalised first.

    Every artifact hashed here is a text format, and the Windows build writes CRLF where Linux
    writes LF. Hashing the raw bytes made a byte-for-byte identical prepared receptor look like a
    divergence at the very first stage -- the one place a false positive costs the most, because
    it sends the reader to inspect a preparation step that in fact agreed.
    """
    return hashlib.sha256(Path(path).read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def _folder(proj: Path, name: str) -> Path:
    """The artifact folder, honouring the legacy Spanish names an older project still uses."""
    try:
        return lay.artifact(proj, name)
    except Exception:
        return proj / name


def entries(project, stages: Iterable = STAGES) -> list:
    """[(stage, name, sha256)] for everything the docking consumed or produced, sorted.

    Sorted by stage order and then by name so two machines produce byte-comparable output and an
    ordinary `diff` is enough. The file name is used rather than the full path, because the two
    machines being compared do not share one.
    """
    proj = Path(project)
    out = []
    for stage, folders, patterns in stages:
        found = {}
        for folder in folders:
            d = _folder(proj, folder)
            if not d.is_dir():
                continue
            for pattern in patterns:
                for f in d.glob(pattern):
                    if f.is_file():
                        found[f.name] = f
        for name in sorted(found):
            out.append((stage, name, sha256(found[name])))
    return out


def versions() -> dict:
    """Versions of everything that can move a number, reused from the Methods block."""
    from . import report as rp
    try:
        return rp._versions()
    except Exception as e:                       # never let the diagnostic be the thing that fails
        return {"error": f"{type(e).__name__}: {e}"}


def render(project, stages: Iterable = STAGES) -> str:
    """The report, as text meant to be diffed rather than read."""
    proj = Path(project)
    lines = [f"# poliscreen fingerprint", f"# project: {proj.name}", "", "## versions"]
    for key, value in sorted(versions().items()):
        lines.append(f"{key:<12} {value}")

    rows = entries(proj, stages)
    lines += ["", "## files (stage, name, sha256)"]
    if not rows:
        lines.append("(nothing found: is this a project folder, and has it been run?)")
    for stage, name, digest in rows:
        lines.append(f"{stage:<15} {name:<44} {digest}")

    counts = {}
    for stage, _name, _digest in rows:
        counts[stage] = counts.get(stage, 0) + 1
    lines += ["", "## totals"] + [f"{s:<15} {counts.get(s, 0)}" for s, _f, _p in stages]
    return "\n".join(lines) + "\n"


def protonation_report(project, n: int = 20, seed: int = 11) -> str:
    """Contacts as PLIP finds them today, beside contacts with the protonation pinned down.

    Today PLIP protonates the fused complex itself, with openbabel. The alternative rebuilds the
    receptor's hydrogen network with PDB2PQR, takes the ligand's hydrogens from the SMILES it was
    built from, and tells PLIP to leave them alone.

    Run on Windows and on Linux over the same 20 complexes, both columns came back identical, so
    on the evidence so far neither path is where a platform difference lives. The report stays
    because that is a claim worth being able to re-check on another target rather than remember:
    sorted and free of paths so an ordinary diff answers it.
    """
    import csv
    import random
    import re
    import shutil
    import subprocess
    import tempfile

    from . import interactions as ix
    from . import ligands as lg
    from . import receptor as rc
    from . import screening as sc

    proj = Path(project)
    poses_dir = proj / "poses"
    comps_dir = next((d for d in (_folder(proj, lay.COMPLEXES), proj / "fused_complexes")
                      if d.is_dir()), None)
    if not poses_dir.is_dir() or comps_dir is None:
        return "nothing found: this project has no poses and complexes to compare.\n"

    receptors = (sorted(_folder(proj, lay.RECEPTORS).glob("*_ready.pdb"))
                 or sorted(_folder(proj, lay.RECEPTORS).glob("*_listo.pdb")))
    if not receptors:
        return "nothing found: no prepared receptor in this project.\n"
    if not rc.pdb2pqr_exe():
        return "PDB2PQR is not installed here. pip install pdb2pqr\n"

    smiles = {}
    for meta in ("ligands_meta.csv", "tablas/ligands_meta.csv"):
        p = proj / meta
        if p.exists():
            with open(p, newline="", encoding="utf-8-sig") as f:
                for row in csv.DictReader(f):
                    if row.get("smiles"):
                        smiles[sc.normalize_key(row["name"])] = row["smiles"]

    def profile(complex_pdb, work, nohydro):
        out = work / ("nohydro" if nohydro else "plip") / complex_pdb.stem
        shutil.rmtree(out, ignore_errors=True)
        out.mkdir(parents=True, exist_ok=True)
        cmd = ["plip", "-f", str(complex_pdb), "-x", "--nopdb", "-o", str(out)]
        if nohydro:
            cmd.append("--nohydro")
        subprocess.run(cmd, capture_output=True)
        found = list(out.glob("*.xml"))
        return {k for k, v in ix.parse_plip_xml(found[0]).items() if v} if found else None

    poses = sorted(poses_dir.glob("*-model1.pdb"))
    random.seed(seed)
    sample = sorted(random.sample(poses, min(n, len(poses))), key=lambda p: p.name)

    L = ["# protonation comparison", f"# receptor: {receptors[0].name}",
         f"# sampled: {len(sample)} of {len(poses)} poses, seed {seed}", ""]
    agree = only_current = only_optimised = 0
    with tempfile.TemporaryDirectory() as td:
        work = Path(td)
        optimised = rc.optimize_hydrogens(receptors[0], work / "receptor_optimised.pdb")
        if optimised is None:
            return "PDB2PQR could not process this receptor; nothing to compare.\n"
        for pose in sample:
            key = re.sub(r"-model\d+$", "", pose.stem.split("compounds_a_")[-1])
            smi = smiles.get(sc.normalize_key(key))
            current = next((c for c in (comps_dir / f"{lay.COMPLEX_PREFIX}{pose.stem}.pdb",
                                        comps_dir / f"Complejo_{pose.stem}.pdb") if c.exists()), None)
            if not smi or current is None:
                L.append(f"{pose.stem}\tSKIPPED\tno smiles or no complex")
                continue
            ligand = lg.with_hydrogens(pose, smi, work / f"lig_{pose.stem}.pdb")
            if ligand is None:
                L.append(f"{pose.stem}\tSKIPPED\tthe template does not match the pose")
                continue
            fused = ix.fuse(optimised, ligand, work / f"opt_{pose.stem}.pdb")
            a, b = profile(current, work, False), profile(fused, work, True)
            if a is None or b is None:
                L.append(f"{pose.stem}\tSKIPPED\tPLIP produced no report")
                continue
            agree += len(a & b); only_current += len(a - b); only_optimised += len(b - a)
            L.append(pose.stem)
            L.append("  current  \t" + " ".join(sorted(a)))
            L.append("  optimised\t" + " ".join(sorted(b)))

    L += ["", "# totals", f"# agreeing\t{agree}", f"# only current\t{only_current}",
          f"# only optimised\t{only_optimised}"]
    return "\n".join(L) + "\n"
