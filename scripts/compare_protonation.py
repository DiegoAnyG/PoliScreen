"""Measure what deterministic protonation costs, and whether it buys portability.

Two machines running the same commit disagreed on 18.6% of the detected contacts, from files that
hashed the same. Every one of the disagreements was a hydrogen bond or a salt bridge and no
hydrophobic contact moved: PLIP protonates the complex with openbabel before profiling it, and
openbabel is compiled per platform.

This compares, complex by complex:

  current      the fused complex as it is, profiled by PLIP, which protonates it itself
  optimised    the receptor's hydrogen network rebuilt by PDB2PQR/PROPKA and the ligand's taken
               from its recorded SMILES, then PLIP told not to protonate anything (--nohydro)

PDB2PQR is noarch -- the same Python file on every platform -- so the second column should agree
between machines where the first does not. That is the whole question, and this script is how it
gets answered rather than argued.

Run it on both machines with the same project and diff the two reports:

    python scripts/compare_protonation.py <project> --out linux.txt

On Windows, inside the installed environment (PDB2PQR is a pure-Python wheel, nothing compiles):

    C:\PoliScreen\python.exe -m pip install pdb2pqr
    C:\PoliScreen\python.exe scripts\compare_protonation.py <project> --out windows.txt

The report is sorted and carries no paths, so an ordinary diff is enough.
"""
from __future__ import annotations

import argparse
import csv
import random
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from poliscreen.core import interactions as ix   # noqa: E402
from poliscreen.core import layout as lay        # noqa: E402
from poliscreen.core import ligands as lg        # noqa: E402
from poliscreen.core import receptor as rc       # noqa: E402
from poliscreen.core import screening as sc      # noqa: E402


def _smiles_by_name(project: Path) -> dict:
    """The SMILES each ligand was built from, which is what gives its pose hydrogens back."""
    out = {}
    for name in ("ligands_meta.csv", "tablas/ligands_meta.csv"):
        p = project / name
        if not p.exists():
            continue
        with open(p, newline="", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                if row.get("smiles"):
                    out[sc.normalize_key(row["name"])] = row["smiles"]
    return out


def _profile(complex_pdb: Path, work: Path, nohydro: bool):
    """The contacts PLIP reports, through the same parser the pipeline uses."""
    out = work / ("nohydro" if nohydro else "plip") / complex_pdb.stem
    shutil.rmtree(out, ignore_errors=True)
    out.mkdir(parents=True, exist_ok=True)
    cmd = ["plip", "-f", str(complex_pdb), "-x", "--nopdb", "-o", str(out)]
    if nohydro:
        cmd.append("--nohydro")
    subprocess.run(cmd, capture_output=True)
    found = list(out.glob("*.xml"))
    return {k for k, v in ix.parse_plip_xml(found[0]).items() if v} if found else None


def _dirs(project: Path):
    """Where the poses and complexes live, in this layout or an exported one."""
    for poses, comps in ((project / "poses", project / lay.COMPLEXES),
                         (project / "poses", project / "fused_complexes"),
                         (project / "poses", project / lay.LEGACY[lay.COMPLEXES])):
        if poses.is_dir() and comps.is_dir():
            return poses, comps
    raise SystemExit(f"{project} has no poses/ and complexes/ to compare.")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("project")
    ap.add_argument("--out", help="write the report here instead of to the screen")
    ap.add_argument("--n", type=int, default=20, help="complexes to sample (default 20)")
    ap.add_argument("--seed", type=int, default=11, help="which ones, held fixed across machines")
    args = ap.parse_args(argv)

    project = Path(args.project)
    poses_dir, comps_dir = _dirs(project)
    receptors = sorted((project / lay.RECEPTORS).glob("*_ready.pdb")) or \
        sorted((project / "receptors").glob("*_ready.pdb")) or \
        sorted((project / lay.LEGACY[lay.RECEPTORS]).glob("*_listo.pdb"))
    if not receptors:
        raise SystemExit(f"No prepared receptor found under {project}.")

    if not rc.pdb2pqr_exe():
        raise SystemExit("PDB2PQR is not installed here. pip install pdb2pqr")

    smiles = _smiles_by_name(project)
    poses = sorted(p for p in poses_dir.glob("*-model1.pdb"))
    random.seed(args.seed)
    sample = sorted(random.sample(poses, min(args.n, len(poses))), key=lambda p: p.name)

    lines = ["# protonation comparison", f"# receptor: {receptors[0].name}",
             f"# sampled: {len(sample)} of {len(poses)} poses, seed {args.seed}", ""]
    agree = lost = gained = 0

    with tempfile.TemporaryDirectory() as td:
        work = Path(td)
        optimised = rc.optimize_hydrogens(receptors[0], work / "receptor_optimised.pdb")
        if optimised is None:
            raise SystemExit("PDB2PQR could not process the receptor; nothing to compare.")

        for pose in sample:
            key = re.sub(r"-model\d+$", "", pose.stem.split("compounds_a_")[-1])
            smi = smiles.get(sc.normalize_key(key))
            current = comps_dir / f"{lay.COMPLEX_PREFIX}{pose.stem}.pdb"
            if not current.exists():
                current = comps_dir / f"Complejo_{pose.stem}.pdb"
            if not smi or not current.exists():
                lines.append(f"{pose.stem}\tSKIPPED\tno smiles or no complex")
                continue
            ligand = lg.with_hydrogens(pose, smi, work / f"lig_{pose.stem}.pdb")
            if ligand is None:
                lines.append(f"{pose.stem}\tSKIPPED\tthe template does not match the pose")
                continue
            fused = ix.fuse(optimised, ligand, work / f"opt_{pose.stem}.pdb")
            a = _profile(current, work, nohydro=False)
            b = _profile(fused, work, nohydro=True)
            if a is None or b is None:
                lines.append(f"{pose.stem}\tSKIPPED\tPLIP produced no report")
                continue
            agree += len(a & b); lost += len(a - b); gained += len(b - a)
            lines.append(f"{pose.stem}")
            lines.append("  current  \t" + " ".join(sorted(a)))
            lines.append("  optimised\t" + " ".join(sorted(b)))

    lines += ["", "# totals",
              f"# agreeing\t{agree}", f"# only current\t{lost}", f"# only optimised\t{gained}"]
    text = "\n".join(lines) + "\n"
    if args.out:
        Path(args.out).write_text(text, encoding="utf-8")
        print(f"Written to {args.out}. Run this on the other machine too and diff the two files.")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
