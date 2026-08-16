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
