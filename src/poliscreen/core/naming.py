"""Check that a IUPAC name really describes the molecule that will be docked.

A name that parses is not necessarily the right name. For benzofuroxans the N-oxide
can sit on either ring nitrogen, and the two are different molecules with the same
formula: naming a 1-oxide as a 3-oxide produces a name OPSIN accepts and a structure
that is not the one being screened. Docking uses the SMILES, so the numbers are
unaffected, but the label in a table or a paper would be wrong.

So every name is round-tripped through OPSIN and compared by InChIKey against its
own product. When it does not match, the known variants are tried and the one that
matches wins; if none does, the name is dropped rather than shown as verified.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Optional, Sequence

TOOLS = Path.home() / "adme" / "tools"


def _tool_dirs() -> list:
    """Where OPSIN and its Java may live, most specific first.

    `~/adme/tools` is where a machine that installed admelab by hand keeps them, and was the only
    place looked at -- which is why a packaged install, having no `~/adme` at all, reported every
    name unparseable instead of missing. An installer puts them beside the environment instead.
    """
    env_root = Path(sys.prefix)
    return [TOOLS, env_root / "tools", env_root / "share" / "opsin", env_root / "Library" / "tools"]


def _opsin_jar() -> Optional[Path]:
    env = os.environ.get("POLISCREEN_OPSIN")
    if env and Path(env).exists():
        return Path(env)
    for d in _tool_dirs():
        jar = d / "opsin.jar"
        if jar.exists():
            return jar
    return None


def _java() -> Optional[str]:
    # java.exe on Windows, so the name is globbed rather than spelled out.
    for d in _tool_dirs():
        for pattern in ("*jre*/bin/java", "*jre*/bin/java.exe", "*jdk*/bin/java", "*jdk*/bin/java.exe"):
            for cand in sorted(d.glob(pattern)):
                if cand.exists():
                    return str(cand)
    for exe in (Path(sys.prefix) / "bin" / "java", Path(sys.prefix) / "Library" / "bin" / "java.exe"):
        if exe.exists():
            return str(exe)
    return shutil.which("java")


def available() -> bool:
    return _opsin_jar() is not None and _java() is not None


def names_to_smiles(names: Sequence[str], timeout: int = 120) -> list:
    """SMILES of each name, or None where OPSIN cannot parse it. Order is preserved."""
    jar, java = _opsin_jar(), _java()
    if not names or jar is None or java is None:
        return [None] * len(names)
    with tempfile.TemporaryDirectory() as td:
        src = Path(td) / "names.txt"
        src.write_text("\n".join(n or "" for n in names) + "\n", encoding="utf-8")
        try:
            r = subprocess.run([java, "-jar", str(jar), "-o", "smi", str(src)],
                               capture_output=True, text=True, timeout=timeout)
        except Exception:
            return [None] * len(names)
    out = [ln.strip() for ln in (r.stdout or "").splitlines()]
    out += [""] * (len(names) - len(out))
    return [(s or None) for s in out[:len(names)]]


def _key(smiles: Optional[str]) -> Optional[str]:
    if not smiles:
        return None
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import inchi
        RDLogger.DisableLog("rdApp.*")
        m = Chem.MolFromSmiles(smiles)
        return inchi.MolToInchiKey(m) if m is not None else None
    except Exception:
        return None


# The N-oxide sits on one of two nitrogens and the composed name does not always pick the right one.
_OXIDE = re.compile(r"\b(?P<loc>[13])-oxido-(?P<core>[\w,\[\]\-]*?)-(?P=loc)-ium\b")

# Trivial radical names arriving with the alcohol; those IUPAC retains are absent on purpose.
_PREFERRED_RADICAL = {
    "furfuryl": "furan-2-ylmethyl", "isoamyl": "3-methylbutyl", "isopentyl": "3-methylbutyl",
    "isobutyl": "2-methylpropyl", "isopropyl": "propan-2-yl", "sec-butyl": "butan-2-yl",
    "neopentyl": "2,2-dimethylpropyl", "phenethyl": "2-phenylethyl",
    "allyl": "prop-2-en-1-yl", "propargyl": "prop-2-yn-1-yl", "amyl": "pentyl",
    "lauryl": "dodecyl", "cetyl": "hexadecyl", "myristyl": "tetradecyl",
    "stearyl": "octadecyl", "oleyl": "octadec-9-en-1-yl",
}
_STRAIGHT = ("methyl", "ethyl", "propyl", "butyl", "pentyl", "hexyl", "heptyl", "octyl",
             "nonyl", "decyl", "undecyl", "dodecyl")
_REDUNDANT_LOCANT = re.compile(r"^1-(" + "|".join(_STRAIGHT) + r")\b", re.I)


def preferred_radical(name: str) -> Optional[str]:
    """The same name with the leading radical in its preferred IUPAC form, or None."""
    if not name:
        return None
    head = name.split(" ", 1)[0]
    rest = name[len(head):]
    low = head.lower()
    if low in _PREFERRED_RADICAL:
        return _PREFERRED_RADICAL[low] + rest
    fixed = _REDUNDANT_LOCANT.sub(r"\1", head)
    return fixed + rest if fixed != head else None


def _swap_oxide(name: str) -> Optional[str]:
    m = _OXIDE.search(name or "")
    if not m:
        return None
    other = "3" if m.group("loc") == "1" else "1"
    return name[:m.start()] + f"{other}-oxido-{m.group('core')}-{other}-ium" + name[m.end():]


def variants(name: str) -> list:
    """Names worth testing, best first: preferred radical before the given one, and for each
    the alternative N-oxide locant."""
    out: list = []
    for form in (preferred_radical(name), name):
        for cand in (form, _swap_oxide(form or "")):
            if cand and cand not in out:
                out.append(cand)
    return out


def verify(names: Sequence[Optional[str]], smiles: Sequence[str]) -> list:
    """[(name, verified)] per product.

    Every candidate for a product is tried in order of preference and the first one that
    round-trips to that exact structure wins. If none does, the name is dropped: a wrong
    name in a table is worse than no name.
    """
    names, smiles = list(names), list(smiles)
    targets = [_key(s) for s in smiles]

    pending: list = []
    for i, name in enumerate(names):
        for cand in variants(name or ""):
            pending.append((i, cand))
    parsed = names_to_smiles([c for _i, c in pending])

    best: dict = {}
    for (i, cand), smi in zip(pending, parsed):
        if i in best or not targets[i]:
            continue
        if _key(smi) == targets[i]:
            best[i] = cand
    return [(best.get(i), i in best) for i in range(len(names))]
