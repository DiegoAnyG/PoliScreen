"""What the installer puts on disk has to be where the code goes looking for it.

OPSIN was searched for only under a hand-installed admelab's tools directory, which exists on the
machine it was written on and on no packaged install at all, so every IUPAC name came back
unparseable rather than unavailable: 0 of 0 named on Windows against 28 of 44 on Linux, same
compounds. The recipe and the lookup are edited in different files, months apart, and nothing at
import time connects them -- this is what connects them.

It is shipped by nobody at the moment: the JVM it needs costs more than the whole rest of the
environment and buys only a label. These assert the two halves stay consistent either way, because
half of it is the state that reads as broken rather than absent.
"""
import re
import sys
from pathlib import Path

from poliscreen.core import naming

ROOT = Path(__file__).resolve().parent.parent
RECIPE = ROOT / "installer" / "construct.yaml"


def _recipe() -> str:
    return RECIPE.read_text(encoding="utf-8")


def _extra_files() -> dict:
    """{source: destination} from the recipe's extra_files, selectors stripped."""
    out, inside = {}, False
    for line in _recipe().splitlines():
        if line.startswith("extra_files:"):
            inside = True
            continue
        if inside:
            if line and not line.startswith((" ", "\t")):
                break
            entry = re.match(r"\s+- (\S+?):\s*(\S+)", re.sub(r"\s+#.*$", "", line))
            if entry:
                out[entry.group(1)] = entry.group(2)
    return out


def _ships_opsin():
    return _extra_files().get("vendor/opsin.jar")


def _ships_a_jvm() -> bool:
    return bool(re.search(r"^\s+- openjdk\b", _recipe(), re.M))


def test_opsin_and_its_jvm_travel_together():
    """A jar with no Java cannot run, and a JVM with no jar is 500 MB of nothing.

    Either half alone is worse than neither: the naming does not degrade when OPSIN is missing,
    it inverts, and every name is dropped as if it were wrong.
    """
    assert bool(_ships_opsin()) == _ships_a_jvm(), (
        "the recipe ships one of opsin.jar / openjdk without the other")


def test_opsin_lands_where_naming_looks_for_it():
    dest = _ships_opsin()
    if not dest:
        return  # not shipped, which is the current decision
    assert Path(dest).name == "opsin.jar", f"naming.py looks for that exact name, not {dest}"
    searched = [d / "opsin.jar" for d in naming._tool_dirs()]
    assert Path(sys.prefix) / dest in searched, (
        f"the installer puts it at {dest} and naming.py never looks there: "
        f"{[str(p) for p in searched]}")


def test_the_screening_does_not_depend_on_opsin():
    """Whatever the recipe decides, docking reads the SMILES. The name is a label."""
    assert not _ships_opsin(), (
        "if OPSIN is shipped again, this test is the place to record why it became worth "
        "the JVM; the screening itself never needed it")
