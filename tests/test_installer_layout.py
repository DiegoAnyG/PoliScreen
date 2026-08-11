"""What the installer puts on disk has to be where the code goes looking for it.

OPSIN was searched for only under ~/adme/tools, which exists on a machine that installed admelab
by hand and on no packaged install at all, so every IUPAC name came back unparseable rather than
unavailable: 0 of 0 named on Windows against 28 of 44 on Linux, same compounds. The recipe and the
lookup are edited in different files, months apart, and nothing at import time connects them --
this is what connects them.
"""
import re
import sys
from pathlib import Path

import pytest

from poliscreen.core import naming

ROOT = Path(__file__).resolve().parent.parent
RECIPE = ROOT / "installer" / "construct.yaml"


def _extra_files() -> dict:
    """{source: destination} from the recipe's extra_files, selectors stripped."""
    out, inside = {}, False
    for line in RECIPE.read_text(encoding="utf-8").splitlines():
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


def test_opsin_lands_where_naming_looks_for_it():
    dest = _extra_files().get("vendor/opsin.jar")
    assert dest, "the recipe stopped shipping opsin.jar"
    assert Path(dest).name == "opsin.jar", f"naming.py looks for that exact name, not {dest}"
    installed = Path(sys.prefix) / dest
    searched = [d / "opsin.jar" for d in naming._tool_dirs()]
    assert installed in searched, (
        f"the installer puts it at {dest} and naming.py never looks there: "
        f"{[str(p) for p in searched]}")


def test_the_jvm_that_runs_it_is_shipped():
    """A jar with no Java is dead weight, and reads exactly like a name that will not parse."""
    specs = RECIPE.read_text(encoding="utf-8")
    assert re.search(r"^\s+- openjdk\b", specs, re.M), "opsin.jar needs a JVM to run"


@pytest.mark.skipif(not naming.available(), reason="needs OPSIN and a JVM")
def test_opsin_still_reads_a_name_it_has_to_get_right():
    """The benzofuroxan N-oxide locant, which is the reason names are verified at all."""
    smiles = naming.names_to_smiles(["ethyl benzoate"])[0]
    assert smiles, "OPSIN parsed nothing"
    assert naming._key(smiles) == naming._key("CCOC(=O)c1ccccc1")
