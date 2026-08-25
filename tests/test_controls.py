"""What counts as a control, and what only looks like one.

`ccd_template` caches the chemical component dictionary's entry beside the control it describes,
because a rerun should need no network and the project should record which definition it used. That
cache is an `.sdf` in the folder controls are read from, and the folder scan took every `.sdf` in
it.

It is not a control. Its coordinates are the dictionary's idealised conformer, generated from the
chemistry alone and sitting nowhere near the crystallographic pose. On 8HTB the two centres are
about 9 A apart. Loaded as a second control it moves the search box off the site, and a run then
docks against empty space with nothing in the interface saying so.
"""
import ast
from pathlib import Path

from poliscreen.core import receptor as rc

APP = Path(__file__).resolve().parent.parent / "src" / "poliscreen" / "ui" / "streamlit_app.py"


def test_the_dictionary_cache_is_recognised():
    assert rc.is_ccd_cache("ZI9_ccd.sdf")
    assert rc.is_ccd_cache(Path("/somewhere/receptors/ZI9_ccd.sdf"))
    assert not rc.is_ccd_cache("control_ZI9.sdf")
    assert not rc.is_ccd_cache("ligand.sdf")


def test_the_cache_is_written_under_the_name_the_filter_looks_for(tmp_path):
    """Two spellings of the same suffix is how this came back. One constant, used both ends."""
    assert rc.CCD_CACHE_SUFFIX == "_ccd.sdf"
    src = Path(rc.__file__).read_text(encoding="utf-8")
    assert '"_ccd.sdf"' not in src.split("CCD_CACHE_SUFFIX = ")[1].split("\n", 1)[1], (
        "the suffix is spelled out again somewhere below its own constant")


def test_the_control_scan_excludes_it():
    """Checked in the source: the scan runs at project load, which no unit test reaches."""
    tree = ast.parse(APP.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        target = node.targets[0]
        if (isinstance(target, ast.Subscript) and isinstance(target.slice, ast.Constant)
                and target.slice.value == "controls"):
            source = ast.dump(node)
            if "iterdir" in source:                       # the folder scan, not a filtered rebuild
                assert "is_ccd_cache" in source, (
                    "the control scan takes every .sdf in the receptors folder again")
                return
    raise AssertionError("the control folder scan was not found; did it move?")
