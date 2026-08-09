"""Everything under assets/ has to travel in the wheel.

package-data listed only data/*.csv, so the installer shipped without the logo, without the
wordmark and without 3Dmol-min.js. None of it fails loudly: the interface falls back to plain text
and the viewer quietly fetches the library from the internet, which is not what a build that claims
to be self-contained should do.
"""
import fnmatch
import sys
from pathlib import Path

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover
    import tomli as tomllib

ROOT = Path(__file__).resolve().parent.parent
PKG = ROOT / "src" / "poliscreen"


def _patterns() -> list:
    data = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    return data["tool"]["setuptools"]["package-data"]["poliscreen"]


def test_every_asset_is_covered_by_package_data():
    pats = _patterns()
    missed = [f.relative_to(PKG).as_posix()
              for f in sorted((PKG / "assets").iterdir()) if f.is_file()
              and not any(fnmatch.fnmatch(f.relative_to(PKG).as_posix(), p) for p in pats)]
    assert not missed, f"these would not reach the wheel: {missed}"


def test_the_files_the_interface_looks_for_are_there():
    """Named individually: they are found by name at runtime, and a rename is silent."""
    from poliscreen.core import viewer as vw
    assert vw.logo_path() is not None, "no logo.* in assets/"
    assert vw.wordmark_path() is not None, "no titulo.* in assets/"
    assert vw._JS_LOCAL.exists(), "without it py3Dmol loads 3Dmol.js from the network"
