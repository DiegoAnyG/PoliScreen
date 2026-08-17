"""One version, in every file that names it.

Two tag builds failed and a third shipped the wrong number, from the same cause: the version was
written out by hand in five places and only some of them were bumped. installer/construct.yaml
asked for a wheel that no longer existed, so the installers never built; __init__.py still said
1.0.0, so the interface, the citation panel and the exported Methods reported 1.0.0 from a 1.0.2
image. Nothing failed loudly enough to be noticed before a release went out.
"""
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def _pyproject() -> str:
    m = re.search(r'^version = "([^"]+)"', (ROOT / "pyproject.toml").read_text(encoding="utf-8"), re.M)
    assert m, "pyproject.toml no longer declares a version"
    return m.group(1)


def test_the_package_reports_the_declared_version():
    """Read from the installed metadata; a literal here is what drifted."""
    import poliscreen
    assert poliscreen.__version__ == _pyproject()


def test_the_citation_matches():
    """Whoever cites the software has to cite the version they ran."""
    assert f"version: {_pyproject()}" in (ROOT / "CITATION.cff").read_text(encoding="utf-8")


def test_the_installer_recipe_matches():
    """It names the wheel by file name: a stale number is a FileNotFoundError halfway through."""
    text = (ROOT / "installer" / "construct.yaml").read_text(encoding="utf-8")
    version = _pyproject()
    assert f"version: {version}" in text, "the installer would be named after the wrong version"
    assert f"poliscreen-{version}-py3-none-any.whl" in text, (
        "the recipe asks for a wheel that will not be built")


def test_the_launcher_banner_matches():
    raw = (ROOT / "scripts" / "PoliScreen-Docker.bat").read_bytes().decode("ascii")
    assert f"v{_pyproject()}" in raw


def test_the_changelog_has_an_entry():
    """A release nobody described is a release nobody can tell apart from the last one."""
    text = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    assert f"## [{_pyproject()}]" in text
