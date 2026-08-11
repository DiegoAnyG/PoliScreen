"""Projects must land where POLISCREEN_PROJECTS says.

The interface used to build the default project from `Path.home()` instead of `default_root()`.
Inside Docker that is the container's own home, not the mounted volume, so a whole screening was
written somewhere that disappeared with the container. The path is checked statically because the
symptom only shows up after the container is gone.
"""
import ast
import os
from pathlib import Path

from poliscreen.core import session as ss

APP = Path(__file__).resolve().parent.parent / "src" / "poliscreen" / "ui" / "streamlit_app.py"


def test_default_root_honours_the_mounted_volume(monkeypatch):
    monkeypatch.setenv("POLISCREEN_PROJECTS", "/data")
    assert ss.default_root() == Path("/data")
    assert str(ss.default_project()).startswith("/data" + os.sep)


def test_default_project_is_one_folder_per_day(monkeypatch):
    monkeypatch.setenv("POLISCREEN_PROJECTS", "/data")
    name = ss.default_project().name
    assert len(name) == 6 and name.isdigit(), name


def test_default_root_keeps_a_legacy_folder_that_holds_projects(monkeypatch, tmp_path):
    monkeypatch.delenv("POLISCREEN_PROJECTS", raising=False)
    monkeypatch.setattr(Path, "home", staticmethod(lambda: tmp_path))
    assert ss.default_root().name == ss.JOBS_DIR
    (tmp_path / ss._LEGACY_JOBS_DIR / "8HTB" / "receptors").mkdir(parents=True)
    (tmp_path / ss._LEGACY_JOBS_DIR / "8HTB" / "receptors" / "ready.pdb").write_text("END\n")
    assert ss.default_root().name == ss._LEGACY_JOBS_DIR


def test_an_empty_legacy_folder_does_not_capture_the_default(monkeypatch, tmp_path):
    """It comes back on its own: an interface left running across the rename recreates it as a
    bare date folder, and the old rule then handed it the default again on the next start."""
    monkeypatch.delenv("POLISCREEN_PROJECTS", raising=False)
    monkeypatch.setattr(Path, "home", staticmethod(lambda: tmp_path))
    (tmp_path / ss._LEGACY_JOBS_DIR / "081126").mkdir(parents=True)
    assert ss.default_root().name == ss.JOBS_DIR


def test_the_interface_never_builds_a_project_path_from_home():
    """`Path.home() / "poliscreen_..."` in the interface means default_root() was bypassed."""
    tree = ast.parse(APP.read_text(encoding="utf-8"))
    offenders = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.BinOp) or not isinstance(node.op, ast.Div):
            continue
        left = ast.dump(node.left)
        if "home" not in left:
            continue
        right = node.right
        if isinstance(right, ast.Constant) and str(right.value).startswith("poliscreen"):
            offenders.append(f"line {node.lineno}: home() / {right.value!r}")
    assert not offenders, "use session.default_root(): " + "; ".join(offenders)
