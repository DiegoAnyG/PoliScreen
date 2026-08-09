"""On native Windows a drive letter is already the real path.

The WSL translation (C:\\... -> /mnt/c/...) is right when PoliScreen runs inside a distribution and
wrong when it is the one-click installer, where it buried the projects in C:\\mnt\\c\\Users\\...
and told the user the app "runs inside Linux (WSL)".
"""
from pathlib import Path

import pytest

from poliscreen.core import session as ss


@pytest.fixture
def on_windows(monkeypatch):
    """Only the flag: os.name itself feeds pathlib, and moving it makes every Path() raise."""
    monkeypatch.setattr(ss, "_ON_WINDOWS", True)


def test_a_drive_letter_survives_untouched(on_windows):
    p, notice = ss.normalize_path(r"C:\Users\Diego\poliscreen_jobs\080826")
    assert "mnt" not in str(p), f"the installer would write to {p}"
    assert not notice, "there is nothing to warn about: the path was already usable"


def test_another_drive_survives_too(on_windows):
    """Only that the drive is not rewritten. Whether `M:\\...` counts as absolute is decided by
    pathlib's flavour, which follows the real platform and cannot be simulated from Linux."""
    p, _ = ss.normalize_path(r"M:\PICS-DOCS\ESCUELA\PoliTest")
    assert "/mnt/m" not in str(p), f"the installer would write to {p}"
    assert str(p).endswith(r"M:\PICS-DOCS\ESCUELA\PoliTest"), p


def test_a_relative_path_hangs_off_the_projects_root(on_windows, monkeypatch, tmp_path):
    monkeypatch.setenv("POLISCREEN_PROJECTS", str(tmp_path))
    p, _ = ss.normalize_path("my_run")
    assert p == tmp_path / "my_run"


def test_wsl_translation_still_happens_off_windows(monkeypatch):
    monkeypatch.setattr(ss, "_ON_WINDOWS", False)
    p, notice = ss.normalize_path(r"C:\Users\Diego\demo")
    assert str(p) == "/mnt/c/Users/Diego/demo"
    assert notice, "reinterpreting a path silently is what loses the results"


def test_an_empty_box_means_todays_project(monkeypatch, tmp_path):
    monkeypatch.setenv("POLISCREEN_PROJECTS", str(tmp_path))
    p, _ = ss.normalize_path("")
    assert p.parent == tmp_path
    assert len(p.name) == 6 and p.name.isdigit(), p.name
