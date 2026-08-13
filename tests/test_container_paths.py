"""In a container, a folder outside the mounted volume is a folder that will be thrown away.

Typing a Windows path into the container does not fail. The translation produces a valid Linux
path, nothing exists at it, so it is created -- inside the container. The screening runs, writes,
and reports success, and every file goes when the container stops. Nothing looks wrong until
everything is gone, and by then there is nothing left to look at.
"""
from pathlib import Path

import pytest

from poliscreen.core import session as ss


@pytest.fixture
def containerised(monkeypatch, tmp_path):
    mounted = tmp_path / "data"
    mounted.mkdir()
    monkeypatch.setenv("POLISCREEN_IN_CONTAINER", "1")
    monkeypatch.setenv("POLISCREEN_PROJECTS", str(mounted))
    return mounted


def test_a_folder_outside_the_mount_is_called_out(containerised, tmp_path):
    notice = ss.warn_if_not_persistent(tmp_path / "somewhere" / "081326")
    assert "disappears" in notice and str(containerised) in notice


def test_a_folder_inside_the_mount_is_silent(containerised):
    assert ss.warn_if_not_persistent(containerised / "081326") == ""


def test_the_mount_itself_is_silent(containerised):
    assert ss.warn_if_not_persistent(containerised) == ""


def test_nothing_is_said_outside_a_container(monkeypatch, tmp_path):
    """The installer writes wherever it is told; only a container throws the rest away."""
    monkeypatch.delenv("POLISCREEN_IN_CONTAINER", raising=False)
    monkeypatch.setattr(ss.Path, "exists", lambda self: False)
    assert ss.warn_if_not_persistent(tmp_path / "anywhere") == ""


def test_the_warning_reaches_the_person_typing_the_path(containerised, monkeypatch):
    """It has to ride on normalize_path's notice, which is what the interface shows."""
    monkeypatch.setattr(ss, "_ON_WINDOWS", False)
    _path, notice = ss.normalize_path(r"C:\Users\someone\Documents\screening")
    assert "disappears" in notice, notice
