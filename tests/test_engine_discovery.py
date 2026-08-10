"""Where the bridge looks for admelab.

admelab lives in its own environment because ADMET-AI drags in torch, which does not coexist with
openbabel/plip/vina. That environment was spelled `.venv/bin/python`, a path that cannot exist on
Windows, so on the one-click installer the engine was unreachable by construction — not merely
absent. Everything in admelab except the ADMET-AI layer is RDKit, which the installer already
ships, so a copy installed beside PoliScreen is worth finding.

A machine that has the venv keeps using it: that is the one with the ADMET predictions.
"""
import sys
from pathlib import Path

import pytest

from poliscreen.core import design as dz


@pytest.fixture
def no_environment(monkeypatch):
    monkeypatch.delenv("POLISCREEN_ADME_PYTHON", raising=False)
    monkeypatch.delenv("POLISCREEN_ADME_ROOT", raising=False)


def _package(root: Path) -> Path:
    (root / "admelab").mkdir(parents=True)
    return root


def test_the_venv_interpreter_is_spelled_for_this_platform():
    assert dz.DEFAULT_PYTHON.parent.name == ("Scripts" if sys.platform == "win32" else "bin"), \
        f"no interpreter on this platform can ever be at {dz.DEFAULT_PYTHON}"


def test_a_copy_installed_alongside_is_found_when_there_is_no_venv(no_environment, monkeypatch,
                                                                  tmp_path):
    site = _package(tmp_path / "site-packages")
    monkeypatch.setattr(dz, "DEFAULT_PYTHON", tmp_path / "absent" / "python")
    monkeypatch.setattr(dz, "_installed", lambda: (Path(sys.executable), site))

    bridge = dz.AdmelabBridge()

    assert bridge.available(), "the installer would report no design engine while shipping one"
    assert bridge.root == site


def test_the_venv_still_wins_where_it_exists(no_environment, monkeypatch, tmp_path):
    """It is the environment with ADMET-AI in it; the copy alongside has only the RDKit half."""
    venv = tmp_path / "adme" / ".venv" / "python"
    venv.parent.mkdir(parents=True)
    venv.write_text("")
    monkeypatch.setattr(dz, "DEFAULT_PYTHON", venv)
    monkeypatch.setattr(dz, "DEFAULT_ROOT", _package(tmp_path / "adme"))
    monkeypatch.setattr(dz, "_installed", lambda: (Path(sys.executable), tmp_path / "elsewhere"))

    assert dz.AdmelabBridge().python == venv


def test_the_environment_variables_still_have_the_last_word(monkeypatch, tmp_path):
    monkeypatch.setenv("POLISCREEN_ADME_PYTHON", str(tmp_path / "chosen"))
    monkeypatch.setenv("POLISCREEN_ADME_ROOT", str(tmp_path / "root"))
    monkeypatch.setattr(dz, "_installed", lambda: (Path(sys.executable), tmp_path / "elsewhere"))

    bridge = dz.AdmelabBridge()

    assert bridge.python == tmp_path / "chosen" and bridge.root == tmp_path / "root"


def test_looking_for_admelab_does_not_import_it(monkeypatch):
    """The docking process must stay clean: find_spec locates the package, nothing runs."""
    monkeypatch.delitem(sys.modules, "admelab", raising=False)
    dz._installed()
    assert "admelab" not in sys.modules
