"""Finding the engines, and refusing to run without them.

Neither CAVER nor CaverDock ships with PoliScreen -- one needs a JVM, the other is a Linux-only
image under an academic licence -- so everything here is about discovery and about the two flags
that decide whether the answer is right. The runs themselves are not exercised: they take minutes
and need the engines installed, which no CI machine has.

Both workarounds below were measured, not assumed, and each one silently changes the answer:
cd-analysis docks through the tunnel extension alone, and more than two MPI processes make the
seed meaningless.
"""
from pathlib import Path

import pytest

from poliscreen.core import caver


def test_the_engines_are_looked_for_and_not_assumed(tmp_path, monkeypatch):
    """A path that exists only on the machine the code was written on reports 'absent' forever."""
    monkeypatch.delenv("POLISCREEN_CAVER", raising=False)
    monkeypatch.delenv("POLISCREEN_CAVERDOCK", raising=False)
    monkeypatch.setattr(caver.shutil, "which", lambda _n: None)
    monkeypatch.setattr(Path, "home", classmethod(lambda _cls: tmp_path))
    assert caver.caver_jar() is None
    assert caver.caverdock_image() is None


def test_the_environment_variable_wins(tmp_path, monkeypatch):
    jar = tmp_path / "caver" / "caver.jar"
    jar.parent.mkdir(parents=True)
    jar.write_bytes(b"")
    sif = tmp_path / "caverdock-1.2.sif"
    sif.write_bytes(b"")

    monkeypatch.setenv("POLISCREEN_CAVER", str(jar))
    monkeypatch.setenv("POLISCREEN_CAVERDOCK", str(sif))
    assert caver.caver_jar() == jar
    assert caver.caverdock_image() == sif


def test_a_folder_of_caver_is_accepted_as_well_as_the_jar(tmp_path, monkeypatch):
    """Nobody remembers whether the variable wants the jar or the folder holding it."""
    home = tmp_path / "caver"
    (home / "lib").mkdir(parents=True)
    (home / "caver.jar").write_bytes(b"")
    monkeypatch.setenv("POLISCREEN_CAVER", str(home))
    assert caver.caver_jar() == home / "caver.jar"


def test_the_config_puts_the_search_box_centre_at_the_starting_point(tmp_path):
    """CAVER measures from the active site outwards, and PoliScreen already knows where that is."""
    from poliscreen.core.docking import Box

    path = caver.write_config(tmp_path, Box(-15.1, -13.7, 18.8, 24, 24, 24))
    text = path.read_text()
    assert "starting_point_coordinates -15.1 -13.7 18.8" in text
    assert "seed 1" in text                       # every stage that touches coordinates is seeded


def test_the_config_can_be_overridden_without_editing_the_file(tmp_path):
    from poliscreen.core.docking import Box

    path = caver.write_config(tmp_path, Box(0, 0, 0, 24, 24, 24), probe_radius=1.2)
    assert "probe_radius 1.2" in path.read_text()
    assert "probe_radius 0.9" not in path.read_text()


def test_the_tunnel_command_keeps_the_whole_tunnel(tmp_path):
    """cd-analysis replaces the tunnel with its own extension: 10 discs of constant radius
    instead of 68, and it finishes in seconds, which is the only sign anything went wrong."""
    cmd = caver.transport_command(tmp_path / "img.sif", Path("/work"), "r.pdbqt", "l.pdbqt",
                                  "t.pdb", direction="in", bound="lb", cpus=2, seed=42)
    assert "--skip-tunnel-extension" in cmd


def test_more_than_two_processes_is_reported_as_losing_the_seed():
    """CaverDock says so once, in the middle of a run log nobody reads."""
    assert caver.reproducible(2) is True
    assert caver.reproducible(4) is False


def test_running_without_the_engine_says_which_one_and_how(tmp_path, monkeypatch):
    monkeypatch.delenv("POLISCREEN_CAVER", raising=False)
    monkeypatch.setattr(caver.shutil, "which", lambda _n: None)
    monkeypatch.setattr(Path, "home", classmethod(lambda _cls: tmp_path))
    with pytest.raises(caver.CaverError) as e:
        caver.find_tunnels(tmp_path / "r.pdb", None, tmp_path / "out")
    assert "POLISCREEN_CAVER" in str(e.value)
