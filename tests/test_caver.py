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


def test_the_config_asks_for_the_directory_everything_downstream_names(tmp_path):
    """Without save_dynamics_visualization, CAVER writes only clusters_timeless and the run looks
    like it found nothing. Caught by running it in the container, not by reading the manual."""
    from poliscreen.core.docking import Box

    assert "save_dynamics_visualization yes" in caver.write_config(
        tmp_path, Box(0, 0, 0, 24, 24, 24)).read_text()


def test_tunnels_are_found_in_either_directory(tmp_path):
    """Which one CAVER wrote depends on a config flag, not on the run. A config that came from
    somewhere else is not ours to assume anything about."""
    data = tmp_path / "data"
    (data / "clusters_timeless").mkdir(parents=True)
    (data / "clusters_timeless" / "tun_cl_001_1.pdb").write_text("")
    assert [p.name for p in caver.clusters(tmp_path)] == ["tun_cl_001_1.pdb"]

    # When both exist, the plain one wins: it is what the run folder names are built from.
    (data / "clusters").mkdir()
    (data / "clusters" / "tun_cl_001.pdb").write_text("")
    assert [p.name for p in caver.clusters(tmp_path)] == ["tun_cl_001.pdb"]


def test_an_empty_cluster_directory_is_not_a_success(tmp_path):
    (tmp_path / "data" / "clusters").mkdir(parents=True)
    assert caver.clusters(tmp_path) == []


PDB = """\
ATOM      1  N   ASP A 199      -1.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ASP A 199       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  H   ASP A 199       0.000   5.000   0.000  1.00  0.00           H
ATOM      4  CA  LEU A 209       0.000   4.000   0.000  1.00  0.00           C
HETATM    5  PA  GDP A 401      10.000  10.000  10.000  1.00  0.00           P
HETATM    6 CA    CA A 402      20.000  20.000  20.000  1.00  0.00          CA
HETATM    7  O   HOH A 501      30.000  30.000  30.000  1.00  0.00           O
END
"""


def test_hydrogens_go_because_they_close_the_narrow_tunnels(tmp_path):
    """CAVER measures the space between van der Waals spheres. On 8HTB the docking-ready file,
    carrying 2237 added hydrogens, finds three tunnels where the bare protein finds six."""
    src = tmp_path / "in.pdb"
    src.write_text(PDB)
    out = caver.prepare_for_caver(src, tmp_path / "out.pdb")
    text = out.read_text()
    assert " H   ASP" not in text
    assert " CA  ASP" in text                     # a calcium-named carbon alpha is not a hydrogen


def test_heterogroups_are_chosen_and_waters_go_by_default(tmp_path):
    src = tmp_path / "in.pdb"
    src.write_text(PDB)

    bare = caver.prepare_for_caver(src, tmp_path / "bare.pdb").read_text()
    assert "GDP" not in bare and "HOH" not in bare and " CA A 402" not in bare

    kept = caver.prepare_for_caver(src, tmp_path / "kept.pdb", keep_hetero=("GDP",)).read_text()
    assert "GDP" in kept and " CA A 402" not in kept

    wet = caver.prepare_for_caver(src, tmp_path / "wet.pdb", keep_waters=True).read_text()
    assert "HOH" in wet


def test_the_starting_point_can_come_from_three_places(tmp_path):
    """The box centre is a cube's middle. A ligand or the catalytic residues say it precisely, and
    on 8HTB the residues land within 2 A of the point CaverWeb used."""
    from poliscreen.core.docking import Box

    assert caver.start_point(Box(1, 2, 3, 24, 24, 24)) == (1, 2, 3)
    assert caver.start_point((1.0, 2.0, 3.0)) == (1.0, 2.0, 3.0)
    assert caver.start_point([(0, 0, 0), (2, 4, 6)]) == (1.0, 2.0, 3.0)

    src = tmp_path / "in.pdb"
    src.write_text(PDB)
    # ASP199's two heavy atoms sit either side of the origin; the hydrogen is off at y=5.
    assert caver.start_point(caver.atoms_of(src, ["ASP199"])) == (0.0, 1.667, 0.0)


def test_a_starting_point_from_nothing_is_refused(tmp_path):
    src = tmp_path / "in.pdb"
    src.write_text(PDB)
    with pytest.raises(caver.CaverError):
        caver.start_point(caver.atoms_of(src, ["TRP999"]))


SDF = """control
  PoliScreen

  3  2  0  0  0  0  0  0  0  0999 V2000
  -14.9700  -13.7270   18.6430 C   0  0
  -15.9700  -13.7270   18.6430 O   0  0
  -13.9700  -13.7270   18.6430 N   0  0
  2  1  1  0
  3  2  2  0
M  END
$$$$
"""


def test_a_bond_is_not_an_atom(tmp_path):
    """A bond in an SDF is written `  2  1  1  0`: four numeric fields, which a permissive reader
    takes for a coordinate. On the real 33-atom control that made 69 atoms and moved the centre
    from inside the site to outside the protein, and CAVER reported the outside as one tunnel."""
    f = tmp_path / "control.sdf"
    f.write_text(SDF)
    atoms = caver.ligand_atoms(f)
    assert len(atoms) == 3
    assert caver.start_point(atoms) == (-14.97, -13.727, 18.643)


def test_a_format_that_cannot_be_read_says_so(tmp_path):
    """Silently returning nothing here is how a wrong centre gets computed from no atoms."""
    f = tmp_path / "control.xyz"
    f.write_text("3\n\nC 0 0 0\n")
    with pytest.raises(caver.CaverError):
        caver.ligand_atoms(f)


def test_a_point_in_open_space_is_caught_before_caver_runs(tmp_path):
    f = tmp_path / "in.pdb"
    f.write_text(PDB)
    assert not caver.inside_structure((500.0, 500.0, 500.0), f)
    # The fixture has only a handful of atoms, so nowhere in it counts as enclosed either.
    assert not caver.inside_structure((0.0, 0.0, 0.0), f)


def test_a_tunnel_is_read_as_spheres_the_viewer_can_draw(tmp_path):
    """CAVER writes a tunnel as a chain of spheres with the radius where the occupancy would be,
    which is the shape fpocket's alpha spheres already arrive in."""
    cluster = tmp_path / "tun_cl_001.pdb"
    cluster.write_text(
        "MODEL        0\n"
        "ATOM      1  H   FIL T   1     -13.771 -12.804  17.384        2.67\n"
        "ATOM      2  H   FIL T   1     -13.709 -13.277  17.523        2.43\n")
    spheres = caver.tunnel_spheres(cluster)
    assert len(spheres) == 2
    assert spheres[0] == (-13.771, -12.804, 17.384, 2.67)


def test_several_structures_are_numbered_so_the_snapshots_keep_their_order(tmp_path, monkeypatch):
    """A trajectory read out of order clusters tunnels across time steps that never followed one
    another. CAVER sorts the folder by file name, so the names have to carry the order."""
    monkeypatch.setattr(caver, "caver_jar", lambda: tmp_path / "caver.jar")
    monkeypatch.setattr(caver, "java_exe", lambda: "/usr/bin/java")
    (tmp_path / "caver.jar").write_bytes(b"")
    snaps = []
    for name in ("frame_b", "frame_a", "frame_c"):
        p = tmp_path / f"{name}.pdb"
        p.write_text(PDB)
        snaps.append(p)

    monkeypatch.setattr(caver.subprocess, "run",
                        lambda *_a, **_k: type("R", (), {"stdout": "", "stderr": ""})())
    with pytest.raises(caver.CaverError):          # no clusters, because nothing really ran
        caver.find_tunnels(snaps, (0, 0, 0), tmp_path / "out")
    written = sorted(p.name for p in (tmp_path / "out" / "structures").iterdir())
    assert written == ["0_frame_b.pdb", "1_frame_a.pdb", "2_frame_c.pdb"]


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
