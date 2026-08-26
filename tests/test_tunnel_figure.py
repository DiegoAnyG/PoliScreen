"""Showing a transport: which poses, and a script that stands on its own.

Three poses are the point of the picture -- the mouth, the barrier, the site -- and the barrier is
the one an evenly spaced sample misses, which is the whole reason the choice is made from the
energy profile rather than by counting.
"""
import pytest

from poliscreen.core import tunnels as tn

pytestmark = pytest.mark.skipif(not tn.available(),
                                reason="caver-translate is not installed in this environment")

# One MODEL per disc. The barrier is disc 2, deliberately not in the middle.
DISCS = [(0, -3.2, 1.5), (1, -2.7, 1.6), (2, -0.5, 1.7), (3, -4.0, 2.2), (4, -6.7, 2.7)]
DSD = "\n".join(f"{i}.0 0.0 0.0 1.0 0.0 0.0 {1.5 + i * 0.1}" for i in range(5)) + "\n"


def model(disc, energy, radius):
    return "\n".join([
        f"MODEL {disc + 1}",
        f"REMARK CAVERDOCK TUNNEL: {disc}     {energy}      {radius}      1.3",
        f"ATOM      1  C   LIG A   1      {disc}.000   0.000   0.000  1.00  0.00           C",
        "ENDMDL"])


@pytest.fixture()
def run(tmp_path):
    folder = tmp_path / "r8HTB-lbenzo-ttun_cl_003-din-upperbound"
    folder.mkdir(parents=True)
    (folder / "analysis-lb.pdbqt").write_text(
        "\n".join(model(d, e, r) for d, e, r in DISCS) + "\n")
    (folder / "tunnel.dsd").write_text(DSD)
    (folder / "8HTB_ready.pdb").write_text("ATOM      1  CA  ALA A   1       0.0   0.0   0.0\nEND\n")
    (folder / "tun_cl_003.pdb").write_text(
        "ATOM      1  H   FIL T   3       0.0   0.0   0.0        2.00\n")
    return folder


def test_the_barrier_is_always_among_the_poses(run):
    """Even spacing lands on it by luck. Disc 2 is the highest energy and is not the middle."""
    poses = tn.chosen_poses(tn.profile_of(run), bound=tn.orientation_of(run))
    tags = [tag for _s, tag, _l, _r in poses]
    assert "barrier" in tags
    barrier = next(s for s, tag, _l, _r in poses if tag == "barrier")
    assert barrier == 3, "state is one-based, so the disc-2 barrier is state 3"


def test_three_poses_without_asking_for_more(run):
    poses = tn.chosen_poses(tn.profile_of(run), bound=tn.orientation_of(run))
    assert 3 <= len(poses) <= 4      # the deepest point joins only if it is deeper than the site


def test_a_longer_tunnel_is_offered_more_poses():
    """Three say everything about a short route; a long one has room before they overlap."""
    from caver_translate.parse import Point

    short = [Point(distance=float(i), disc=i, radius=1.5, energy_lb=-1.0) for i in range(11)]
    long_ = [Point(distance=float(i), disc=i, radius=1.5, energy_lb=-1.0) for i in range(31)]
    assert tn.suggested_extra(short) == 0
    assert tn.suggested_extra(long_) >= 1


def test_a_pose_comes_back_as_the_atoms_of_that_state(run):
    """State N is model N, one-based. Off by one here draws the wrong point of the route."""
    blocks = tn.pose_blocks(run / "analysis-lb.pdbqt", [1, 3])
    assert len(blocks) == 2
    assert "      0.000   0.000   0.000" in blocks[0]     # disc 0
    assert "      2.000   0.000   0.000" in blocks[1]     # disc 2, the barrier


def test_the_script_carries_everything_it_needs(run, tmp_path):
    """caver-translate's script draws into a session someone must already have open. This one is
    for a person who has the folder and nothing else."""
    poses = tn.chosen_poses(tn.profile_of(run), bound=tn.orientation_of(run))
    out = tmp_path / "figure.pml"
    text = tn.pymol_script(run, run / "8HTB_ready.pdb", run / "tun_cl_003.pdb", poses, out)

    assert out.exists()
    assert "load" in text and "reinitialize" in text
    assert text.count("create pose_") == len(poses)
    assert "show mesh, tunnel" in text              # a solid tube hides the poses inside it
    assert "set cartoon_transparency, 0.8, receptor" in text


def test_the_script_is_ascii_only(run, tmp_path):
    """A .pml read with the console code page turns a stray accent into a syntax error."""
    poses = tn.chosen_poses(tn.profile_of(run), bound=tn.orientation_of(run))
    out = tmp_path / "figure.pml"
    tn.pymol_script(run, run / "8HTB_ready.pdb", run / "tun_cl_003.pdb", poses, out)
    out.read_bytes().decode("ascii")                # raises if anything is not


def test_the_script_points_at_files_beside_it(run, tmp_path):
    """So the folder can be moved, or sent to someone, whole."""
    poses = tn.chosen_poses(tn.profile_of(run), bound=tn.orientation_of(run))
    text = tn.pymol_script(run, run / "8HTB_ready.pdb", run / "tun_cl_003.pdb", poses,
                           run / "figure.pml")
    assert "load 8HTB_ready.pdb, receptor" in text
    assert str(tmp_path) not in text, "an absolute path was written where a relative one fits"


def test_runs_are_found_by_their_trajectory(run, tmp_path):
    assert tn.runs_in(tmp_path) == [run]


def test_the_three_landmarks_are_where_the_numbers_come_from(run):
    """Ea and dE_BS are differences between these points. A reader who cannot see which points
    they came from has to take them on trust."""
    profile = tn.profile_of(run)
    marks = tn.landmarks(profile, bound=tn.orientation_of(run))
    assert set(marks) == {"surface", "barrier", "site"}

    _at, e_surface = marks["surface"]
    _at, e_barrier = marks["barrier"]
    site_at, e_site = marks["site"]
    assert e_barrier == -0.5                     # the highest point of the fixture
    assert site_at == 0.0, "distance is measured from the site, so the site is at zero"
    assert pytest.approx(e_barrier - e_surface) == 2.7      # Ea
    assert pytest.approx(e_site - e_surface) == -3.5        # dE_BS


def test_the_plot_draws_and_marks_three_points(run):
    profile = tn.profile_of(run)
    fig = tn.draw_profile(profile, bound=tn.orientation_of(run), title="test")
    ax = fig.axes[0]
    marked = sum(len(c.get_offsets()) for c in ax.collections)
    assert marked == 3
    assert "active site" in ax.get_xlabel()


def test_an_empty_profile_still_returns_a_figure():
    """The panel asks before it knows there is anything to draw."""
    fig = tn.draw_profile([], bound="last")
    assert fig.axes
