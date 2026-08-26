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


def points(with_ub_until=None):
    from caver_translate.parse import Point

    out = []
    for i, (_d, e, r) in enumerate(DISCS):
        ub = e + 0.5 if (with_ub_until is None or i < with_ub_until) else None
        out.append(Point(distance=float(i), disc=i, radius=r, energy_lb=e,
                         energy_ub_min=ub, energy_ub_max=ub))
    return out


def test_an_upper_bound_that_stops_partway_does_not_break_the_plot():
    """The case that prompted this: a refused upper bound leaves the later discs empty, and
    plotting the pair without filtering draws a line to nowhere or raises."""
    fig = tn.draw_profile(points(with_ub_until=2), bound="last")
    ax = fig.axes[0]
    assert sum(len(c.get_offsets()) for c in ax.collections) == 3
    drawn = [ln for ln in ax.lines if len(ln.get_xdata())]
    assert any(len(ln.get_xdata()) == 2 for ln in drawn), "the partial upper bound is not drawn"


def test_a_profile_with_no_upper_bound_at_all_draws_one_line():
    fig = tn.draw_profile(points(with_ub_until=0), bound="last")
    assert len([ln for ln in fig.axes[0].lines if len(ln.get_xdata())]) == 1


def test_both_bounds_share_one_plot():
    """Two charts for one route is two things to line up by eye."""
    fig = tn.draw_profile(points(), bound="last")
    labels = [ln.get_label() for ln in fig.axes[0].lines]
    assert "upper bound" in labels and "lower bound" in labels


def test_a_pose_is_coloured_like_its_point_on_the_curve():
    """A dot on the plot and a molecule on screen should be one thing said twice."""
    assert tn.pose_color("start") == tn.LANDMARKS["surface"]
    assert tn.pose_color("barrier") == tn.LANDMARKS["barrier"]
    assert tn.pose_color("end") == tn.LANDMARKS["site"]
    # Anything else takes a colour that is none of theirs.
    others = {tn.pose_color("step", i) for i in range(6)} | {tn.pose_color("lowest")}
    assert not others & set(tn.LANDMARKS.values())


def test_a_run_is_named_for_reading():
    """`r8HTB_ready-lBenzofuroxan-ttun_cl_003-din-lowerbound` is not a chart title."""
    assert tn.short_name(
        "r8HTB_ready-lBenzofuroxan-ttun_cl_003-din-lowerbound"
    ) == "Benzofuroxan · Tunnel 3 · 8HTB (in)"
    assert tn.short_name("something_else") == "something_else"


def test_the_two_bounds_of_one_route_are_one_entry(tmp_path):
    """Six folders for three routes is a menu nobody can read."""
    for bound in ("lowerbound", "upperbound"):
        f = tmp_path / f"r8HTB-lbenzo-ttun_cl_003-din-{bound}"
        f.mkdir()
        (f / "analysis-lb.pdbqt").write_text(model(0, -1.0, 1.5) + "\n")
        if bound == "upperbound":
            (f / "analysis-ub.dat").write_text("0.0 0 -1.0 -1.0 1.5 -1.0\n")
    runs = tn.runs_in(tmp_path)
    assert len(runs) == 1
    assert runs[0].name.endswith("upperbound"), "the folder with both profiles is the one to keep"


def test_how_many_poses_fit_comes_from_the_tunnel():
    from caver_translate.parse import Point

    short = [Point(distance=float(i), disc=i, radius=1.5, energy_lb=-1.0) for i in range(11)]
    long_ = [Point(distance=float(i) * 2, disc=i, radius=1.5, energy_lb=-1.0) for i in range(21)]
    assert tn.most_extra(short) == 0
    assert tn.most_extra(long_) >= 4
