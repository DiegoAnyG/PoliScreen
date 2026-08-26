"""Labels for several poses of one route must not land on each other.

Every callout was offset the same way, which works until two poses are close -- and at the end of
a tunnel they always are. On the reference route the deepest point and the active site sit 1.3 A
apart, and their labels overlapped into one unreadable word: "active site int RDEST STEP".
"""
import math

from poliscreen.core.viewer import label_positions


def closest(points):
    return min(math.dist(a, b) for i, a in enumerate(points) for b in points[i + 1:])


def test_labels_end_up_further_apart_than_the_poses():
    """The two at the end of the route are the pair that overlapped."""
    poses = [(-9.3, -24.6, 13.4),      # entrance
             (-9.9, -22.6, 15.2),      # barrier
             (-16.2, -11.2, 19.6),     # deepest point
             (-15.2, -10.6, 19.0)]     # active site, 1.3 A from the one above
    labels = label_positions(poses)
    assert closest(poses) < 2.0, "the fixture is meant to have a crowded pair"
    assert closest(labels) > 4 * closest(poses)


def test_two_poses_in_the_same_place_still_get_separate_labels():
    """The degenerate case: identical centres, which a fixed offset maps onto one point."""
    labels = label_positions([(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)])
    assert closest(labels) > 8.0


def test_a_label_stays_near_the_thing_it_points_at():
    """Spread them too far and the arrow crosses the protein to reach its pose."""
    poses = [(0.0, 0.0, 0.0), (10.0, 0.0, 0.0), (0.0, 10.0, 0.0)]
    for pose, label in zip(poses, label_positions(poses)):
        assert math.dist(pose, label) < 20.0


def test_no_points_no_labels():
    assert label_positions([]) == []
