"""Crossing midnight must not look like losing the analysis.

Projects are named after the day they start, deliberately: a fixed name meant yesterday's
receptors and ligands were still sitting there and got silently reused. The cost is that reopening
the interface the next morning resolves to a brand-new empty folder, and the work finished last
night is one directory over with nothing on screen saying so.
"""
import os
import time
from pathlib import Path

from poliscreen.core import session as ss


def _project(root: Path, name: str, with_content: bool = True) -> Path:
    p = root / name
    (p / "receptors").mkdir(parents=True) if with_content else p.mkdir(parents=True)
    if with_content:
        (p / "receptors" / "8HTB_ready.pdb").write_text("END\n")
    return p


def test_it_points_at_yesterdays_work(tmp_path):
    old = _project(tmp_path, "080926")
    today = _project(tmp_path, "081026", with_content=False)
    assert ss.previous_project(today) == old


def test_the_newest_one_wins(tmp_path):
    older = _project(tmp_path, "080826")
    newer = _project(tmp_path, "080926")
    # mtime is what orders them, and both were just written.
    os.utime(older, (time.time() - 7200, time.time() - 7200))
    today = _project(tmp_path, "081026", with_content=False)
    assert ss.previous_project(today) == newer


def test_an_empty_neighbour_is_not_offered(tmp_path):
    """Pointing at another empty folder would be worse than saying nothing."""
    _project(tmp_path, "080926", with_content=False)
    today = _project(tmp_path, "081026", with_content=False)
    assert ss.previous_project(today) is None


def test_the_first_project_ever_has_nothing_to_point_at(tmp_path):
    assert ss.previous_project(_project(tmp_path, "081026", with_content=False)) is None
