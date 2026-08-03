"""A session is meant to be shared: it must not carry the folder layout of the machine.

Absolute paths travelled inside run.json and the saved widget values, so sending a
.poliscreen file to a colleague or attaching it to a paper disclosed the author's
directory tree.
"""
import json
import re
import zipfile

from poliscreen.core import session as ss

PERSONAL = re.compile(r"(/home/[^/\s\"']+|/Users/[^/\s\"']+|/mnt/[a-z]/|[A-Za-z]:\\\\)")


def test_paths_inside_the_project_become_relative():
    proj = "/home/someone/projects/demo"
    state = {"vis_rec_sel": f"{proj}/receptors/8HTB_ready.pdb",
             "pockets": {f"{proj}/receptors/8HTB_ready.pdb": [1, 2]}}
    clean = ss.strip_paths(state, proj)
    assert not PERSONAL.search(json.dumps(clean))
    assert ss.restore_paths(clean, proj) == state


def test_a_path_outside_the_project_keeps_only_its_name():
    clean = ss.strip_paths({"src_pdb": "/home/someone/Downloads/4D44.pdb"}, "/tmp/proj")
    assert clean == {"src_pdb": "4D44.pdb"}


def test_restoring_into_another_folder_rewrites_the_paths():
    """The same session opened on another machine must point at the new folder."""
    original = {"vis_rec_sel": "/home/a/proj/receptors/x.pdb"}
    clean = ss.strip_paths(original, "/home/a/proj")
    assert ss.restore_paths(clean, "/data/b") == {"vis_rec_sel": "/data/b/receptors/x.pdb"}


def test_a_saved_session_discloses_no_local_path(tmp_path):
    proj = tmp_path / "proj"
    (proj / "receptors").mkdir(parents=True)
    (proj / "receptors" / "8HTB_ready.pdb").write_text("END\n")
    (proj / "ranking.csv").write_text("compound,best_dock\na,-8.0\n")
    state = ss.strip_paths({"widgets": {"vis_rec_sel": str(proj / "receptors" / "8HTB_ready.pdb")}}, proj)
    dest = ss.save_session(proj, tmp_path / "s", state_=state)
    with zipfile.ZipFile(dest) as zf:
        for member in zf.namelist():
            if member.endswith((".json", ".txt")):
                text = zf.read(member).decode("utf-8", "ignore")
                assert not PERSONAL.search(text), f"{member} leaks a local path"
