"""fpocket must not be handed a path at all.

It works out where its output goes by taking apart the path it is given, and on Windows every step
of that goes wrong. It splits on '/', so backslashes leave it with nothing and the cross-compiled
binary dies with an access violation. Handed forward slashes it gets past that and then creates
each component of "C:/.../<name>_out" in turn, starting with the bare drive letter "C:", which
fails whenever the process has no current directory on that drive; it then returns without writing
anything and still exits 0. Measured with the binary on 4D44: the same file, same run, produced its
pockets from the folder that held it and nothing at all from an absolute path.

Both symptoms were silent: detect() returned an empty list, which is also what a protein with no
cavities returns.
"""
from pathlib import Path

from poliscreen.core import pockets as pk


class _Run:
    """Captures how fpocket would be started."""

    def __init__(self):
        self.argv = None
        self.kw = None
        self.found_its_input = False

    def __call__(self, cmd, **kw):
        # Checked here: detect() removes the folder before returning, so afterwards there is
        # nothing left to look at.
        self.argv, self.kw = cmd, kw
        self.found_its_input = (Path(kw["cwd"]) / cmd[cmd.index("-f") + 1]).exists()
        return type("P", (), {"returncode": 0, "stdout": "", "stderr": ""})()


def test_fpocket_is_given_a_bare_name_inside_its_own_folder(monkeypatch, tmp_path):
    spy = _Run()
    monkeypatch.setattr(pk, "fpocket_available", lambda: True)
    monkeypatch.setattr(pk.subprocess, "run", spy)
    pdb = tmp_path / "4D44.pdb"
    pdb.write_text("END\n")

    pk.detect(pdb)

    handed = spy.argv[spy.argv.index("-f") + 1]
    assert handed == "4D44.pdb", f"a path is what fpocket takes apart wrongly: {handed}"
    assert spy.kw.get("cwd"), "without a cwd the output lands wherever the app was started"
    assert spy.found_its_input, "the copy is not in the folder fpocket runs from"


def test_a_run_that_wrote_nothing_says_so(monkeypatch, tmp_path):
    """Silence here is indistinguishable from a protein with no cavities."""
    def _run(cmd, **kw):
        return type("P", (), {"returncode": -1073741819, "stdout": "",
                              "stderr": "***** POCKET HUNTING BEGINS *****"})()

    monkeypatch.setattr(pk, "fpocket_available", lambda: True)
    monkeypatch.setattr(pk.subprocess, "run", _run)
    pdb = tmp_path / "4D44.pdb"
    pdb.write_text("END\n")

    said = []
    assert pk.detect(pdb, on_notice=said.append) == []
    assert said, "the crash was swallowed"
    assert "4D44.pdb" in said[0] and "POCKET HUNTING BEGINS" in said[0], said[0]
