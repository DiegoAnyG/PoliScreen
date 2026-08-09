"""fpocket must be handed a POSIX-style path.

It splits the path on '/' to decide where its output goes. Given backslashes it finds none and the
Windows build dies with an access violation (0xC0000005) before writing anything — measured on 4D44
with the cross-compiled binary: "C:\\fp\\test.pdb" crashes, "C:/fp/test.pdb" returns its pockets.

The symptom was silent: detect() returned an empty list, which is also what a protein with no
cavities returns.
"""
from pathlib import Path

from poliscreen.core import pockets as pk


class _Run:
    """Captures the argv fpocket would be started with."""

    def __init__(self):
        self.argv = None

    def __call__(self, cmd, **kw):
        self.argv = cmd
        return type("P", (), {"returncode": 0, "stdout": "", "stderr": ""})()


def test_the_path_reaches_fpocket_with_forward_slashes(monkeypatch, tmp_path):
    spy = _Run()
    monkeypatch.setattr(pk, "fpocket_available", lambda: True)
    monkeypatch.setattr(pk.subprocess, "run", spy)
    pdb = tmp_path / "4D44.pdb"
    pdb.write_text("END\n")

    pk.detect(pdb)

    handed = spy.argv[spy.argv.index("-f") + 1]
    assert "\\" not in handed, f"backslashes crash the Windows build: {handed}"
    assert handed.endswith("4D44.pdb"), handed


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
