"""Closing the window has to end the screening, not hide it.

A run spawns vina, obabel, plip and the design engine as grandchildren of the Streamlit server.
Terminating the server alone leaves them holding CPU and the project files while the user, having
closed the window, believes they are gone -- and the next start then finds the port taken by the
copy that is supposedly not running.

So the test is deliberately about the *grandchild*: killing the process you launched is easy, and
is not the thing that was going wrong.
"""
import os
import subprocess
import sys
import time

import pytest

from poliscreen.ui import desktop

pytestmark = pytest.mark.skipif(os.name == "nt", reason="POSIX process groups; Windows uses taskkill")

# A stand-in for the real shape: Streamlit, with a docking tool underneath it that ignores the
# polite signal. Only a group kill reaches the second one. The PID goes through a file because
# spawn_group deliberately leaves the streams alone -- the server's output belongs on the console.
PARENT = """
import subprocess, sys, time
child = subprocess.Popen([sys.executable, "-c",
    "import signal, time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(120)"])
open(sys.argv[1], "w").write(str(child.pid))
time.sleep(120)
"""


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _wait_for(path, timeout=10.0):
    deadline = time.time() + timeout
    while time.time() < deadline:
        if path.exists() and path.read_text().strip():
            return int(path.read_text().strip())
        time.sleep(0.05)
    raise AssertionError("the stand-in never reported its child")


def test_the_grandchild_dies_with_the_window(tmp_path):
    pidfile = tmp_path / "grandchild.pid"
    proc = desktop.spawn_group([sys.executable, "-c", PARENT, str(pidfile)])
    try:
        grandchild = _wait_for(pidfile)
        assert _alive(grandchild), "the stand-in did not start properly"
        desktop.kill_group(proc, grace=3.0)
        assert proc.poll() is not None, "the server itself survived"
        time.sleep(0.3)
        assert not _alive(grandchild), "the docking tool outlived the window that started it"
    finally:
        desktop.kill_group(proc, grace=1.0)


def test_spawning_puts_the_server_in_its_own_group():
    """The group is what makes the kill reach a grandchild whose parent already exited."""
    proc = desktop.spawn_group([sys.executable, "-c", "import time; time.sleep(30)"])
    try:
        assert os.getpgid(proc.pid) != os.getpgid(os.getpid())
    finally:
        desktop.kill_group(proc, grace=3.0)


def test_the_group_is_still_reachable_after_the_server_was_waited_on(tmp_path):
    """The browser-tab route waits on the server, then cleans up in its `finally`.

    Waiting reaps the PID, and the PID is the only way to look the group up -- which is precisely
    when the tools that outlived the server still need ending. So the group is recorded at spawn.
    """
    pidfile = tmp_path / "grandchild.pid"
    proc = desktop.spawn_group([sys.executable, "-c", PARENT, str(pidfile)])
    grandchild = _wait_for(pidfile)
    proc.terminate()
    proc.wait()                                   # the PID is gone from here on
    desktop.kill_group(proc, grace=3.0)
    time.sleep(0.3)
    assert not _alive(grandchild), "nothing could reach the group once the server had been reaped"


def test_killing_twice_is_harmless():
    """It runs from a `finally`, and may run again after the process is already gone."""
    proc = desktop.spawn_group([sys.executable, "-c", "import time; time.sleep(30)"])
    desktop.kill_group(proc, grace=3.0)
    desktop.kill_group(proc, grace=3.0)
    desktop.kill_group(None)


def test_the_browser_window_gets_a_profile_of_its_own():
    """Without one, an already-running browser takes the URL and exits at once.

    Waiting on that would return immediately and shut the server down under a window the user has
    only just seen open.
    """
    command = desktop._browser_app_command("http://localhost:8501", "/tmp/profile")
    if command is None:
        pytest.skip("no chromium-family browser on this machine")
    assert any(a.startswith("--user-data-dir=") for a in command), command
    assert any(a.startswith("--app=") for a in command), command


def test_waiting_gives_up_when_the_server_dies():
    """Otherwise a broken install is two minutes of blank window before anything is said."""
    proc = subprocess.Popen([sys.executable, "-c", "raise SystemExit(1)"])
    proc.wait()
    started = time.time()
    assert desktop.wait_until_serving(1, tries=100, pause=0.05, proc=proc) is False
    assert time.time() - started < 2, "it waited out the whole timeout on a dead server"
