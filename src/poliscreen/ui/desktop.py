"""Run the interface in a window, and make sure closing it stops everything it started.

Two separate problems, and the second is the one that bites.

**The window.** Streamlit only speaks HTTP, so something has to render it. pywebview draws a real
native window; on Windows it embeds the Edge WebView2 runtime that ships with the system. It is
not on conda-forge and pulls pythonnet on Windows, so it may simply be absent -- and a browser in
app mode (`--app=`) gives the same chrome-less window with nothing to install, because it is the
same engine either way. Both are tried, then an ordinary browser tab, so this never ends in a
blank screen.

**The processes.** A screening spawns vina, obabel, plip and the design engine as grandchildren of
Streamlit. Killing the Streamlit process alone leaves them running, holding CPU and the project
files, and the user has closed the window and believes they are gone. So the server is started in
its own process group and the whole group is taken down at once, which does not depend on the
intermediate process still being alive to walk the tree.
"""
from __future__ import annotations

import os
import shutil
import signal
import subprocess
import sys
import time
from typing import Optional, Sequence

WINDOWS = os.name == "nt"


def spawn_group(cmd: Sequence[str]) -> subprocess.Popen:
    """Starts the server so that it and every child can be terminated together."""
    if WINDOWS:
        # A new process group is what lets taskkill /T find the whole tree from this one PID.
        return subprocess.Popen(list(cmd), creationflags=subprocess.CREATE_NEW_PROCESS_GROUP)
    # setsid: the group survives whatever happens to the intermediate process, so a docking tool
    # that outlives Streamlit is still reachable.
    proc = subprocess.Popen(list(cmd), start_new_session=True)
    # Recorded now, while the process certainly exists. Once it has been waited on, its PID is
    # gone and with it the only route to the group -- which is exactly when the children that
    # outlived it still need ending.
    proc.poliscreen_pgid = os.getpgid(proc.pid)
    return proc


def _group_alive(pgid: int) -> bool:
    """True while any process in the group is still running."""
    try:
        os.killpg(pgid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def kill_group(proc: Optional[subprocess.Popen], grace: float = 5.0) -> None:
    """Ends the server and everything it started. Safe to call twice, and on an already-dead one.

    Asks the whole group first and insists on the whole group afterwards. Waiting on the server
    alone is what let a docking tool survive: Streamlit exits politely and immediately, so watching
    only its exit meant leaving before the tools underneath it had gone -- the very processes this
    is here to end. A tool is given the grace period to finish its write, because a pose file
    truncated mid-write later reads as a bad result rather than an interrupted one.
    """
    if proc is None:
        return
    if WINDOWS:
        # ponytail: taskkill walks the parent-PID chain, so a grandchild already orphaned when this
        # runs is missed. A Job Object would guarantee it; needs ctypes, and has not been worth it
        # while Streamlit is alive for the whole window's life and the tree is one level deep.
        subprocess.run(["taskkill", "/F", "/T", "/PID", str(proc.pid)],
                       capture_output=True, check=False)
        proc.poll()
        return

    pgid = getattr(proc, "poliscreen_pgid", None)
    if pgid is None:
        try:
            pgid = os.getpgid(proc.pid)
        except (ProcessLookupError, PermissionError):
            proc.poll()
            return
    try:
        os.killpg(pgid, signal.SIGTERM)
    except (ProcessLookupError, PermissionError):
        proc.poll()
        return
    deadline = time.time() + grace
    while time.time() < deadline and _group_alive(pgid):
        time.sleep(0.1)
    try:
        os.killpg(pgid, signal.SIGKILL)     # harmless once the group is empty
    except (ProcessLookupError, PermissionError):
        pass
    proc.poll()


def wait_until_serving(port: int, tries: int = 240, pause: float = 0.5,
                       proc: Optional[subprocess.Popen] = None) -> bool:
    """True once the port answers. False if it never does, or the server died trying.

    Watching the process as well as the port matters: a server that exits immediately -- a bad
    port, a broken install -- would otherwise be waited on for two minutes behind a blank window.
    """
    from ..cli import port_in_use
    for _ in range(tries):
        if port_in_use(port):
            return True
        if proc is not None and proc.poll() is not None:
            return False
        time.sleep(pause)
    return False


def _browser_app_command(url: str, profile: str) -> Optional[list]:
    """A browser told to open one page as its own window: no tabs, no address bar.

    The same engine pywebview embeds, reached without installing anything. Edge is on every
    supported Windows; the others are there for a Linux desktop.

    The throwaway profile is not a detail. Launched without one, a browser that is already running
    hands the URL to the existing process and exits immediately -- so waiting on it would return at
    once, and the server would be shut down under a window that had only just appeared. A separate
    profile forces a process of its own, which is the thing there is to wait on.
    """
    candidates = ("msedge", "chrome", "google-chrome", "chromium", "chromium-browser", "brave")
    for name in candidates:
        found = shutil.which(name)
        if found:
            return [found, f"--app={url}", f"--user-data-dir={profile}", "--no-first-run",
                    "--no-default-browser-check"]
    return None


def open_window(url: str, title: str = "PoliScreen", width: int = 1280, height: int = 860) -> str:
    """Shows the interface and returns only when the user closes it.

    Returns which route was taken, so the caller can say so. The browser-tab route returns
    immediately -- there is no window to wait on -- and the caller keeps the server running.
    """
    try:
        import webview
    except ImportError:
        webview = None

    if webview is not None:
        webview.create_window(title, url, width=width, height=height)
        webview.start()                      # returns when the last window is closed
        return "pywebview"

    import tempfile
    with tempfile.TemporaryDirectory(prefix="poliscreen-window-") as profile:
        app = _browser_app_command(url, profile)
        if app:
            # Waited on, not fired and forgotten: this call returning is what tells the caller the
            # user closed the window, which is what triggers the shutdown.
            subprocess.run(app, check=False)
            return "browser-app"

    import webbrowser
    webbrowser.open(url)
    return "browser-tab"
