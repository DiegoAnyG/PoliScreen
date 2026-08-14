"""Ctrl+C must stop the container and leave nothing running.

The Windows launcher already meets this: it asks "Terminate batch job (Y/N)?" and leaves no
orphans. The Docker route has to match, and two settings were in the way -- a restart policy that
brings the container back, and PID 1 being exempt from the default signal handlers so the SIGTERM
was ignored and the container was SIGKILLed with the run half-written.

Neither fails at build time, or on a short run, or anywhere a test would normally look: the
container comes back minutes later, or on the next Docker Desktop start. Parsed as text rather
than YAML on purpose -- the recipe test does the same, and neither is worth a dependency.
"""
import re
from pathlib import Path

DOCKER = Path(__file__).resolve().parent.parent / "docker"
COMPOSE = DOCKER / "docker-compose.yml"


def _settings() -> str:
    """The file with its comments removed, so prose about a setting is not read as the setting."""
    return "\n".join(re.sub(r"#.*$", "", line) for line in COMPOSE.read_text(encoding="utf-8").splitlines())


def test_nothing_restarts_the_container_behind_the_user():
    """A screening is something you start, not a service that comes back on its own."""
    policy = re.search(r"^\s*restart:\s*(\S+)", _settings(), re.M)
    assert policy is None or policy.group(1) == "no", (
        f"restart: {policy.group(1) if policy else ''} outlives Ctrl+C and starts again with "
        "Docker Desktop, which is the orphan this route is supposed not to leave")


def test_the_signal_reaches_streamlit():
    """PID 1 ignores signals it has no handler for, and Python installs none for SIGTERM."""
    assert re.search(r"^\s*init:\s*true\s*$", _settings(), re.M), (
        "without init, tini is not PID 1, the SIGTERM from Ctrl+C is ignored and the container "
        "is killed ten seconds later mid-write")


def test_the_inner_tini_can_reap_what_the_docking_leaves():
    """Two tinis stack: Docker's takes PID 1, the base image's lands behind it and needs telling."""
    assert "TINI_SUBREAPER=1" in (DOCKER / "Dockerfile").read_text(encoding="utf-8"), (
        "the image's own tini is not PID 1 under init: true, so without this it reaps nothing "
        "and warns about it at every start")
