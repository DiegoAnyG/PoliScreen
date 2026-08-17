"""The double-click route has to point at the image that actually gets published.

The Windows installer existed because "download it and run it" is the only setup most people will
do. The container gives the same thing without the native build's hydrogen-bond divergence, but
only if launching it is one double-click and the image is already built. Two names have to agree
for that: the one the workflow pushes and the one the launcher pulls. Nothing else checks them.
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LAUNCHER = ROOT / "scripts" / "PoliScreen-Docker.bat"
WORKFLOW = ROOT / ".github" / "workflows" / "installer.yml"


def _repo() -> str:
    """owner/name as GitHub knows it, from the file that has to be right for the citation."""
    cff = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    m = re.search(r"github\.com/([\w.-]+)/([\w.-]+)", cff)
    assert m, "CITATION.cff no longer names the repository"
    return f"{m.group(1)}/{m.group(2)}".lower()


def test_the_launcher_pulls_the_image_the_workflow_pushes():
    m = re.search(r"set \"IMAGE=ghcr\.io/([^:\"]+):", LAUNCHER.read_text(encoding="utf-8"))
    assert m, "the launcher no longer names a ghcr image"
    assert m.group(1) == _repo(), (
        f"launcher pulls ghcr.io/{m.group(1)}, workflow pushes ghcr.io/{_repo()}")


def test_the_workflow_lower_cases_the_name():
    """ghcr rejects capitals, and this account has one in it."""
    assert "tr '[:upper:]' '[:lower:]'" in WORKFLOW.read_text(encoding="utf-8")


def test_the_launcher_never_starts_docker_itself():
    """Starting a background service on someone's machine is not a launcher's business."""
    text = LAUNCHER.read_text(encoding="utf-8").lower()
    assert "start \"\" \"docker desktop" not in text
    assert "dockercli" not in text and "-switchdaemon" not in text


def test_the_port_stays_on_the_loopback():
    """The interface has no authentication, so it must not be published to the network."""
    assert "127.0.0.1:8501:8501" in LAUNCHER.read_text(encoding="utf-8")


def test_it_checks_for_a_newer_image_when_one_is_already_there():
    """Otherwise somebody runs a version fixed months ago and is never told."""
    text = LAUNCHER.read_text(encoding="utf-8")
    after_first_run = text.split("else (", 1)
    assert len(after_first_run) == 2, "no branch for the image already being present"
    assert "docker pull" in after_first_run[1], "an existing image is never refreshed"


def test_being_offline_is_not_treated_as_a_failure():
    """The screening itself needs no network; refusing to start without one would be wrong."""
    text = LAUNCHER.read_text(encoding="utf-8").lower()
    assert "no network" in text


def test_the_release_carries_the_launcher():
    """A double-click needs a direct download, not four clicks through a source tree."""
    wf = WORKFLOW.read_text(encoding="utf-8")
    assert "files: scripts/PoliScreen-Docker.bat" in wf
