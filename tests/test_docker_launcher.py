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


def test_the_banner_states_the_version_it_ships():
    """A launcher that does not say which version it is leaves support guessing."""
    import re
    version = re.search(r'^version = "([^"]+)"',
                        (ROOT / "pyproject.toml").read_text(encoding="utf-8"), re.M).group(1)
    assert f"v{version}" in LAUNCHER.read_text(encoding="utf-8"), (
        f"the banner does not name v{version}; it drifted from pyproject")


def test_it_is_seven_bit_ascii():
    """Block-drawing characters needed chcp 65001, and changing the code page part-way through a
    batch file shifts the parser's byte offset: lines split mid-command and every fragment came
    back as "is not recognized as an internal or external command"."""
    raw = LAUNCHER.read_bytes()
    assert not raw.startswith(bytes([0xEF, 0xBB, 0xBF])), "a BOM breaks the first line"
    try:
        raw.decode("ascii")
    except UnicodeDecodeError as e:
        raise AssertionError(f"non-ASCII byte at offset {e.start}; cmd needs a code page for it")


def test_it_uses_crlf():
    """cmd is unreliable with LF-only batch files, in the same way and with the same message."""
    raw = LAUNCHER.read_bytes()
    assert raw.count(b"\r\n") == raw.count(b"\n"), "some lines end in a bare LF"


def test_line_endings_are_pinned_for_checkout():
    """A .gitattributes is what stops the next clone handing out LF again."""
    ga = (ROOT / ".gitattributes").read_text(encoding="utf-8")
    assert "*.bat text eol=crlf" in ga


def test_the_banner_fits_a_default_console():
    """80 columns is the stock cmd width and a line that fills it wraps, doubling the banner."""
    import base64
    import re

    m = re.search(r"FromBase64String\('([A-Za-z0-9+/=]+)'\)", LAUNCHER.read_text(encoding="utf-8"))
    assert m, "the banner is no longer carried as base64"
    art = base64.b64decode(m.group(1)).decode("utf-8")
    for line in art.splitlines():
        drawn = re.sub(r"\[[0-9;]*m", "", line)
        assert len(drawn) < 80, f"{len(drawn)} columns wraps in an 80-column window"


def test_the_banner_never_reaches_cmd_as_text():
    """cmd reads this file with the console code page; the characters must not be in it at all."""
    LAUNCHER.read_bytes().decode("ascii")


def test_a_container_already_running_is_noticed():
    """One left from another window keeps port 8501, and the browser shows that build instead."""
    text = LAUNCHER.read_text(encoding="utf-8")
    assert "--name poliscreen" in text, "an unnamed container cannot be looked for"
    assert "choice /C YN" in text, "stopping someone's running screening must be asked, not assumed"


def test_only_our_own_stale_images_are_removed():
    """A blanket prune would take other projects' images on a machine that is not ours to tidy."""
    text = LAUNCHER.read_text(encoding="utf-8")
    assert "docker images ghcr.io/diegoanyg/poliscreen --filter \"dangling=true\"" in text
    assert "docker image prune" not in text, "too broad: that removes every project's leftovers"
