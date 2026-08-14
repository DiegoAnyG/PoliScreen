"""Choosing an engine has to set every switch that engine needs.

gnina takes two: the build argument and the overlay that reserves the GPU. An image built with
one of them reports gnina as not installed, and the report arrives at the end of a 4.5 GB build.
"""
import os
import subprocess
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "install-docker.sh"


def _command(choice: str) -> str:
    env = {**os.environ, "POLISCREEN_PRINT_ONLY": "1"}
    out = subprocess.run(["bash", str(SCRIPT), choice], capture_output=True, text=True, env=env)
    assert out.returncode == 0, out.stderr
    return out.stdout.strip().splitlines()[-1]


def test_gnina_is_asked_for_with_the_overlay_that_gives_it_the_gpu():
    for choice in ("3", "4"):
        assert "docker-compose.gpu.yml" in _command(choice)


def test_the_base_image_does_not_drag_in_the_gpu_overlay():
    for choice in ("1", "2"):
        assert "docker-compose.gpu.yml" not in _command(choice)


def test_an_answer_outside_the_menu_stops_instead_of_guessing():
    out = subprocess.run(["bash", str(SCRIPT), "9"], capture_output=True, text=True)
    assert out.returncode == 2
