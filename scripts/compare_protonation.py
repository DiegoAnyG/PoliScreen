"""Thin wrapper: the comparison lives in the package, where the installed launcher can reach it.

    PoliScreen.bat protonation-check <project> --out windows.txt

is the same thing without a checkout, and is what to use on a machine that only has the installer.
This file stays because a repository checkout is the other way people arrive at it.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from poliscreen.cli import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main(["protonation-check"] + sys.argv[1:]))
