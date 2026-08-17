"""PoliScreen - reproducible virtual screening with analogue design and pharmacophoric quality control.

Provisional name. Target workflow (closed loop):
    design analogues -> filter by synthesizability -> prepare receptor -> dock
    -> profile interactions (PLIP) -> score quality vs. control -> ADMET/toxicity
"""

# Read from the installed distribution rather than written here twice. A literal in this file is
# how the interface, the citation panel and the exported Methods all reported 1.0.0 from a 1.0.2
# image: pyproject.toml was bumped and this was not, and nothing failed to say so.
try:                                              # noqa: SIM105
    from importlib.metadata import version as _version

    __version__ = _version("poliscreen")
except Exception:                                 # a source tree with nothing installed
    import re as _re
    from pathlib import Path as _Path

    _pp = _Path(__file__).resolve().parents[2] / "pyproject.toml"
    _m = _re.search(r'^version = "([^"]+)"', _pp.read_text(encoding="utf-8"), _re.M) if _pp.exists() else None
    __version__ = _m.group(1) if _m else "0+unknown"
