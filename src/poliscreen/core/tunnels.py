"""Transport tunnels: reading what CAVER and CaverDock produced.

A docking score says how well a compound sits in the site. It says nothing about whether the
compound can *reach* it. CAVER finds the routes through the protein and CaverDock costs one of
them, and the two together answer a question the rest of this pipeline cannot ask.

**PoliScreen does not run either of them.** It reads their output, which is the half that costs
nothing: no engine to install, no licence to accept, no gigabyte in the image, and the same
behaviour in the container, in a development checkout and in the one-click installer. Running them
is a separate, opt-in step -- the pattern already used for ADCP -- and a folder of results reads
identically whether it was produced here or downloaded from CaverWeb.

The reading itself lives in caver-translate, a separate GPL-3 package with no dependencies of its
own. It is imported rather than run out-of-process, because unlike admelab it drags nothing in that
could conflict: the whole of it is the standard library. When it is absent this module reports so
and the interface says how to install it, exactly as the ADMET bridge does.

What the numbers mean, and what quietly ruins them, is in caver-translate's own documentation. The
short version, because it decides how the table is read:

    Ea      = E_max - E_surface     what entering costs -- the number that compares tunnels
    dE_BS   = E_bound - E_surface   how much better the site is than the outside

Neither is a binding free energy. They compare; they do not measure.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

# Why a row cannot be read at face value. Kept here rather than imported so the interface can show
# it without caver-translate installed, and so the wording is translated with the rest of the UI.
FLAG_TEXT = {
    "short_tunnel": "Shorter than 2 A: this is the mouth of the pocket, not a route through the "
                    "protein. It scores well because there is no distance to cross.",
    "positive_surface": "The energy at the mouth is positive, so the ligand already clashes there. "
                        "dE_BS looks favourable only because a positive number was subtracted.",
    "direction_mismatch": "The tunnel widens the opposite way to what the file name claims. One of "
                          "the two is wrong; the radius was trusted.",
    "orientation_from_name": "Constant radius, so which end is the binding site was taken from the "
                             "file name rather than from the geometry.",
    "lower_bound_only": "No upper-bound trajectory. The lower bound can pass through "
                        "discontinuities and understate a barrier.",
    "failed": "No profile in this calculation. A combination that fails leaves no log behind.",
}

INSTALL_HINT = "pip install git+https://github.com/DiegoAnyG/caver-translate.git"


def available() -> bool:
    """Whether the reader is installed in this interpreter's environment.

    find_spec locates the package without importing it, so a Results tab that never opens the
    tunnels sub-tab pays nothing for it.
    """
    from importlib.util import find_spec
    try:
        return find_spec("caver_translate") is not None
    except (ImportError, ValueError):
        return False


def read(folder) -> tuple:
    """Every tunnel calculation under this folder, as a table and a coverage count.

    Accepts what CaverWeb hands back -- one sub-folder per receptor, each holding hash-named
    archives -- and what a local run leaves behind, without being told which it is.

    The coverage count is not decoration. A combination that fails leaves nothing behind at all,
    so the gap in the ligand x tunnel x direction grid is the only evidence it was attempted, and
    a table that silently lists what succeeded reads as a complete study.
    """
    from caver_translate.metrics import coverage
    from caver_translate.parse import scan
    from caver_translate.report import COLUMNS, rows

    tunnels_, jobs = scan(Path(folder))
    records = rows(tunnels_, jobs) if jobs else []
    table = pd.DataFrame(records, columns=COLUMNS)
    cov = coverage(jobs) if jobs else {"missing": [], "duplicated": [], "expected": 0, "present": 0}
    return table, cov


def export(folder, out_dir) -> Path:
    """The two tables and the one page, written where the caller asks.

    The page is a single self-contained file -- the profiles are SVG written by hand, no plotting
    library -- so it can be handed to someone who has none of this installed.
    """
    from caver_translate.parse import scan
    from caver_translate.report import COLUMNS, TUNNEL_COLUMNS, rows, tunnel_rows, write_csv, write_html

    out = Path(out_dir)
    tunnels_, jobs = scan(Path(folder))
    write_csv(out / "transport.csv", rows(tunnels_, jobs), COLUMNS)
    if tunnels_:
        write_csv(out / "tunnels.csv", tunnel_rows(tunnels_), TUNNEL_COLUMNS)
    write_html(out / "report.html", tunnels_, jobs)
    return out


def flags_in(table: pd.DataFrame) -> list:
    """The distinct caveats present in a table, worst first is not meaningful here -- all matter."""
    if table.empty or "flags" not in table:
        return []
    seen = set()
    for value in table["flags"].fillna(""):
        seen.update(str(value).split())
    return sorted(seen)


def version() -> Optional[str]:
    """Which caver-translate produced these numbers, for the methods section."""
    if not available():
        return None
    from caver_translate import __version__
    return __version__
