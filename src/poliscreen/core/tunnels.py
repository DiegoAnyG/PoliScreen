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

import re
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
    "upper_bound_failed": "The upper bound was calculated and did not converge: with its rotation "
                          "constrained the ligand does not get past one of the discs. That is a "
                          "result, not a gap.",
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


def runs_in(folder, one_per_route: bool = True) -> list:
    """The finished transport runs under a folder.

    A run is a folder with a lower-bound trajectory in it; that is what everything downstream
    needs, and it is present whether or not the upper bound was asked for or succeeded.

    Asking for both bounds used to write two folders for one route, and the upper-bound one holds
    both profiles. One entry per route is what a person is choosing between, so the richer folder
    is the one kept.
    """
    root = Path(folder)
    if not root.is_dir():
        return []
    found = sorted({p.parent for p in root.rglob("*-lb.pdbqt")})
    if not one_per_route:
        return found
    best = {}
    for run in found:
        key = route_key(run)
        richer = bool(next(iter(run.glob("*-ub.dat")), None))
        if key not in best or (richer and not best[key][0]):
            best[key] = (richer, run)
    return [run for _richer, run in best.values()]


def profile_of(run_folder):
    """The energy profile of one run, or an empty list if it produced none."""
    if not available():
        return []
    from caver_translate.local import parse_local_job

    return parse_local_job(Path(run_folder)).profile


def orientation_of(run_folder) -> str:
    """Which end of the profile is the binding site, judged the way the tables judge it."""
    if not available():
        return "last"
    from caver_translate.local import parse_local_job
    from caver_translate.metrics import evaluate

    return evaluate(parse_local_job(Path(run_folder))).orientation or "last"


def chosen_poses(profile, bound: str = "last", extra: int = 0) -> list:
    """The states worth showing: (state, tag, label, reason).

    Three are always there -- the mouth, the barrier, the site -- because those are what a figure
    of a transport exists to show, and the barrier is the one an evenly spaced sample misses.
    `extra` adds context poses between them.
    """
    if not available() or not profile:
        return []
    from caver_translate.figures import choose_states

    return choose_states(profile, bound=bound, extra=extra)


def most_extra(profile, per_angstrom: float = 5.0) -> int:
    """The most context poses this tunnel can carry before they start overlapping.

    A short route says everything with three. Room for more is room in angstroms, so the ceiling
    comes from the tunnel rather than from one number applied to every route.
    """
    if not profile:
        return 0
    span = abs(profile[-1].distance - profile[0].distance)
    return max(0, min(6, int(span // per_angstrom) - 2))


def centroid_of(pdb_block: str) -> tuple:
    """The middle of a pose, which is where a label pointing at it belongs."""
    points = []
    for line in (pdb_block or "").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                points.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                continue
    if not points:
        return (0.0, 0.0, 0.0)
    n = len(points)
    return tuple(sum(p[i] for p in points) / n for i in range(3))


def suggested_extra(profile, per_angstrom: float = 5.0) -> int:
    """How many context poses a tunnel of this length can carry without becoming a crowd.

    Three poses say everything about a short route. A long one has room for more before they
    start overlapping: roughly one more for every five angstrom past fifteen.
    """
    if not profile:
        return 0
    span = abs(profile[-1].distance - profile[0].distance)
    return max(0, min(4, int((span - 15.0) // per_angstrom) + 1)) if span > 15.0 else 0


def pose_blocks(trajectory_pdbqt, states) -> list:
    """The named MODEL of a trajectory, as text a viewer can load.

    State N is model N, one-based, which is how PyMOL counts them and how choose_states reports
    them. Everything else in the file is the same molecule at another point of the route.
    """
    wanted = {int(s) for s in states}
    out, current, keeping, index = {}, [], False, 0
    for line in Path(trajectory_pdbqt).read_text(errors="ignore").splitlines():
        if line.startswith("MODEL"):
            index += 1
            keeping = index in wanted
            current = []
        elif line.startswith("ENDMDL"):
            if keeping:
                out[index] = "\n".join(current) + "\n"
            keeping = False
        elif keeping and line.startswith(("ATOM", "HETATM")):
            current.append(line)
    return [out.get(int(s), "") for s in states]


# The order poses appear in, matching the viewer so a figure and the screen agree.
POSE_COLORS = ["marine", "cyan", "green", "yellow", "orange", "magenta", "salmon", "purple"]

# The three points the profile is read for, and the colour each is marked with. The poses in the
# 3D view take the same colours, so a dot on the plot and a molecule on screen are the same thing
# said twice rather than two things to match up by reading.
LANDMARKS = {"surface": "#1971C2", "barrier": "#E03131", "site": "#2F9E44"}

# choose_states' tag for each of those three, and the colours for everything else it may add.
TAG_COLOR = {"start": LANDMARKS["surface"], "barrier": LANDMARKS["barrier"],
             "end": LANDMARKS["site"], "lowest": "#9C36B5"}
CONTEXT_COLORS = ["#E8A33D", "#00A5A5", "#C25E8C", "#7A8B3D", "#B5651D", "#5C7CBA"]


def pose_color(tag: str, index: int = 0) -> str:
    """The colour a pose is drawn in: its landmark's, or one that is nobody else's."""
    if tag in TAG_COLOR:
        return TAG_COLOR[tag]
    return CONTEXT_COLORS[index % len(CONTEXT_COLORS)]


def short_name(run_folder) -> str:
    """A run named for reading: `Benzofuroxan - Tunnel 3 - 8HTB`.

    The folder is named r<receptor>-l<ligand>-t<tunnel>-d<direction>-<bound>, which is what makes
    it machine-readable and what makes it unreadable on a chart title or in a menu.
    """
    name = Path(run_folder).name
    m = re.match(r"^r(?P<receptor>.*?)-l(?P<ligand>.*?)-t(?P<tunnel>.*?)"
                 r"-d(?P<direction>in|out)-(?P<bound>lowerbound|upperbound)$", name)
    if not m:
        return name
    receptor = re.sub(r"[_-]?(ready|listo|prepared)$", "", m.group("receptor"), flags=re.I)
    digits = re.search(r"\d+", m.group("tunnel"))
    tunnel = f"Tunnel {int(digits.group())}" if digits else m.group("tunnel")
    arrow = "in" if m.group("direction") == "in" else "out"
    return f"{m.group('ligand')} · {tunnel} · {receptor} ({arrow})"


def route_key(run_folder) -> tuple:
    """What makes two runs the same route: everything but which bound was asked for."""
    name = Path(run_folder).name
    m = re.match(r"^r(?P<receptor>.*?)-l(?P<ligand>.*?)-t(?P<tunnel>.*?)"
                 r"-d(?P<direction>in|out)-(?P<bound>lowerbound|upperbound)$", name)
    return (m.group("receptor"), m.group("ligand"), m.group("tunnel"),
            m.group("direction")) if m else (name,)


def landmarks(profile, bound: str = "last") -> dict:
    """Where the three numbers come from, as (distance from the site, energy) for each.

    Distance is measured from the active site outwards, which is the axis the profile is read
    along. A run stored the other way round has the same three points at the same energies; only
    the direction of travel differs, and taking the file's own order would flip the plot.
    """
    if not profile:
        return {}
    n = len(profile)
    surface_i = 0 if bound == "last" else n - 1
    site_i = n - 1 if bound == "last" else 0
    barrier_i = max(range(n), key=lambda i: profile[i].energy_lb)
    site_at = profile[site_i].distance
    return {name: (abs(site_at - profile[i].distance), profile[i].energy_lb)
            for name, i in (("surface", surface_i), ("barrier", barrier_i), ("site", site_i))}


def draw_profile(profile, bound: str = "last", title: str = "", figsize=(5.6, 3.4)):
    """The energy along the tunnel, with the three points that produce the numbers marked.

    The plot is the argument: Ea and dE_BS are differences between points on this curve, and a
    reader who cannot see which points they came from has to take them on trust.
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=figsize)
    if not profile:
        ax.text(0.5, 0.5, "no profile", ha="center", va="center", transform=ax.transAxes)
        return fig

    n = len(profile)
    site_i = n - 1 if bound == "last" else 0
    site_at = profile[site_i].distance
    x = [abs(site_at - p.distance) for p in profile]
    y = [p.energy_lb for p in profile]
    order = sorted(range(n), key=lambda i: x[i])
    ax.plot([x[i] for i in order], [y[i] for i in order], color="#495057", lw=1.6, zorder=2)

    # The upper bound, where there is one, on the same axes. It is the same route with rotation
    # constrained, and the gap between the two is what "the lower bound can understate" means.
    #
    # Only the points that have a value: an upper bound that did not converge stops partway, and
    # the discs past that point carry None. Plotting the pair without filtering draws a line to
    # nowhere, or raises, depending on the version.
    ub = [(x[i], profile[i].energy_ub_min) for i in order
          if profile[i].energy_ub_min is not None]
    if ub:
        ax.plot([a for a, _b in ub], [b for _a, b in ub], color="#ADB5BD", lw=1.4,
                ls="--", zorder=1, label="upper bound")
        ax.plot([], [], color="#495057", lw=1.6, label="lower bound")
        ax.legend(frameon=False, fontsize=7, loc="best")
        if len(ub) < len(profile):
            # Said on the plot, because a line that simply stops looks like the end of the tunnel.
            stopped = max(a for a, _b in ub)
            ax.axvline(stopped, color="#ADB5BD", lw=0.8, ls=":", zorder=1)
            ax.annotate("upper bound stops here", (stopped, max(y)), fontsize=7,
                        color="#868E96", ha="right", va="top",
                        textcoords="offset points", xytext=(-4, -2))

    # The three points the reported numbers come from. They are lower-bound quantities, which is
    # why they sit on the lower-bound line and not between the two.
    for name, (at, energy) in landmarks(profile, bound).items():
        ax.scatter([at], [energy], s=52, color=LANDMARKS[name], zorder=3,
                   edgecolors="white", linewidths=0.8)
        ax.annotate(name, (at, energy), textcoords="offset points", xytext=(0, 9),
                    ha="center", fontsize=8, color=LANDMARKS[name], weight="bold")

    ax.set_xlabel("distance from the active site (A)", fontsize=8)
    ax.set_ylabel("energy (kcal/mol)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.25, lw=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    if title:
        ax.set_title(title, fontsize=9)
    fig.tight_layout()
    return fig


def pymol_script(run_folder, receptor, tunnel, states, out_path=None) -> str:
    """A self-contained .pml: loads the receptor, the tunnel and the chosen poses.

    Not the same thing as caver-translate's script, which draws into a CaverWeb session that is
    already open. This one carries its own `load` lines, because everything it needs is a file in
    the project and a figure should not depend on a session someone has to reconstruct.

    Paths are written relative to where the script is saved when that is possible, so the folder
    can be moved or sent to someone else whole.
    """
    run = Path(run_folder)
    base = Path(out_path).parent if out_path else run

    def ref(path):
        path = Path(path)
        try:
            return path.resolve().relative_to(base.resolve()).as_posix()
        except ValueError:
            return path.resolve().as_posix()

    trajectory = next(iter(sorted(run.glob("*-lb.pdbqt"))), None)
    lines = [
        "# Transport of one compound through one tunnel, written by PoliScreen.",
        "#",
        "# Each pose below is a point of the energy profile and the comment says why it is here.",
        "# To drop one, delete its lines. Everything it needs is loaded by this script, so it can",
        "# be run on its own: pymol this_file.pml",
        "",
        "reinitialize",
        f"load {ref(receptor)}, receptor",
        "hide everything, receptor",
        "show cartoon, receptor",
        "color grey80, receptor",
        # The route is the subject; the protein is context and gets out of its way.
        "set cartoon_transparency, 0.8, receptor",
        "remove receptor and solvent",
        "",
    ]
    if tunnel:
        lines += [
            f"load {ref(tunnel)}, tunnel",
            "hide everything, tunnel",
            # A mesh rather than a surface: a solid tube hides the poses it is drawn around.
            "show mesh, tunnel",
            "color grey60, tunnel",
            "set mesh_width, 0.4",
            "",
        ]
    if trajectory is not None:
        lines += [f"load {ref(trajectory)}, trajectory", "hide everything, trajectory", ""]
        for i, (state, tag, label, reason) in enumerate(states):
            colour = POSE_COLORS[i % len(POSE_COLORS)]
            name = f"pose_{i + 1}_{tag}"
            lines += [
                f"# {reason}",
                f"create {name}, trajectory, {state}, 1",
                f"show sticks, {name}",
                f"color {colour}, {name} and elem C",
                f'pseudoatom label_{i + 1}, selection=({name}), label="{label}"',
                "",
            ]
        lines.append("delete trajectory")
    lines += ["", "orient receptor", "bg_color white", "set ray_opaque_background, 0"]
    text = "\n".join(lines) + "\n"

    if out_path is not None:
        out = Path(out_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        # ASCII only: a .pml read with the console code page turns a stray accent into mojibake
        # and the line it is on into a syntax error.
        out.write_text(text.encode("ascii", "replace").decode("ascii"), encoding="ascii")
    return text


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
