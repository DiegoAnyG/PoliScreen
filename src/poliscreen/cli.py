"""PoliScreen CLI.

Everything the web interface does must also be possible here: the command line is what makes the
workflow scriptable, reproducible and citable in an article.
"""
from __future__ import annotations

import argparse
import json
import sys

from . import __version__
from .core import AdmelabBridge, AdmelabError


def cmd_info(args) -> int:
    import shutil
    import subprocess

    print(f"PoliScreen {__version__}")

    tools = [
        ("vina", ["--version"], "docking", True),
        ("obabel", ["-V"], "conversion and protonation", True),
        ("obrms", [], "RMSD between poses (confidence stability)", False),
        ("plip", [], "protein-ligand interactions", True),
        ("fpocket", [], "cavity detection", False),
    ]
    print("external tools:")
    faltan = []
    for exe, flags, para, critica in tools:
        if not shutil.which(exe):
            print(f"  {exe:8s} NOT FOUND  ({para})")
            if critica:
                faltan.append(exe)
            continue
        ver = ""
        if flags:
            try:
                r = subprocess.run([exe] + flags, capture_output=True, text=True, timeout=15)
                ver = next((l.strip() for l in ((r.stdout or "") + (r.stderr or "")).splitlines() if l.strip()), "")
            except Exception:
                ver = ""
        print(f"  {exe:8s} OK  {ver[:60]}")

    b = AdmelabBridge()
    print(f"design/ADMET engine (admelab): {'available' if b.available() else 'NOT available (optional)'}")
    print(f"  python: {b.python}")
    print(f"  root  : {b.root}")
    if b.available():
        i = b.info()
        print(f"  isolated environment: python {i.get('python')} | torch {i.get('torch')} | cuda: {i.get('cuda')}")
        print(f"  modules: {', '.join(i.get('modules', []))}")
        if not i.get("torch"):
            # What the installer ships. Design and the RDKit descriptors work; the predicted
            # endpoints do not, and a table missing columns does not say why on its own.
            print("  no ADMET-AI there: analogue design and the RDKit descriptors work, the "
                  "predicted endpoints do not. See docs/INSTALL.md.")
    else:
        print("  without it, screening works; only analogue design and ADMET are disabled.")
        print("  set POLISCREEN_ADME_PYTHON and POLISCREEN_ADME_ROOT if the paths differ.")

    # gnina was missing from this report entirely, so a machine with an NVIDIA GPU that had
    # just installed PoliScreen had no way to find out why the second scoring never appeared.
    from .core import docking as dk
    g = dk.gnina_exe()
    if g:
        print(f"second scoring (gnina): {g}")
    elif sys.platform == "win32":
        print("second scoring (gnina): not on Windows -- it is only built for Linux. A GPU on "
              "this machine can still run it through Docker, docker-compose.gpu.yml.")
    else:
        print("second scoring (gnina): not installed (optional). scripts/get_gnina.sh, NVIDIA GPU.")

    if sys.platform == "win32":
        print("platform: native Windows build. Hydrogen-bond detection differs from the reference")
        print("  container -- 18.6% of contacts on identical input, all hydrogen-dependent, fewer")
        print("  found here. Fine to explore and triage; use the container for anything published.")

    if faltan:
        print(f"\nMISSING required tools: {', '.join(faltan)}. See docs/INSTALL.md.")
        return 1
    return 0


def cmd_design(args) -> int:
    b = AdmelabBridge()
    r = b.design(
        args.lead,
        use_ml=not args.no_ml,
        n_substitutions=args.n_sub,
        positions=args.positions,
        max_decor=args.max_decor,
        max_rows=args.top,
    )
    print(f"generated: {r.n_generated} | scored: {r.n_scored}")
    if args.out:
        r.to_dataframe().to_csv(args.out, index=False)
        print("written:", args.out)
    else:
        keys = [k for k in ("SMILES", "score", "MW", "LogP", "LD50_mg_per_kg", "GHS_category") if k in r.columns]
        for row in r.rows[: args.top]:
            print(json.dumps({k: row.get(k) for k in keys}, ensure_ascii=False))
    return 0


def port_in_use(port: int) -> bool:
    """Whether something is already listening there.

    Worth checking before starting: Streamlit exits when the port is taken, and in the Windows
    installer that closes the console before anyone can read why. What was serving the page then
    was another PoliScreen — a container Docker Desktop had brought back on its own — so the
    interface answered, from a different machine's filesystem, and the paths made no sense.
    """
    import socket
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.settimeout(0.5)
        return s.connect_ex(("127.0.0.1", port)) == 0


def _open_when_serving(url: str, port: int, tries: int = 120):
    """Opens the page in the browser once Streamlit answers, in the background.

    Double-clicking the launcher has to end in a browser: the address printed in a console is
    clickable only with Ctrl, which nobody outside a terminal knows. Streamlit's own auto-open is
    not usable here — it comes with --server.headless false, and that brings back its first-run
    e-mail prompt, which would sit there waiting on a console the user is not reading. Waiting for
    the port instead of guessing a delay means a slow first start opens a page, not an error.
    """
    import threading
    import time
    import webbrowser

    def wait():
        for _ in range(tries):
            if port_in_use(port):
                webbrowser.open(url)
                return
            time.sleep(0.5)

    t = threading.Thread(target=wait, daemon=True)
    t.start()
    return t


def cmd_ui(args) -> int:
    import subprocess
    from pathlib import Path

    if port_in_use(args.port):
        print(f"ERROR: something is already listening on port {args.port}, so this interface "
              f"cannot start.", file=sys.stderr)
        print("  Whatever answers at that address is NOT the copy you just launched.", file=sys.stderr)
        print("  A container left running is the usual cause; docker/docker-compose.yml sets "
              "`restart: unless-stopped`, so Docker Desktop brings it back by itself.", file=sys.stderr)
        print("  Check it with:   docker ps", file=sys.stderr)
        print("  Stop it with:    docker compose -f docker/docker-compose.yml down", file=sys.stderr)
        print(f"  Or use another:  poliscreen ui --port {args.port + 1}", file=sys.stderr)
        return 1

    address = "0.0.0.0" if getattr(args, "expose", False) else "127.0.0.1"
    app = Path(__file__).parent / "ui" / "streamlit_app.py"
    cmd = [sys.executable, "-m", "streamlit", "run", str(app),
           "--server.headless", "true",
           "--server.address", address,
           "--server.port", str(args.port),
           "--browser.gatherUsageStats", "false"]
    url = f"http://localhost:{args.port}"
    if getattr(args, "window", False):
        return _ui_in_a_window(cmd, url, args.port)

    if address == "127.0.0.1":
        print(f"Local interface at {url} (this machine only). Ctrl+C to close.")
        _open_when_serving(url, args.port)
    else:
        print(f"WARNING: exposed on the local network, no authentication. http://<your-ip>:{args.port}")
    return subprocess.call(cmd)


def _ui_in_a_window(cmd, url: str, port: int) -> int:
    """The interface as a desktop window, with the server tied to that window's lifetime.

    Kept apart from the ordinary path above, which is untouched: that one hands the terminal to
    Streamlit so Ctrl+C reaches it, and it is what the Windows launcher relies on for its
    "Terminate batch job" prompt. Here there is no terminal to hand over, so the server is a
    background process this function owns and has to clean up -- when the window is closed, on
    Ctrl+C, and on the way out of an exception. The `finally` is the whole point: a return path
    that skips it is a screening still running with nothing on screen.
    """
    from .ui import desktop

    proc = None
    try:
        proc = desktop.spawn_group(cmd)
        if not desktop.wait_until_serving(port, proc=proc):
            print("ERROR: the interface did not start. Run `poliscreen ui` without --window to "
                  "see what it says.", file=sys.stderr)
            return proc.returncode or 1
        if desktop.open_window(url) == "browser-tab":
            # Only a tab was opened, and a tab cannot be waited on. Falling back to the ordinary
            # behaviour beats exiting and killing the server the user is about to read.
            print(f"No window backend available; opened {url} in the browser. Ctrl+C to close.")
            return proc.wait()
    except KeyboardInterrupt:
        return 0
    finally:
        desktop.kill_group(proc)
    return 0


def cmd_protonation_check(args) -> int:
    """Whether pinning the protonation down costs contacts, and whether it buys portability.

    Reachable from the installed launcher on purpose: the machine that has to run it is the one
    with the divergence, and asking someone to clone a repository to answer that is how the
    question stays unanswered.
    """
    from pathlib import Path

    from .core import fingerprint as fp

    proj = Path(args.project)
    if not proj.is_dir():
        print(f"ERROR: {proj} is not a folder.", file=sys.stderr)
        return 1
    text = fp.protonation_report(proj, n=args.n, seed=args.seed)
    if args.out:
        Path(args.out).write_text(text, encoding="utf-8")
        print(f"Written to {args.out}. Run this on the other machine too and diff the two files.")
    else:
        print(text, end="")
    return 0


def cmd_fingerprint(args) -> int:
    """Print what the docking consumed and produced, hashed, so two machines can be compared.

    Run it on both and diff the output. The first line that differs names the stage that caused
    the divergence, which is the question that repeated guessing could not answer.
    """
    from pathlib import Path

    from .core import fingerprint as fp

    proj = Path(args.project)
    if not proj.is_dir():
        print(f"ERROR: {proj} is not a folder.", file=sys.stderr)
        return 1
    text = fp.render(proj)
    if args.out:
        Path(args.out).write_text(text, encoding="utf-8")
        print(f"Written to {args.out}. Run this on the other machine too and diff the two files.")
    else:
        print(text, end="")
    return 0


def cmd_prep(args) -> int:
    from pathlib import Path

    from .core import receptor as rc
    from .core import screening as sc

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    src = rc.fetch_pdb(args.pdb_id, out) if args.pdb_id else Path(args.pdb)
    st = rc.inspect(src)
    print(st.summary())
    if args.list:
        print("\nUse --chains, --keep-het and --extract with the keys above "
              "(format RESNAME|CHAIN|NUMBER).")
        return 0
    dest = out / f"{src.stem}{sc.READY_SUFFIX}.pdb"
    rc.prepare(src, dest, keep_chains=args.chains, keep_het=args.keep_het or (),
               ph=args.ph, add_missing_residues=args.add_loops)
    print(f"\nreceptor ready: {dest}")
    for key in (args.extract or []):
        het = st.find(key)
        if het is None:
            print(f"  notice: group {key} does not exist")
            continue
        p = rc.extract_ligand(src, het, out / f"control_{het.resname}.sdf", ph=args.ph, smiles=args.smiles)
        print(f"  control extracted: {p}")
    return 0


def cmd_run(args) -> int:
    from pathlib import Path

    from .core import pipeline as pl

    cat = {}
    for spec in (args.catalytic or []):
        if ":" in spec:
            k, v = spec.split(":", 1)
            cat[k] = [x.strip() for x in v.split(",") if x.strip()]

    receptors = [Path(p) for p in args.receptor]
    boxes = {}
    for spec in (args.box or []):
        stem, _, nums = spec.partition(":")
        vals = [float(v) for v in nums.split(",")]
        target = next((r for r in receptors if r.stem == stem or stem in r.stem), None)
        if target is None or len(vals) != 6:
            raise SystemExit(f"--box invalid: {spec}. Format: STEM:cx,cy,cz,sx,sy,sz")
        boxes[str(target)] = pl.dk.Box(*vals)

    cfg = pl.RunConfig(
        receptors=receptors,
        boxes=boxes,
        out_dir=Path(args.out),
        lead=args.lead,
        ligands=[Path(p) for p in (args.ligands or [])],
        controls=[Path(p) for p in (args.control or [])],
        catalytic=cat,
        n_analogs=args.n_analogs,
        n_substitutions=args.n_sub,
        use_ml=not args.no_ml,
        seed=args.seed,
        exhaustiveness=args.exhaustiveness,
        n_poses=args.poses,
        cpu=args.cpu,
        workers=args.workers,
    )
    r = pl.run(cfg, on_step=lambda n, d: print(f"[{n}] {d}"))
    if r.ranking is not None and not r.ranking.empty:
        cols = [c for c in ("receptor", "compound", "best_dock", "best_inter", "effectiveness_pct", "type")
                if c in r.ranking.columns]
        print()
        print(r.ranking[cols].head(15).to_string(index=False))
    print(f"\nresults in: {r.out_dir}")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="poliscreen", description="Reproducible virtual screening")
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="cmd", required=True)

    pi = sub.add_parser("info", help="environment and engine status")
    pi.set_defaults(func=cmd_info)

    pdsg = sub.add_parser("design", help="generate analogues of a lead molecule with ADME and toxicity")
    pdsg.add_argument("lead", help="SMILES of the lead molecule")
    pdsg.add_argument("--no-ml", action="store_true", help="RDKit descriptors only (fast, no ADMET-AI)")
    pdsg.add_argument("--n-sub", type=int, nargs="+", default=[1], help="number of substitutions, e.g. 1 2")
    pdsg.add_argument("--positions", type=int, nargs="*", default=None, help="growth points (atom indices)")
    pdsg.add_argument("--max-decor", type=int, default=300, help="cap on analogues per decoration")
    pdsg.add_argument("--top", type=int, default=10, help="how many rows to show/save")
    pdsg.add_argument("--out", help="output CSV")
    pdsg.set_defaults(func=cmd_design)

    pui = sub.add_parser("ui", help="open the graphical interface (this machine only)")
    pui.add_argument("--port", type=int, default=8501)
    pui.add_argument("--expose", action="store_true",
                     help="listen on the local network (no authentication); by default only 127.0.0.1")
    pui.add_argument("--window", action="store_true",
                     help="open in a desktop window instead of a browser tab; closing the window "
                          "stops the interface and everything it started")
    pui.set_defaults(func=cmd_ui)

    pfp = sub.add_parser("fingerprint",
                         help="hashes of what went into the docking, to compare two machines")
    pfp.add_argument("project", help="project folder")
    pfp.add_argument("--out", help="write to this file instead of the screen")
    pfp.set_defaults(func=cmd_fingerprint)

    ppc = sub.add_parser("protonation-check",
                         help="does pinning the protonation down cost contacts, and does it agree "
                              "across machines")
    ppc.add_argument("project", help="project folder with poses and fused complexes")
    ppc.add_argument("--out", help="write to this file instead of the screen")
    ppc.add_argument("--n", type=int, default=20, help="complexes to sample")
    ppc.add_argument("--seed", type=int, default=11, help="which ones; keep it equal on both machines")
    ppc.set_defaults(func=cmd_protonation_check)

    pprep = sub.add_parser("prep", help="prepare a receptor: download, inspect and clean")
    pprep.add_argument("--pdb-id", help="PDB identifier, e.g. 4D44")
    pprep.add_argument("--pdb", help="local .pdb file (alternative to --pdb-id)")
    pprep.add_argument("--out", required=True, help="output folder")
    pprep.add_argument("--list", action="store_true", help="inspect only and exit")
    pprep.add_argument("--chains", nargs="*", help="chains to keep (all by default)")
    pprep.add_argument("--keep-het", nargs="*", help="hetero groups to keep, e.g. NAP|A|400")
    pprep.add_argument("--extract", nargs="*", help="hetero groups to extract as control, e.g. JA3|A|1259")
    pprep.add_argument("--smiles", help="SMILES template of the extracted ligand, fixes the bonds")
    pprep.add_argument("--ph", type=float, default=7.4)
    pprep.add_argument("--add-loops", action="store_true", help="rebuild missing loops")
    pprep.set_defaults(func=cmd_prep)

    prun = sub.add_parser("run", help="full cycle: design, docking, interactions and ranking")
    prun.add_argument("--lead", help="SMILES of the lead molecule (omit it if you pass --ligands)")
    prun.add_argument("--ligands", nargs="*", help="already prepared ligands")
    prun.add_argument("--receptor", nargs="+", required=True, help="receptors (.pdb)")
    prun.add_argument("--control", nargs="*", help="co-crystallized controls: define the reference fingerprint")
    prun.add_argument("--out", required=True, help="results folder")
    prun.add_argument("--catalytic", nargs="*", help="anchor residues per receptor, e.g. 4D44:Tyr157,Tyr147")
    prun.add_argument("--box", nargs="*", help="box per receptor: STEM:cx,cy,cz,sx,sy,sz (automatic by default)")
    prun.add_argument("--n-analogs", type=int, default=20)
    prun.add_argument("--n-sub", type=int, nargs="+", default=[1])
    prun.add_argument("--no-ml", action="store_true", help="no ADMET-AI (fast)")
    prun.add_argument("--seed", type=int, default=42)
    prun.add_argument("--exhaustiveness", type=int, default=24)
    prun.add_argument("--poses", type=int, default=10)
    prun.add_argument("--cpu", type=int, default=1, help="threads per docking; 1 = reproducible")
    prun.add_argument("--workers", type=int, default=0, help="dockings in parallel; 0 = automatic")
    prun.set_defaults(func=cmd_run)

    args = p.parse_args(argv)
    try:
        return args.func(args)
    except AdmelabError as e:
        print("ERROR:", e, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
