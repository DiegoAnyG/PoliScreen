"""CLI de PoliScreen.

Todo lo que haga la interfaz web debe poder hacerse también aquí: la línea de comandos
es lo que hace el flujo scriptable, reproducible y citable en un articulo.
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

    # Herramientas externas del nucleo: se listan con su versión para verificar la instalación
    # y para poder citarlas en la sección de Metodos.
    herramientas = [
        ("vina", ["--version"], "docking", True),
        ("obabel", ["-V"], "conversion y protonacion", True),
        ("obrms", [], "RMSD entre poses (estabilidad de la confianza)", False),
        ("plip", [], "interacciones proteina-ligando", True),   # sin --versión limpio
        ("fpocket", [], "deteccion de cavidades", False),
    ]
    print("herramientas externas:")
    faltan = []
    for exe, flags, para, critica in herramientas:
        if not shutil.which(exe):
            print(f"  {exe:8s} NO ENCONTRADO  ({para})")
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
    print(f"motor de diseno/ADMET (admelab): {'disponible' if b.available() else 'NO disponible (opcional)'}")
    print(f"  python: {b.python}")
    print(f"  raiz  : {b.root}")
    if b.available():
        i = b.info()
        print(f"  entorno aislado: python {i.get('python')} | torch {i.get('torch')} | cuda: {i.get('cuda')}")
        print(f"  modulos: {', '.join(i.get('modules', []))}")
    else:
        print("  sin el, el cribado funciona; solo se desactivan diseno de analogos y ADMET.")
        print("  define POLISCREEN_ADME_PYTHON y POLISCREEN_ADME_ROOT si las rutas son otras.")

    if faltan:
        print(f"\nFALTAN herramientas necesarias: {', '.join(faltan)}. Ver docs/INSTALL.md.")
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
    print(f"generados: {r.n_generated} | puntuados: {r.n_scored}")
    if args.out:
        r.to_dataframe().to_csv(args.out, index=False)
        print("escrito:", args.out)
    else:
        keys = [k for k in ("SMILES", "score", "MW", "LogP", "LD50_mg_per_kg", "GHS_category") if k in r.columns]
        for row in r.rows[: args.top]:
            print(json.dumps({k: row.get(k) for k in keys}, ensure_ascii=False))
    return 0


def cmd_ui(args) -> int:
    import subprocess
    from pathlib import Path

    # Por defecto el servidor escucha SOLO en la maquina local (127.0.0.1): no es alcanzable
    # desde la red ni desde internet. --expose lo abre a la red local, que carece de
    # autenticación y solo debería usarse en una red de confianza.
    address = "0.0.0.0" if getattr(args, "expose", False) else "127.0.0.1"
    app = Path(__file__).parent / "ui" / "streamlit_app.py"
    cmd = [sys.executable, "-m", "streamlit", "run", str(app),
           "--server.headless", "true",
           "--server.address", address,
           "--server.port", str(args.port),
           "--browser.gatherUsageStats", "false"]
    if address == "127.0.0.1":
        print(f"Interfaz local en http://localhost:{args.port} (solo este equipo). Ctrl+C para cerrar.")
    else:
        print(f"AVISO: expuesto en la red local, sin autenticacion. http://<tu-ip>:{args.port}")
    return subprocess.call(cmd)


def cmd_prep(args) -> int:
    from pathlib import Path

    from .core import receptor as rc

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    src = rc.fetch_pdb(args.pdb_id, out) if args.pdb_id else Path(args.pdb)
    st = rc.inspect(src)
    print(st.summary())
    if args.list:
        print("\nUsa --chains, --keep-het y --extract con las claves de arriba "
              "(formato RESNAME|CADENA|NUMERO).")
        return 0
    dest = out / f"{src.stem}_listo.pdb"
    rc.prepare(src, dest, keep_chains=args.chains, keep_het=args.keep_het or (),
               ph=args.ph, add_missing_residues=args.add_loops)
    print(f"\nreceptor listo: {dest}")
    for key in (args.extract or []):
        het = st.find(key)
        if het is None:
            print(f"  aviso: no existe el grupo {key}")
            continue
        p = rc.extract_ligand(src, het, out / f"control_{het.resname}.sdf", ph=args.ph, smiles=args.smiles)
        print(f"  control extraido: {p}")
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
            raise SystemExit(f"--box invalido: {spec}. Formato: STEM:cx,cy,cz,sx,sy,sz")
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
        cols = [c for c in ("receptor", "compound", "best_dock", "best_inter", "efectividad_pct", "tipo")
                if c in r.ranking.columns]
        print()
        print(r.ranking[cols].head(15).to_string(index=False))
    print(f"\nresultados en: {r.out_dir}")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="poliscreen", description="Cribado virtual reproducible")
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="cmd", required=True)

    pi = sub.add_parser("info", help="estado del entorno y de los motores")
    pi.set_defaults(func=cmd_info)

    pdsg = sub.add_parser("design", help="genera analogos de una molecula lider con ADME y toxicidad")
    pdsg.add_argument("lead", help="SMILES de la molecula lider")
    pdsg.add_argument("--no-ml", action="store_true", help="solo descriptores RDKit (rapido, sin ADMET-AI)")
    pdsg.add_argument("--n-sub", type=int, nargs="+", default=[1], help="numero de sustituciones, p. ej. 1 2")
    pdsg.add_argument("--positions", type=int, nargs="*", default=None, help="puntos de crecimiento (indices de atomo)")
    pdsg.add_argument("--max-decor", type=int, default=300, help="tope de analogos por decoracion")
    pdsg.add_argument("--top", type=int, default=10, help="cuantas filas mostrar/guardar")
    pdsg.add_argument("--out", help="CSV de salida")
    pdsg.set_defaults(func=cmd_design)

    pui = sub.add_parser("ui", help="abre la interfaz grafica (solo en este equipo)")
    pui.add_argument("--port", type=int, default=8501)
    pui.add_argument("--expose", action="store_true",
                     help="escucha en la red local (sin autenticacion); por defecto solo 127.0.0.1")
    pui.set_defaults(func=cmd_ui)

    pprep = sub.add_parser("prep", help="prepara un receptor: descarga, inspecciona y limpia")
    pprep.add_argument("--pdb-id", help="identificador del PDB, p. ej. 4D44")
    pprep.add_argument("--pdb", help="archivo .pdb local (alternativa a --pdb-id)")
    pprep.add_argument("--out", required=True, help="carpeta de salida")
    pprep.add_argument("--list", action="store_true", help="solo inspeccionar y salir")
    pprep.add_argument("--chains", nargs="*", help="cadenas a conservar (por defecto todas)")
    pprep.add_argument("--keep-het", nargs="*", help="heterogrupos a conservar, p. ej. NAP|A|400")
    pprep.add_argument("--extract", nargs="*", help="heterogrupos a extraer como control, p. ej. JA3|A|1259")
    pprep.add_argument("--smiles", help="plantilla SMILES del ligando extraido, corrige los enlaces")
    pprep.add_argument("--ph", type=float, default=7.4)
    pprep.add_argument("--add-loops", action="store_true", help="reconstruir lazos ausentes")
    pprep.set_defaults(func=cmd_prep)

    prun = sub.add_parser("run", help="ciclo completo: diseno, docking, interacciones y ranking")
    prun.add_argument("--lead", help="SMILES de la molecula lider (omitelo si das --ligands)")
    prun.add_argument("--ligands", nargs="*", help="ligandos ya preparados")
    prun.add_argument("--receptor", nargs="+", required=True, help="receptores (.pdb)")
    prun.add_argument("--control", nargs="*", help="controles co-cristalizados: definen la huella de referencia")
    prun.add_argument("--out", required=True, help="carpeta de resultados")
    prun.add_argument("--catalytic", nargs="*", help="residuos ancla por receptor, p. ej. 4D44:Tyr157,Tyr147")
    prun.add_argument("--box", nargs="*", help="caja por receptor: STEM:cx,cy,cz,sx,sy,sz (por defecto automatica)")
    prun.add_argument("--n-analogs", type=int, default=20)
    prun.add_argument("--n-sub", type=int, nargs="+", default=[1])
    prun.add_argument("--no-ml", action="store_true", help="sin ADMET-AI (rapido)")
    prun.add_argument("--seed", type=int, default=42)
    prun.add_argument("--exhaustiveness", type=int, default=24)
    prun.add_argument("--poses", type=int, default=10)
    prun.add_argument("--cpu", type=int, default=1, help="hilos por acoplamiento; 1 = reproducible")
    prun.add_argument("--workers", type=int, default=0, help="acoplamientos en paralelo; 0 = automatico")
    prun.set_defaults(func=cmd_run)

    args = p.parse_args(argv)
    try:
        return args.func(args)
    except AdmelabError as e:
        print("ERROR:", e, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
