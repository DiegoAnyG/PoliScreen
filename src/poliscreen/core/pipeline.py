"""Orquestador: encadena diseño, docking, interacciones y puntuación.

Escribe todo en carpetas, sin base de datos. Cada etapa es reanudable: si el resultado
ya esta en disco no se recalcula, así que volver a lanzar una corrida es barato.
"""
from __future__ import annotations

import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional, Sequence

import pandas as pd

from . import adcp
from . import docking as dk
from . import interactions as ix
from . import ligands as lg
from . import peptides as pep
from . import screening as sc
from . import validation as vl
from .design import AdmelabBridge

DEFAULT_WEIGHTS = dict(sc.DEFAULT_WEIGHTS)


@dataclass
class RunConfig:
    receptors: Sequence[Path]
    out_dir: Path
    lead: Optional[str] = None
    ligands: Sequence[Path] = field(default_factory=list)
    controls: Sequence[Path] = field(default_factory=list)
    control_map: dict = field(default_factory=dict)      # ckey -> stem del receptor
    boxes: dict = field(default_factory=dict)            # ruta receptor -> Box
    site_boxes: dict = field(default_factory=dict)       # ruta receptor -> [(etiqueta_sitio, Box)] (docking hibrido)
    catalytic: dict = field(default_factory=dict)        # stem receptor -> [residuos]
    weights: dict = field(default_factory=lambda: dict(DEFAULT_WEIGHTS))
    n_analogs: int = 20
    n_substitutions: Sequence[int] = (1,)
    use_ml: bool = True
    seed: int = 42
    exhaustiveness: int = 24
    n_poses: int = 10
    energy_range: float = 3.0
    ph: float = 7.4          # protonación de receptor y ligando al pasar a pdbqt
    cpu: int = 1
    workers: int = 0
    reuse: bool = False     # por defecto se recalcula todo: reutilizar poses de otra caja falsea el resultado
    rescoring_cnn: bool = False   # re-puntuar las poses con la red neuronal de gnina (segunda opinion)
    # Muestreo de ADCP: su equivalente a la exhaustividad de Vina, pero no se deduce de ella (son
    # algoritmos distintos). Un péptido largo necesita más; sin poder subirlo, la única salida ante
    # una pose mal convergida era repetir la corrida idéntica.
    adcp_pasos: int = 250_000
    adcp_replicas: int = 20


DERIVADOS = ("poses", "Complejos_Fusionados", "xml_plip", "san", "prep", "ligands", "xtal")
ARCHIVOS_DERIVADOS = ("resultados_docking.csv", "interacciones.csv", "ranking.csv",
                      "resumen.csv", "analogos.csv", "run.json")


def clean(out_dir) -> None:
    """Borra lo calculado, no lo aportado por el usuario (receptores y controles se conservan)."""
    out = Path(out_dir)
    for d in DERIVADOS:
        shutil.rmtree(out / d, ignore_errors=True)
    for f in ARCHIVOS_DERIVADOS:
        (out / f).unlink(missing_ok=True)


@dataclass
class RunResult:
    out_dir: Path
    analogs: Optional[pd.DataFrame] = None
    docking: Optional[pd.DataFrame] = None
    interactions: Optional[pd.DataFrame] = None
    ranking: Optional[pd.DataFrame] = None
    validation: Optional[pd.DataFrame] = None
    ref_info: dict = field(default_factory=dict)
    errores_docking: list = field(default_factory=list)


def _assign_controls(controls: Sequence[Path], receptors: Sequence[Path], given: dict) -> dict:
    """Empareja cada control con su receptor. Prioridad: lo indicado a mano > geometría > nombre.

    La geometría es el criterio fiable: un control cocristalizado comparte el sistema de
    coordenadas de SU receptor, así que sus átomos se solapan con los de esa estructura y quedan
    lejos de las demas. El nombre no sirve, porque el control suele llamarse por su ligando
    (control_ZI9) y el receptor por su PDB (8HTB), sin ninguna palabra en comun.
    """
    out = dict(given)
    stems = [Path(r).stem for r in receptors]
    for c in controls:
        ck = sc.normalize_key(Path(c).stem)
        if ck in out:
            continue
        if len(stems) == 1:
            out[ck] = stems[0]
            continue
        hit = _receptor_por_geometria(c, receptors)
        if hit is None:
            hit = next((Path(r).stem for r in receptors
                        if any(len(t) >= 4 and t.lower() in Path(c).stem.lower()
                               for t in re.split(r"[_\-.]", Path(r).stem))), None)
        if hit:
            out[ck] = hit
    return out


def _receptor_por_geometria(control, receptors) -> Optional[str]:
    """Receptor en cuyo espacio cae el control, o None. Compara el centroide del control con los
    átomos de cada receptor: pertenece al más cercano, y solo si está dentro (bajo un umbral), para
    no forzar la asignación de un control suelto."""
    try:
        cc = dk.coords_from_file(control)
    except Exception:
        cc = []
    if not cc:
        return None
    n = len(cc)
    centro = [sum(p[i] for p in cc) / n for i in range(3)]
    mejor, mejor_d = None, 1e18
    for r in receptors:
        try:
            pts = dk._coords(r)
        except Exception:
            continue
        if not pts:
            continue
        d = min((centro[0] - p[0]) ** 2 + (centro[1] - p[1]) ** 2 + (centro[2] - p[2]) ** 2
                for p in pts)
        if d < mejor_d:
            mejor, mejor_d = Path(r).stem, d
    # 8 Å al cuadrado: el centroide de un ligando cocristalizado está a pocos ángstrom de su
    # receptor; por encima, no pertenece a ninguno de los cargados.
    return mejor if mejor_d <= 64.0 else None


def _separar_peptidos(proj, ligandos) -> tuple:
    """({nombre: (secuencia, ciclico)} para ADCP, [rutas] para Vina).

    La procedencia se lee de ligands_meta.csv (lo escribe la interfaz al construir la biblioteca):
    sobrevive a reiniciar y distingue un péptido de cualquier ligando sin deducirlo de la estructura.
    Solo van a ADCP los que entran en su intervalo; por debajo del mínimo, Vina es practicable.
    """
    meta_p = Path(proj) / "ligands_meta.csv"
    if not (adcp.available() and meta_p.exists()):
        return {}, list(ligandos)
    try:
        meta = pd.read_csv(meta_p)
    except Exception:
        return {}, list(ligandos)
    if "fuente" not in meta.columns or "producto" not in meta.columns:
        return {}, list(ligandos)

    seqs = {}
    for _, r in meta[meta["fuente"].astype(str).str.startswith("pép")].iterrows():
        s = "".join(ch for ch in str(r.get("producto") or "").upper() if ch.isalpha())
        if adcp.MIN_RESIDUOS <= len(s) <= adcp.MAX_RESIDUOS:
            seqs[str(r["name"])] = (s, str(r["name"]).lower().startswith("ciclo"))

    a_adcp, a_vina = {}, []
    for l in ligandos:
        nombre = Path(l).stem
        if nombre in seqs:
            a_adcp[nombre] = seqs[nombre]
            continue
        # Los controles no figuran en la tabla de ligandos: su naturaleza se deduce de la estructura.
        # Un control peptídico debe ir también a ADCP; en Vina, la fila de referencia compararía
        # energías de dos funciones distintas frente al resto de la tabla.
        est = pep.secuencia_de_estructura(l)
        if est and adcp.MIN_RESIDUOS <= len(est[0]) <= adcp.MAX_RESIDUOS:
            a_adcp[nombre] = est
        else:
            a_vina.append(l)
    return a_adcp, a_vina


def run(cfg: RunConfig, on_step: Optional[Callable[[str, str], None]] = None) -> RunResult:
    """Ejecuta el ciclo completo. on_step(etapa, detalle) informa del avance."""
    def step(name, detail=""):
        if on_step:
            on_step(name, detail)

    out = Path(cfg.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if not cfg.reuse:
        clean(out)
        step("limpieza", "corrida desde cero; no se reutiliza nada previo")
    res = RunResult(out_dir=out)

    # 1. Diseño de análogos (opcional: si no hay lider, se usan los ligandos dados)
    lig_files = [Path(p) for p in cfg.ligands]
    if cfg.lead:
        step("diseno", "generando analogos y prediciendo ADMET")
        d = AdmelabBridge().design(cfg.lead, use_ml=cfg.use_ml,
                                   n_substitutions=list(cfg.n_substitutions),
                                   max_rows=cfg.n_analogs)
        res.analogs = d.to_dataframe()
        smis = d.smiles()[: cfg.n_analogs]
        step("diseno", f"{d.n_generated} generados, {len(smis)} seleccionados")
        step("3d", "generando confomeros")
        made = lg.materialize(smis, out / "ligands", seed=cfg.seed)
        if res.analogs is not None and not res.analogs.empty:
            names = {s: n for n, _, s in made}
            res.analogs.insert(0, "name", res.analogs["SMILES"].map(names) if "SMILES" in res.analogs else None)
            res.analogs.to_csv(out / "analogos.csv", index=False)
        lig_files += [p for _, p, _ in made]
        step("3d", f"{len(made)} estructuras 3D")

    dock_set = list(lig_files) + [Path(c) for c in cfg.controls]
    if not dock_set:
        raise ValueError("No hay nada que acoplar: da una molecula lider o ligandos.")
    if not cfg.receptors:
        raise ValueError("No hay receptores.")

    # 2. Caja. Prioridad: la indicada > el control cocristalizado de esa diana > automática. El
    # control marca el sitio real; el centro geométrico o un cofactor apuntarían a otro.
    control_keys = {sc.normalize_key(Path(c).stem) for c in cfg.controls}
    control_assign = _assign_controls(cfg.controls, cfg.receptors, cfg.control_map)
    boxes = dict(cfg.boxes)
    for r in cfg.receptors:
        if str(r) in boxes:
            step("caja", f"{Path(r).name}: indicada por el usuario")
            continue
        ctrl = next((c for c in cfg.controls
                     if control_assign.get(sc.normalize_key(Path(c).stem)) == Path(r).stem), None)
        if ctrl:
            boxes[str(r)] = dk.box_from_file(ctrl)
            step("caja", f"{Path(r).name}: centrada en el control {Path(ctrl).name}")
        else:
            boxes[str(r)] = dk.auto_box(r)
            step("caja", f"{Path(r).name}: automatica, sin control asignado")

    # 2b. Sitios. Por defecto uno por receptor. Con cfg.site_boxes se activa el híbrido: varios
    # bolsillos del mismo receptor, cada uno un sitio 'receptor~Etiqueta' que el scoring trata como
    # un receptor, separando el ranking por sitio.
    import re as _re
    def _sid(stem, label):
        return f"{stem}~{_re.sub(r'[^A-Za-z0-9]+', '', str(label))[:24] or 'sitio'}"
    def _dedup(pares):
        # Colapsa sitios de caja idéntica: evita acoplar dos veces el mismo bolsillo cuando ya es la
        # caja principal. Conserva el primero.
        vistos, out = set(), []
        for lab, bx in pares:
            k = tuple(round(v, 1) for v in (bx.cx, bx.cy, bx.cz, bx.sx, bx.sy, bx.sz))
            if k not in vistos:
                vistos.add(k); out.append((lab, bx))
        return out
    targets, stem_to_file, site_box, sites_of = [], {}, {}, {}
    for r in cfg.receptors:
        sb = cfg.site_boxes.get(str(r)) if cfg.site_boxes else None
        pares = _dedup(sb) if sb else []
        # Con un solo sitio, el id es el nombre del receptor SIN sufijo: el sufijo entra en el nombre
        # de las poses y la validación del redocking dejaría de encontrarlas.
        if len(pares) <= 1:
            lst = [(Path(r).stem, pares[0][1] if pares else boxes[str(r)])]
        else:
            lst = [(_sid(Path(r).stem, lab), bx) for lab, bx in pares]
        sites_of[str(r)] = [sid for sid, _ in lst]
        for sid, bx in lst:
            targets.append((r, sid, bx)); stem_to_file[sid] = str(r); site_box[sid] = bx
    hibrido = any(len(v) > 1 for v in sites_of.values())

    # 3. Docking
    step("docking", f"{len(dock_set)} ligandos x {len(targets)} sitio(s)"
                    + (" (hibrido)" if hibrido else ""))
    engine = dk.VinaEngine(cpu=cfg.cpu, seed=cfg.seed, exhaustiveness=cfg.exhaustiveness,
                           n_poses=cfg.n_poses, energy_range=cfg.energy_range)
    # Enrutado por motor, decidido por el ligando y no por el usuario: elegir mal aquí no da un
    # resultado peor sino uno sin sentido. Los péptidos van a ADCP (Vina no cubre su espacio
    # conformacional) sin configurar nada.
    peptidos, restantes = _separar_peptidos(out, dock_set)
    rows, errors = [], []
    if peptidos:
        step("docking", f"{len(peptidos)} peptido(s) -> ADCP; "
                        f"{len(restantes)} ligando(s) -> Vina")
        # Los hilos de ADCP son independientes de los de Vina: aquel paraleliza sus réplicas
        # internamente y sigue siendo reproducible con la semilla, mientras que en Vina más de un
        # hilo rompe el determinismo. Se le dan todos los núcleos para que la tanda no tarde de más.
        import os as _os
        _nuc = max(1, (_os.cpu_count() or 2) - 2)
        f_ad, e_ad = adcp.dock_sitios(targets, peptidos, out,
                                      {sid: r for r, sid, _b in targets},
                                      n_poses=cfg.n_poses, nucleos=_nuc, semilla=cfg.seed,
                                      n_pasos=cfg.adcp_pasos,
                                      n_repeticiones=max(cfg.adcp_replicas, cfg.n_poses))
        rows += f_ad
        errors += e_ad
    if restantes:
        f_vi, e_vi = dk.dock_batch(cfg.receptors, restantes, boxes, out, engine=engine,
                                   workers=cfg.workers, ph=cfg.ph, targets=targets)
        for f in f_vi:
            f.setdefault("motor", "vina")
        rows += f_vi
        errors += e_vi
    res.docking = pd.DataFrame(rows)
    if not res.docking.empty:
        res.docking.to_csv(out / "resultados_docking.csv", index=False)
    step("docking", f"{len(rows)} poses" + (f", {len(errors)} fallidos" if errors else ""))
    # Se detallan los fallos: un compuesto que no se acopla desaparece del ranking de ese sitio, y
    # sin aviso la ausencia pasa inadvertida y el cribado parece completo.
    if errors:
        res.errores_docking = [(b, e) for b, e in errors]
        for base, err in errors[:8]:
            step("docking: fallo", f"{base}: {err}")
    # Cobertura real por sitio: avisa si algun sitio quedo con menos compuestos que otro.
    if rows:
        import collections
        por_sitio = collections.Counter(r["receptor"] for r in rows)
        esperados = len(dock_set)
        for sid, _n in por_sitio.items():
            distintos = len({r["compound_name"] for r in rows if r["receptor"] == sid})
            if distintos < esperados:
                step("docking: aviso",
                     f"{sid}: {distintos} de {esperados} compuestos con poses")

    # 4. Complejos + interacciones
    step("complejos", "fusionando receptor y pose")
    ix.fuse_batch(cfg.receptors, out / "poses", out / "Complejos_Fusionados", cache_dir=out / "prep",
                  stem_to_file=stem_to_file)
    complexes = sorted((out / "Complejos_Fusionados").glob("*.pdb"))
    step("plip", f"{len(complexes)} complejos")
    res.interactions = ix.plip_batch(complexes, out, workers=cfg.workers, force=not cfg.reuse)

    # 5. Puntuación objetiva por pocket. Referencia = huella del ligando cristalográfico (no el
    # docking del control); pocket = residuos dentro de la caja.
    step("puntuacion", "referencia cristalografica + interacciones del pocket")
    base_crystal = ix.crystal_fingerprints(cfg.receptors, cfg.controls, control_assign,
                                           out / "xtal", cache_dir=out / "prep") if cfg.controls else {}
    # Referencia por sitio (clave en híbrido): el del control usa su huella cristalográfica; uno con
    # cofactor dentro de la caja usa la de ESE cofactor; uno sin ligando se marca para puntuar por
    # sus residuos catalíticos.
    crystal_feats, ref_src = {}, {}
    for r in cfg.receptors:
        stem = Path(r).stem
        bf = base_crystal.get(stem)
        ctrl = next((c for c in cfg.controls
                     if control_assign.get(sc.normalize_key(Path(c).stem)) == stem), None)
        cc = None
        if ctrl:
            try:
                pts = dk.coords_from_file(ctrl)
                cc = [sum(a) / len(a) for a in zip(*pts)]
            except Exception:
                cc = None
        for sid in sites_of[str(r)]:
            bx = site_box[sid]
            en_caja = cc is not None and (abs(cc[0] - bx.cx) <= bx.sx / 2
                                          and abs(cc[1] - bx.cy) <= bx.sy / 2
                                          and abs(cc[2] - bx.cz) <= bx.sz / 2)
            if bf and en_caja:
                crystal_feats[sid] = bf; ref_src[sid] = f"cristalografica ({Path(ctrl).stem})"
            else:
                feats, etq = ix.hetero_fingerprint(r, bx, out / "xtal")
                if feats:
                    crystal_feats[sid] = feats; ref_src[sid] = f"cofactor {etq}"
                else:
                    ref_src[sid] = "residuos del pocket"
    if crystal_feats:
        step("referencia", f"huella de ligando en {len(crystal_feats)} de {len(site_box)} sitio(s)")
    pocket_res = {}
    for sid, bx in site_box.items():
        pocket_res[sid] = sorted(dk.residues_in_box(stem_to_file[sid], bx))
    inter = sc.prepare_interactions(res.interactions, control_keys)
    dc = res.docking.copy()
    if not dc.empty:
        dc["ckey"] = dc["compound_name"].apply(sc.normalize_key)
    step("confianza", "estabilidad geometrica de poses (obrms)")
    pose_stab = sc.pose_stability_map(dc, out / "poses")

    # Segunda opinión opcional: gnina re-puntúa la mejor pose de cada compuesto en cada sitio. No
    # se re-acopla (el muestreo sigue siendo el de Vina); una pose por compuesto basta para el
    # consenso y evita que el coste crezca con el número de modos.
    cnn_map = {}
    if cfg.rescoring_cnn and not dc.empty and dk.gnina_available():
        objetivos, receptores_pose = [], {}
        for (sid, ck), sub in dc.groupby(["receptor", "ckey"]):
            mejor = sub.sort_values("docking_score").iloc[0]["pose_name"]
            p = out / "poses" / f"{mejor}.pdb"
            if p.exists():
                objetivos.append((ck, p))
                receptores_pose[ck] = stem_to_file.get(sid, cfg.receptors[0])
        step("rescoring", f"red neuronal (gnina) sobre {len(objetivos)} poses")
        for ck, p in objetivos:
            r = dk.rescore_poses(receptores_pose[ck], [p])
            for _nombre, vals in r.items():
                cnn_map[ck] = vals
        step("rescoring", f"{len(cnn_map)} compuestos re-puntuados")
    elif cfg.rescoring_cnn and not dk.gnina_available():
        step("rescoring", "gnina no esta disponible; se omite la segunda puntuacion")
    ref_info, icols, dscore = sc.build_ref_info(inter, dc, control_keys, control_assign,
                                                crystal_feats=crystal_feats)
    cat = {r: (cfg.catalytic.get(r) or cfg.catalytic.get(sc.base_of(r))
               or ref_info.get(r, {}).get("autocat", [])) for r in ref_info}
    res.ref_info = ref_info

    # SMILES para ADME: los análogos ya los traen del diseño; los ligandos de archivo se leen del disco.
    smiles_map = {}
    if res.analogs is not None and not res.analogs.empty and {"name", "SMILES"} <= set(res.analogs.columns):
        for _, r in res.analogs.iterrows():
            if r.get("name"):
                smiles_map[sc.normalize_key(r["name"])] = r["SMILES"]
    for carpeta in {Path(p).parent for p in list(lig_files) + [Path(c) for c in cfg.controls]}:
        for k, v in sc.build_smiles_map(carpeta).items():
            smiles_map.setdefault(k, v)

    res.ranking = sc.compute_ranking(inter, dc, control_keys, control_assign,
                                     ref_info, icols, dscore, cat, cfg.weights,
                                     smiles_map=smiles_map, pocket_res_map=pocket_res,
                                     pose_stability=pose_stab, cnn_map=cnn_map)

    # Se arrastran las columnas ADMET del diseño para no perder lo que admelab ya calculó.
    if res.analogs is not None and not res.analogs.empty and "name" in res.analogs.columns:
        adm = res.analogs.copy()
        adm["ckey"] = adm["name"].map(lambda n: sc.normalize_key(n) if n else None)
        cols_admet = [c for c in ("MW", "LogP", "TPSA", "QED", "Lipinski_violations",
                                  "LD50_mg_per_kg", "GHS_category") if c in adm.columns]
        if cols_admet:
            res.ranking = res.ranking.merge(adm[["ckey"] + cols_admet].drop_duplicates("ckey"),
                                            on="ckey", how="left")
    res.ranking.to_csv(out / "ranking.csv", index=False)

    # metadatos de la corrida: permiten reabrir y repuntuar los resultados sin repetir nada
    import json
    (out / "run.json").write_text(json.dumps({
        "lead": cfg.lead,
        "receptors": [str(r) for r in cfg.receptors],
        "controls": [str(c) for c in cfg.controls],
        "control_assign": control_assign,
        "control_keys": sorted(control_keys),
        "catalytic": cat,
        "reference": {r: ref_info.get(r, {}).get("src") for r in ref_info},
        "crystal_feats": {r: sorted(f) for r, f in crystal_feats.items()},
        "pocket_residues": pocket_res,
        "sites": {Path(r).name: sids for r, sids in sites_of.items()},
        "site_reference": ref_src,
        "hibrido": hibrido,
        "pose_stability": pose_stab,
        "weights": cfg.weights,
        "seed": cfg.seed,
        "exhaustiveness": cfg.exhaustiveness,
        "n_poses": cfg.n_poses,
        "energy_range": cfg.energy_range,
        "ph": cfg.ph,
        "rescoring_cnn": bool(cnn_map),
    }, indent=2))

    cols = [c for c in ("receptor", "compound", "best_dock", "pred_ki_M", "best_inter",
                        "cat_coverage", "efectividad_pct", "key_interaction", "tipo")
            if c in res.ranking.columns]
    res.ranking[cols].to_csv(out / "resumen.csv", index=False)

    # 6. Validación: el control debe recuperar su propia postura o el montaje no es fiable
    if cfg.controls:
        step("validacion", "redocking de los controles")
        res.validation = vl.redock_validation(cfg.controls, control_assign, out / "poses")
        res.validation.to_csv(out / "validacion_redocking.csv", index=False)
        step("validacion", vl.resumen(res.validation))

    step("listo", str(out))
    return res
