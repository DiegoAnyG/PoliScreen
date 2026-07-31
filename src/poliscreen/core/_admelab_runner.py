"""Runner que se ejecuta DENTRO del entorno de admelab (venv Python 3.12 con torch).

Deliberadamente NO importa poliscreen: solo admelab + stdlib + pandas. Lee sus parámetros
de un JSON y escribe el resultado en otro JSON. Así el diseño de análogos + ADMET vive
aislado del entorno de docking (openbabel/plip/vina), que es incompatible, y ese mismo
límite se convierte manana en un contenedor separado SIN cambiar la interfaz.

Uso:  python _admelab_runner.py parámetros.json salida.json
"""
import json
import math
import sys


def _clean(v):
    """NaN/Inf -> None para que el JSON sea valido y portable."""
    if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
        return None
    return v


def _do_design(p):
    from admelab import pipeline
    kw = dict(
        methods=tuple(p.get("methods") or ["decoration"]),
        use_ml=bool(p.get("use_ml", True)),
        scope=p.get("scope", "aromatic_ch"),
        positions=p.get("positions"),
        n_substitutions=tuple(p.get("n_substitutions") or [1]),
        max_decor=int(p.get("max_decor", 300)),
        max_brics=int(p.get("max_brics", 60)),
        include_lead=bool(p.get("include_lead", True)),
    )
    if p.get("substituents"):
        kw["substituents"] = p["substituents"]
    if p.get("filters"):
        kw["filters"] = p["filters"]
    r = pipeline.run_pipeline(p["lead_smiles"], **kw)
    df = r.ranked
    if p.get("max_rows"):
        df = df.head(int(p["max_rows"]))
    df = df.astype(object).where(df.notna(), None)
    return {
        "ok": True,
        "action": "design",
        "n_generated": int(r.n_generated),
        "n_scored": int(r.n_scored),
        "columns": [str(c) for c in df.columns],
        "rows": [{str(k): _clean(v) for k, v in rec.items()} for rec in df.to_dict(orient="records")],
    }


def _do_reaction_sites(p):
    """Sitios reactivos de una molécula: OH clasificados y si tiene acido carboxilico."""
    from rdkit import Chem
    from admelab import esterification as est
    m = Chem.MolFromSmiles(p["smiles"])
    if m is None:
        return {"ok": False, "error": "SMILES no valido: %s" % p["smiles"]}
    return {"ok": True, "action": "reaction_sites",
            "sites": est.classify_hydroxyls(m),
            "has_cooh": bool(m.HasSubstructMatch(est._COOH_PATTERN)),
            "viability_rules": dict(est.FISCHER_VIABILITY)}


def _do_esterify(p):
    """Esterifica el acido con cada alcohol. policy: 'preferred' (OH más favorable) o 'all'."""
    from admelab import esterification as est
    out = []
    for a in p["alcohols"]:
        smi = a["smiles"] if isinstance(a, dict) else a
        nombre = a.get("name") if isinstance(a, dict) else None
        for pr in est.esterify(p["acid"], smi, policy=p.get("policy", "preferred")):
            r = dict(pr)
            r["alcohol"] = nombre or smi
            r["alcohol_smiles"] = smi
            # Se normaliza aquí para que el resto del código no dependa de como nombre admelab
            # sus claves. 'sintetizable' es lo que usa el filtro antes de dockear.
            v = pr.get("viabilidad_fischer") or pr.get("viabilidad") or "desconocida"
            r["viabilidad"] = v
            r["sintetizable"] = not any(w in str(v).lower() for w in ("desfavorable", "dificil"))
            out.append({k: _clean(v2) for k, v2 in r.items()})
    return {"ok": True, "action": "esterify", "products": out}


def _do_predict(p):
    """ADMET completo (descriptores RDKit + reglas + ADMET-AI + toxicidad) para una lista de SMILES."""
    from admelab import predict, toxicity
    df = predict.predict_batch(list(p["smiles"]), use_ml=bool(p.get("use_ml", True)))
    try:
        df = toxicity.annotate_toxicity(df)
    except Exception:
        pass
    df = df.astype(object).where(df.notna(), None)
    return {"ok": True, "action": "predict", "columns": [str(c) for c in df.columns],
            "rows": [{str(k): _clean(v) for k, v in r.items()} for r in df.to_dict(orient="records")]}


def _do_name_esters(p):
    """Nombre IUPAC verificado (OPSIN) de cada ester. Offline con los nombres de alcohol dados;
    cae a web solo si use_web y OPSIN no cuadro. Limpia los nombres para que OPSIN los parsee."""
    import re
    from admelab import naming_smart as ns

    def clean(n):
        n = re.sub(r"\s*\(.*?\)", "", str(n or "")).strip()   # quita parentesis
        n = re.sub(r"^\(.*?\)-?", "", n).strip()               # quita estereo al inicio
        return n or None

    esters = list(p["ester_smiles"])
    alcs = list(p["alcohol_smiles"])
    anames = [clean(a) for a in (p.get("alcohol_names") or [None] * len(esters))]
    df = ns.name_esters_batch(esters, alcs, acid_smiles=p.get("acid_smiles"),
                              alcohol_names=anames, use_web=bool(p.get("use_web", True)), progress=False)
    out = []
    for _, r in df.iterrows():
        nm = r["iupac_name"]
        out.append({"smiles": r["SMILES"], "iupac_name": (None if (nm is None or str(nm) == "nan") else str(nm)),
                    "verified": bool(r["iupac_verified"])})
    return {"ok": True, "action": "name_esters", "names": out, "opsin": ns.opsin_available()}


def _do_info(p):
    import admelab  # noqa: F401
    mods = []
    for m in ["generation", "esterification", "predict", "toxicity", "ranking",
              "naming", "naming_smart", "viz", "pipeline"]:
        try:
            __import__("admelab." + m)
            mods.append(m)
        except Exception:
            pass
    info = {"ok": True, "action": "info", "modules": mods, "python": sys.version.split()[0]}
    try:
        import torch
        info["torch"] = torch.__version__
        info["cuda"] = bool(torch.cuda.is_available())
    except Exception:
        info["torch"] = None
        info["cuda"] = False
    return info


def main():
    if len(sys.argv) < 3:
        print("uso: _admelab_runner.py parametros.json salida.json", file=sys.stderr)
        return 2
    params_path, out_path = sys.argv[1], sys.argv[2]
    try:
        p = json.load(open(params_path))
        action = p.get("action", "design")
        if action == "design":
            out = _do_design(p)
        elif action == "reaction_sites":
            out = _do_reaction_sites(p)
        elif action == "esterify":
            out = _do_esterify(p)
        elif action == "name_esters":
            out = _do_name_esters(p)
        elif action == "predict":
            out = _do_predict(p)
        elif action == "info":
            out = _do_info(p)
        else:
            out = {"ok": False, "error": "accion desconocida: %s" % action}
    except Exception as e:
        import traceback
        out = {"ok": False, "error": str(e), "traceback": traceback.format_exc()[-2000:]}
    with open(out_path, "w") as f:
        json.dump(out, f)
    return 0 if out.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
