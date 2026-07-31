"""Reporte ADMET legible a partir de la salida de admelab.predict (estilo ADMETlab/ProTox).

Sin dependencias de UI: devuelve secciones estructuradas y una figura de radar. La interfaz decide
como pintarlas. Cada propiedad trae un veredicto (bueno/medio/malo/info) para colorear.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

# veredicto -> se traduce a color en la UI
GOOD, MID, BAD, INFO = "bueno", "medio", "malo", "info"


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _rng(v, lo, hi, margen=0.15):
    """Verde dentro de [lo,hi]; ambar en un margen; rojo fuera."""
    x = _num(v)
    if x is None:
        return INFO
    span = (hi - lo) * margen
    if lo <= x <= hi:
        return GOOD
    if lo - span <= x <= hi + span:
        return MID
    return BAD


def _max(v, t, margen=0.2):
    x = _num(v)
    if x is None:
        return INFO
    return GOOD if x <= t else (MID if x <= t * (1 + margen) else BAD)


def _min(v, t):
    x = _num(v)
    if x is None:
        return INFO
    return GOOD if x >= t else (MID if x >= t - 1 else BAD)


def _plow(v):   # probabilidad de un evento MALO (toxicidad, inhibicion): bajo = bueno
    x = _num(v)
    if x is None:
        return INFO
    return GOOD if x < 0.3 else (MID if x < 0.7 else BAD)


def _phigh(v):  # probabilidad de un evento BUENO (absorcion): alto = bueno
    x = _num(v)
    if x is None:
        return INFO
    return GOOD if x > 0.7 else (MID if x > 0.3 else BAD)


def _flag(v):   # 0 = sin alerta (bueno)
    x = _num(v)
    if x is None:
        return INFO
    return GOOD if x == 0 else BAD


def _pass(v):
    return GOOD if v in (True, "True", 1, 1.0) else (BAD if v in (False, "False", 0, 0.0) else INFO)


def _sa(v):     # SAscore: 1 fácil ... 10 difícil de sintetizar
    x = _num(v)
    if x is None:
        return INFO
    return GOOD if x <= 4 else (MID if x <= 6 else BAD)


AVISO_LD50 = (
    "**LD50 is a single-model estimate and tends to be optimistic.** Calibration against experimental rat values: acetaminophen 2274 vs 1944 mg/kg and ibuprofen 756 vs 636 (~20 % error), but **caffeine 738 vs 192 mg/kg: almost 4x less toxic than real**. The unit conversion is verified; the bias comes from the model. For chemotypes underrepresented in training (such as benzofuroxans) the discrepancy can reach two orders of magnitude vs ProTox. **Always contrast with a second predictor and do not publish «low toxicity» based on a single model.**"
)


def tox_index(row) -> Optional[float]:
    """Índice de toxicidad 0..1 (mayor = peor): promedio de las probabilidades de los eventos
    toxicos clave que prediga ADMET-AI. Sirve de eje 'tox' único en el ranking."""
    g = row.get if hasattr(row, "get") else (lambda k, d=None: None)
    vals = [_num(g(k)) for k in ("DILI", "AMES", "Carcinogens_Lagunin", "hERG")]
    vals = [v for v in vals if v is not None]
    return round(sum(vals) / len(vals), 3) if vals else None


def _pct(v):
    x = _num(v)
    return f"{x * 100:.0f}%" if x is not None else "-"


def _f(v, n=2):
    x = _num(v)
    return f"{x:.{n}f}" if x is not None else (str(v) if v not in (None, "") else "-")


def _fclamp(v, hi=100.0, n=1):
    """Formatea acotando a [0,hi]. Para magnitudes fisicamente acotadas (p.ej. % de unión a proteína)
    que el regresor ADMET-AI puede predecir por encima del límite. Marca el recorte con '~'."""
    x = _num(v)
    if x is None:
        return "-"
    return (f"~{hi:.{n}f}" if x > hi else (f"~0" if x < 0 else f"{x:.{n}f}"))


def sections(row: dict) -> list:
    """Devuelve [(titulo, [(propiedad, valor, veredicto), ...]), ...]."""
    g = row.get
    fis = [
        ("Molecular weight", _f(g("MW"), 1), _rng(g("MW"), 150, 500)),
        ("LogP", _f(g("LogP")), _rng(g("LogP"), 0, 5)),
        ("TPSA (A^2)", _f(g("TPSA"), 1), _rng(g("TPSA"), 20, 140)),
        ("H-bond donors", _f(g("HBD"), 0), _max(g("HBD"), 5)),
        ("H-bond acceptors", _f(g("HBA"), 0), _max(g("HBA"), 10)),
        ("Rotatable bonds", _f(g("RotB"), 0), _max(g("RotB"), 10)),
        ("Fraction sp3", _f(g("FractionCSP3")), INFO),
        ("QED (drug-likeness)", _f(g("QED")), _min(g("QED"), 0.5)),
        ("ESOL solubility (logS)", _f(g("ESOL_logS")), _min(g("ESOL_logS"), -4)),
        ("Synthesizability (SAscore 1-10)", _f(g("sa_score") if g("sa_score") is not None else g("SAscore")),
         _sa(g("sa_score") if g("sa_score") is not None else g("SAscore"))),
    ]
    reglas = [(nm, "pass" if _pass(g(k)) == GOOD else "no", _pass(g(k)))
              for nm, k in [("Lipinski", "Lipinski_pass"), ("Veber", "Veber_pass"),
                            ("Egan", "Egan_pass"), ("Ghose", "Ghose_pass"), ("Lead-like", "LeadLike_pass")]]
    alertas = [("PAINS", _f(g("PAINS_alert"), 0), _flag(g("PAINS_alert"))),
               ("BRENK", _f(g("BRENK_alert"), 0), _flag(g("BRENK_alert"))),
               ("Structural alerts", str(g("alerts") or "none"), _flag(g("n_alerts")))]
    absor = [
        ("Intestinal absorption (HIA)", _pct(g("HIA_Hou")), _phigh(g("HIA_Hou"))),
        ("Oral bioavailability", _pct(g("Bioavailability_Ma")), _phigh(g("Bioavailability_Ma"))),
        ("PAMPA permeability", _pct(g("PAMPA_NCATS")), _phigh(g("PAMPA_NCATS"))),
        ("Caco-2 (log cm/s)", _f(g("Caco2_Wang")), INFO),
        ("P-gp substrate", _pct(g("Pgp_Broccatelli")), _plow(g("Pgp_Broccatelli"))),
    ]
    distr = [
        ("Blood-brain barrier", _pct(g("BBB_Martins")), INFO),
        ("Plasma protein binding (%)", _fclamp(g("PPBR_AZ"), 100.0, 1), INFO),
        ("Volume of distribution (log)", _f(g("VDss_Lombardo")), INFO),
    ]
    metab = [(f"{cyp} inhibitor", _pct(g(k)), _plow(g(k))) for cyp, k in
             [("CYP1A2", "CYP1A2_Veith"), ("CYP2C9", "CYP2C9_Veith"), ("CYP2C19", "CYP2C19_Veith"),
              ("CYP2D6", "CYP2D6_Veith"), ("CYP3A4", "CYP3A4_Veith")]]
    excr = [("Hepatocyte clearance", _f(g("Clearance_Hepatocyte_AZ"), 1), INFO),
            ("Half-life (Obach)", _f(g("Half_Life_Obach")), INFO)]
    tox = [
        ("Mutagenicity (AMES)", _pct(g("AMES")), _plow(g("AMES"))),
        ("Cardiotoxicity (hERG)", _pct(g("hERG")), _plow(g("hERG"))),
        ("Hepatotoxicity (DILI)", _pct(g("DILI")), _plow(g("DILI"))),
        ("Carcinogenicity", _pct(g("Carcinogens_Lagunin")), _plow(g("Carcinogens_Lagunin"))),
        ("Clinical toxicity (ClinTox)", _pct(g("ClinTox")), _plow(g("ClinTox"))),
        ("Skin irritation", _pct(g("Skin_Reaction")), _plow(g("Skin_Reaction"))),
        ("Acute oral LD50 (mg/kg)*", _f(g("LD50_mg_per_kg"), 0), _min(g("LD50_mg_per_kg"), 500)),
        ("GHS category*", str(g("GHS_category") or "-"), INFO),
    ]
    return [("Physicochemical", fis), ("Druggability rules", reglas), ("Structural alerts", alertas),
            ("Absorption", absor), ("Distribution", distr), ("Metabolism (CYP450)", metab),
            ("Excretion", excr), ("Toxicity", tox)]


def _versions() -> dict:
    """Versiones del software para el bloque de Metodos (reproducibilidad exigible en publicación)."""
    import sys
    import subprocess
    v = {"python": sys.version.split()[0]}
    try:
        import rdkit
        v["rdkit"] = rdkit.__version__
    except Exception:
        pass
    for tool, args in (("vina", ["--version"]), ("obabel", ["-V"])):
        try:
            r = subprocess.run([tool] + args, capture_output=True, text=True, timeout=15)
            line = next((l for l in ((r.stdout or "") + (r.stderr or "")).splitlines() if l.strip()), "")
            if line:
                v[tool] = line.strip()[:90]
        except Exception:
            pass
    # PLIP y fpocket no exponen --versión limpio: PLIP se lee de los metadatos del paquete y fpocket
    # se reporta como presente. Volcar su texto de uso en la sección de Metodos sería ruido.
    # La versión se consulta con importlib.metadata, que lee los metadatos SIN importar el modulo.
    # Es deliberado: PLIP se distribuye bajo GPL y PoliScreen lo ejecuta como proceso aparte, de
    # modo que no incorpora código suyo; importarlo aquí, aunque fuera solo para leer un número,
    # enturbiaria esa separacion.
    import shutil as _sh
    try:
        import importlib.metadata as _md
        v["plip"] = _md.version("plip")
    except Exception:
        if _sh.which("plip"):
            v["plip"] = "instalado"
    if _sh.which("fpocket"):
        v["fpocket"] = "instalado"
    return v


def methods_text(meta: dict, weights: Optional[dict] = None, catalytic: Optional[dict] = None,
                 secondary: Optional[dict] = None) -> str:
    """Bloque de Metodos en Markdown a partir del run.json: parámetros, caja, pesos, referencia y
    versiones. Pensado para pegar en el paper y garantizar reproducibilidad."""
    meta = meta or {}
    weights = weights or meta.get("weights", {})
    catalytic = catalytic or meta.get("catalytic", {})
    secondary = secondary or meta.get("secondary", {})
    def _l(x):
        return ", ".join(Path(str(p)).name for p in x) if isinstance(x, (list, tuple)) else str(x)
    L = []
    L.append("## Methods (in silico) — generated by PoliScreen\n")
    L.append(f"- Receptors: {_l(meta.get('receptors', []))}")
    L.append(f"- Controls (co-crystallized): {_l(meta.get('controls', []))}")
    ref = meta.get("reference", {})
    L.append(f"- Interaction reference: {', '.join(f'{k}={v}' for k, v in ref.items()) or 'n/a'} "
             "(crystallographic = PLIP fingerprint of the co-crystallized ligand in its real pose)")
    L.append(f"- Seed: {meta.get('seed')} · CPU per docking: 1 (deterministic)")
    L.append(f"- Vina: exhaustiveness={meta.get('exhaustiveness')}, num_modes={meta.get('n_poses')}, "
             f"energy_range={meta.get('energy_range')} kcal/mol")
    L.append(f"- Protonation (OpenBabel): pH {meta.get('ph')}")
    pr = meta.get("pocket_residues", {})
    if pr:
        L.append("- Pocket residues (inside the box): "
                 + "; ".join(f"{k}: {len(v)}" for k, v in pr.items()))
    if catalytic:
        L.append("- Catalytic residues (gate, mandatory): "
                 + "; ".join(f"{k}: {', '.join(v)}" for k, v in catalytic.items() if v))
    if secondary:
        L.append("- Secondary residues (bonus, not mandatory): "
                 + "; ".join(f"{k}: {', '.join(v)}" for k, v in secondary.items() if v))
    if weights:
        ejes = {k: weights.get(k) for k in ("dock", "inter", "adme", "ki", "tox")}
        L.append("- Axis weights (auto-normalized): "
                 + ", ".join(f"{k}={v}" for k, v in ejes.items()))
        L.append("- Affinity-axis metric: "
                 + ("ligand efficiency LE = -dG/heavy atoms (corrects size bias)"
                    if str(weights.get("dock_metric", "dock")).lower() == "le" else "raw docking score"))
        L.append("- Ki: shown but NOT scored (deterministic transform of the affinity; scoring it would count docking twice)")
        L.append(f"- Objective-model knobs: w_cat={weights.get('w_cat')}, "
                 f"w_sec={weights.get('w_sec')}, w_out={weights.get('w_out')}, cat_gate={weights.get('cat_gate')}")
        tw = weights.get("type_weights")
        if tw:
            L.append("- Weights by interaction type: " + ", ".join(f"{k}={v}" for k, v in tw.items()))
    L.append("- Scoring: objective per-pocket model; quality = sum of interactions weighted by type and residue role, normalized by the crystallographic ligand (=100%); PLIP for interactions.")
    if meta.get("rescoring_cnn"):
        L.append("- Second scoring function: the poses generated by Vina were re-evaluated with gnina's convolutional neural network, without repeating the conformational search. The agreement between both scores enters as consensus evidence in the confidence metric.")
    L.append("- Confidence (0-1): geometric mean of binding-mode convergence (Tanimoto between the PLIP fingerprints of the best poses) and affinity-interaction agreement; reduced if the target fails its redocking. The geometric spread of the poses (obrms) is reported as a diagnostic but does not enter the number: with AutoDock Vina it is nearly constant and does not discriminate. Orthogonal to effectiveness: it quantifies the reliability of the result, not its magnitude.")
    L.append("\n### Versions")
    for k, val in _versions().items():
        L.append(f"- {k}: {val}")
    L.append("\n*Limitations: rigid docking, non-covalent; Ki estimated from Vina's affinity (dG approximation), informative, not measured.*")
    return "\n".join(L)


def radar_fig(row: dict, title: str = "", fig=None):
    """Radar de drogabilidad: cada eje 0..1 donde 1 = dentro del rango optimo."""
    import matplotlib.pyplot as plt
    import numpy as np

    def score_rng(v, lo, hi):
        x = _num(v)
        if x is None:
            return 0.0
        if lo <= x <= hi:
            return 1.0
        span = (hi - lo)
        d = (lo - x) if x < lo else (x - hi)
        return max(0.0, 1.0 - d / span)

    ejes = [
        ("MW", score_rng(row.get("MW"), 150, 500)),
        ("LogP", score_rng(row.get("LogP"), 0, 5)),
        ("TPSA", score_rng(row.get("TPSA"), 20, 140)),
        ("HBD", 1.0 - min(1.0, (_num(row.get("HBD")) or 0) / 8)),
        ("HBA", 1.0 - min(1.0, (_num(row.get("HBA")) or 0) / 14)),
        ("RotB", 1.0 - min(1.0, (_num(row.get("RotB")) or 0) / 12)),
        ("Sol.", score_rng(row.get("ESOL_logS"), -5, 0)),
        ("QED", min(1.0, (_num(row.get("QED")) or 0))),
    ]
    labels = [e[0] for e in ejes]
    vals = [e[1] for e in ejes]
    ang = np.linspace(0, 2 * np.pi, len(labels), endpoint=False).tolist()
    vals += vals[:1]; ang += ang[:1]
    if fig is None:
        fig = plt.figure(figsize=(4.6, 4.6))
    ax = fig.add_subplot(111, polar=True)
    ax.plot(ang, vals, color="#1b9e77", lw=2)
    ax.fill(ang, vals, color="#1b9e77", alpha=0.25)
    ax.set_xticks(ang[:-1]); ax.set_xticklabels(labels, fontsize=9)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0]); ax.set_yticklabels(["", "", "", ""])
    ax.set_ylim(0, 1)
    if title:
        ax.set_title(title, fontsize=10, pad=14)
    return fig
