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
    "**La LD50 es una estimación de un solo modelo y tiende a ser optimista.** Calibración frente a "
    "valores experimentales en rata: paracetamol 2274 vs 1944 mg/kg y ibuprofeno 756 vs 636 (error ~20 %), "
    "pero **cafeína 738 vs 192 mg/kg: casi 4 veces menos tóxica de lo real**. La conversión de unidades está "
    "verificada; el sesgo viene del modelo. Para quimiotipos poco representados en el entrenamiento (como los "
    "benzofuroxanos) la discrepancia puede llegar a dos órdenes de magnitud frente a ProTox. "
    "**Contrasta siempre con un segundo predictor y no publiques «baja toxicidad» apoyándote en un solo modelo.**"
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
        ("Peso molecular", _f(g("MW"), 1), _rng(g("MW"), 150, 500)),
        ("LogP", _f(g("LogP")), _rng(g("LogP"), 0, 5)),
        ("TPSA (A^2)", _f(g("TPSA"), 1), _rng(g("TPSA"), 20, 140)),
        ("Donadores de H", _f(g("HBD"), 0), _max(g("HBD"), 5)),
        ("Aceptores de H", _f(g("HBA"), 0), _max(g("HBA"), 10)),
        ("Enlaces rotables", _f(g("RotB"), 0), _max(g("RotB"), 10)),
        ("Fraccion sp3", _f(g("FractionCSP3")), INFO),
        ("QED (calidad de farmaco)", _f(g("QED")), _min(g("QED"), 0.5)),
        ("Solubilidad ESOL (logS)", _f(g("ESOL_logS")), _min(g("ESOL_logS"), -4)),
        ("Sintetizabilidad (SAscore 1-10)", _f(g("sa_score") if g("sa_score") is not None else g("SAscore")),
         _sa(g("sa_score") if g("sa_score") is not None else g("SAscore"))),
    ]
    reglas = [(nm, "cumple" if _pass(g(k)) == GOOD else "no", _pass(g(k)))
              for nm, k in [("Lipinski", "Lipinski_pass"), ("Veber", "Veber_pass"),
                            ("Egan", "Egan_pass"), ("Ghose", "Ghose_pass"), ("Lead-like", "LeadLike_pass")]]
    alertas = [("PAINS", _f(g("PAINS_alert"), 0), _flag(g("PAINS_alert"))),
               ("BRENK", _f(g("BRENK_alert"), 0), _flag(g("BRENK_alert"))),
               ("Alertas estructurales", str(g("alerts") or "ninguna"), _flag(g("n_alerts")))]
    absor = [
        ("Absorcion intestinal (HIA)", _pct(g("HIA_Hou")), _phigh(g("HIA_Hou"))),
        ("Biodisponibilidad oral", _pct(g("Bioavailability_Ma")), _phigh(g("Bioavailability_Ma"))),
        ("Permeabilidad PAMPA", _pct(g("PAMPA_NCATS")), _phigh(g("PAMPA_NCATS"))),
        ("Caco-2 (log cm/s)", _f(g("Caco2_Wang")), INFO),
        ("Sustrato de P-gp", _pct(g("Pgp_Broccatelli")), _plow(g("Pgp_Broccatelli"))),
    ]
    distr = [
        ("Barrera hematoencefalica", _pct(g("BBB_Martins")), INFO),
        ("Union a proteinas plasma (%)", _fclamp(g("PPBR_AZ"), 100.0, 1), INFO),
        ("Volumen de distribucion (log)", _f(g("VDss_Lombardo")), INFO),
    ]
    metab = [(f"Inhibidor {cyp}", _pct(g(k)), _plow(g(k))) for cyp, k in
             [("CYP1A2", "CYP1A2_Veith"), ("CYP2C9", "CYP2C9_Veith"), ("CYP2C19", "CYP2C19_Veith"),
              ("CYP2D6", "CYP2D6_Veith"), ("CYP3A4", "CYP3A4_Veith")]]
    excr = [("Aclaramiento hepatocito", _f(g("Clearance_Hepatocyte_AZ"), 1), INFO),
            ("Vida media (Obach)", _f(g("Half_Life_Obach")), INFO)]
    tox = [
        ("Mutagenicidad (AMES)", _pct(g("AMES")), _plow(g("AMES"))),
        ("Cardiotoxicidad (hERG)", _pct(g("hERG")), _plow(g("hERG"))),
        ("Hepatotoxicidad (DILI)", _pct(g("DILI")), _plow(g("DILI"))),
        ("Carcinogenicidad", _pct(g("Carcinogens_Lagunin")), _plow(g("Carcinogens_Lagunin"))),
        ("Toxicidad clinica (ClinTox)", _pct(g("ClinTox")), _plow(g("ClinTox"))),
        ("Irritacion de piel", _pct(g("Skin_Reaction")), _plow(g("Skin_Reaction"))),
        ("LD50 oral aguda (mg/kg)*", _f(g("LD50_mg_per_kg"), 0), _min(g("LD50_mg_per_kg"), 500)),
        ("Categoria GHS*", str(g("GHS_category") or "-"), INFO),
    ]
    return [("Fisicoquimica", fis), ("Reglas de drogabilidad", reglas), ("Alertas estructurales", alertas),
            ("Absorcion", absor), ("Distribucion", distr), ("Metabolismo (CYP450)", metab),
            ("Excrecion", excr), ("Toxicidad", tox)]


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
    L.append("## Metodos (in silico) — generado por PoliScreen\n")
    L.append(f"- Receptores: {_l(meta.get('receptors', []))}")
    L.append(f"- Controles (co-cristalizados): {_l(meta.get('controls', []))}")
    ref = meta.get("reference", {})
    L.append(f"- Referencia de interacciones: {', '.join(f'{k}={v}' for k, v in ref.items()) or 'n/d'} "
             "(cristalografica = huella PLIP del ligando co-cristalizado en su pose real)")
    L.append(f"- Semilla: {meta.get('seed')} · CPU por acoplamiento: 1 (determinista)")
    L.append(f"- Vina: exhaustiveness={meta.get('exhaustiveness')}, num_modes={meta.get('n_poses')}, "
             f"energy_range={meta.get('energy_range')} kcal/mol")
    L.append(f"- Protonacion (OpenBabel): pH {meta.get('ph')}")
    pr = meta.get("pocket_residues", {})
    if pr:
        L.append("- Residuos del pocket (dentro de la caja): "
                 + "; ".join(f"{k}: {len(v)}" for k, v in pr.items()))
    if catalytic:
        L.append("- Residuos cataliticos (gate, obligatorios): "
                 + "; ".join(f"{k}: {', '.join(v)}" for k, v in catalytic.items() if v))
    if secondary:
        L.append("- Residuos secundarios (bonificacion, no obligatorios): "
                 + "; ".join(f"{k}: {', '.join(v)}" for k, v in secondary.items() if v))
    if weights:
        ejes = {k: weights.get(k) for k in ("dock", "inter", "adme", "ki", "tox")}
        L.append("- Pesos de eje (auto-normalizados): "
                 + ", ".join(f"{k}={v}" for k, v in ejes.items()))
        L.append("- Metrica del eje de afinidad: "
                 + ("eficiencia de ligando LE = -dG/atomos pesados (corrige sesgo de tamano)"
                    if str(weights.get("dock_metric", "dock")).lower() == "le" else "score de docking crudo"))
        L.append("- Ki: se muestra pero NO se puntua (es transformacion determinista de la afinidad; "
                 "puntuarla contaria el docking dos veces)")
        L.append(f"- Perillas del modelo objetivo: w_cat={weights.get('w_cat')}, "
                 f"w_sec={weights.get('w_sec')}, w_out={weights.get('w_out')}, cat_gate={weights.get('cat_gate')}")
        tw = weights.get("type_weights")
        if tw:
            L.append("- Pesos por tipo de interaccion: " + ", ".join(f"{k}={v}" for k, v in tw.items()))
    L.append("- Puntuacion: modelo objetivo por pocket; calidad = suma de interacciones ponderadas por tipo "
             "y por rol del residuo, normalizada por el ligando cristalografico (=100%); PLIP para interacciones.")
    if meta.get("rescoring_cnn"):
        L.append("- Segunda funcion de puntuacion: las poses generadas por Vina se re-evaluaron con la red "
                 "neuronal convolucional de gnina, sin repetir la busqueda conformacional. El acuerdo entre "
                 "ambas puntuaciones se incorpora como evidencia de consenso en la metrica de confianza.")
    L.append("- Confianza (0-1): media geometrica de la convergencia del modo de union (Tanimoto entre las "
             "huellas PLIP de las mejores poses) y la concordancia afinidad-interaccion; reducida si la diana no "
             "valida su redocking. La dispersion geometrica de las poses (obrms) se reporta como diagnostico pero "
             "no entra en el numero: con AutoDock Vina es casi constante y no discrimina. Ortogonal a la "
             "efectividad: cuantifica la fiabilidad del resultado, no su magnitud.")
    L.append("\n### Versiones")
    for k, val in _versions().items():
        L.append(f"- {k}: {val}")
    L.append("\n*Limitaciones: docking rigido, no covalente; Ki estimada desde la afinidad de Vina "
             "(aproximacion de dG), informativa, no medida.*")
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
