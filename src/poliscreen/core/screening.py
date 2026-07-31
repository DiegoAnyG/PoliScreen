"""Cribado: interacciones (PLIP) y puntuación de calidad frente al control.

Sin UI. Motor único de interacciones: PLIP, con numeración de autor.
Puntuación: Tversky ponderada contra la huella del control; los contactos de más no suman.
Sin dianas ni residuos cableados: todo se deriva del control y de lo que elija el usuario.
"""
import os, re, json, subprocess, shutil
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np, pandas as pd

REC_EXT = (".pdb", ".pdbqt")
LIG_EXT = (".mol2", ".sdf", ".smi", ".mol", ".pdbqt")
# patrones de archivos GENERADOS (no son entradas): complejos, poses por modelo, intermedios de prep/PLIP
GEN_PAT = re.compile(r"(^Complejo_|_compounds_a_|-model\d|_protH|_protonly|^fx_|^_cof_|_conv\.|_rec\.|_only\.)", re.I)
PLIP_MAP = {"hbond": ["hydrogen_bond"], "hydrophobic": ["hydrophobic_interaction"], "saltbridge": ["salt_bridge"],
            "water": ["water_bridge"], "halogen": ["halogen_bond"], "pistack": ["pi_stack"],
            "pication": ["pi_cation_interaction"], "metal": ["metal_complex"]}
# color, etiqueta, estilo de línea por tipo (para el diagrama 2D)
TYPE_STYLE = {"hbond": ("#1f77b4", "H-bond", "-"), "hydrophobic": ("#7f7f7f", "Hidrofobica", "--"),
              "saltbridge": ("#d62728", "Puente salino", "-"), "pistack": ("#2ca02c", "Pi-stacking", "-"),
              "pication": ("#9467bd", "Pi-cation", "-"), "halogen": ("#17becf", "Halogeno", "-"),
              "water": ("#8c564b", "Puente de agua", ":"), "metal": ("#e377c2", "Metal", "-")}

# ---------------------------------------------------------------- parsers
# Interacciones específicas y orientadas. Añadir una donde el control no la tiene es mérito real;
# acumular hidrofóbicas o puentes de agua es promiscuidad y no cuenta.
DIRECCIONALES = {"hbond", "saltbridge", "pistack", "pication", "halogen", "metal"}

# Mérito relativo por tipo de interacción (heurístico, guiado por literatura, editable por el
# usuario en weights["type_weights"]). Iónico y H-bond, fuertes y específicos; pi/halógeno,
# intermedios; hidrofóbica, débil pero numerosa; agua, incierta. No son energías: son prioridades.
TYPE_WEIGHTS = {"saltbridge": 1.00, "metal": 0.95, "hbond": 0.85, "pication": 0.80,
                "pistack": 0.65, "halogen": 0.55, "hydrophobic": 0.35, "water": 0.30}


def feat_type(feat):
    """Tipo de interacción de una columna 'Residuo#_tipo' (p.ej. 'Tyr157_hbond' -> 'hbond')."""
    return str(feat).rsplit("_", 1)[-1]


def normalize_key(v): return re.sub(r"[^a-z0-9]+", "", str(v).lower())
def base_of(receptor):
    """Receptor base de un id de sitio. En docking hibrido el id es 'receptor~PocketN'; así el control
    y los cataliticos (asignados al receptor) se reconocen en todos sus sitios. Sin sufijo, identidad."""
    return str(receptor).split("~", 1)[0]
def receptor_from_name(n):
    m = re.search(r"docking_(.+?)_compounds_a_", str(n)); return m.group(1) if m else "receptor"
def compound_from_pose_name(v):
    t = str(v).strip()
    m = re.search(r"(?i)compounds?[_\-]?a[_\-](.+?)[_\-]model\d+\s*$", t)
    if m: return m.group(1).strip(" :_-")
    m = re.search(r"(?i)([^\\/]+?)[_\-:]?model\d+\s*$", t)
    return m.group(1).strip(" :_-") if m else t
def model_of(n):
    m = re.search(r"model(\d+)", str(n)); return int(m.group(1)) if m else 1
def pose_key(n):
    """Clave comun receptor|compuesto|modelo, tolerante a prefijos (Complejo_) y sufijos, para
    casar filas de interacciones con filas de docking aunque el nombre no sea idéntico."""
    return f"{normalize_key(receptor_from_name(n))}|{normalize_key(compound_from_pose_name(n))}|{model_of(n)}"
def resname(feat): return feat.split("_")[0]
def resnum(res):
    d = "".join(ch for ch in str(res) if ch.isdigit()); return int(d) if d else 0

# ---------------------------------------------------------------- escaneo de entradas
def scan_inputs(data_dir):
    """Receptores (.pdb) y ligandos (.mol2/.sdf/.smi/.mol) de la carpeta de datos,
    EXCLUYENDO archivos generados (complejos, poses, intermedios). Evita el problema de los 2000+."""
    root = Path(data_dir); recs, ligs = [], []
    if not root.exists(): return recs, ligs
    for p in sorted(root.rglob("*"), key=lambda x: str(x).lower()):
        if not p.is_file() or GEN_PAT.search(p.name): continue
        suf = p.suffix.lower()
        if suf == ".pdb": recs.append(p)
        elif suf in (".mol2", ".sdf", ".smi", ".mol"): ligs.append(p)
    return recs, ligs

# ---------------------------------------------------------------- carga de resultados existentes
def load_results(work_dir):
    """Lee interacciones.csv, resultados_docking.csv y selección.json de la carpeta de trabajo."""
    W = Path(work_dir)
    inter = pd.read_csv(W / "interacciones.csv") if (W / "interacciones.csv").exists() else pd.DataFrame()
    dock = pd.read_csv(W / "resultados_docking.csv") if (W / "resultados_docking.csv").exists() else pd.DataFrame()
    sel = json.loads((W / "seleccion.json").read_text()) if (W / "seleccion.json").exists() else {}
    return inter, dock, sel

# ---------------------------------------------------------------- ADME / Ki
def adme_score_from_smiles(smi):
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski, Crippen, rdMolDescriptors, QED
    except Exception:
        return (np.nan, np.nan)
    if not smi: return (np.nan, np.nan)
    m = Chem.MolFromSmiles(smi)
    if m is None: return (np.nan, np.nan)
    mw, logp, tpsa = Descriptors.MolWt(m), Crippen.MolLogP(m), rdMolDescriptors.CalcTPSA(m)
    hbd, hba = Lipinski.NumHDonors(m), Lipinski.NumHAcceptors(m)
    viol = int(mw > 500) + int(logp > 5) + int(hbd > 5) + int(hba > 10)
    def itv(x, lo, hi, sl, sh):
        if lo <= x <= hi: return 1.
        if sl < x < lo: return max(0, (x - sl) / (lo - sl))
        if hi < x < sh: return max(0, (sh - x) / (sh - hi))
        return 0.
    comp = {"MW": itv(mw, 150, 500, 100, 650), "LogP": itv(logp, 0, 5, -1, 7), "TPSA": itv(tpsa, 20, 140, 0, 180),
            "HBD": itv(hbd, 0, 5, 0, 8), "HBA": itv(hba, 0, 10, 0, 14), "QED": QED.qed(m)}
    wts = {"MW": .18, "LogP": .22, "TPSA": .18, "HBD": .1, "HBA": .1, "QED": .22}
    return (100 * sum(comp[k] * wts[k] for k in wts), viol)

def _sascore(m):
    """Synthetic Accessibility score (Ertl & Schuffenhauer 2009): 1 fácil - 10 difícil.
    Usa el sascorer que RDKit trae en Contrib; si no esta, devuelve NaN."""
    try:
        import os, sys
        from rdkit.Chem import RDConfig
        sap = os.path.join(RDConfig.RDContribDir, "SA_Score")
        if sap not in sys.path:
            sys.path.append(sap)
        import sascorer
        return round(float(sascorer.calculateScore(m)), 2)
    except Exception:
        return np.nan


def extra_props_from_smiles(smi):
    """Sintetizabilidad (SAscore), alerta PAINS y booleanos de Lipinski. Todo offline con RDKit.
    Booleanos por regla: 1 cumple, 0 incumple. pains: 1 = tiene subestructura promiscua (malo)."""
    out = dict(sa_score=np.nan, pains=np.nan, ro5_mw=np.nan, ro5_logp=np.nan,
               ro5_hbd=np.nan, ro5_hba=np.nan, ro5_pass=np.nan, logp=np.nan, n_heavy=np.nan)
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski, Crippen
    except Exception:
        return out
    if not smi:
        return out
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return out
    mw, logp = Descriptors.MolWt(m), Crippen.MolLogP(m)
    hbd, hba = Lipinski.NumHDonors(m), Lipinski.NumHAcceptors(m)
    out["logp"] = round(float(logp), 2); out["n_heavy"] = int(m.GetNumHeavyAtoms())
    out["ro5_mw"] = int(mw <= 500); out["ro5_logp"] = int(logp <= 5)
    out["ro5_hbd"] = int(hbd <= 5); out["ro5_hba"] = int(hba <= 10)
    out["ro5_pass"] = int((out["ro5_mw"] + out["ro5_logp"] + out["ro5_hbd"] + out["ro5_hba"]) >= 3)  # Ro5: <=1 violacion
    try:
        from rdkit.Chem import FilterCatalog
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        out["pains"] = int(FilterCatalog.FilterCatalog(params).HasMatch(m))
    except Exception:
        pass
    out["sa_score"] = _sascore(m)
    return out


def ki_from_dg(dg, temp_k=298.15):
    RT = 0.0019872 * temp_k
    return np.exp(pd.to_numeric(dg, errors="coerce") / RT)


def _fp_convergence(sub, icols, dscore, topk: int = 5):
    """Convergencia del modo de unión: Tanimoto medio entre las huellas PLIP de las MEJORES poses de un
    compuesto. Alto = el docking encuentra el mismo modo una y otra vez (senal de confianza). NaN si <2 poses."""
    rows = sub.copy()
    rows["_sc"] = rows["name"].map(lambda n: dscore.get(pose_key(n), np.nan))
    rows = rows.sort_values("_sc").head(topk)
    fps = [frozenset(c for c in icols if r[c] > 0) for _, r in rows.iterrows()]
    fps = [f for f in fps if f]
    if len(fps) < 2:
        return np.nan
    sims = []
    for i in range(len(fps)):
        for j in range(i + 1, len(fps)):
            u = fps[i] | fps[j]
            sims.append(len(fps[i] & fps[j]) / len(u) if u else 0.0)
    return round(float(np.mean(sims)), 3) if sims else np.nan


def _geomean(vals):
    """Media geométrica de las componentes disponibles (ignora NaN). Exige que TODAS concuerden: una
    componente baja hunde el resultado, a diferencia de la media aritmetica. Suelo 1e-3 para evitar log(0)."""
    xs = [float(v) for v in vals if v is not None and not pd.isna(v) and float(v) >= 0]
    if not xs:
        return np.nan
    return round(float(np.exp(np.mean([np.log(max(v, 1e-3)) for v in xs]))), 3)


def fp_recovery(ref_feats, pose_feats) -> dict:
    """Coincidencia entre dos huellas de interacción (conjuntos de features Res#_tipo).
    recovery = |intersección|/|referencia| (cuanto del cristalográfico reproduce la pose dockeada);
    tanimoto = |intersección|/|unión|. Es la Validación de que el docking recupera el modo real."""
    ref, pose = set(ref_feats or []), set(pose_feats or [])
    if not ref:
        return dict(recovery=np.nan, tanimoto=np.nan, shared=0, ref_n=0, extra=len(pose))
    shared = ref & pose
    union = ref | pose
    return dict(recovery=round(len(shared) / len(ref), 3),
                tanimoto=round(len(shared) / len(union), 3) if union else np.nan,
                shared=len(shared), ref_n=len(ref), extra=len(pose - ref))

# ---------------------------------------------------------------- referencia (huella del control)
def build_ref_info(inter, dc, control_keys, control_assign, crystal_feats=None):
    """Por receptor: huella de REFERENCIA. Si se pasa crystal_feats[R] (huella PLIP del ligando
    cristalográfico en su pose real), esa es la referencia objetiva. Si no, se cae a la mejor pose
    DOCKEADA del control (compatibilidad). Devuelve también autocat: los residuos con los que la
    referencia hace interacciones DIRECCIONALES, sugerencia de residuos clave cuando el usuario no
    los conoce."""
    crystal_feats = crystal_feats or {}
    icols = [c for c in inter.columns if "_" in c and c not in ("name", "compound", "ckey", "receptor", "is_control")]
    dscore = {}
    if not dc.empty:
        for _, r in dc.iterrows(): dscore[pose_key(r["pose_name"])] = r["docking_score"]
    def best_row(sub):
        if sub is None or sub.empty: return None
        sc = sub["name"].map(lambda n: dscore.get(pose_key(n), np.nan))
        return sub.loc[sc.idxmin()] if sc.notna().any() else sub.iloc[0]
    info = {}
    for R in sorted(inter["receptor"].unique()):
        rin = inter[inter["receptor"] == R]
        cks = {ck for ck, rc in control_assign.items() if rc == base_of(R)} or set(rin[rin["is_control"]]["ckey"].unique())
        refck0 = None
        if cks:
            if (not dc.empty) and dc["ckey"].isin(cks).any():
                refck0 = dc[dc["ckey"].isin(cks)].sort_values("docking_score").iloc[0]["ckey"]
            else:
                refck0 = sorted(cks)[0]
        cfeats = crystal_feats.get(R)
        if cfeats:                                   # referencia cristalográfica del sitio (ligando/cofactor)
            feats = sorted(set(cfeats))
            freq = {c: 1.0 for c in feats}
            src = "cristalografica"
        elif crystal_feats:                          # modo sitio: este bolsillo no tiene ligando de referencia
            info[R] = dict(ckey=refck0, feats=[], freq={}, residues=[], autocat=[], src="residuos/relativa")
            continue
        else:                                        # compat (sin huellas): control dockeado
            cpre = rin[rin["ckey"].isin(cks)]
            if cpre.empty:
                info[R] = dict(ckey=refck0, feats=[], freq={}, residues=[], autocat=[], src="ninguna"); continue
            cref = cpre[cpre["ckey"] == refck0]; best = best_row(cref)
            feats = [c for c in icols if best is not None and best[c] > 0]
            freq = {c: float((cref[c] > 0).mean()) for c in feats}
            src = "control dockeado"
        residues = sorted({resname(c) for c in feats}, key=lambda r: (resnum(r), r))
        autocat = sorted({resname(c) for c in feats if feat_type(c) in DIRECCIONALES},
                         key=lambda r: (resnum(r), r))
        info[R] = dict(ckey=refck0, feats=feats, freq=freq, residues=residues, autocat=autocat, src=src)
    return info, icols, dscore

# ---------------------------------------------------------------- ranking (Tversky calidad)
def compute_ranking(inter, dc, control_keys, control_assign, ref_info, icols, dscore,
                    cat_map, weights, smiles_map=None, protox_map=None, temp_k=298.15,
                    pocket_res_map=None, sec_map=None, pose_stability=None, reliable_map=None,
                    cnn_map=None):
    """Ranking OBJETIVO por pocket. La calidad de interacción no mide PARECIDO al control: suma el
    valor de cada interacción (peso por TIPO) según el ROL del residuo -catalitico (gate, x w_cat),
    secundario (ancla conocida, x w_sec), de pocket (x1) o externo (x w_out)- y se normaliza por la
    calidad del ligando Cristalográfico (=100%). Así un compuesto supera al control haciendo Más/mejores
    contactos productivos, no copiandolo. Los residuos los designa el usuario (cat_map / sec_map); si no,
    los cataliticos se sugieren desde la referencia (ref_info)."""
    smiles_map = smiles_map or {}; protox_map = protox_map or {}; pocket_res_map = pocket_res_map or {}
    sec_map = sec_map or {}
    TW = dict(TYPE_WEIGHTS); TW.update(weights.get("type_weights") or {})
    W_CAT = float(weights.get("w_cat", weights.get("key", 3.0)))    # catalitico (gate): necesario y premiado
    W_SEC = float(weights.get("w_sec", 1.5))                         # secundario: ancla conocida, no obligatoria
    W_OUT = float(weights.get("w_out", 0.15))                        # contacto fuera del pocket: poco merito
    CAT_GATE = float(weights.get("cat_gate", 0.5))                   # cuanto penaliza faltar a un catalitico (0..1)
    def best_row(sub):
        if sub is None or sub.empty: return None
        sc = sub["name"].map(lambda n: dscore.get(pose_key(n), np.nan))
        return sub.loc[sc.idxmin()] if sc.notna().any() else sub.iloc[0]
    _adme, _extra = {}, {}
    def adme(k):
        if k not in _adme: _adme[k] = adme_score_from_smiles(smiles_map.get(k))
        return _adme[k]
    def extra(k):
        if k not in _extra: _extra[k] = extra_props_from_smiles(smiles_map.get(k))
        return _extra[k]
    bloques = []
    for R in sorted(inter["receptor"].unique()):
        rin = inter[inter["receptor"] == R].copy()
        info = ref_info.get(R, {}); feats = info.get("feats", [])
        cats = {s.lower() for s in cat_map.get(R, [])}
        secs = {s.lower() for s in sec_map.get(R, [])} - cats   # el gate manda sobre lo secundario
        # residuos del pocket: por geometría de la caja si se pasan; si no, los que toca cualquiera
        pocket = {s.lower() for s in pocket_res_map.get(R, [])} or {resname(c).lower() for c in icols}
        def role_mult(res):
            r = res.lower()
            if r in cats: return W_CAT
            if r in secs: return W_SEC
            return 1.0 if r in pocket else W_OUT
        def quality(fs):
            return float(sum(TW.get(feat_type(c), 0.3) * role_mult(resname(c)) for c in fs))
        ctrl_cks = {ck for ck, rc in control_assign.items() if rc == base_of(R)} or set(rin[rin["is_control"]]["ckey"].unique())
        Q0 = quality(feats)   # referencia cristalográfica (ligando/cofactor del sitio): la línea del 100%
        # Sitio sin ligando de referencia (bolsillo secundario del híbrido): el 100% se define por
        # sus residuos catalíticos (una interacción direccional con cada uno). Sin catalíticos en el
        # pocket, se normaliza al mejor compuesto del sitio.
        cats_en_pocket = cats & pocket
        if Q0 <= 0 and cats_en_pocket:
            Q0 = float(len(cats_en_pocket)) * TW.get("hbond", 0.85) * W_CAT
        # Primera pasada: huella y calidad de cada compuesto (el control usa su huella cristalográfica
        # solo si el sitio la tiene; en un sitio sin referencia usa su pose dockeada, como los demas).
        comp = []
        for ck, sub in rin.groupby("ckey"):
            r0 = best_row(sub)
            fs = set(feats) if (ck in ctrl_cks and feats) else {c for c in icols if r0 is not None and r0[c] > 0}
            comp.append((ck, sub, r0, fs, quality(fs)))
        if Q0 <= 0:
            qmax = max((q for *_, q in comp), default=0.0)
            Q0 = qmax if qmax > 0 else 1.0     # normalización relativa: el mejor del sitio = 100%
        recs = []
        for ck, sub, r0, fs, Q in comp:
            comp_res = {resname(c).lower() for c in fs}
            cov = (sum(1 for cr in cats if cr in comp_res) / len(cats)) if cats else np.nan
            gate = 1.0 if not cats else ((1.0 - CAT_GATE) + CAT_GATE * cov)   # faltar a un catalitico penaliza
            T = (Q / Q0 * gate) if Q0 > 0 else np.nan
            destacar = sorted([c for c in fs if resname(c).lower() in cats or resname(c).lower() in pocket],
                              key=lambda c: TW.get(feat_type(c), 0.0) * role_mult(resname(c)), reverse=True)
            recs.append(dict(ckey=ck, compound=sub["compound"].iloc[0], is_control=int(sub["is_control"].max()),
                             inter_quality=T, cat_coverage=cov, key_interaction="; ".join(destacar[:12]),
                             pose=(model_of(r0["name"]) if r0 is not None else None),
                             conv=_fp_convergence(sub, icols, dscore)))   # convergencia del modo de unión
        isum = pd.DataFrame(recs)
        dR = dc[dc["receptor"] == R].groupby("ckey").agg(best_dock=("docking_score", "min")).reset_index() if not dc.empty else pd.DataFrame(columns=["ckey", "best_dock"])
        mR = isum.merge(dR, on="ckey", how="outer"); mR["receptor"] = R; mR["compound"] = mR["compound"].fillna(mR["ckey"])
        mR["is_control"] = mR["is_control"].fillna(0).astype(int)
        mR["pred_ki_M"] = ki_from_dg(mR["best_dock"], temp_k)
        mR["ld50_mgkg"] = mR["ckey"].map(lambda k: protox_map.get(k, {}).get("ld50", np.nan))
        mR["tox_class"] = mR["ckey"].map(lambda k: protox_map.get(k, {}).get("tox_class", np.nan))
        adv = mR["ckey"].apply(lambda k: pd.Series(adme(k), index=["adme", "lip_viol"])); mR = pd.concat([mR, adv], axis=1)
        ext = mR["ckey"].apply(lambda k: pd.Series(extra(k))); mR = pd.concat([mR, ext], axis=1)
        # pKi (numérico, sortable, evita el lio de unidades), y eficiencia: LE = -dG/átomos pesados,
        # LLE = pKi - LogP. Premian unión por átomo, no molécula grande: guarda anti-"greaseball".
        bd = pd.to_numeric(mR["best_dock"], errors="coerce"); nh = pd.to_numeric(mR["n_heavy"], errors="coerce")
        mR["pKi"] = (-np.log10(pd.to_numeric(mR["pred_ki_M"], errors="coerce"))).round(2)
        mR["LE"] = (-bd / nh.where(nh > 0)).round(3)
        mR["LLE"] = (mR["pKi"] - pd.to_numeric(mR["logp"], errors="coerce")).round(2)
        mR["es_control_del_target"] = mR.apply(lambda r: (r["is_control"] == 1) and (control_assign.get(r["ckey"]) == base_of(R)), axis=1)
        refck = info.get("ckey")
        def refval(col, mode):
            v = pd.to_numeric(mR[mR["ckey"] == refck][col], errors="coerce").dropna() if refck else pd.Series(dtype=float)
            if v.empty: v = pd.to_numeric(mR[mR["es_control_del_target"]][col], errors="coerce").dropna()
            if v.empty: v = pd.to_numeric(mR[col], errors="coerce").dropna()
            return np.nan if v.empty else (v.min() if mode == "low" else v.max())
        ref_dock = refval("best_dock", "low"); ref_adme = refval("adme", "high"); ref_ki = refval("pred_ki_M", "low"); ref_tox = refval("ld50_mgkg", "high")
        def eff(col, ref, invert=False):
            v = pd.to_numeric(mR[col], errors="coerce")
            if ref is None or pd.isna(ref) or ref == 0: return pd.Series(np.nan, index=mR.index)
            return (ref / v) if invert else (v / ref)
        # Eje de afinidad: 'dock' (score crudo, sesgado hacia moléculas grandes) o 'le' (eficiencia
        # de ligando, -dG/átomos pesados, corrige el sesgo). Ambos normalizados contra el control.
        if str(weights.get("dock_metric", "dock")).lower() == "le":
            dock_axis = eff("LE", refval("LE", "high"))
        else:
            dock_axis = eff("best_dock", ref_dock)
        axes = {"dock": dock_axis, "sim": pd.to_numeric(mR["inter_quality"], errors="coerce"),
                "adme": eff("adme", ref_adme), "ki": eff("pred_ki_M", ref_ki, invert=True), "tox": eff("ld50_mgkg", ref_tox)}
        ws = pd.Series(0., index=mR.index); wa = pd.Series(0., index=mR.index)
        for key, wt in [("dock", weights.get("dock", 0.0)), ("sim", weights.get("inter", 0.0)),
                        ("adme", weights.get("adme", 0.0)), ("ki", weights.get("ki", 0.0)), ("tox", weights.get("tox", 0.0))]:
            if wt <= 0: continue
            v = pd.to_numeric(axes[key], errors="coerce"); ok = v.notna(); ws[ok] += v[ok] * wt; wa[ok] += wt
        mR["efectividad_pct"] = np.where(wa > 0, 100 * ws / wa, np.nan)
        # Percentil dentro de la biblioteca (por receptor): comparable entre dianas, al contrario que
        # el % vs control, que varía con lo fuerte que sea el control de cada diana.
        mR["percentil"] = (pd.to_numeric(mR["efectividad_pct"], errors="coerce").rank(pct=True) * 100).round(0)
        mR["best_inter"] = pd.to_numeric(mR["inter_quality"], errors="coerce").round(3)
        mR["cat_coverage"] = pd.to_numeric(mR["cat_coverage"], errors="coerce").round(3)
        # --- Confianza (consenso multi-evidencia, 0-1). Ortogonal a la efectividad: no mide cuán
        # bueno es el compuesto sino cuánto fiarse del número. Se reduce a la mitad si la diana no
        # valida (su control no re-dockea).
        aff_pct = (-pd.to_numeric(mR["best_dock"], errors="coerce")).rank(pct=True)
        int_pct = pd.to_numeric(mR["inter_quality"], errors="coerce").rank(pct=True)
        mR["conc"] = (1.0 - (aff_pct - int_pct).abs()).round(3)
        ps = pose_stability or {}
        mR["geom"] = pd.to_numeric(mR["ckey"].map(lambda k: ps.get(k, np.nan)), errors="coerce")
        # Confianza = media geométrica de conv (reproducibilidad del modo de unión vía huella PLIP) y
        # conc (concordancia afinidad-interacción). geom (dispersión geométrica) se reporta aparte
        # pero NO entra: con Vina es casi constante (~0.08, num_modes devuelve poses diversas) y no
        # varía con la caja ni la exhaustividad, así que hundía la confianza de todos por igual.
        # Tercera evidencia opcional (si se re-puntuó con gnina): consenso entre dos funciones
        # independientes, la empírica de Vina y la red de gnina. Que dos métodos con supuestos
        # distintos ordenen igual un compuesto es una señal que ninguno da por sí solo.
        cm = cnn_map or {}
        componentes = ["conv", "conc"]
        if cm:
            mR["cnn_score"] = pd.to_numeric(mR["ckey"].map(lambda k: cm.get(k, {}).get("cnn_score")),
                                            errors="coerce")
            mR["cnn_affinity"] = pd.to_numeric(mR["ckey"].map(lambda k: cm.get(k, {}).get("cnn_affinity")),
                                               errors="coerce")
            if mR["cnn_affinity"].notna().any():
                cnn_pct = mR["cnn_affinity"].rank(pct=True)
                mR["consenso"] = (1.0 - (aff_pct - cnn_pct).abs()).round(3)
                componentes.append("consenso")
        mR["confidence"] = mR.apply(lambda r: _geomean([r.get(c) for c in componentes]), axis=1)
        if reliable_map and not reliable_map.get(R, True):
            mR["confidence"] = (pd.to_numeric(mR["confidence"], errors="coerce") * 0.5).round(3)
        mR = mR[(mR["is_control"] == 0) | (mR["es_control_del_target"])].copy()
        bloques.append(mR)
    rk = pd.concat(bloques, ignore_index=True)
    rk["tipo"] = np.where(rk["is_control"] == 1, "Control",
                 np.where(rk["efectividad_pct"] >= 105, "Supera control",
                 np.where(rk["efectividad_pct"] >= 95, "Comparable", "Inferior")))
    return rk.sort_values(["receptor", "efectividad_pct"], ascending=[True, False])

def prepare_interactions(inter, control_keys):
    """Anade columnas compound/ckey/receptor/is_control a interacciones.csv."""
    inter = inter.copy()
    inter["compound"] = inter["name"].apply(compound_from_pose_name)
    inter["ckey"] = inter["compound"].apply(normalize_key)
    inter["receptor"] = inter["name"].apply(receptor_from_name)
    inter["is_control"] = inter["ckey"].isin(control_keys)
    return inter

# ---------------------------------------------------------------- diagrama 2D (PLIP, esquema radial)
def draw_2d(row, title, fig=None, reference=None, figsize=(5.0, 5.0)):
    """Esquema radial de las interacciones de una pose (fila de interacciones.csv).

    Si se pasa `reference` (huella del control: conjunto de features residuo_tipo), colorea en VERDE lo
    que la reproduce y en GRIS lo que no cuenta (contacto de más o mismo residuo con otro enlace). El
    estilo de línea sigue codificando el tipo de interacción. Mismo motor y numeración que la tabla.
    """
    import matplotlib.pyplot as plt
    feats = [(c, int(row[c])) for c in row.index
             if "_" in c and c.rsplit("_", 1)[-1] in TYPE_STYLE and pd.notna(row[c]) and row[c] > 0]
    fig, ax = plt.subplots(figsize=figsize) if fig is None else (fig, fig.subplots())
    ax.set_aspect("equal"); ax.axis("off")
    if not feats:
        ax.text(0, 0, "Sin interacciones PLIP en esta pose", ha="center", va="center"); return fig
    ref = set(reference) if reference else None
    MATCH, EXTRA = "#2ca02c", "#9aa0a6"
    byres = {}
    for feat, n in feats:
        r, t = feat[:feat.rfind("_")], feat[feat.rfind("_") + 1:]
        byres.setdefault(r, []).append((t, n, feat))
    res_list = list(byres.keys()); m = len(res_list)
    ang = (np.linspace(90, 90 + 360, m, endpoint=False)) * np.pi / 180.0
    ax.scatter([0], [0], s=1700, c="#ffe08a", edgecolors="#b8860b", linewidths=1.5, zorder=5)
    ax.text(0, 0, "LIG", ha="center", va="center", fontweight="bold", zorder=6)
    used = set()
    for a, res in zip(ang, res_list):
        x, y = np.cos(a), np.sin(a); types = byres[res]
        res_match = ref is not None and any(f in ref for _, _, f in types)
        for j, (typ, n, feat) in enumerate(types):
            _, _, ls = TYPE_STYLE.get(typ, ("#333333", typ, "-"))
            col = TYPE_STYLE.get(typ, ("#333333",))[0] if ref is None else (MATCH if feat in ref else EXTRA)
            off = (j - (len(types) - 1) / 2.0) * 0.06
            ax.plot([-np.sin(a) * off, x * 0.80 - np.sin(a) * off], [np.cos(a) * off, y * 0.80 + np.cos(a) * off],
                    color=col, ls=ls, lw=2.4, zorder=2); used.add(typ)
        node_c = "#cfe8ff" if ref is None else ("#c8f0cf" if res_match else "#e8e8e8")
        edge_c = "#33618f" if ref is None else (MATCH if res_match else "#8a8f94")
        ax.scatter([x], [y], s=1050, c=node_c, edgecolors=edge_c, linewidths=1.3, zorder=4)
        ax.text(x, y, res, ha="center", va="center", fontsize=9, zorder=5)
        ntot = sum(n for _, n, _ in types)
        if ntot > 1: ax.text(x * 1.26, y * 1.26, f"x{ntot}", ha="center", va="center", fontsize=8, color="dimgray")
    ax.set_xlim(-1.55, 1.55); ax.set_ylim(-1.55, 1.55)
    import matplotlib.pyplot as _plt
    if ref is None:
        handles = [_plt.Line2D([0], [0], color=TYPE_STYLE[t][0], ls=TYPE_STYLE[t][2], lw=2.2, label=TYPE_STYLE[t][1])
                   for t in TYPE_STYLE if t in used]
    else:
        handles = [_plt.Line2D([0], [0], color=MATCH, lw=2.4, label="Reproduce el control"),
                   _plt.Line2D([0], [0], color=EXTRA, lw=2.4, label="Extra / otro enlace")]
        handles += [_plt.Line2D([0], [0], color="#555", ls=TYPE_STYLE[t][2], lw=1.6, label=TYPE_STYLE[t][1])
                    for t in TYPE_STYLE if t in used]
    if handles: ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(-0.02, 1.03), fontsize=8, frameon=False)
    ax.set_title(title, fontsize=11)
    return fig

# ---------------------------------------------------------------- redocking RMSD (obrms)
def redocking_rmsd(ref_file, pose_file):
    try:
        r = subprocess.run(["obrms", str(ref_file), str(pose_file)], capture_output=True, text=True)
        nums = re.findall(r"[-+]?\d+\.\d+", r.stdout)
        if nums: return float(nums[-1])
    except Exception:
        pass
    return np.nan


def pose_stability_map(dc, poses_dir, soft: float = 4.0, topk: int = 20) -> dict:
    """{ckey: estabilidad geométrica de las poses} (obrms vs la mejor pose). Componente de la confianza.

    Crédito GRADUADO por pose: 1 si coincide con la mejor y decae linealmente a 0 en `soft` A, en vez
    de un corte binario a 2 A que castigaba de más. Usa TODAS las poses generadas (hasta topk), no solo
    unas pocas, para aprovechar el muestreo que pidio el usuario. Se precalcula (obrms es costoso).
    Un valor bajo es informativo: significa que el docking no converge en un modo de unión (ligando
    promiscuo, caja demasiado grande o busqueda poco exhaustiva)."""
    out = {}
    if dc is None or getattr(dc, "empty", True):
        return out
    pdir = Path(poses_dir)
    for ck, sub in dc.groupby("ckey"):
        s = sub.sort_values("docking_score").head(topk)
        files = [pdir / f"{n}.pdb" for n in s["pose_name"]]
        files = [f for f in files if f.exists()]
        if len(files) < 2:
            out[ck] = np.nan
            continue
        rmsds = [redocking_rmsd(files[0], f) for f in files[1:]]
        rmsds = [r for r in rmsds if not np.isnan(r)]
        if not rmsds:
            out[ck] = np.nan
            continue
        cred = [max(0.0, 1.0 - r / soft) for r in rmsds]
        out[ck] = round(sum(cred) / len(cred), 3)
    return out

# ---------------------------------------------------------------- SMILES para ADME
def build_smiles_map(data_dir):
    """SMILES por compuesto: primero de csv (name,smiles), luego de los propios archivos de ligando."""
    smiles = {}
    root = Path(data_dir)
    if not root.exists(): return smiles
    for csvf in root.rglob("*.csv"):
        try: t = pd.read_csv(csvf)
        except Exception: continue
        cn = [c for c in t.columns if c.lower() in ("name", "nombre", "compound", "compound_name")]
        cs = [c for c in t.columns if c.lower() in ("smiles", "smile")]
        if cn and cs:
            for _, r in t.iterrows():
                k = normalize_key(r[cn[0]])
                if k and k not in smiles: smiles[k] = str(r[cs[0]])
    try:
        from rdkit import Chem
        for lp in root.rglob("*"):
            if lp.suffix.lower() in (".mol2", ".sdf", ".mol", ".smi"):
                k = normalize_key(lp.stem)
                if not k or k in smiles: continue
                try:
                    if lp.suffix.lower() == ".smi":
                        smiles[k] = lp.read_text().split()[0].strip()
                        continue
                    if lp.suffix.lower() == ".mol2":
                        m = Chem.MolFromMol2File(str(lp))
                    elif lp.suffix.lower() == ".sdf":
                        m = next(iter(Chem.SDMolSupplier(str(lp))), None)
                    else:
                        m = Chem.MolFromMolFile(str(lp))
                    s = Chem.MolToSmiles(Chem.RemoveHs(m)) if m is not None else None
                    # RDKit no siempre percibe bien los heterociclos con carga formal (el N-oxido
                    # del benzofuroxano hace que no pueda kekulizar). OpenBabel si los interpreta,
                    # así que se usa como respaldo antes de dar el compuesto por ilegible.
                    if not s:
                        s = smiles_via_obabel(lp)
                    if s:
                        smiles[k] = s
                except Exception:
                    s = smiles_via_obabel(lp)
                    if s:
                        smiles[k] = s
    except Exception:
        pass
    return smiles


def smiles_via_obabel(archivo):
    """SMILES leido con OpenBabel y validado con RDKit. Devuelve None si no sale nada utilizable."""
    try:
        r = subprocess.run(["obabel", str(archivo), "-osmi"], capture_output=True, text=True, timeout=60)
        bruto = (r.stdout or "").strip().split()
        if not bruto:
            return None
        from rdkit import Chem
        m = Chem.MolFromSmiles(bruto[0])
        return Chem.MolToSmiles(m) if m is not None else None
    except Exception:
        return None

def find_data_dir(work_dir):
    """Busca la carpeta de entrada (la que tiene más .mol2 / *clean*.pdb) bajo work/entrada o work."""
    W = Path(work_dir); cands = []
    for base in [W / "entrada", W]:
        if not base.exists(): continue
        dirs = [base] + [p for p in base.rglob("*") if p.is_dir()]
        for d in dirs:
            n = len(list(d.glob("*.mol2"))) + len(list(d.glob("*clean*.pdb")))
            if n: cands.append((n, str(d)))
    return max(cands)[1] if cands else ""

def controls_and_assign(sel):
    """control_keys y control->receptor desde selección.json."""
    return set(sel.get("control_keys", [])), dict(sel.get("control_target_map", {}))

# ki a 0 a propósito: deriva del score de docking, así que puntuarla sería contarlo dos veces. Se
# muestra en las tablas pero no entra en el ranking.
# Ejes (dock/inter/adme/ki/tox): media ponderada auto-normalizada por su suma, así que no tienen que
# sumar 1 y ponerlos todos a 1.0 = promedio simple. Perillas del modelo: w_cat (peso del catalítico),
# w_out (contacto fuera del pocket), cat_gate (penalización por faltar a un catalítico).
DEFAULT_WEIGHTS = dict(dock=0.50, inter=0.50, adme=0.0, ki=0.0, tox=0.0,
                       dock_metric="dock",   # 'dock' = score crudo | 'le' = eficiencia de ligando (corrige sesgo de tamaño)
                       w_cat=3.0, w_sec=1.5, w_out=0.15, cat_gate=0.5,
                       # compat con el modelo anterior (no intervienen en el modelo objetivo):
                       key=3.0, alpha=0.3, beta=1.0, bonus=0.25, type_rigor=1.0)
