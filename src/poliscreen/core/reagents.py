"""Gestion de reactivos (p. ej. alcoholes para esterificar).

Reproducibilidad primero: la fuente fiable es un archivo que el usuario controla (biblioteca
interna o su propia carga). PubChem es un complemento de descubrimiento, best-effort: si no
responde, el cribado continua con lo que ya hay. Todo se deduplica por InChIKey y conserva su
procedencia, para poder distinguir en la tabla lo que aporto el usuario.
"""
from __future__ import annotations

import csv
import json
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

from . import reactions as rx

# procedencia, de mayor a menor prioridad al deduplicar (lo del usuario manda)
FUENTES = ("tuyo", "interno", "pubchem")
PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


@dataclass
class Reagent:
    name: str
    smiles: str
    inchikey: str
    source: str


def _inchikey(smiles: str) -> Optional[str]:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles)
    return Chem.MolToInchiKey(m) if m is not None else None


def inchikey(smiles: str) -> Optional[str]:
    """InChIKey de un SMILES; None si no es valido. Clave de deduplicacion."""
    return _inchikey(smiles)


def _mk(name, smiles, source) -> Optional[Reagent]:
    ik = _inchikey(smiles)
    if not ik:
        return None
    return Reagent(str(name or smiles), smiles, ik, source)


def matches_group(reagent: Reagent, reaction: rx.Reaction) -> bool:
    """El reactivo tiene el grupo que la reacción necesita (p. ej. un OH para esterificar)."""
    from rdkit import Chem
    m = Chem.MolFromSmiles(reagent.smiles)
    patt = Chem.MolFromSmarts(reaction.partner_smarts)
    return m is not None and patt is not None and m.HasSubstructMatch(patt)


# ------------------------------------------------------------------ fuentes
def load_internal(reaction: rx.Reaction) -> list:
    out = []
    for row in rx.load_library(reaction):
        r = _mk(row.get("name"), row.get("smiles"), "interno")
        if r:
            out.append(r)
    return out


def load_user_files(paths: Sequence) -> list:
    """Alcoholes del usuario desde csv (name,smiles), .smi, .sdf, .mol o .mol2. Acepta una carpeta."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    archivos = []
    for p in paths:
        p = Path(p)
        if p.is_dir():
            archivos += [q for q in p.rglob("*") if q.suffix.lower() in (".csv", ".smi", ".sdf", ".mol", ".mol2")]
        elif p.exists():
            archivos.append(p)
    out = []
    for f in archivos:
        suf = f.suffix.lower()
        try:
            if suf in (".csv", ".xlsx", ".xls"):
                import pandas as pd
                t = pd.read_csv(f) if suf == ".csv" else pd.read_excel(f)
                cols = {str(c).lower().strip(): c for c in t.columns}
                # se aceptan cabeceras en espanol/ingles y las del inventario bfx_esters_v2
                ncol = next((cols[k] for k in ("name", "nombre", "compound", "compuesto",
                                               "alcohol origen", "nombre clave") if k in cols), None)
                scol = next((cols[k] for k in ("smiles", "smile", "smiles alcohol") if k in cols), None)
                if scol is None:
                    continue
                for _, rr in t.iterrows():
                    s = rr[scol]
                    if pd.isna(s):
                        continue
                    r = _mk(rr[ncol] if ncol is not None and pd.notna(rr[ncol]) else s, str(s), "tuyo")
                    if r:
                        out.append(r)
            elif suf == ".smi":
                for line in f.read_text(errors="ignore").splitlines():
                    parts = line.split()
                    if parts:
                        r = _mk(parts[1] if len(parts) > 1 else parts[0], parts[0], "tuyo")
                        if r:
                            out.append(r)
            elif suf in (".sdf", ".mol"):
                for m in (Chem.SDMolSupplier(str(f)) if suf == ".sdf" else [Chem.MolFromMolFile(str(f))]):
                    if m is not None:
                        r = _mk(m.GetProp("_Name") if m.HasProp("_Name") else f.stem, Chem.MolToSmiles(m), "tuyo")
                        if r:
                            out.append(r)
            elif suf == ".mol2":
                m = Chem.MolFromMol2File(str(f))
                if m is not None:
                    r = _mk(f.stem, Chem.MolToSmiles(m), "tuyo")
                    if r:
                        out.append(r)
        except Exception:
            continue
    return out


def from_pubchem(reaction: rx.Reaction, max_records: int = 25, timeout: int = 30,
                 retries: int = 2) -> tuple:
    """Busca reactivos en PubChem por el grupo de la reacción. Best-effort.

    Devuelve (lista, aviso). Si PubChem no responde, la lista va vacia y el aviso explica por que;
    el cribado no se detiene por ello.
    """
    smarts = reaction.partner_smarts
    url = (f"{PUG}/compound/fastsubstructure/smarts/{urllib.parse.quote(smarts)}/property/"
           f"IsomericSMILES,Title/JSON?MaxRecords={int(max_records)}")
    ultimo = ""
    for intento in range(retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                props = json.loads(resp.read())["PropertyTable"]["Properties"]
            out = []
            for pr in props:
                smi = pr.get("IsomericSMILES") or pr.get("CanonicalSMILES")
                if not smi:
                    continue
                r = _mk(pr.get("Title") or f"CID sin nombre", smi, "pubchem")
                if r and matches_group(r, reaction):
                    out.append(r)
            return out, (None if out else "PubChem respondio pero sin resultados utiles.")
        except Exception as e:
            ultimo = f"{type(e).__name__}: {str(e)[:80]}"
            time.sleep(1.5 * (intento + 1))
    return [], f"PubChem no respondio ({ultimo}). Se continua con las demas fuentes."


# ------------------------------------------------------------------ combinar
def dedup(reagents: Sequence) -> list:
    """Une por InChIKey. Ante duplicados, gana la fuente de mayor prioridad (tuyo > interno > pubchem)."""
    prio = {s: i for i, s in enumerate(FUENTES)}
    best = {}
    for r in reagents:
        if not r or not r.inchikey:
            continue
        prev = best.get(r.inchikey)
        if prev is None or prio.get(r.source, 9) < prio.get(prev.source, 9):
            best[r.inchikey] = r
    return list(best.values())


def build(reaction: rx.Reaction, use_internal: bool = True, user_paths: Sequence = (),
          use_pubchem: bool = False, pubchem_max: int = 25) -> tuple:
    """Reune reactivos de las fuentes elegidas, ya deduplicados y filtrados al grupo de la reacción.

    Devuelve (lista, info) donde info trae los conteos por fuente y el aviso de PubChem si lo hubo.
    """
    todos, aviso = [], None
    if user_paths:
        todos += load_user_files(user_paths)
    if use_internal:
        todos += load_internal(reaction)
    if use_pubchem:
        pc, aviso = from_pubchem(reaction, max_records=pubchem_max)
        todos += pc
    todos = [r for r in todos if matches_group(r, reaction)]
    combinados = dedup(todos)
    from collections import Counter
    conteo = Counter(r.source for r in combinados)
    return combinados, {"total": len(combinados), "por_fuente": dict(conteo), "aviso_pubchem": aviso}
