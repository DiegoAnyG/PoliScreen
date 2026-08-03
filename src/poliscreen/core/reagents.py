"""Reagent management (e.g. alcohols to esterify).

Reproducibility first: the reliable source is a file the user controls (internal library or their own
upload). PubChem is a best-effort discovery complement: if it does not respond, the screening
continues with what is already there. Everything is deduplicated by InChIKey and keeps its
provenance, so the table can tell apart what the user provided.
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

SOURCES = ("yours", "internal", "pubchem")
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
    """InChIKey of a SMILES; None if invalid. Deduplication key."""
    return _inchikey(smiles)


def _mk(name, smiles, source) -> Optional[Reagent]:
    ik = _inchikey(smiles)
    if not ik:
        return None
    return Reagent(str(name or smiles), smiles, ik, source)


def matches_group(reagent: Reagent, reaction: rx.Reaction) -> bool:
    """The reagent has the group the reaction needs (e.g. an OH to esterify)."""
    from rdkit import Chem
    m = Chem.MolFromSmiles(reagent.smiles)
    patt = Chem.MolFromSmarts(reaction.partner_smarts)
    return m is not None and patt is not None and m.HasSubstructMatch(patt)


def load_internal(reaction: rx.Reaction) -> list:
    out = []
    for row in rx.load_library(reaction):
        r = _mk(row.get("name"), row.get("smiles"), "internal")
        if r:
            out.append(r)
    return out


def load_user_files(paths: Sequence) -> list:
    """User alcohols from csv (name,smiles), .smi, .sdf, .mol or .mol2. Accepts a folder."""
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
                ncol = next((cols[k] for k in ("name", "nombre", "compound", "compuesto",
                                               "alcohol origen", "nombre clave") if k in cols), None)
                scol = next((cols[k] for k in ("smiles", "smile", "smiles alcohol") if k in cols), None)
                if scol is None:
                    continue
                for _, rr in t.iterrows():
                    s = rr[scol]
                    if pd.isna(s):
                        continue
                    r = _mk(rr[ncol] if ncol is not None and pd.notna(rr[ncol]) else s, str(s), "yours")
                    if r:
                        out.append(r)
            elif suf == ".smi":
                for line in f.read_text(errors="ignore").splitlines():
                    parts = line.split()
                    if parts:
                        r = _mk(parts[1] if len(parts) > 1 else parts[0], parts[0], "yours")
                        if r:
                            out.append(r)
            elif suf in (".sdf", ".mol"):
                for m in (Chem.SDMolSupplier(str(f)) if suf == ".sdf" else [Chem.MolFromMolFile(str(f))]):
                    if m is not None:
                        r = _mk(m.GetProp("_Name") if m.HasProp("_Name") else f.stem, Chem.MolToSmiles(m), "yours")
                        if r:
                            out.append(r)
            elif suf == ".mol2":
                m = Chem.MolFromMol2File(str(f))
                if m is not None:
                    r = _mk(f.stem, Chem.MolToSmiles(m), "yours")
                    if r:
                        out.append(r)
        except Exception:
            continue
    return out


def from_pubchem(reaction: rx.Reaction, max_records: int = 25, timeout: int = 30,
                 retries: int = 2) -> tuple:
    """Searches PubChem for reagents by the reaction group. Best-effort.

    Returns (list, notice). If PubChem does not respond, the list is empty and the notice explains
    why; the screening does not stop because of it.
    """
    smarts = reaction.partner_smarts
    url = (f"{PUG}/compound/fastsubstructure/smarts/{urllib.parse.quote(smarts)}/property/"
           f"IsomericSMILES,Title/JSON?MaxRecords={int(max_records)}")
    last_ = ""
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                props = json.loads(resp.read())["PropertyTable"]["Properties"]
            out = []
            for pr in props:
                smi = pr.get("IsomericSMILES") or pr.get("CanonicalSMILES")
                if not smi:
                    continue
                r = _mk(pr.get("Title") or "unnamed CID", smi, "pubchem")
                if r and matches_group(r, reaction):
                    out.append(r)
            return out, (None if out else "PubChem responded but with no useful results.")
        except Exception as e:
            last_ = f"{type(e).__name__}: {str(e)[:80]}"
            time.sleep(1.5 * (attempt + 1))
    return [], f"PubChem did not respond ({last_}). Continuing with the other sources."


def dedup(reagents: Sequence) -> list:
    """Merges by InChIKey. On duplicates, the highest-priority source wins (tuyo > interno > pubchem)."""
    prio = {s: i for i, s in enumerate(SOURCES)}
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
    """Gathers reagents from the chosen sources, already deduplicated and filtered to the reaction group.

    Returns (list, info) where info carries the per-source counts and the PubChem notice if any.
    """
    all_items, notice = [], None
    if user_paths:
        all_items += load_user_files(user_paths)
    if use_internal:
        all_items += load_internal(reaction)
    if use_pubchem:
        pc, notice = from_pubchem(reaction, max_records=pubchem_max)
        all_items += pc
    all_items = [r for r in all_items if matches_group(r, reaction)]
    combinados = dedup(all_items)
    from collections import Counter
    count_ = Counter(r.source for r in combinados)
    return combinados, {"total": len(combinados), "por_fuente": dict(count_), "aviso_pubchem": notice}
