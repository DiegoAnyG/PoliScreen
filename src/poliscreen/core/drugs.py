"""Approved drugs as a ligand source, with property filters.

Screening a library of drugs that already exist is a different question from building analogues of
a lead: nothing here is designed or made, so there is no synthesizability to judge and no reaction
to be feasible under. The two sources coexist in one run on purpose -- five compounds from the
reaction builder beside five from here is a normal, and informative, experiment.

ChEMBL rather than DrugBank: DrugBank's full dataset cannot be redistributed or fetched on a
user's behalf, which would put a licence problem inside a published tool. ChEMBL is CC0 and its
`max_phase = 4` set is precisely "approved somewhere", about 4200 compounds.

The descriptors are computed here with RDKit rather than taken from the fields ChEMBL already
provides. ChEMBL's alogp is not RDKit's Crippen logP, and a filter that admits a compound on one
definition while the ranking table reports the other is a contradiction the user has to
reconcile. One definition, used everywhere, and an uploaded compound is filtered by the same rule
as a fetched one.
"""
from __future__ import annotations

import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Optional, Sequence

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"

# The properties a filter can act on, with the range each one is scored against elsewhere in
# PoliScreen. Keeping the names identical to the ranking table's columns is deliberate: the user
# filters on MW and later reads MW, and they are the same number.
PROPERTIES = ("MW", "LogP", "TPSA", "HBD", "HBA", "RotB", "QED")

# Lipinski's rule of five, as bounds rather than as a violation count, because a slider is easier
# to reason about than "at most one violation of four". Left open where the rule says nothing.
LIPINSKI = {"MW": (0.0, 500.0), "LogP": (None, 5.0), "HBD": (None, 5.0), "HBA": (None, 10.0)}

# Veber, the other rule most often wanted alongside it: oral bioavailability by flexibility and
# polar surface, which Lipinski does not look at.
VEBER = {"RotB": (None, 10.0), "TPSA": (None, 140.0)}


def descriptors(smiles: str) -> dict:
    """The properties a filter acts on, or an empty dict if the SMILES does not parse."""
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
        RDLogger.DisableLog("rdApp.*")
    except Exception:
        return {}
    m = Chem.MolFromSmiles(smiles) if smiles else None
    if m is None:
        return {}
    return {
        "MW": round(float(Descriptors.MolWt(m)), 2),
        "LogP": round(float(Crippen.MolLogP(m)), 2),
        "TPSA": round(float(rdMolDescriptors.CalcTPSA(m)), 2),
        "HBD": int(Lipinski.NumHDonors(m)),
        "HBA": int(Lipinski.NumHAcceptors(m)),
        "RotB": int(rdMolDescriptors.CalcNumRotatableBonds(m)),
        "QED": round(float(QED.qed(m)), 3),
    }


def passes(row: dict, limits: dict) -> bool:
    """True when every bound in `limits` holds. A bound of None is "no limit on that side".

    A property the molecule has no value for fails rather than passes: an unparseable structure
    must not slip through a filter it was never measured against.
    """
    for prop, bounds in (limits or {}).items():
        lo, hi = bounds if isinstance(bounds, (tuple, list)) else (None, bounds)
        value = row.get(prop)
        if value is None:
            return False
        if lo is not None and value < lo:
            return False
        if hi is not None and value > hi:
            return False
    return True


def apply_filters(rows: Sequence[dict], limits: Optional[dict] = None) -> list:
    """Rows with their descriptors merged in, keeping only those the limits admit."""
    out = []
    for r in rows:
        merged = dict(r)
        merged.update(descriptors(r.get("smiles", "")))
        if not merged.get("MW"):          # did not parse; nothing to filter it on
            continue
        if passes(merged, limits or {}):
            out.append(merged)
    return out


def _rows_from_page(payload: dict) -> list:
    """The usable small molecules on one page of the ChEMBL response.

    Biotherapeutics -- antibodies, peptides above what a docking box holds -- carry no SMILES at
    all, so they are dropped here rather than failing later at 3D generation.
    """
    out = []
    for mol in payload.get("molecules") or []:
        smiles = ((mol.get("molecule_structures") or {}).get("canonical_smiles") or "").strip()
        if not smiles:
            continue
        out.append({"name": (mol.get("pref_name") or mol.get("molecule_chembl_id") or "").strip(),
                    "smiles": smiles,
                    "chembl_id": mol.get("molecule_chembl_id") or "",
                    "source": "chembl"})
    return out


def fetch_approved(cache: Optional[Path] = None, max_records: int = 5000, page: int = 1000,
                   timeout: int = 60, retries: int = 2) -> tuple:
    """Approved drugs (ChEMBL max_phase 4). Returns (rows, notice), best-effort like PubChem.

    Written to `cache` as CSV on the first successful fetch and read from there afterwards, so a
    machine that has run this once keeps working offline. Delete the file to refresh it.
    """
    if cache is not None and Path(cache).exists():
        rows = read_csv(cache)
        if rows:
            return rows, None

    rows, offset, last = [], 0, ""
    while len(rows) < max_records:
        query = urllib.parse.urlencode({"max_phase": 4, "limit": min(page, max_records - len(rows)),
                                        "offset": offset})
        for attempt in range(retries + 1):
            try:
                with urllib.request.urlopen(f"{CHEMBL}?{query}", timeout=timeout) as resp:
                    payload = json.loads(resp.read())
                break
            except Exception as e:
                last = f"{type(e).__name__}: {str(e)[:80]}"
                payload = None
                time.sleep(1.5 * (attempt + 1))
        if payload is None:
            if rows:
                break                      # keep what did arrive; the notice says it is partial
            return [], f"ChEMBL did not respond ({last}). Upload a file of your own instead."
        batch = _rows_from_page(payload)
        rows.extend(batch)
        offset += page
        if not (payload.get("page_meta") or {}).get("next"):
            break

    if not rows:
        return [], "ChEMBL responded but returned no usable structures."
    if cache is not None:
        write_csv(cache, rows)
    return rows, (None if not last else
                  f"Fetched {len(rows)} before ChEMBL stopped responding ({last}). Partial list.")


def read_csv(path) -> list:
    try:
        with Path(path).open(newline="", encoding="utf-8") as fh:
            return [r for r in csv.DictReader(fh) if r.get("smiles")]
    except OSError:
        return []


def write_csv(path, rows: Sequence[dict]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["name", "smiles", "chembl_id", "source"])
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in writer.fieldnames})
