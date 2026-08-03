"""Reaction registry.

Deliberately declarative: adding a new reaction means adding an entry here, not touching the engine.
Each reaction describes which group the lead molecule must have and which group the coupling reagent
contributes. The chemistry is run by admelab; this module only decides what is applicable and where.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


@dataclass(frozen=True)
class Reaction:
    key: str
    name_: str
    lead_smarts: str
    lead_grupo: str
    partner_smarts: str
    partner_grupo: str
    description_: str
    kind: str = "coupling"
    library_: Optional[str] = None


REACTIONS = {
    "esterificacion": Reaction(
        key="esterificacion",
        name_="Fischer esterification",
        lead_smarts="[CX3](=[OX1])[OX2H1]",
        lead_grupo="carboxylic acid",
        # Alcohol or phenol OH; excludes the carboxylic-acid OH, which is not a reagent.
        partner_smarts="[OX2H1][CX4,c]",
        partner_grupo="alcohol",
        description_=("Carboxylic acid + alcohol gives ester. Feasibility depends on the OH type: "
                     "primary and methanol work well, secondary moderate, phenolic and tertiary "
                     "need another route (acyl chloride or Steglich)."),
        kind="coupling",
        library_="alcoholes.csv",
    ),
    "decoracion": Reaction(
        key="decoracion",
        name_="Aromatic decoration (R-group)",
        lead_smarts="[cH]",
        lead_grupo="aromatic carbon with H",
        partner_smarts="",
        partner_grupo="substituent",
        description_=("Replaces H with small groups (F, Cl, CN, OMe...) on aromatic carbons. "
                     "Exploratory: not a single concrete reaction, it sweeps the chemical space."),
        kind="decoration",
        library_=None,
    ),
}


def get(key: str) -> Reaction:
    if key not in REACTIONS:
        raise KeyError(f"Unknown reaction: {key}. Available: {', '.join(REACTIONS)}")
    return REACTIONS[key]


def applicable(smiles: str) -> list:
    """Reactions the molecule can undergo, by having the required reactive group."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return []
    out = []
    for r in REACTIONS.values():
        patt = Chem.MolFromSmarts(r.lead_smarts)
        if patt is not None and m.HasSubstructMatch(patt):
            out.append(r)
    return out


def lead_sites(smiles: str, reaction: Reaction) -> list:
    """Positions on the lead molecule where the reaction can occur."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles)
    patt = Chem.MolFromSmarts(reaction.lead_smarts)
    if m is None or patt is None:
        return []
    return [{"atomos": list(match)} for match in m.GetSubstructMatches(patt)]


def load_library(reaction: Reaction) -> list:
    """Reagents bundled in the package for that reaction."""
    if not reaction.library_:
        return []
    path = DATA_DIR / reaction.library_
    if not path.exists():
        return []
    with path.open(encoding="utf-8") as fh:
        return [dict(row) for row in csv.DictReader(fh)]


def library_from_csv(path, name_col: str = "name", smiles_col: str = "smiles") -> list:
    """Reagents from a user csv. Accepts Spanish or English headers."""
    alias_n = {name_col, "name", "nombre", "compound", "compuesto"}
    alias_s = {smiles_col, "smiles", "smile"}
    out = []
    with Path(path).open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            keys_ = {k.lower().strip(): v for k, v in row.items() if k}
            n = next((keys_[k] for k in keys_ if k in alias_n), None)
            s = next((keys_[k] for k in keys_ if k in alias_s), None)
            if s:
                out.append({"name": n or s, "smiles": s})
    if not out:
        raise ValueError("The csv has no recognizable columns. 'name' and 'smiles' are expected.")
    return out
