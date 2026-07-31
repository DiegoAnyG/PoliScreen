"""Registro de reacciones.

Declarativo a proposito: anadir una reacción nueva es anadir una entrada aquí, no tocar el motor.
Cada reacción describe que grupo debe tener la molécula lider y que grupo aporta el reactivo que
se le une. La química la ejecuta admelab; este modulo solo decide que es aplicable y donde.
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
    nombre: str
    lead_smarts: str        # grupo que debe tener la molécula lider
    lead_grupo: str
    partner_smarts: str     # grupo que aporta el reactivo que se une
    partner_grupo: str
    descripcion: str
    kind: str = "coupling"               # 'coupling' une dos moléculas; 'decoration' sustituye in situ
    biblioteca: Optional[str] = None     # csv de reactivos incluido en el paquete


REACTIONS = {
    "esterificacion": Reaction(
        key="esterificacion",
        nombre="Esterificacion de Fischer",
        lead_smarts="[CX3](=[OX1])[OX2H1]",
        lead_grupo="acido carboxilico",
        # OH sobre carbono sp3 o aromatico: alcohol o fenol. Excluye el OH del acido carboxilico
        # (unido a un carbonilo), que no es un alcohol y no debe entrar como reactivo.
        partner_smarts="[OX2H1][CX4,c]",
        partner_grupo="alcohol",
        descripcion=("Acido carboxilico + alcohol da ester. La viabilidad depende del tipo de OH: "
                     "primario y metanol van bien, secundario moderado, fenolico y terciario "
                     "requieren otra ruta (cloruro de acilo o Steglich)."),
        kind="coupling",
        biblioteca="alcoholes.csv",
    ),
    "decoracion": Reaction(
        key="decoracion",
        nombre="Decoracion aromatica (R-group)",
        lead_smarts="[cH]",                 # carbono aromatico con H, posición decorable
        lead_grupo="carbono aromatico con H",
        partner_smarts="",                  # no hay reactivo externo: usa sustituyentes pequeños internos
        partner_grupo="sustituyente",
        descripcion=("Sustituye H por grupos pequenos (F, Cl, CN, OMe...) en carbonos aromaticos. "
                     "Exploratorio: no es una reaccion unica concreta, sirve para barrer el espacio quimico."),
        kind="decoration",
        biblioteca=None,
    ),
}


def get(key: str) -> Reaction:
    if key not in REACTIONS:
        raise KeyError(f"Reaccion desconocida: {key}. Disponibles: {', '.join(REACTIONS)}")
    return REACTIONS[key]


def applicable(smiles: str) -> list:
    """Reacciones que la molécula puede sufrir, por tener el grupo reactivo necesario."""
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
    """Posiciones de la molécula lider donde la reacción puede ocurrir."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles)
    patt = Chem.MolFromSmarts(reaction.lead_smarts)
    if m is None or patt is None:
        return []
    return [{"atomos": list(match)} for match in m.GetSubstructMatches(patt)]


def load_library(reaction: Reaction) -> list:
    """Reactivos incluidos en el paquete para esa reacción."""
    if not reaction.biblioteca:
        return []
    path = DATA_DIR / reaction.biblioteca
    if not path.exists():
        return []
    with path.open(encoding="utf-8") as fh:
        return [dict(row) for row in csv.DictReader(fh)]


def library_from_csv(path, name_col: str = "name", smiles_col: str = "smiles") -> list:
    """Reactivos desde un csv del usuario. Acepta cabeceras en espanol o ingles."""
    alias_n = {name_col, "name", "nombre", "compound", "compuesto"}
    alias_s = {smiles_col, "smiles", "smile"}
    out = []
    with Path(path).open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            claves = {k.lower().strip(): v for k, v in row.items() if k}
            n = next((claves[k] for k in claves if k in alias_n), None)
            s = next((claves[k] for k in claves if k in alias_s), None)
            if s:
                out.append({"name": n or s, "smiles": s})
    if not out:
        raise ValueError("El csv no tiene columnas reconocibles. Se esperan 'name' y 'smiles'.")
    return out
