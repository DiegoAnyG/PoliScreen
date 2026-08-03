"""Names of the files and folders a project is made of.

The names are English. Projects created by earlier versions used Spanish ones, so
`artifact()` hands back the legacy path when that is the one on disk. Nothing is
renamed behind the user's back: an existing project keeps working exactly as it is
and a new one is written in English.
"""
from pathlib import Path

RECEPTORS = "receptors"
INPUT_LIGANDS = "input_ligands"
COMPLEXES = "fused_complexes"
DOCKING_CSV = "docking_results.csv"
INTERACTIONS_CSV = "interactions.csv"
VALIDATION_CSV = "redocking_validation.csv"
SUMMARY_CSV = "summary.csv"
ANALOGUES_CSV = "analogues.csv"
SELECTION_JSON = "selection.json"
CONTENTS_TXT = "CONTENTS.txt"

COMPLEX_PREFIX = "Complex_"
COMPLEX_PREFIXES = (COMPLEX_PREFIX, "Complejo_")

LEGACY = {
    RECEPTORS: "receptores",
    INPUT_LIGANDS: "ligandos_entrada",
    COMPLEXES: "Complejos_Fusionados",
    DOCKING_CSV: "resultados_docking.csv",
    INTERACTIONS_CSV: "interacciones.csv",
    VALIDATION_CSV: "validacion_redocking.csv",
    SUMMARY_CSV: "resumen.csv",
    ANALOGUES_CSV: "analogos.csv",
    SELECTION_JSON: "seleccion.json",
    CONTENTS_TXT: "CONTENIDO.txt",
}


def artifact(proj, name) -> Path:
    """Path of a project artifact, preferring the legacy Spanish name if it exists."""
    p = Path(proj) / name
    if p.exists():
        return p
    legacy = LEGACY.get(name)
    if legacy:
        old = Path(proj) / legacy
        if old.exists():
            return old
    return p


def strip_complex_prefix(name: str) -> str:
    """Complex file name without its prefix, in either spelling."""
    for pref in COMPLEX_PREFIXES:
        if str(name).startswith(pref):
            return str(name)[len(pref):]
    return str(name)
