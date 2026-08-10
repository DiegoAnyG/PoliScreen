"""Analogue design + ADMET: bridge to the `admelab` engine, isolated by environment.

Why a subprocess and not a direct import: `admelab` needs torch/ADMET-AI (Python 3.12 venv) and the
docking engine needs openbabel/plip/vina (conda 3.11). They do not coexist well. Isolating it behind
this interface allows:
  - it to work TODAY without Docker (calling the existing venv),
  - it to become two containers TOMORROW, without touching the code that uses it.

It only depends on the standard library (pandas is optional, for `to_dataframe`).
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Sequence

_RUNNER = Path(__file__).with_name("_admelab_runner.py")
# A venv keeps its interpreter in bin/ on POSIX and in Scripts/ on Windows. Spelling only the
# first left the default unreachable on Windows, where no path can ever match it.
_VENV_PYTHON = ("Scripts", "python.exe") if os.name == "nt" else ("bin", "python")
DEFAULT_PYTHON = Path.home() / "adme" / ".venv" / _VENV_PYTHON[0] / _VENV_PYTHON[1]
DEFAULT_ROOT = Path.home() / "adme"


def _installed() -> Optional[tuple]:
    """(python, root) if admelab is installed in THIS interpreter's environment.

    The separate environment exists because ADMET-AI drags in torch, which does not sit well with
    openbabel/plip/vina. Everything else in admelab — the reactions, the descriptors, the domain —
    is RDKit, which the docking environment already has, so admelab installed alongside works and
    only the ADMET-AI layer is missing (admelab degrades to its RDKit descriptors on its own).
    That is what makes the engine reachable from the one-click installer, which ships no venv.

    find_spec locates the package without importing it: nothing of admelab is loaded in this
    process, the runner still goes out to its own.
    """
    from importlib.util import find_spec
    try:
        spec = find_spec("admelab")
    except (ImportError, ValueError):
        return None
    where = list(getattr(spec, "submodule_search_locations", None) or []) if spec else []
    return (Path(sys.executable), Path(where[0]).parent) if where else None


_EXTERNAL_VALUES = {
    "primario": "primary", "secundario": "secondary", "terciario": "tertiary",
    "buena": "good", "bueno": "good", "moderada": "moderate", "moderado": "moderate",
    "desfavorable": "unfavorable", "dificil": "difficult", "difícil": "difficult",
    "desconocida": "unknown", "desconocido": "unknown",
    "fenol": "phenol", "fenolico": "phenolic", "fenólico": "phenolic",
    "alilico": "allylic", "alílico": "allylic",
    "bencilico": "benzylic", "bencílico": "benzylic",
}
# admelab qualifies a verdict with an explanation in Spanish, e.g.
# 'dificil (fenol poco nucleofilo; usar cloruro de acilo/DMAP)'. The verdict is what the tables
# read and sort by, so it is translated even when the explanation is not one we know.
_EXTERNAL_NOTES = {
    "fenol poco nucleofilo": "poorly nucleophilic phenol",
    "usar cloruro de acilo": "use acyl chloride",
    "impedimento esterico": "steric hindrance",
    "impedimento estérico": "steric hindrance",
    "alcohol terciario": "tertiary alcohol",
}
_GHS_WORDS = {
    "Muy alta toxicidad": "Very high toxicity", "Alta toxicidad": "High toxicity",
    "Toxicidad moderada": "Moderate toxicity", "Moderada toxicidad": "Moderate toxicity",
    "Baja toxicidad": "Low toxicity", "Muy baja toxicidad": "Very low toxicity",
    "No clasificado": "Not classified",
}


def english_value(v):
    """One admelab label in English, or the value unchanged if it is not one."""
    if not isinstance(v, str):
        return v
    hit = _EXTERNAL_VALUES.get(v.strip().lower())
    if hit:
        return hit
    for es, en in _GHS_WORDS.items():
        if es in v:
            return v.replace(es, en)
    # A verdict followed by its explanation: the verdict is what the tables read and sort by, so
    # it is translated even when the note that follows is not one we know.
    head, sep, rest = v.partition(" (")
    verdict = _EXTERNAL_VALUES.get(head.strip().lower())
    if verdict and sep:
        for es, en in _EXTERNAL_NOTES.items():
            rest = rest.replace(es, en)
        return f"{verdict} ({rest}"
    return v


def english_values(row: dict) -> dict:
    return {k: english_value(v) for k, v in row.items()} if isinstance(row, dict) else row


class AdmelabError(RuntimeError):
    """Failure invoking admelab; the message includes the probable cause."""


@dataclass
class DesignResult:
    """Generated analogues, already scored by ADME/toxicity and ranked."""
    rows: list = field(default_factory=list)
    columns: list = field(default_factory=list)
    n_generated: int = 0
    n_scored: int = 0

    def __len__(self) -> int:
        return len(self.rows)

    def to_dataframe(self):
        import pandas as pd
        return pd.DataFrame(self.rows, columns=self.columns or None)

    def smiles(self) -> list:
        """SMILES in ranking order (what is sent to dock)."""
        return [r.get("SMILES") for r in self.rows if r.get("SMILES")]


class AdmelabBridge:
    """Invokes admelab in its own environment.

    Environment parameters (useful in Docker):
      POLISCREEN_ADME_PYTHON  path to admelab's venv python
      POLISCREEN_ADME_ROOT    folder containing the admelab/ package
    """

    def __init__(self, python: Optional[os.PathLike] = None,
                 root: Optional[os.PathLike] = None, timeout: int = 3600):
        # The venv comes first, so a machine that has one keeps using it, ADMET-AI and all; the
        # copy installed alongside is the fallback for the machines that have no venv to find.
        fallback = (DEFAULT_PYTHON, DEFAULT_ROOT)
        if not DEFAULT_PYTHON.exists():
            fallback = _installed() or fallback
        self.python = Path(python or os.environ.get("POLISCREEN_ADME_PYTHON", fallback[0]))
        self.root = Path(root or os.environ.get("POLISCREEN_ADME_ROOT", fallback[1]))
        self.timeout = timeout

    def available(self) -> bool:
        return self.python.exists() and (self.root / "admelab").is_dir()

    def _call(self, params: dict) -> dict:
        if not self.available():
            raise AdmelabError(
                f"Cannot find admelab. Expected python: {self.python} ; root: {self.root}. "
                "Probable cause: different paths on this machine. "
                "Solution: set POLISCREEN_ADME_PYTHON and POLISCREEN_ADME_ROOT."
            )
        env = dict(os.environ)
        env["PYTHONPATH"] = os.pathsep.join([str(self.root)] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else []))
        with tempfile.TemporaryDirectory() as td:
            pin, pout = Path(td) / "in.json", Path(td) / "out.json"
            pin.write_text(json.dumps(params))
            try:
                r = subprocess.run(
                    [str(self.python), str(_RUNNER), str(pin), str(pout)],
                    cwd=str(self.root), env=env, capture_output=True, text=True, timeout=self.timeout,
                )
            except subprocess.TimeoutExpired:
                raise AdmelabError(
                    f"admelab exceeded the time limit ({self.timeout}s). "
                    "Probable cause: too many analogues or the first download of the ADMET-AI model."
                )
            if not pout.exists():
                raise AdmelabError(
                    f"admelab produced no output (code {r.returncode}). "
                    f"stderr: {(r.stderr or '')[-800:]}"
                )
            out = json.loads(pout.read_text())
        if not out.get("ok"):
            raise AdmelabError((out.get("error") or "unknown error") + "\n" + (out.get("traceback") or ""))
        return out

    def info(self) -> dict:
        """Checks the bridge: available modules, python version, torch and CUDA."""
        return self._call({"action": "info"})

    def predict(self, smiles: Sequence, use_ml: bool = True) -> "DesignResult":
        """Full ADMET of a list of SMILES. Returns a DesignResult (rows/columns)."""
        out = self._call({"action": "predict", "smiles": list(smiles), "use_ml": bool(use_ml)})
        rows = [english_values(r) for r in out["rows"]]
        return DesignResult(rows, out["columns"], len(rows), len(rows))

    def applicability(self, smiles: Sequence, endpoints: Optional[Sequence] = None,
                      percentile: float = 5.0, detail: bool = False) -> "DesignResult":
        """Applicability domain of the ADMET-AI predictions: one row per molecule.

        Says how close each structure sits to the closest compound the model was trained on, and on
        how many endpoints that distance is inside the training set's own spread. It qualifies the
        ADMET numbers rather than replacing them: outside the domain the prediction is unsupported,
        not wrong. Needs admelab >= 0.3; older ones have no `domain` module.
        """
        out = self._call({"action": "applicability", "smiles": list(smiles),
                          "endpoints": list(endpoints) if endpoints else None,
                          "percentile": float(percentile), "detail": bool(detail)})
        rows = out["rows"]
        return DesignResult(rows, out["columns"], len(rows), len(rows))

    def has_applicability(self) -> bool:
        """Whether the installed admelab carries the domain module. An older one does not, and a
        machine may have no admelab at all, so every caller has to cope with its absence anyway."""
        try:
            return "domain" in (self.info().get("modules") or [])
        except AdmelabError:
            return False

    def reaction_sites(self, smiles: str) -> dict:
        """Reactive sites of a molecule: OHs classified with their feasibility and whether it has -COOH."""
        return self._call({"action": "reaction_sites", "smiles": smiles})

    def esterify(self, acid: str, alcohols: Sequence, policy: str = "preferred") -> list:
        """Esterifies the acid with each alcohol. policy 'preferred' uses only the most favorable OH."""
        products = self._call({"action": "esterify", "acid": acid,
                               "alcohols": list(alcohols), "policy": policy})["products"]
        return [english_values(p) for p in products]

    def name_esters(self, ester_smiles: Sequence, alcohol_smiles: Sequence,
                    acid_smiles: Optional[str] = None, alcohol_names: Optional[Sequence] = None,
                    use_web: bool = True) -> list:
        """Verified IUPAC name (OPSIN) of each ester. Returns [{smiles, iupac_name, verified}]."""
        return self._call({"action": "name_esters", "ester_smiles": list(ester_smiles),
                           "alcohol_smiles": list(alcohol_smiles), "acid_smiles": acid_smiles,
                           "alcohol_names": list(alcohol_names) if alcohol_names else None,
                           "use_web": bool(use_web)})["names"]

    def design(self, lead_smiles: str,
               methods: Sequence[str] = ("decoration",),
               use_ml: bool = True,
               positions: Optional[Sequence[int]] = None,
               n_substitutions: Sequence[int] = (1,),
               scope: str = "aromatic_ch",
               max_decor: int = 300,
               max_brics: int = 60,
               substituents: Optional[dict] = None,
               filters: Optional[dict] = None,
               include_lead: bool = True,
               max_rows: Optional[int] = None) -> DesignResult:
        """Generates analogues of a lead molecule and returns them with ADME/toxicity and ranking.

        positions        growth point(s) (atom indices); None = automatic
        n_substitutions  number of substitutions (e.g. [1, 2])
        use_ml           True uses ADMET-AI (GPU if available); False only RDKit descriptors (fast)
        """
        out = self._call({
            "action": "design",
            "lead_smiles": lead_smiles,
            "methods": list(methods),
            "use_ml": bool(use_ml),
            "positions": list(positions) if positions else None,
            "n_substitutions": list(n_substitutions),
            "scope": scope,
            "max_decor": int(max_decor),
            "max_brics": int(max_brics),
            "substituents": substituents,
            "filters": filters,
            "include_lead": bool(include_lead),
            "max_rows": max_rows,
        })
        return DesignResult(out["rows"], out["columns"], out["n_generated"], out["n_scored"])
