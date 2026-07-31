"""Diseño de análogos + ADMET: puente al motor `admelab`, aislado por entorno.

Por que un subproceso y no un import directo: `admelab` necesita torch/ADMET-AI
(venv Python 3.12) y el motor de docking necesita openbabel/plip/vina (conda 3.11).
No conviven bien. Aislarlo detras de esta interfaz permite:
  - que HOY funcione sin Docker (llamando al venv existente),
  - que MANANA sean dos contenedores, sin tocar el código que lo usa.

Solo depende de la librería estandar (pandas es opcional, para `to_dataframe`).
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Sequence

_RUNNER = Path(__file__).with_name("_admelab_runner.py")
DEFAULT_PYTHON = Path.home() / "adme" / ".venv" / "bin" / "python"
DEFAULT_ROOT = Path.home() / "adme"


class AdmelabError(RuntimeError):
    """Fallo al invocar admelab; el mensaje incluye la causa probable."""


@dataclass
class DesignResult:
    """Análogos generados, ya puntuados por ADME/toxicidad y rankeados."""
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
        """SMILES en el orden del ranking (lo que se manda a dockear)."""
        return [r.get("SMILES") for r in self.rows if r.get("SMILES")]


class AdmelabBridge:
    """Invoca admelab en su propio entorno.

    Parámetros por entorno (útiles en Docker):
      POLISCREEN_ADME_PYTHON  ruta al python del venv de admelab
      POLISCREEN_ADME_ROOT    carpeta que contiene el paquete admelab/
    """

    def __init__(self, python: Optional[os.PathLike] = None,
                 root: Optional[os.PathLike] = None, timeout: int = 3600):
        self.python = Path(python or os.environ.get("POLISCREEN_ADME_PYTHON", DEFAULT_PYTHON))
        self.root = Path(root or os.environ.get("POLISCREEN_ADME_ROOT", DEFAULT_ROOT))
        self.timeout = timeout

    def available(self) -> bool:
        return self.python.exists() and (self.root / "admelab").is_dir()

    def _call(self, params: dict) -> dict:
        if not self.available():
            raise AdmelabError(
                f"No encuentro admelab. Python esperado: {self.python} ; raiz: {self.root}. "
                "Causa probable: rutas distintas en esta maquina. "
                "Solucion: define POLISCREEN_ADME_PYTHON y POLISCREEN_ADME_ROOT."
            )
        # `cwd` NO basta: python anade a sys.path el directorio del SCRIPT, no el de trabajo.
        # Por eso la raiz de admelab se inyecta explicitamente vía PYTHONPATH.
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
                    f"admelab excedio el tiempo limite ({self.timeout}s). "
                    "Causa probable: demasiados analogos o la primera descarga del modelo ADMET-AI."
                )
            if not pout.exists():
                raise AdmelabError(
                    f"admelab no produjo salida (codigo {r.returncode}). "
                    f"stderr: {(r.stderr or '')[-800:]}"
                )
            out = json.loads(pout.read_text())
        if not out.get("ok"):
            raise AdmelabError((out.get("error") or "error desconocido") + "\n" + (out.get("traceback") or ""))
        return out

    def info(self) -> dict:
        """Comprueba el puente: modulos disponibles, versión de python, torch y CUDA."""
        return self._call({"action": "info"})

    def predict(self, smiles: Sequence, use_ml: bool = True) -> "DesignResult":
        """ADMET completo de una lista de SMILES. Devuelve un DesignResult (rows/columns)."""
        out = self._call({"action": "predict", "smiles": list(smiles), "use_ml": bool(use_ml)})
        return DesignResult(out["rows"], out["columns"], len(out["rows"]), len(out["rows"]))

    def reaction_sites(self, smiles: str) -> dict:
        """Sitios reactivos de una molécula: OH clasificados con su viabilidad y si tiene -COOH."""
        return self._call({"action": "reaction_sites", "smiles": smiles})

    def esterify(self, acid: str, alcohols: Sequence, policy: str = "preferred") -> list:
        """Esterifica el acido con cada alcohol. policy 'preferred' usa solo el OH más favorable."""
        return self._call({"action": "esterify", "acid": acid,
                           "alcohols": list(alcohols), "policy": policy})["products"]

    def name_esters(self, ester_smiles: Sequence, alcohol_smiles: Sequence,
                    acid_smiles: Optional[str] = None, alcohol_names: Optional[Sequence] = None,
                    use_web: bool = True) -> list:
        """Nombre IUPAC verificado (OPSIN) de cada ester. Devuelve [{smiles, iupac_name, verified}]."""
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
        """Genera análogos de una molécula lider y los devuelve con ADME/toxicidad y ranking.

        positions        punto(s) de crecimiento (indices de átomo); None = automático
        n_substitutions  número de sustituciones (p. ej. [1, 2])
        use_ml           True usa ADMET-AI (GPU si hay); False solo descriptores RDKit (rápido)
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
