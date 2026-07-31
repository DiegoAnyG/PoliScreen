"""Validación por redocking.

Si el control co-cristalizado no recupera su postura al reacoplarlo, el montaje de esa diana está
mal (caja fuera del sitio, receptor mal preparado, ligando de otro cristal) y ningún número
posterior es fiable. Se ejecuta en cada corrida para avisar antes, no después.
"""
from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from .screening import (base_of, compound_from_pose_name, model_of, normalize_key,
                        receptor_from_name, redocking_rmsd)

RMSD_ACEPTABLE = 2.0


def _formula(pdb) -> dict:
    """Recuento de átomos pesados por elemento, para comparar dos estructuras."""
    out = {}
    try:
        for l in Path(pdb).read_text(errors="ignore").splitlines():
            if l.startswith(("ATOM", "HETATM")):
                e = (l[76:78].strip() or l[12:16].strip()[:1]).upper()
                if e and e != "H":
                    out[e] = out.get(e, 0) + 1
    except Exception:
        pass
    return out


def _motivo_sin_rmsd(ref, pose) -> str:
    """Por qué no hay RMSD. obrms devuelve infinito cuando no empareja los átomos, señal de que se
    comparan dos moléculas distintas."""
    fr, fp = _formula(ref), _formula(pose)
    if fr and fp and fr != fp:
        sobran = {e: n - fp.get(e, 0) for e, n in fr.items() if n > fp.get(e, 0)}
        detalle = ", ".join(f"{e}×{n}" for e, n in sorted(sobran.items())) or "composicion distinta"
        return (f"referencia y pose no son la misma molécula "
                f"({sum(fr.values())} vs {sum(fp.values())} átomos pesados; sobran {detalle} "
                f"en la referencia)")
    return "no calculable"


def redock_validation(controls: Sequence, control_assign: dict, poses_dir,
                      rmsd_ok: float = RMSD_ACEPTABLE, max_poses: int = 10) -> pd.DataFrame:
    """RMSD de cada control frente a su archivo co-cristalizado.

    rmsd_pose1 es la pose de mejor energía; rmsd_min la más parecida al cristal entre las primeras.
    Valido si el RMSD baja de rmsd_ok (2 A por convencion).
    """
    poses_dir = Path(poses_dir)
    filas = []
    for ref in controls:
        ref = Path(ref)
        ck = normalize_key(ref.stem)
        diana = control_assign.get(ck)
        # En híbrido el control se acopla en varios sitios ('receptor~sitio'). Se valida contra el
        # que mejor recupera la postura —donde el control realmente une—; hacerlo contra otro
        # bolsillo daría un fallo esperable que no dice nada del montaje.
        todas = [p for p in poses_dir.glob("*-model*.pdb")
                 if normalize_key(compound_from_pose_name(p.stem)) == ck
                 and (not diana or base_of(receptor_from_name(p.stem)) == base_of(diana))]
        por_sitio = {}
        for p in todas:
            por_sitio.setdefault(receptor_from_name(p.stem), []).append(p)
        if not por_sitio:
            filas.append({"control": ref.stem, "diana": diana or "sin asignar", "rmsd_pose1_A": np.nan,
                          "rmsd_min_A": np.nan, "validado": "faltan poses"})
            continue
        mejor = None
        for sitio, ps in por_sitio.items():
            ps = sorted(ps, key=lambda p: model_of(p.stem))[:max_poses]
            rs = [redocking_rmsd(ref, p) for p in ps]
            v = [r for r in rs if pd.notna(r)]
            clave = min(v) if v else float("inf")
            if mejor is None or clave < mejor[0]:
                mejor = (clave, sitio, rs)
        diana = mejor[1]
        rmsds = mejor[2]
        validos = [r for r in rmsds if pd.notna(r)]
        r1 = rmsds[0] if rmsds else np.nan
        rmin = min(validos) if validos else np.nan
        if pd.notna(r1) and r1 <= rmsd_ok:
            veredicto = "SI"
        elif pd.notna(rmin) and rmin <= rmsd_ok:
            veredicto = "SI (en otra pose)"
        elif validos:
            veredicto = "NO"
        else:
            # Casi siempre la causa es que referencia y pose no son la misma molécula; se comprueba.
            veredicto = _motivo_sin_rmsd(ref, por_sitio[diana][0])
        filas.append({"control": ref.stem, "diana": diana or "sin asignar",
                      "rmsd_pose1_A": round(r1, 3) if pd.notna(r1) else np.nan,
                      "rmsd_min_A": round(rmin, 3) if pd.notna(rmin) else np.nan,
                      "validado": veredicto})
    return pd.DataFrame(filas)


def resumen(val: pd.DataFrame) -> str:
    """Frase corta para avisar al usuario en la interfaz o en la consola."""
    if val is None or val.empty:
        return "Sin controles: no se puede validar el montaje."
    malos = val[~val["validado"].astype(str).str.startswith("SI")]
    n, m = len(val), len(malos)
    sujeto = "El control recupera" if n == 1 else f"Los {n} controles recuperan"
    if malos.empty:
        return f"{sujeto} su postura cristalografica: el montaje es fiable."
    nombres = ", ".join(malos["diana"].astype(str))
    cuantos = "El control no recupera" if n == 1 else f"{m} de {n} controles NO recuperan"
    return (f"ATENCION: {cuantos} su postura ({nombres}). "
            "Revisa la caja o la preparacion de esa diana antes de fiarte del ranking.")
