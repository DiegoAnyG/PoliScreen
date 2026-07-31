"""De SMILES a estructura 3D lista para acoplar.

Determinista: semilla fija en la generación de confomeros y elección del de menor energía.
Sin esto la reproducibilidad se rompe antes de llegar al docking.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional, Sequence


def safe_name(text: str, fallback: str = "lig") -> str:
    s = re.sub(r"[^A-Za-z0-9._-]+", "-", str(text)).strip("-")
    return s or fallback


def smiles_to_3d(smiles: str, name: str, out_dir, seed: int = 42, n_confs: int = 50,
                 max_iters: int = 2000) -> Optional[Path]:
    """SMILES -> 3D optimizado. Anade hidrogenos, genera n_confs confomeros (ETKDGv3, semilla fija),
    los optimiza con MMFF94s y escribe el de MENOR energía en SDF. 50 confomeros = misma receta que
    ligand_prep.py del notebook (mejor geometría). La protonación a pH 7.4 se aplica después, al pasar
    a pdbqt con obabel -p 7.4. Determinista con la semilla."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")

    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    m = Chem.AddHs(m)
    # En moléculas muy flexibles (péptidos) generar decenas de confomeros cuesta minutos y aporta
    # poco: Vina vuelve a muestrear cada torsión durante el acoplamiento, así que el confomero de
    # partida apenas influye. Se reduce el número según la flexibilidad, sin perder determinismo.
    from rdkit.Chem import rdMolDescriptors
    if rdMolDescriptors.CalcNumRotatableBonds(m) > 15 or m.GetNumHeavyAtoms() > 45:
        n_confs = min(n_confs, 8)
        max_iters = min(max_iters, 600)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    cids = list(AllChem.EmbedMultipleConfs(m, numConfs=n_confs, params=params))
    if not cids:
        # El algoritmo de distancias falla en moléculas grandes y muy flexibles: no encuentra
        # coordenadas que satisfagan todas las restricciones dentro del número de intentos por
        # defecto. Partir de coordenadas aleatorias y ampliar los intentos lo resuelve casi siempre,
        # y no altera lo anterior porque solo se recurre a ello cuando la vía normal no da nada.
        # Sin este segundo intento, un péptido largo desaparecía del cribado sin mensaje alguno.
        params.useRandomCoords = True
        params.maxIterations = 2000
        params.enforceChirality = False
        cids = list(AllChem.EmbedMultipleConfs(m, numConfs=max(4, n_confs // 2), params=params))
    if not cids:
        return None
    try:
        res = AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=max_iters, mmffVariant="MMFF94s")
        ok = [(e, cid) for cid, (conv, e) in zip(cids, res) if conv == 0]
        best = min(ok)[1] if ok else min(zip([e for _, e in res], cids))[1]
    except Exception:
        best = cids[0]

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{safe_name(name)}.sdf"
    w = Chem.SDWriter(str(path))
    m.SetProp("_Name", str(name))
    w.write(m, confId=best)
    w.close()
    return path if path.exists() and path.stat().st_size > 0 else None


def materialize(smiles_list: Sequence[str], out_dir, names: Optional[Sequence[str]] = None,
                seed: int = 42, prefix: str = "ana", on_progress=None) -> list:
    """Convierte una lista de SMILES en archivos 3D. Devuelve [(nombre, ruta, smiles)]."""
    made = []
    for i, smi in enumerate(smiles_list):
        name = names[i] if names and i < len(names) else f"{prefix}{i + 1:03d}"
        p = smiles_to_3d(smi, name, out_dir, seed=seed)
        if p:
            made.append((name, p, smi))
        if on_progress:
            on_progress(i + 1, len(smiles_list), name, p is not None)
    return made
