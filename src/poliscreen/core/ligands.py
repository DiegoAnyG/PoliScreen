"""From SMILES to a 3D structure ready to dock.

Deterministic: fixed seed for conformer generation and lowest-energy selection.
Without this, reproducibility breaks before reaching the docking.
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
    """SMILES -> optimized 3D. Adds hydrogens, generates n_confs conformers (ETKDGv3, fixed seed),
    optimizes them with MMFF94s and writes the LOWEST-energy one to SDF. 50 conformers = same recipe
    as the notebook's ligand_prep.py (better geometry). Protonation at pH 7.4 is applied later, when
    converting to pdbqt with obabel -p 7.4. Deterministic given the seed."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")

    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    m = Chem.AddHs(m)
    # The starting conformer barely matters: Vina re-samples every torsion while docking.
    from rdkit.Chem import rdMolDescriptors
    if rdMolDescriptors.CalcNumRotatableBonds(m) > 15 or m.GetNumHeavyAtoms() > 45:
        n_confs = min(n_confs, 8)
        max_iters = min(max_iters, 600)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    cids = list(AllChem.EmbedMultipleConfs(m, numConfs=n_confs, params=params))
    if not cids:
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
    """Converts a list of SMILES into 3D files. Returns [(name, path, smiles)]."""
    made = []
    for i, smi in enumerate(smiles_list):
        name = names[i] if names and i < len(names) else f"{prefix}{i + 1:03d}"
        p = smiles_to_3d(smi, name, out_dir, seed=seed)
        if p:
            made.append((name, p, smi))
        if on_progress:
            on_progress(i + 1, len(smiles_list), name, p is not None)
    return made
