"""Smoke test of the admelab bridge. Runs with ANY python3 (stdlib only).

    python3 tests/smoke_bridge.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from poliscreen.core import AdmelabBridge, AdmelabError  # noqa: E402

b = AdmelabBridge()
print("bridge available:", b.available())
print("  admelab python:", b.python)
print("  admelab root  :", b.root)
if not b.available():
    sys.exit("Cannot find admelab; adjust POLISCREEN_ADME_PYTHON / POLISCREEN_ADME_ROOT.")

info = b.info()
print("\nisolated environment info:")
print("  python :", info.get("python"))
print("  torch  :", info.get("torch"), "| cuda:", info.get("cuda"))
print("  modules:", ", ".join(info.get("modules", [])))

# Quick design (no ML) to validate the round trip: ibuprofen, monosubstitution.
print("\nanalogue design (use_ml=False, fast)...")
r = b.design("CC(C)Cc1ccc(C(C)C(=O)O)cc1", use_ml=False, n_substitutions=[1], max_decor=25, max_rows=5)
print(f"  generated: {r.n_generated} | scored: {r.n_scored} | rows returned: {len(r)}")
cols = [c for c in ("SMILES", "is_lead", "MW", "LogP", "TPSA", "QED", "score") if c in r.columns]
print("  sample columns:", cols)
for row in r.rows[:3]:
    print("   ", {k: row.get(k) for k in cols})
smis = r.smiles()
assert smis, "the bridge returned no SMILES"
print(f"\nSMILES ready to dock: {len(smis)} (e.g. {smis[0]})")
print("\nBRIDGE OK: design+ADMET runs isolated and returns data usable by the docking engine.")
