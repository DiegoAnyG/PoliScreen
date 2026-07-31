"""Prueba de humo del puente a admelab. Se ejecuta con CUALQUIER python3 (solo stdlib).

    python3 tests/smoke_bridge.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from poliscreen.core import AdmelabBridge, AdmelabError  # noqa: E402

b = AdmelabBridge()
print("puente disponible:", b.available())
print("  python admelab:", b.python)
print("  raiz admelab  :", b.root)
if not b.available():
    sys.exit("No encuentro admelab; ajusta POLISCREEN_ADME_PYTHON / POLISCREEN_ADME_ROOT.")

info = b.info()
print("\ninfo del entorno aislado:")
print("  python :", info.get("python"))
print("  torch  :", info.get("torch"), "| cuda:", info.get("cuda"))
print("  modulos:", ", ".join(info.get("modules", [])))

# Diseño rápido (sin ML) para validar el ida y vuelta: ibuprofeno, monosustitucion.
print("\ndiseno de analogos (use_ml=False, rapido)...")
r = b.design("CC(C)Cc1ccc(C(C)C(=O)O)cc1", use_ml=False, n_substitutions=[1], max_decor=25, max_rows=5)
print(f"  generados: {r.n_generated} | puntuados: {r.n_scored} | filas devueltas: {len(r)}")
cols = [c for c in ("SMILES", "is_lead", "MW", "LogP", "TPSA", "QED", "score") if c in r.columns]
print("  columnas de muestra:", cols)
for row in r.rows[:3]:
    print("   ", {k: row.get(k) for k in cols})
smis = r.smiles()
assert smis, "el puente no devolvio SMILES"
print(f"\nSMILES listos para dockear: {len(smis)} (ej. {smis[0]})")
print("\nPUENTE OK: diseno+ADMET corre aislado y devuelve datos usables por el motor de docking.")
