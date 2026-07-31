#!/usr/bin/env python3
"""Corrige un defecto numerico de ADFRsuite que impide preparar ciertos receptores.

PyBabel.util.bond_angle calcula el coseno del angulo dividiendo productos escalares y se lo pasa
a math.acos. Por redondeo en coma flotante ese coseno puede salir ligeramente fuera de [-1, 1]
cuando tres atomos estan casi alineados, y entonces acos lanza 'math domain error'. El codigo
original protege el extremo -1 pero no el +1.

Se observo con 1GAG (dominio quinasa del receptor de insulina, con fosfotirosinas y dos Mg): la
preparacion abortaba y, en consecuencia, ADCP no acoplaba ningun peptido contra esa diana.

Acotar el coseno al intervalo valido es la formulacion numericamente estable habitual, no un
apano: para tres atomos alineados el angulo es exactamente 0 o 180 grados.

Uso:  python scripts/parche_adfr.py [ruta_a_ADFRsuite]
Es idempotente: aplicarlo dos veces no cambia nada.
"""
import os
import sys
from pathlib import Path

ORIGINAL = "    else: angle = (math.acos(cos_theta)) * RAD_TO_DEG"
PARCHE = ("    # PoliScreen: el coseno puede salir de [-1, 1] por redondeo cuando los tres atomos\n"
          "    # estan casi alineados, y acos lanza entonces 'math domain error'.\n"
          "    else: angle = (math.acos(max(-1.0, min(1.0, cos_theta)))) * RAD_TO_DEG")


def raiz() -> Path:
    if len(sys.argv) > 1:
        return Path(sys.argv[1])
    env = os.environ.get("POLISCREEN_ADCP")
    if env:
        return Path(env)
    base = Path.home() / "poliscreen_tools"
    cand = sorted(base.glob("adfrsuite*"))
    if not cand:
        sys.exit("No encuentro ADFRsuite. Pasa su ruta como argumento.")
    return cand[-1]


def main() -> int:
    r = raiz()
    objetivo = next(iter(r.rglob("PyBabel/util.py")), None)
    if objetivo is None:
        print(f"No encuentro PyBabel/util.py bajo {r}")
        return 1
    txt = objetivo.read_text(errors="ignore")
    if "PoliScreen: el coseno puede salir" in txt:
        print(f"Ya estaba parcheado: {objetivo}")
        return 0
    if ORIGINAL not in txt:
        print(f"El codigo no coincide con lo esperado; reviselo a mano: {objetivo}")
        return 1
    copia = objetivo.with_suffix(".py.orig")
    if not copia.exists():
        copia.write_text(txt)
    objetivo.write_text(txt.replace(ORIGINAL, PARCHE, 1))
    print(f"Parcheado {objetivo}\nCopia del original en {copia}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
