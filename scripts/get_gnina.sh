#!/usr/bin/env bash
# Instala gnina (segunda función de puntuación, red neuronal, GPU NVIDIA).
#   bash scripts/get_gnina.sh
# El binario "static" no es autocontenido: se genera un envoltorio que fija LD_LIBRARY_PATH.
set -euo pipefail

DEST="${1:-$HOME/poliscreen_tools}"
URL="https://github.com/gnina/gnina/releases/latest/download/gnina.cuda12.8.static"

if ! command -v nvidia-smi >/dev/null 2>&1; then
  echo "AVISO: no detecto una GPU NVIDIA. gnina la necesita; la instalacion seguira, pero no funcionara."
fi

mkdir -p "${DEST}"
echo "Descargando gnina (unos 2 GB) -> ${DEST}/gnina"
curl -fL --progress-bar "${URL}" -o "${DEST}/gnina.parcial"
mv "${DEST}/gnina.parcial" "${DEST}/gnina"
chmod +x "${DEST}/gnina"

echo "Instalando las librerias de CUDA 12 y cuDNN 9 que el binario necesita..."
python -m pip install --quiet --upgrade --target "${DEST}/cuda12" \
  nvidia-cuda-runtime-cu12 nvidia-cudnn-cu12 nvidia-cublas-cu12 nvidia-cusparse-cu12 \
  nvidia-cufft-cu12 nvidia-curand-cu12 nvidia-cusolver-cu12 nvidia-nvjitlink-cu12

# Lanzador: fija LD_LIBRARY_PATH con las rutas reales de las librerias recien instaladas.
python - "$DEST" <<'PY'
import stat, sys
from pathlib import Path
dest = Path(sys.argv[1])
dirs = sorted({str(p.parent) for p in (dest / "cuda12").rglob("*.so*")})
w = dest / "gnina-run"
w.write_text("\n".join([
    "#!/usr/bin/env bash",
    "# Lanzador de gnina: el binario oficial necesita CUDA 12 y cuDNN 9 junto a el.",
    'export LD_LIBRARY_PATH="' + ":".join(dirs) + ':${LD_LIBRARY_PATH}"',
    'exec "' + str(dest / "gnina") + '" "$@"',
    "",
]))
w.chmod(w.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
print("lanzador:", w)
PY

echo -n "Version instalada: "
"${DEST}/gnina-run" --version | head -1

cat <<EOF

Para que PoliScreen lo encuentre, exporta la ruta del lanzador:

    export POLISCREEN_GNINA=${DEST}/gnina-run

Anadelo a tu ~/.bashrc si quieres que persista. Despues activa la casilla
"Re-puntuar las poses con gnina" en los ajustes de docking.
EOF
