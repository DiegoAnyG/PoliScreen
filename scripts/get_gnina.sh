#!/usr/bin/env bash
# Installs gnina: a second scoring function (neural network, NVIDIA GPU).
#   bash scripts/get_gnina.sh
# The "static" binary is not self-contained, so a wrapper setting LD_LIBRARY_PATH is generated.
set -euo pipefail

DEST="${1:-$HOME/poliscreen_tools}"
URL="https://github.com/gnina/gnina/releases/latest/download/gnina.cuda12.8.static"

if ! command -v nvidia-smi >/dev/null 2>&1; then
  echo "WARNING: no NVIDIA GPU detected. gnina needs one; installation continues but it will not run."
fi

mkdir -p "${DEST}"
echo "Downloading gnina (about 2 GB) -> ${DEST}/gnina"
curl -fL --progress-bar "${URL}" -o "${DEST}/gnina.parcial"
mv "${DEST}/gnina.parcial" "${DEST}/gnina"
chmod +x "${DEST}/gnina"

echo "Installing the CUDA 12 and cuDNN 9 libraries the binary needs..."
python -m pip install --quiet --upgrade --target "${DEST}/cuda12" \
  nvidia-cuda-runtime-cu12 nvidia-cudnn-cu12 nvidia-cublas-cu12 nvidia-cusparse-cu12 \
  nvidia-cufft-cu12 nvidia-curand-cu12 nvidia-cusolver-cu12 nvidia-nvjitlink-cu12

# Wrapper: sets LD_LIBRARY_PATH to the real paths of the libraries just installed.
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

echo -n "Installed version: "
"${DEST}/gnina-run" --version | head -1

cat <<EOF

So that PoliScreen finds it, export the path of the wrapper:

    export POLISCREEN_GNINA=${DEST}/gnina-run

Add it to your ~/.bashrc to make it persist, then tick "Re-score the poses with
gnina" in the docking settings.
EOF
