#!/usr/bin/env bash
# Instala AutoDock CrankPep (ADCP), en la suite ADFR de Scripps, para el acoplamiento de péptidos.
#   bash scripts/get_adcp.sh
# Ocupa ~900 MB. Aplica después scripts/parche_adfr.py (corrige un defecto numérico de ADFRsuite).
set -euo pipefail

DEST="${1:-$HOME/poliscreen_tools}"
URL="https://sourceforge.net/projects/adfrsuite/files/ADFRsuite-1.1dev/ADFRsuite_Linux-x86_64_1.1dev_Install/download"

mkdir -p "${DEST}"
cd "${DEST}"

echo "Descargando la suite ADFR (unos 110 MB)..."
curl -fL --progress-bar "${URL}" -o adfr_install
chmod +x adfr_install

echo "Extrayendo..."
./adfr_install --mode silent --prefix "${DEST}/adfr" >/dev/null 2>&1
cd "${DEST}/adfr"
tar xzf ADFRsuite_x86_64Linux_*.tar.gz
cd ADFRsuite_x86_64Linux_*/

echo "Instalando (acepta la licencia academica; para uso comercial, ver LICENSE.txt)..."
yes Y | ./install.sh -d "${DEST}/adfrsuite" -c 0 >/dev/null 2>&1

BIN="$(echo "${DEST}"/adfrsuite/ADFRsuite-*/bin)"
if [ -x "${BIN}/adcp" ]; then
  echo "Instalado en ${BIN}"
  cat <<EOF

PoliScreen lo detecta solo si esta en esa ruta. Si lo instalaste en otro sitio:

    export POLISCREEN_ADCP=${BIN}

Nota: los binarios necesitan libgomp (OpenMP). PoliScreen usa la del entorno conda activo, asi que
no hace falta instalarla en el sistema.
EOF
else
  echo "ERROR: no encuentro el ejecutable adcp tras la instalacion." >&2
  exit 1
fi
