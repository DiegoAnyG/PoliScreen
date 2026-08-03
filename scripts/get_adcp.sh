#!/usr/bin/env bash
# Installs AutoDock CrankPep (ADCP), part of the Scripps ADFR suite, for peptide docking.
#   bash scripts/get_adcp.sh
# About 900 MB. Then apply scripts/parche_adfr.py, which fixes a numerical defect in ADFRsuite.
set -euo pipefail

DEST="${1:-$HOME/poliscreen_tools}"
URL="https://sourceforge.net/projects/adfrsuite/files/ADFRsuite-1.1dev/ADFRsuite_Linux-x86_64_1.1dev_Install/download"

mkdir -p "${DEST}"
cd "${DEST}"

echo "Downloading the ADFR suite (about 110 MB)..."
curl -fL --progress-bar "${URL}" -o adfr_install
chmod +x adfr_install

echo "Extrayendo..."
./adfr_install --mode silent --prefix "${DEST}/adfr" >/dev/null 2>&1
cd "${DEST}/adfr"
tar xzf ADFRsuite_x86_64Linux_*.tar.gz
cd ADFRsuite_x86_64Linux_*/

echo "Installing (accepts the academic licence; for commercial use see LICENSE.txt)..."
yes Y | ./install.sh -d "${DEST}/adfrsuite" -c 0 >/dev/null 2>&1

BIN="$(echo "${DEST}"/adfrsuite/ADFRsuite-*/bin)"
if [ -x "${BIN}/adcp" ]; then
  echo "Installed in ${BIN}"
  cat <<EOF

PoliScreen finds it there automatically. If you installed it elsewhere:

    export POLISCREEN_ADCP=${BIN}

The binaries need libgomp (OpenMP). PoliScreen uses the one from the active conda
environment, so it does not have to be installed system-wide.
EOF
else
  echo "ERROR: the adcp executable was not found after installation." >&2
  exit 1
fi
