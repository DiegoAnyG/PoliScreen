#!/usr/bin/env bash
# Installs the official AutoDock Vina binary.
#
# Not via conda: bioconda only ships autodock-vina 1.1.2 (2011), whose scoring function
# differs from the 1.2.x series. PoliScreen was validated with 1.2.5, so the version is
# pinned explicitly to keep results reproducible.
#
#   bash scripts/get_vina.sh              # installs 1.2.5 into the active conda env
#   bash scripts/get_vina.sh 1.2.7        # another version
#   bash scripts/get_vina.sh 1.2.5 /usr/local/bin
set -euo pipefail

VERSION="${1:-1.2.5}"
DEST="${2:-${CONDA_PREFIX:-/usr/local}/bin}"
URL="https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v${VERSION}/vina_${VERSION}_linux_x86_64"

# SHA256 of the verified versions. Pins the exact binary: if the release asset changed or the
# download were corrupted, installation fails instead of producing
# resultados distintos en silencio.
case "${VERSION}" in
  1.2.5) SHA256="fa0126a28a9ea9162d1b161dfa92bc76e632416db28ca246278ea4b2dc6860cb" ;;
  *)     SHA256="" ;;
esac

mkdir -p "${DEST}"
echo "Downloading AutoDock Vina ${VERSION} -> ${DEST}/vina"
curl -fsSL "${URL}" -o "${DEST}/vina"

if [ -n "${SHA256}" ]; then
  echo "${SHA256}  ${DEST}/vina" | sha256sum -c - >/dev/null 2>&1 || {
    echo "ERROR: SHA256 mismatch. Corrupted download or altered binary; aborting."
    rm -f "${DEST}/vina"
    exit 1
  }
  echo "Integrity verified (SHA256)."
else
  echo "WARNING: no known sum for ${VERSION}; integrity was not verified."
fi

chmod +x "${DEST}/vina"

echo -n "Installed: "
"${DEST}/vina" --version | head -1

# obrms ships with OpenBabel and is checked here because the confidence metric relies on it
# for the geometric stability of the poses.
if ! command -v obrms >/dev/null 2>&1; then
  echo "WARNING: obrms (OpenBabel) not found. The geometric stability of the confidence will be empty."
fi
