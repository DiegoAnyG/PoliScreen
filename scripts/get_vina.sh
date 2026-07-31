#!/usr/bin/env bash
# Instala el binario oficial de AutoDock Vina.
#
# Por que no via conda: bioconda solo distribuye autodock-vina 1.1.2 (2011), cuya funcion de
# puntuacion difiere de la serie 1.2.x. PoliScreen se valido con 1.2.5, asi que la version se
# fija explicitamente para que los resultados sean reproducibles.
#
#   bash scripts/get_vina.sh              # instala 1.2.5 en el entorno conda activo
#   bash scripts/get_vina.sh 1.2.7        # otra version
#   bash scripts/get_vina.sh 1.2.5 /usr/local/bin
set -euo pipefail

VERSION="${1:-1.2.5}"
DEST="${2:-${CONDA_PREFIX:-/usr/local}/bin}"
URL="https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v${VERSION}/vina_${VERSION}_linux_x86_64"

# Suma SHA256 de las versiones verificadas. Fija la identidad EXACTA del binario: si el asset de
# la release cambiara o la descarga se corrompiera, la instalacion falla en vez de producir
# resultados distintos en silencio.
case "${VERSION}" in
  1.2.5) SHA256="fa0126a28a9ea9162d1b161dfa92bc76e632416db28ca246278ea4b2dc6860cb" ;;
  *)     SHA256="" ;;
esac

mkdir -p "${DEST}"
echo "Descargando AutoDock Vina ${VERSION} -> ${DEST}/vina"
curl -fsSL "${URL}" -o "${DEST}/vina"

if [ -n "${SHA256}" ]; then
  echo "${SHA256}  ${DEST}/vina" | sha256sum -c - >/dev/null 2>&1 || {
    echo "ERROR: la suma SHA256 no coincide. Descarga corrupta o binario alterado; se aborta."
    rm -f "${DEST}/vina"
    exit 1
  }
  echo "Integridad verificada (SHA256)."
else
  echo "AVISO: no hay suma conocida para ${VERSION}; no se verifico la integridad."
fi

chmod +x "${DEST}/vina"

echo -n "Instalado: "
"${DEST}/vina" --version | head -1

# obrms viene con OpenBabel (incluido en environment.yml); se comprueba porque la metrica de
# confianza depende de el para la estabilidad geometrica de las poses.
if ! command -v obrms >/dev/null 2>&1; then
  echo "AVISO: no encuentro 'obrms' (OpenBabel). La estabilidad geometrica de la confianza quedara vacia."
fi
