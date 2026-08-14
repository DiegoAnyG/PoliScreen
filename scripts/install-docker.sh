#!/usr/bin/env bash
# The Docker route with a choice of engines:
#
#   bash scripts/install-docker.sh        # asks
#   bash scripts/install-docker.sh 3      # answers up front, for a rerun or a script
#
# Everything here is something `docker compose` already accepts. What it saves is knowing which
# switch each engine needs, and gnina needs two: the build argument and the overlay that hands
# the container the GPU. Setting only the first produced an image that reported gnina missing on
# the one machine assembled to run it.
set -euo pipefail

cd "$(dirname "$0")/.."

choice="${1:-}"
if [ -z "$choice" ]; then
    cat <<'MENU'
Select what to install:

  1. Base           docking, scoring and the ADMET report
  2. Base + ADCP    adds peptide docking      (+900 MB, accepts the Scripps academic licence)
  3. Base + gnina   adds the second score     (+4.5 GB, needs an NVIDIA GPU)
  4. Full           both

Base already carries the ADMET engine (+1.5 GB); without it the analogue builder and the ADMET
report are disabled. POLISCREEN_WITH_ADMET=0 leaves it out.

MENU
    read -r -p "Choice [1]: " choice
fi

case "${choice:-1}" in
    1) adcp=0; gnina=0 ;;
    2) adcp=1; gnina=0 ;;
    3) adcp=0; gnina=1 ;;
    4) adcp=1; gnina=1 ;;
    *) echo "Not one of 1-4: ${choice}" >&2; exit 2 ;;
esac

files=(-f docker/docker-compose.yml)
if [ "$gnina" = 1 ]; then
    # The overlay sets WITH_GNINA itself and reserves the device; both halves live in that file.
    files+=(-f docker/docker-compose.gpu.yml)
    if ! command -v nvidia-smi >/dev/null 2>&1; then
        echo "Warning: no nvidia-smi here. gnina will build and will not run without an NVIDIA GPU." >&2
    fi
fi

export POLISCREEN_WITH_ADMET="${POLISCREEN_WITH_ADMET:-1}"
export POLISCREEN_WITH_ADCP="$adcp"

cmd=(docker compose "${files[@]}" up --build)
printf '%s\n' "${cmd[*]}"
[ -n "${POLISCREEN_PRINT_ONLY:-}" ] || exec "${cmd[@]}"
