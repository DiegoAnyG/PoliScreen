#!/bin/bash
# Runs inside the freshly installed environment. PREFIX is set by constructor.
set -euo pipefail

# The installed prefix is a bare environment: it has no conda and therefore no bin/activate.
# Everything here, and the launcher below, works by putting its bin first on PATH.
export PATH="$PREFIX/bin:$PATH"

"$PREFIX/bin/python" -m pip install --no-deps --no-index --find-links "$PREFIX" poliscreen

# The design engine, in this same environment. --no-deps because everything it needs is already in
# the specs; what is missing is ADMET-AI, and admelab reports that itself. Not fatal: the screening
# runs without it, only the analogue design would go.
if ! "$PREFIX/bin/python" -m pip install --no-deps --no-index --find-links "$PREFIX" admelab; then
    echo "WARNING: analogue design was not installed."
fi

# Vina is not a conda package and its 1.2.x series is not on any channel, so the official binary
# is fetched with its SHA256 verified, exactly as scripts/get_vina.sh does.
if ! bash "$PREFIX/scripts/get_vina.sh" 1.2.5 "$PREFIX/bin"; then
    echo "WARNING: Vina was not installed. Run scripts/get_vina.sh later."
fi

# The launcher travels as extra_files; only its executable bit has to be restored here.
chmod +x "$PREFIX/PoliScreen"
