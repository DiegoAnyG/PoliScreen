#!/bin/bash
# Entry point of the installed environment. Without arguments it opens the interface; with them it
# forwards to the CLI, so `./PoliScreen info` checks the installation without activating anything.
here="$(cd "$(dirname "$0")" && pwd)"
export PATH="$here/bin:$PATH"

if [ "$#" -eq 0 ]; then
    exec "$here/bin/poliscreen" ui
fi
exec "$here/bin/poliscreen" "$@"
