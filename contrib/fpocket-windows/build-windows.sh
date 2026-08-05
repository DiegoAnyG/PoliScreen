#!/bin/bash
# Cross-compiles fpocket for Windows from Linux with mingw-w64, applying the two patches
# next to this script. Produces static .exe files whose only imports are KERNEL32 and
# msvcrt, both part of Windows, so they need no redistributable runtime.
#
#   sudo apt-get install mingw-w64
#   bash build-windows.sh /path/to/fpocket-4.2.3 ./out
#
# Builds fpocket, tpocket and dpocket. mdpocket is left out: it analyses trajectories
# through the VMD molfile plugins, which is exactly what 0002 removes.
set -euo pipefail

SRC="${1:?usage: build-windows.sh <fpocket source dir> <output dir>}"
OUT="${2:?usage: build-windows.sh <fpocket source dir> <output dir>}"
CC="${CC:-x86_64-w64-mingw32-gcc}"

Q="$SRC/src/qhull/src"
INC="-I$SRC/headers -I$SRC/plugins/include -I$SRC/plugins/WIN64/molfile"
CFLAGS="-DM_NO_MOLFILE -DMD_NOT_USE_GSL -DMNO_MEM_DEBUG -std=gnu99 -O2"

# The object lists of the makefile, minus read_mmcif and topology, which are the two
# molfile-dependent translation units.
FPOBJ="fpmain psorting pscoring utils pertable memhandler voronoi sort calc writepdb rpdb
       tparams fparams pocket refine descriptors aa fpocket write_visu fpout atom
       writepocket voronoi_lst asa clusterlib energy"

TPOBJ="tpmain psorting pscoring utils pertable memhandler voronoi sort calc writepdb rpdb
       tparams fparams pocket refine tpocket descriptors aa fpocket write_visu fpout atom
       writepocket voronoi_lst neighbor asa clusterlib energy"

DPOBJ="dpmain psorting pscoring dpocket dparams voronoi sort rpdb descriptors neighbor atom
       aa pertable calc utils writepdb memhandler pocket refine fparams fpocket fpout
       writepocket write_visu asa voronoi_lst clusterlib energy"

QHULL="geom2 geom global io libqhull mem merge poly2 poly qset random rboxlib stat user
       usermem userprintf userprintf_rbox"

rm -rf "$OUT" && mkdir -p "$OUT/.qhull"

# qhull is shared by the three binaries, so it is compiled once.
cd "$OUT/.qhull"
for f in $QHULL; do
    "$CC" -c "$Q/libqhull/$f.c" -I"$Q" -I"$Q/libqhull" -O2 -o "q_$f.o"
done
for f in qvoronoi qconvex; do
    "$CC" -c "$Q/$f/$f.c" -I"$Q" -I"$Q/libqhull" -O2 -o "q_$f.o"
done

build() {
    local target="$1" list="$2" d
    d="$OUT/.$target"
    mkdir -p "$d" && cd "$d"
    for s in $list; do
        "$CC" -c "$SRC/src/$s.c" $INC $CFLAGS
    done
    "$CC" ./*.o "$OUT"/.qhull/*.o -o "$OUT/$target.exe" -lm -static
    cd "$OUT" && rm -rf "$d"
}

build fpocket "$FPOBJ"
build tpocket "$TPOBJ"
build dpocket "$DPOBJ"
rm -rf "$OUT/.qhull"

cd "$OUT" && ls -l ./*.exe
