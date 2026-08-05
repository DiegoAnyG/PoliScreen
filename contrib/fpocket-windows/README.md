# fpocket for Windows

conda-forge builds fpocket for `linux-64`, `osx-64` and `osx-arm64` only, so the Windows
installer has no cavity detection and PoliScreen falls back to centring the search box on
the co-crystallized control. That is the better choice when a control exists, but it leaves
a novel or predicted target (SrrAB, for instance, which is an AlphaFold model) with no way
to find the site.

This directory holds what is needed to close that gap: two patches against fpocket 4.2.3 and
a cross-compilation script. Both were built and checked, and the results are below.

## The patches

**`0001`** replaces three things fpocket handed to a shell with the C library equivalents:
`system("mkdir -p …")`, `system("chmod +x …")` and a hardcoded `/tmp` for the qhull scratch
files. The last one is what actually stopped a Windows build: `TMPDIR` is unset there,
`/tmp` does not exist, `fopen()` returned `NULL`, and `load_vvertices()` wrote through it,
an access violation inside the tessellation, with no message. The two `fopen()` calls are
now checked. This patch is platform-neutral and worth having on Linux too, if only because a
filename no longer reaches a shell unquoted.

**`0002`** adds `-DM_NO_MOLFILE`, which leaves out mmCIF input and Amber topologies. Both
come from the vendored VMD `libmolfile_plugin`, and the copy fpocket ships for Windows is
unusable twice over: it is compiled with MSVC, so its objects need `__security_cookie` and
`__GSHandlerCheck`, which MinGW cannot resolve; and it has no `pdbx` plugin at all, so mmCIF
would not link even under MSVC. PDB input is unaffected.

## What was verified

Both builds were run on the same 2833-atom structure (4D44, saFabI).

| | pockets | descriptors |
|---|---|---|
| Linux, upstream 4.2.3 | 140 | reference |
| Linux, patched | 140 | identical except Volume |
| Windows, patched | 140 | identical except Volume |

`Volume` is a Monte Carlo estimate seeded from `time(NULL)`, so it moves on its own. Two runs
of the *same* Linux binary differ by 1.41% on average (5.44% worst case); Windows against
Linux differs by 1.51% (6.79%). **The port differs from Linux by no more than Linux differs
from itself**, and every other descriptor (score, druggability, alpha spheres, SASA,
hydrophobic density) matches exactly.

Feeding both output directories to PoliScreen's own parser (`core/pockets.py`) gives the
same 140 cavities with **no differing pocket** in druggability, score, box centre, box size
or lining residues. The top cavity is pocket 1 on both, centred at (13.37, -20.54, -12.78)
with the same 41 residues. Even under `sizing='volume'`, the only mode that consumes Volume,
the box edge moves at most 0.30 A between platforms against 0.20 A between two Linux runs,
below Vina's 0.375 A grid spacing.

`fpocket.exe`, `tpocket.exe` and `dpocket.exe` are about 1.0 MB each, statically linked, and
import only `KERNEL32.dll` and `msvcrt.dll`.

## Rebuilding

```bash
sudo apt-get install mingw-w64
git clone --depth 1 --branch 4.2.3 https://github.com/Discngine/fpocket.git
cd fpocket
git apply ../0001-*.patch ../0002-*.patch
bash ../build-windows.sh . ../out
```

`check.c` and `cluster.c` are excluded: they belong to no binary and do not compile as they
stand, on Linux either.

## What the Windows build does not do

| | Windows | Why |
|---|---|---|
| PDB input | yes | |
| mmCIF (`.cif`) input | **no** | `read_mmcif.c` needs the molfile `pdbx` plugin, absent from the shipped WIN64 library |
| Amber topology (`--topology_file`) | **no** | `topology.c` needs the molfile `parm7` plugin |
| `fpocket`, `tpocket`, `dpocket` | yes | |
| `mdpocket` | **no** | it reads trajectories through molfile; that is what it is for |

Nothing else changes: every descriptor, the tessellation, the pocket files and the VMD and
PyMOL scripts are produced as on Linux.

For PoliScreen none of this is reachable: it works in PDB throughout, and `core/pockets.py`
calls `fpocket -f` on a PDB copy and uses no other binary.

## Getting it upstream

Both patches are open as
[Discngine/fpocket#184](https://github.com/Discngine/fpocket/pull/184), branched off `master`,
which is at the 4.2.3 tag. An issue is filed on the
[conda-forge feedstock](https://github.com/conda-forge/fpocket-feedstock) asking for `win-64`.
Once that lands, `fpocket 4.2.*` resolves on Windows and the installer picks it up: the
`vendor/` lines and the `fpocket-windows` job come out, and the recipe goes back to a single
`specs` entry.

The open question for the feedstock is the compiler. conda-forge's `{{ compiler('c') }}` is
MSVC on Windows, and fpocket will not build with it without considerably more work: it uses
`getopt_long`, `strtok_r` and at least one variable-length array, none of which MSVC has.
The realistic option is `m2w64-gcc`, which conda-forge publishes for `win-64` and keeps
current (GCC 15). That needs the maintainers' agreement before the recipe is written.
