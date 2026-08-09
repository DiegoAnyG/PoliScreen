# Installing PoliScreen

**No Linux or WSL knowledge required.** On Windows the recommended route is
Docker: install Docker Desktop and launch everything with one command.

| System | Recommended route | Cavity detection (fpocket) |
|---|---|---|
| Windows | `scripts\install-windows.ps1`, or the one-click installer | Yes |
| Windows (native, conda) | `environment-windows.yml` | **No**, conda-forge has no Windows build |
| macOS | conda or Docker | Yes |
| Linux | conda or Docker | Yes |

Receptor preparation, docking, PLIP, scoring, confidence and redocking validation
work the same on all of them.

## What every route includes, and what it does not

The four optional engines are **not** part of any default install. This is the same
for Docker, conda and the one-click installer, so the choice of route does not
change what you get:

| | Docker | one-click installer | conda |
|---|---|---|---|
| fpocket | **yes** | **yes** | **yes** |
| ADMET engine (analogue design) | **yes** | no | no, separate environment |
| ADCP (peptide docking) | `-WithAdcp` | never | `scripts/get_adcp.sh` |
| gnina (second scoring function) | `-WithGnina` | never | `scripts/get_gnina.sh` |

ADCP and gnina are Linux-only builds, so on Windows they are reachable **only**
through Docker, never through the native installer.

The ADMET engine is [admelab](https://github.com/DiegoAnyG/admelab), pinned to a
tag so that two people building the image months apart predict with the same
code. It goes in its own environment inside the image, because admet-ai pulls
torch and chemprop whose versions conflict with the docking dependencies; that
conflict is why PoliScreen talks to it over a subprocess. The CPU build of torch
is deliberate: these are small models and the CUDA build would add about 4 GB for
no gain. Leave it out with `-NoAdmet`.

Without it the screening cycle is complete, but the analogue builder and the
ADMET report are disabled. **This is why Docker is the recommended route on every
platform**, and the one-click installer is the fallback for machines where
virtualization cannot be enabled.

---

## Option A — Docker (one command, nothing to configure)

**On Windows, run this instead and skip the rest of this section:**

```powershell
.\scripts\install-windows.ps1
```

From an Administrator PowerShell. It checks the machine, enables the Windows
features WSL 2 needs, updates the WSL kernel, installs Docker Desktop if it is
missing, waits for the engine, builds the image and opens the interface. Add
`-WithAdcp` or `-WithGnina` to build the optional engines into the image.

Run it again after any reboot: it looks at the state of the machine on each run
and continues from wherever it stopped. It saves no state and schedules nothing.

Two things it cannot do, and says so instead of failing quietly:

* **Enable virtualization.** VT-x and AMD-SVM live in the BIOS/UEFI and no program
  running under Windows can change them. The script detects it and tells you which
  setting to look for.
* **Avoid the reboot** Windows requires after enabling WSL.

On macOS and Linux, or to do it by hand:

```bash
git clone https://github.com/DiegoAnyG/PoliScreen.git
cd PoliScreen
docker compose -f docker/docker-compose.yml up --build
```

Open <http://localhost:8501>. Projects are created under `/data`, which is mounted
to `./proyectos` on your machine, so results persist even if you delete the image;
leaving the project field empty already puts them there. To use another location,
set `POLISCREEN_PROJECTS`.

The image pins AutoDock Vina 1.2.5 and every environment version, so two people on
different machines obtain the same result.

---

## Option B — conda (development)

Requires [Miniforge](https://github.com/conda-forge/miniforge) (or Miniconda). On
Windows, inside WSL2.

```bash
git clone https://github.com/DiegoAnyG/PoliScreen.git
cd PoliScreen
conda env create -f environment.yml
conda activate poliscreen
pip install -e .
bash scripts/get_vina.sh
poliscreen ui
```

**Why Vina is installed separately.** The **1.2.x series is not on conda**:
bioconda only ships `autodock-vina 1.1.2` (2011), whose scoring function differs.
`scripts/get_vina.sh` downloads the official release binary and pins the version
(1.2.5 by default), the one PoliScreen was validated with.

---

## Option D — one-click installer

The releases carry an installer per platform, built from `installer/construct.yaml`. It ships the
whole environment: no conda, no Python, nothing to install first. Run it, then launch
`PoliScreen` from the folder it creates.

> **On Windows, choose a destination folder with no accented characters, on an SSD.** Accept the
> one offered, or use something like `C:\PoliScreen`. A path such as
> `…\Investigación\BIOINFORMÁTICA\…` gets through the installer's own check and then breaks conda
> halfway, leaving nothing installed; why, and how to recover, is in *Common problems* below. The
> install unpacks about 48,000 files, so a mechanical drive — with the antivirus reading every one
> of them — turns minutes into a long wait.

It contains PoliScreen, its dependencies and AutoDock Vina. The optional engines are not in it and
are added afterwards with `scripts/get_adcp.sh` and `scripts/get_gnina.sh`, which the installer
copies alongside: ADCP is downloaded from Scripps under an academic licence that cannot be
redistributed and is only built for Linux, and gnina is 2 GB plus the CUDA runtime and does
nothing without an NVIDIA GPU.

It also runs the command line, so the install can be checked without a shell that knows the
environment: `PoliScreen info` on Windows, `./PoliScreen info` elsewhere. Without arguments
it opens the interface.

**Cavity detection is included on all three systems.** conda-forge has no Windows fpocket, so
the release workflow cross-compiles one from `contrib/fpocket-windows/`, where the patches,
the build script and the verification against the Linux build are kept. That binary reads PDB,
which is what PoliScreen uses throughout; it does not read mmCIF, and `mdpocket` is not part
of it.

---

## Option C — native Windows (no WSL, no Docker)

**First install [Miniforge for Windows](https://github.com/conda-forge/miniforge/releases/latest)**
and open the **Miniforge Prompt** it adds to the Start menu. Windows PowerShell does
not know `conda` until you run `conda init powershell` once from that prompt and
open a new window; without this every command below fails with
`The term 'conda' is not recognized`.

Then:

```powershell
git clone https://github.com/DiegoAnyG/PoliScreen.git
cd PoliScreen
conda env create -f environment-windows.yml
conda activate poliscreen
pip install -e .
powershell -ExecutionPolicy Bypass -File scripts\get_vina.ps1
poliscreen ui
```

Vina is downloaded from the official release (`vina_1.2.5_win.exe`) with its
SHA256 verified, as on Linux. **Only limitation: no automatic cavity detection**
(fpocket does not exist for Windows). The app detects and reports this; center the
box on the co-crystallized control.

---

## Optional — second scoring function (gnina, GPU)

`gnina` re-evaluates with a neural network the poses Vina already generated. It
**does not replace Vina or speed up docking**: it provides an independent score,
and the agreement between both feeds the confidence metric.

```bash
bash scripts/get_gnina.sh
export POLISCREEN_GNINA=$HOME/poliscreen_tools/gnina-run
```

Then enable *Re-score poses with gnina* in the step 3 docking settings.

- **NVIDIA GPU required.**
- **~4.5 GB**: the binary is 2 GB and also needs CUDA 12 and cuDNN 9 libraries
  (~2.5 GB). The official *static* binary is **not** self-contained; the installer
  generates a launcher that sets `LD_LIBRARY_PATH`.
- ~2 seconds per compound (only its best pose is re-scored).
- Optional: without gnina, confidence is computed from the other evidence.

---

## Optional — peptide docking (ADCP)

Vina treats the ligand as a torsion tree and its sampling does not cover a
peptide's flexibility. **ADCP** (AutoDock CrankPep) docks it by generating the
conformation with a rotamer library. PoliScreen routes peptides to ADCP
**automatically** when installed; otherwise peptides fall back to Vina with a
warning.

```bash
bash scripts/get_adcp.sh
python scripts/parche_adfr.py        # fixes a numerical defect in ADFRsuite (below)
```

- Only needed to **screen peptides**. Small molecules do not use it.
- Docks peptides of **5–20 residues**. Below 5, Vina is still practical.
- `scripts/parche_adfr.py` is **required**: without it, certain targets
  (phosphoproteins, structures with nearly collinear atoms) crash ADFRsuite
  preparation with a `math domain error`. The patch clamps an angle computation to
  the valid interval and keeps a copy of the original.
- Cycles of **few residues** (5–6) come out strained and ADCP does not close the
  ring well; the app warns. For reliable cycles use ≥ 7–8 residues.

---

## Check the installation

```bash
poliscreen info
```

It should find `vina`, `obabel`, `obrms`, `plip` and `fpocket`. If `obrms` is
missing, the geometric-stability component of the confidence metric will be empty
(everything else works).

---

## Optional — ADMET engine

ADMET prediction and analogue design use `admelab` + ADMET-AI + PyTorch. They live
in a **separate** environment because their dependencies are incompatible with the
docking ones; PoliScreen invokes them as a subprocess. **The full screening works
without this engine**: only the reaction-based analogue builder and the ADMET
report are disabled.

```bash
export POLISCREEN_ADME_PYTHON=/path/to/venv/bin/python
export POLISCREEN_ADME_ROOT=/path/that/contains/admelab
```

---

## Long-term reproducibility

What is pinned and what risk remains, so the install keeps giving the same results
years from now:

| Component | How it is pinned | Residual risk |
|---|---|---|
| Base image | Exact tag `condaforge/miniforge3:26.3.2-3` | Low. Published tags are not rewritten. |
| AutoDock Vina | Release URL **+ verified SHA256** | Low. If the binary changed, the build fails instead of giving other numbers. |
| conda packages | Versions pinned in `environment.yml` | Medium. conda-forge keeps old versions, but an old solve can get hard over time. |
| streamlit | `==1.59.2` from PyPI | Low. PyPI does not delete versions. |

The weak point is the conda solve, resolved on each build. Two ways to remove it:
**publish the pre-built image** to a registry and cite the digest, or **archive on
Zenodo** for a citable DOI. Meanwhile, `docker compose up --build` rebuilds from
scratch and the `sha256sum` guarantees the docking engine is bit-for-bit identical
to the validated one.

---

## Common problems

| Symptom | Cause and fix |
|---|---|
| `Cannot find the vina executable` | Run `bash scripts/get_vina.sh` with the environment active, or set `POLISCREEN_VINA` to the binary path. |
| `PLIP is not installed` | `conda install -c conda-forge plip`. |
| `fpocket not installed` | `conda install -c conda-forge fpocket`. Without it there is no cavity detection; center the box on the control. |
| `AttributeError` on a method that exists | Streamlit does not re-import loaded modules. Restart the process (Ctrl+C and `poliscreen ui`); reloading the page is not enough. |
| The control produces no poses | Extract the ligand with its SMILES in step 1 to fix bond orders; a fragmented ligand yields a PDBQT that Vina rejects. |
| `The term 'conda' is not recognized` (and then the same for `pip`) | conda is not installed, or not on this shell's PATH. Install Miniforge and use the **Miniforge Prompt**, or run `conda init powershell` once from it and open a new window. `pip` fails for the same reason: it comes with the environment. |
| Docker Desktop shows `Virtualization support not detected` and *Engine stopped* | Docker needs hardware virtualization. Enable **Intel VT-x** or **AMD SVM** in the BIOS/UEFI, then, in an administrator PowerShell, `wsl --install` and reboot. Confirm it under Task Manager → Performance → CPU → *Virtualization: Enabled*. On a machine whose BIOS you cannot change, install with conda instead (Option B or C): it needs no virtualization. |
| `docker compose` fails with `500 Internal Server Error` or `check if the server supports the requested API version` | The engine is not running, so compose cannot reach it. Open Docker Desktop and fix what it reports; the compose error is only the symptom. |
| `docker-credential-desktop.exe: exec format error` on build | Docker Desktop leaves a Windows credential helper in WSL that Linux cannot run. Images are public and do not need it: `cp ~/.docker/config.json ~/.docker/config.json.bak && echo '{}' > ~/.docker/config.json`. |
| `docker` not found inside WSL | Docker Desktop → Settings → Resources → **WSL Integration** → enable your distro → *Apply & Restart*. |
| The Windows installer worked once and now fails in the same folder | The first install *wrote* its path into `%USERPROFILE%\.conda\environments.txt`; the second one *reads* that file, and only reading fails (see the row below). Uninstalling does not remove the line: that file belongs to the user's conda, not to PoliScreen — the same file lists any other conda-based tool, so an uninstaller that deleted it would deregister those too. It should drop its own line and does not. Delete the file; conda writes it again. |
| The Windows install is extremely slow at *Setting up the package cache* | About 48,000 small files are being unpacked. On a mechanical drive, with real-time antivirus reading each one, this is the worst possible shape of work. Install on an SSD, or add the destination folder to the antivirus exclusions. Let it finish rather than cancelling halfway. |
| The Windows installer stops with `UnicodeDecodeError('charmap', …)` and `Failed to link extracted packages to …` | The destination folder contains an accented character. conda records installed environments in `%USERPROFILE%\.conda\environments.txt`, writes that file as UTF-8 and reads it back with the system's ANSI codepage; `Á` is two bytes in UTF-8 (`C3 81`) and `0x81` does not exist in cp1252, so the read raises and the transaction rolls back. **Two steps, both needed:** delete `%USERPROFILE%\.conda\environments.txt` (conda writes it again by itself; the failing line is already in there from the previous attempt, so a clean path alone will not help), then install to a folder with no accents. The installer does test the path, but it accepts anything the ANSI codepage can represent — which `Investigación` can, one byte at a time — so it lets this one through. |
