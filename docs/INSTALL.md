# Installing PoliScreen

**No Linux or WSL knowledge required.** On Windows the recommended route is
Docker: install Docker Desktop and launch everything with one command.

| System | Recommended route | Cavity detection (fpocket) |
|---|---|---|
| Windows | one-click installer, or Docker Desktop | Yes |
| Windows (native, conda) | `environment-windows.yml` | **No** — conda-forge has no Windows build |
| macOS | conda or Docker | Yes |
| Linux | conda or Docker | Yes |

Everything else — receptor preparation, docking, PLIP, scoring, confidence and
ADMET — works the same on all three platforms. Without fpocket you only lose
automatic cavity detection; the search box can be centered on the co-crystallized
ligand, which is the recommended option anyway.

---

## Option A — Docker (one command, nothing to configure)

Requires Docker. On Windows: Docker Desktop with WSL2 integration.

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
