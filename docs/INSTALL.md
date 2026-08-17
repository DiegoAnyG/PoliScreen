# Installing PoliScreen

| Route | For | fpocket |
|---|---|---|
| **Docker** | Windows, Linux. Recommended. | yes |
| One-click installer | Windows, no Docker | yes |
| conda | development | yes |
| conda on native Windows | no WSL, no Docker | **no** — conda-forge has no Windows build |

The screening cycle is identical on all of them. Optional engines are never installed by default.

| | Docker | installer | conda |
|---|---|---|---|
| ADMET engine (analogue design) | yes | yes, without ADMET-AI | separate environment |
| ADCP (peptides) | build flag | no | `scripts/get_adcp.sh` |
| gnina (second score) | GPU overlay | no | `scripts/get_gnina.sh` |

ADCP and gnina are Linux-only builds, so on Windows they exist only through Docker.

---

## Docker

**Windows, no terminal.** Install Docker Desktop, download `scripts/PoliScreen-Docker.bat`,
double-click it. It pulls the published image the first time and starts the interface; projects go
to `%USERPROFILE%\PoliScreen`. Nothing to build and nothing to clone.

**From a checkout**, when you want to choose the engines or run your own code:

```bash
git clone https://github.com/DiegoAnyG/PoliScreen.git
cd PoliScreen
bash scripts/install-docker.sh
```

It asks which engines go in the image, then builds and starts it. Answer up front to skip the
question: `bash scripts/install-docker.sh 4`. The plain equivalent of option 1 is
`docker compose -f docker/docker-compose.yml up --build`.

**The published image**, to pull without a checkout:

```bash
docker run --rm -it --init -p 127.0.0.1:8501:8501 -v "$PWD/proyectos:/data"   ghcr.io/diegoanyg/poliscreen:latest
```

Cite its **digest**, not the tag: the digest pins every binary inside it, which is what makes a
result reproducible. The release workflow prints it in the run summary.

<http://localhost:8501>. Projects go to `/data`, mounted from `./proyectos`, so they survive the
container. Set `POLISCREEN_PROJECTS` to move them.

**Anything written outside `/data` is lost when the container stops.** The interface says so before
a run if the project folder is not inside it.

Later runs: `docker compose -f docker/docker-compose.yml up` — no `--build`. After a code change,
`--build` reuses the cached environment layers and takes seconds.

gnina is option 3 of that menu: ~4.5 GB, an NVIDIA GPU and the Container Toolkit. Both halves,
the build and the device, come from the overlay:

```bash
docker compose -f docker/docker-compose.yml -f docker/docker-compose.gpu.yml up --build
```

On Windows, `scripts\install-windows.ps1` from an Administrator PowerShell does the whole setup:
WSL 2 features, Docker Desktop, image, interface. Re-run it after any reboot; it continues where it
stopped. It cannot enable virtualization (BIOS/UEFI) or skip the reboot, and says so instead of
failing quietly.

## One-click installer

From the releases. Ships the whole environment — no conda, no Python. Run it, then launch
`PoliScreen` from the folder it creates; `PoliScreen info` checks the install.

> **Choose a destination with no accented characters, on an SSD.** `C:\PoliScreen` is fine. An
> accented path passes the installer's own check and then breaks conda halfway — see
> troubleshooting. It unpacks ~48,000 files, so a mechanical drive with live antivirus is slow.

Includes fpocket on Windows, cross-compiled in the release workflow from `contrib/fpocket-windows/`
because conda-forge has no Windows build.

## conda (development)

Needs [Miniforge](https://github.com/conda-forge/miniforge). On Windows, inside WSL 2.

```bash
conda env create -f environment.yml
conda activate poliscreen
pip install -e .
bash scripts/get_vina.sh
poliscreen ui
```

Vina is installed separately because the 1.2.x series is not on conda — bioconda ships 1.1.2
(2011), whose scoring function differs.

## conda on native Windows

Install Miniforge, then use the **Miniforge Prompt**; PowerShell does not know `conda` until
`conda init powershell` and a new window.

```powershell
conda env create -f environment-windows.yml
conda activate poliscreen
pip install -e .
powershell -ExecutionPolicy Bypass -File scripts\get_vina.ps1
poliscreen ui
```

No cavity detection. Centre the box on the co-crystallised control.

---

## Optional engines

**gnina** — re-scores poses Vina already produced, with a neural network. It does not replace Vina
or speed anything up; the agreement between the two feeds the confidence metric. Needs an NVIDIA
GPU and ~4.5 GB, because the "static" binary still needs CUDA 12 and cuDNN 9.

```bash
bash scripts/get_gnina.sh
export POLISCREEN_GNINA=$HOME/poliscreen_tools/gnina-run
```

**ADCP** — peptide docking, 5–20 residues. Peptides are routed to it automatically when present,
and fall back to Vina with a warning when not.

```bash
bash scripts/get_adcp.sh
python scripts/parche_adfr.py     # required: clamps an angle computation that otherwise
                                  # crashes ADFRsuite preparation on some targets
```

**ADMET-AI predictions** — the rest of `admelab` is RDKit and already ships. Only the ML endpoints
need ADMET-AI, which pulls torch and cannot share the docking environment. Install it separately
and point PoliScreen at it; it runs as a subprocess, so nothing is shared:

```bash
python -m venv <env>
<env>/bin/pip install admet-ai "admelab @ git+https://github.com/DiegoAnyG/admelab"
export POLISCREEN_ADME_PYTHON=<env>/bin/python
export POLISCREEN_ADME_ROOT=<env>/lib/python3.*/site-packages
```

**OPSIN** — IUPAC name verification, off in packaged builds because it needs a JVM. Point
`POLISCREEN_OPSIN` at `opsin.jar` with Java on `PATH`, or drop the jar in the environment's
`share/opsin/`.

---

## Verify

```bash
poliscreen info                      # vina, obabel, obrms, plip, fpocket
poliscreen fingerprint <project>     # per-stage hashes + versions, to compare two machines
```

Reproducibility is pinned by the base image tag, the versions in `environment.yml`, and Vina's
verified SHA256. The conda solve, resolved on each build, is the weak point; publishing the built
image and citing its digest removes it.

---

## Troubleshooting

| Symptom | Fix |
|---|---|
| `Cannot find the vina executable` | `bash scripts/get_vina.sh`, or set `POLISCREEN_VINA`. |
| `PLIP is not installed` / `fpocket not installed` | `conda install -c conda-forge plip fpocket`. Without fpocket, centre the box on the control. |
| `AttributeError` on a method that exists | Streamlit does not re-import loaded modules. Restart the process; reloading the page is not enough. |
| The control produces no poses | Extract the ligand with its SMILES in step 1 — a fragmented ligand gives a PDBQT that Vina rejects. |
| `The term 'conda' is not recognized` | Use the Miniforge Prompt, or run `conda init powershell` once and open a new window. |
| `docker` not found inside WSL | Usually **the engine is not running**, not the integration. Check `/var/run/docker.sock`, `/mnt/wsl/docker-desktop`, and that Docker Desktop is started; only then Settings → Resources → WSL Integration. |
| Docker Desktop: `Virtualization support not detected` | Enable Intel VT-x / AMD SVM in BIOS/UEFI, then `wsl --install` and reboot. Otherwise use conda, which needs no virtualization. |
| `docker-credential-desktop.exe: exec format error` | `cp ~/.docker/config.json ~/.docker/config.json.bak && echo '{}' > ~/.docker/config.json`. |
| Windows installer fails in a folder that worked before | Delete `%USERPROFILE%\.conda\environments.txt`; conda rewrites it. The first install left its path there and uninstalling does not remove it. |
| Windows installer: `UnicodeDecodeError('charmap', …)` | The destination path has an accented character: conda writes `environments.txt` as UTF-8 and reads it back as ANSI. Delete that file **and** install to a path with no accents — both steps. |
| Windows install stuck at *Setting up the package cache* | ~48,000 small files against live antivirus. Use an SSD or exclude the folder, and let it finish. |
