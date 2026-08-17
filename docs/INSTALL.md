# Installing PoliScreen

What you need, including the obvious: **[REQUIREMENTS.md](REQUIREMENTS.md)**. The short version is
Docker Desktop, 15 GB free, and internet the first time.

| Route | For | Everything included |
|---|---|---|
| **Container** | Windows, Linux. Recommended. | yes, and it is the reproducible reference |
| One-click installer | Windows without Docker | no ADCP, no gnina, no ADMET-AI |
| conda | development | you install the engines |

The screening cycle is the same on all of them.

---

## Container

### Windows, no terminal

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop).
2. Download **`PoliScreen-Docker.bat`** from the [latest release](https://github.com/DiegoAnyG/PoliScreen/releases/latest).
3. Double-click it.

It downloads the image the first time (~3 GB), checks for a newer one on every later start, and
opens the interface. Projects go to `%USERPROFILE%\PoliScreen`. Nothing to build, nothing to clone.

If something is missing it says which piece, in plain words.

### Any machine with Docker

```bash
docker run --rm -it --init -p 127.0.0.1:8501:8501 -v "$PWD/proyectos:/data" \
  ghcr.io/diegoanyg/poliscreen:latest
```

<http://localhost:8501>. Projects live in `./proyectos`, mounted at `/data`, so they survive the
container. **Anything written outside `/data` is lost when the container stops** — the interface
says so before a run if the project folder is not inside it.

For a published result, pin the **digest** rather than the tag:

```
ghcr.io/diegoanyg/poliscreen@sha256:a22a7ccbac87c5174b81a6080de026a48700f5867b8e92378060270742daa20f
```

A tag can be moved; a digest pins every binary inside the image, which is what makes the number
reproducible by someone else.

### From a checkout

To choose the engines, or to run your own code:

```bash
git clone https://github.com/DiegoAnyG/PoliScreen.git
cd PoliScreen
bash scripts/install-docker.sh
```

It asks what to include — base, +ADCP for peptides, +gnina for the second score, or both — then
builds and starts it. `bash scripts/install-docker.sh 4` answers up front.

gnina needs an NVIDIA GPU and adds ~4.5 GB; its build flag and the device handover both come from
the same overlay:

```bash
docker compose -f docker/docker-compose.yml -f docker/docker-compose.gpu.yml up --build
```

Later runs need no `--build`.

---

## One-click installer

From the [releases](https://github.com/DiegoAnyG/PoliScreen/releases/latest). Ships the whole
environment — no conda, no Python. Run it, then launch `PoliScreen` from the folder it creates;
`PoliScreen info` checks the install.

> **Install to a path with no accented characters, on an SSD.** `C:\PoliScreen` is fine. An accented
> path passes the installer's own check and then breaks conda halfway. It unpacks ~48,000 files, so
> a mechanical drive with live antivirus is slow.

ADCP, gnina and ADMET-AI are not in it: the first two are Linux-only builds, and ADMET-AI conflicts
with the docking dependencies. On Windows they exist only through the container.

Includes fpocket, cross-compiled in the release workflow because conda-forge has no Windows build.

---

## conda (development)

Needs [Miniforge](https://github.com/conda-forge/miniforge). On Windows, inside WSL 2.

```bash
conda env create -f environment.yml
conda activate poliscreen
pip install -e .
bash scripts/get_vina.sh
poliscreen ui
```

Vina is installed separately because 1.2.x is not on conda — bioconda ships 1.1.2 (2011), whose
scoring function differs.

**Native Windows** (no WSL, no Docker): same steps with `environment-windows.yml`, from the
**Miniforge Prompt**, and `scripts\get_vina.ps1`. No cavity detection — centre the box on the
co-crystallised control.

### Optional engines

```bash
bash scripts/get_gnina.sh          # second score, NVIDIA GPU, ~4.5 GB
export POLISCREEN_GNINA=$HOME/poliscreen_tools/gnina-run

bash scripts/get_adcp.sh           # peptide docking, 5-20 residues
python scripts/parche_adfr.py      # required: fixes an angle computation that crashes
                                   # ADFRsuite preparation on some targets
```

**ADMET-AI** needs its own environment, because it pulls torch and cannot share the docking one.
PoliScreen talks to it as a subprocess:

```bash
python -m venv <env>
<env>/bin/pip install admet-ai "admelab @ git+https://github.com/DiegoAnyG/admelab"
export POLISCREEN_ADME_PYTHON=<env>/bin/python
export POLISCREEN_ADME_ROOT=<env>/lib/python3.*/site-packages
```

**OPSIN** verifies IUPAC names and needs a JVM, so it is off in packaged builds. Point
`POLISCREEN_OPSIN` at `opsin.jar` with Java on `PATH`.

---

## Verify

```bash
poliscreen info                    # vina, obabel, obrms, plip, fpocket, ADMET, gnina
poliscreen fingerprint <project>   # per-stage hashes and versions, to compare two machines
```

Inside one image, two runs of the same configuration give the same result. Across different
installs they will not: generating the ligand's 3D coordinates, converting to PDBQT and docking all
do floating-point arithmetic in libraries compiled separately for each platform. That is why the
container is the reference and its digest is what gets cited.

---

## Troubleshooting

| Symptom | Fix |
|---|---|
| `Docker is not responding` | Start Docker Desktop and wait for the whale icon to stop animating. |
| `denied` / `unauthorized` pulling the image | The package is private. Make it public once, in its GitHub package settings. |
| `docker` not found inside WSL | Usually **the engine is not running**, not the integration. Check `/var/run/docker.sock` and that Docker Desktop is started; only then Settings → Resources → WSL Integration. |
| Docker Desktop: `Virtualization support not detected` | Enable Intel VT-x / AMD SVM in BIOS/UEFI, then `wsl --install` and reboot. |
| `docker-credential-desktop.exe: exec format error` | `cp ~/.docker/config.json ~/.docker/config.json.bak && echo '{}' > ~/.docker/config.json` |
| `Cannot find the vina executable` | `bash scripts/get_vina.sh`, or set `POLISCREEN_VINA`. |
| `PLIP is not installed` / `fpocket not installed` | `conda install -c conda-forge plip fpocket`. Without fpocket, centre the box on the control. |
| The control produces no poses | Extract it with its SMILES in step 1 — a fragmented ligand gives a PDBQT that Vina rejects. |
| `AttributeError` on a method that exists | Streamlit does not re-import loaded modules. Restart the process; reloading the page is not enough. |
| `The term 'conda' is not recognized` | Use the Miniforge Prompt, or run `conda init powershell` once and open a new window. |
| Windows installer fails in a folder that worked before | Delete `%USERPROFILE%\.conda\environments.txt`; conda rewrites it. Uninstalling does not remove it. |
| Windows installer: `UnicodeDecodeError('charmap', …)` | The path has an accented character. Delete that file **and** install to a path without accents. |
| Windows install stuck at *Setting up the package cache* | ~48,000 small files against live antivirus. Use an SSD or exclude the folder. |
