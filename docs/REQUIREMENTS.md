# What you need to run PoliScreen

Everything, including the obvious. If a row says **you provide it**, PoliScreen cannot install it
for you.

## The short answer

**A 64-bit computer, Docker Desktop, an internet connection for the first run, and 15 GB free.**
Nothing else. No Python, no conda, no compilers.

---

## 1. Hardware

| | Minimum | Comfortable | Why |
|---|---|---|---|
| CPU | 64-bit, 2 cores | 8 cores | docking is one thread per ligand; more cores means more ligands at once |
| RAM | 8 GB | 16 GB | the ADMET models load into memory beside the docking |
| Free disk | 15 GB | 30 GB | 7 GB image, the rest is poses and complexes, which grow fast |
| GPU | none | NVIDIA, 4 GB | **only** for gnina's second scoring. Everything else is CPU |

An ARM machine (Apple Silicon, Snapdragon) is not supported: AutoDock Vina ships an x86-64 binary.

## 2. Operating system

| | Works | Notes |
|---|---|---|
| Windows 10 (2004+) or 11, 64-bit | yes | needs WSL 2, which Docker Desktop installs |
| Linux, 64-bit | yes | any distribution with a current Docker |
| macOS | not supported | dropped deliberately; nothing is built or tested for it |

Windows needs **virtualization enabled in BIOS/UEFI** (Intel VT-x or AMD SVM). It is usually on
already; if it is not, no installer can switch it on for you — it is a firmware setting.

## 3. Software you provide

**Docker Desktop** (Windows) or **Docker Engine + Compose** (Linux). That is the whole list.

Free for personal use, education and small businesses; a paid licence applies to large companies —
Docker's terms, not ours.

Everything else travels inside the image: Python, RDKit, AutoDock Vina, Open Babel, PLIP, fpocket,
OpenMM, PDBFixer, Streamlit and the ADMET engine. You never install them and you cannot get their
versions wrong.

## 4. Network

| When | What for | Can you avoid it |
|---|---|---|
| First run | downloading the image, about 3 GB compressed | no |
| Downloading a receptor by PDB code | fetching it from RCSB | yes — supply the `.pdb` yourself |
| Extracting a control | its bond orders, from the PDB dictionary | yes — type the SMILES instead |
| Searching reagents or approved drugs | PubChem, ChEMBL | yes — upload your own list |
| Later runs | nothing | the screening itself is entirely offline |

No account, no API key, no registration. Nothing is uploaded: every request is a download.

## 5. A browser

The interface is a web page served on your own machine. Any current browser. Nothing is published
to the network — the port is bound to `127.0.0.1`, so only your computer can reach it.

## 6. What you need to have ready

Not software — inputs:

- **A receptor**: a PDB code, or a `.pdb` file.
- **A control** (strongly recommended): the co-crystallised ligand. It defines the reference every
  compound is measured against, marks the real site, and enables redocking validation.
- **Ligands**, from any of: a SMILES core plus a reagent library, approved drugs from ChEMBL,
  peptide sequences, or your own prepared files.
- **The catalytic residues**, if you know them. If you do not, PoliScreen suggests them from the
  control — as a starting point, not an answer.

## 7. What is optional, and what it costs

| | Adds | Needs |
|---|---|---|
| **gnina**, second scoring | ~4.5 GB | an NVIDIA GPU handed to the container |
| **ADCP**, peptide docking | ~900 MB | accepting the Scripps academic licence |
| **ADMET-AI** | already included | nothing |
| **Tunnel reading** (caver-translate) | already included | nothing |
| **Tunnel search** (CAVER) | ~330 MB | the CAVER zip at build time |
| **Tunnel transport** (CaverDock) | ~490 MB | a Linux install outside the container |

gnina and ADCP are Linux-only builds, which is another reason the container route is the one that
has everything.

### The two tunnel engines, and why only one is in the image

**CAVER** is GPL-3 and written in Java, so it may be redistributed and it behaves the same
everywhere. It is in the image, with the Java runtime. caver.cz serves the zip behind a form rather
than a stable URL, so the build takes it one of two ways:

```bash
docker build --build-arg CAVER_URL=https://.../caver_3.0.2.zip --build-arg CAVER_SHA256=... .
# or put caver_*.zip in installer/vendor/ and build normally
```

Without either, the image still reads tunnel results; it just cannot compute them.

**CaverDock** is not in the image and cannot easily be. It is a native binary built against Ubuntu
20.04's libraries, and the base here is newer: its MPI daemons fail to start. Its Apptainer image
solves that, but opening one inside Docker needs privileges no screening tool should ask for.

So CaverDock runs **outside** the container — on Linux or WSL directly, where PoliScreen finds the
`.sif` by itself or through `POLISCREEN_CAVERDOCK`. Point the Results tab at the folder afterwards.
Finding the tunnels, which is the interactive half at seconds per structure, works in the container
for everyone.

## 8. What you do NOT need

- Python, conda, pip, a compiler, Visual Studio, a C toolchain
- Administrator rights after Docker Desktop is installed
- An account, licence key or registration for PoliScreen
- A cluster, a server, or a GPU (unless you want gnina)
- Any of AutoDock Vina, Open Babel, PLIP or fpocket installed by hand

## 9. Honest limits

- **Two machines will not produce identical numbers unless both run the container.** Ligand 3D
  generation, the PDBQT conversion and the docking engine each do floating-point arithmetic in
  libraries compiled separately per platform. Inside the same image, they agree.
- Rigid docking: the receptor does not change conformation.
- Non-covalent only.
- Ki is estimated from the docking score — informative, not measured.
