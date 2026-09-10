![PoliScreen](docs/banner.png)

# PoliScreen

[Leer en Español (Guía completa)](README.es.md)

Reproducible virtual screening that closes the loop **design → synthesizability
filter → docking → interaction-quality scoring → ADMET**, with an objective
per-cavity scoring function, orthogonal confidence metric, and multi-objective
Pareto optimization.

> **v1.1.0** — Multi-objective Pareto landscape, interaction footprint polygons,
> transport tunnels (CAVER / CaverDock), and one-click container distribution.

---

## Quick Start

The recommended distribution route is **Docker** (or Linux/WSL). Because docking
and interaction scoring involve floating-point operations across multiple C/Fortran
libraries, running inside the container ensures exact reproducibility across any
computer.

### Windows (No terminal needed)

1. Install and start [Docker Desktop](https://www.docker.com/products/docker-desktop).
2. Download `PoliScreen-Docker.bat` and `PoliScreen.ico` from the [latest release](https://github.com/DiegoAnyG/PoliScreen/releases/latest).
3. Double-click `PoliScreen-Docker.bat`. It will pull the official pre-built image,
   create a desktop shortcut with the official icon, and open `http://localhost:8501`
   in your browser.

### Linux / WSL2 (Single command)

```bash
docker run --rm -it --init -p 127.0.0.1:8501:8501 \
  -v "$PWD/proyectos:/data" ghcr.io/diegoanyg/poliscreen:latest
# then open http://localhost:8501
```

For native conda installation for local development: **[docs/INSTALL.md](docs/INSTALL.md)**.
System requirements: **[docs/REQUIREMENTS.md](docs/REQUIREMENTS.md)**.

---

## 4-Step Workflow

PoliScreen features a modular web interface (available in both **English and Spanish**
under Settings → Language):

1. **Receptors** — Download from the PDB (e.g. `8HTB`) or upload a local `.pdb` file.
   Inspect chains, keep cofactors, and extract the co-crystallized ligand as the control.
   *(Tip: click **Example (8HTB)** to load the worked example immediately).*
2. **Ligands** — Upload your molecules (`.sdf`, `.mol2`, `.smi`), build a series by
   reaction (with synthesizability verification), screen approved drugs from ChEMBL,
   or enumerate peptides.
   *(Tip: click **Example ligands (8HTB)** to load 3 active Cruzain inhibitors and 3 decoys).*
3. **Run** — Configure the search box around the cavity and launch docking (AutoDock Vina
   for small molecules, ADCP for peptides) and optional tunnel calculation (CAVER).
4. **Results** — Explore the ranking through:
   - **Interactive Pareto Frontier**: Multi-objective landscape balancing binding affinity,
     interaction quality, and confidence score, featuring live 2D structure tooltips on hover.
   - **TOP Interaction Leaderboard**: Geometric footprint polygons where vertices and edges
     visually represent specific protein contact residues and bond types.
   - **ADMET Profiler**: Radar plots and endpoint estimations.
   - **Transport Tunnels**: Activation barrier ($E_a$) and binding energy profile.

A detailed tutorial with the scientific rationale of each control is in
**[docs/TUTORIAL.md](docs/TUTORIAL.md)**.

---

## Why another docking front-end

Most screening panels rank strictly by *affinity* and *contact count*. That rewards
promiscuity: a molecule touching many irrelevant residues can outrank one that
anchors exactly where it should. PoliScreen is built on four core principles:

1. **Objective per-cavity scoring, not similarity to the control.** Interaction
   quality is the sum of each contact weighted by **bond type** (salt bridge >
   H-bond > pi-stacking > hydrophobic) and by **residue role** (catalytic, secondary,
   cavity, external). It is normalized against the fingerprint of the
   **crystallographic ligand in its real pose** — so a compound beats the
   control by making **more and better productive contacts**, not by copying it.

2. **Confidence metric, orthogonal to the score.** `confidence` (0–1) is the
   **geometric** mean of binding-mode convergence across poses, affinity-interaction
   agreement, and (if enabled) Vina-neural-network consensus, attenuated when the
   control fails to reproduce its crystallographic pose. It quantifies **how much
   to trust** a result, not its magnitude. A high score with low confidence is a
   red flag.

3. **Synthesizability from the design stage.** Analogues are filtered by real
   reaction feasibility (regioselectivity, OH classification, steric hindrance),
   so what gets docked is what a chemist can actually make.

4. **Explicit reproducibility.** Fixed seed, single thread per docking task, and a
   Methods export with every parameter and version. Inside one container image,
   two runs of the same configuration give the same result — and the image digest
   is what a paper cites. `poliscreen fingerprint` hashes each stage so two machines
   can be compared.

The co-crystallised control defines the reference everything is measured against,
so its chemistry is not guessed: bond orders come from the **PDB's own chemical
component dictionary**, looked up from the ligand code. A single interaction engine
(**PLIP**) feeds both the table and the diagrams: same author numbering, same bond
types, no mismatch between what is measured and what is drawn.

---

## Engines

- **AutoDock Vina 1.2.5**: Small molecule docking engine.
- **AutoDock CrankPep (ADCP)**: Peptide docking (5–20 residues).
- **PLIP**: Macromolecular interaction profiler (hydrogen bonds, salt bridges, pi-stacking, hydrophobic contacts).
- **fpocket**: Cavity and binding pocket detection.
- **CAVER 3.0.2**: Transport tunnel detection (included in the container image).
- **CaverDock**: Transport trajectory and energy barrier profiling (via caver-translate).
- **gnina**: Optional deep learning neural-network re-scoring (GPU-accelerated).
- **ADMET-AI**: Molecular property and pharmacokinetic predictions.

---

## Limitations

- **Rigid receptor**: No side-chain flexibility during standard Vina docking.
- **Non-covalent docking**: Vina and compatible engines do not model covalent bond formation.
- **ADMET predictions are estimates**: For prioritization and filtering, not experimental values.

---

## Citation

See **[CITATION.cff](CITATION.cff)** (GitHub renders it under *Cite this repository*).
Cite also the underlying tools PoliScreen invokes; the full bibliography is provided
in the app's *How to cite* panel and in the exported Methods file.

## License

**GNU GPL v3 or later** (see [LICENSE](LICENSE)). The scientific tools PoliScreen
invokes (Vina, ADCP, PLIP, RDKit, Open Babel, fpocket, gnina, CAVER) are independent
programs keeping their own respective licenses; PoliScreen does not incorporate their
source code.

## Author

**Diego Cesar Anaya Guerrero**
