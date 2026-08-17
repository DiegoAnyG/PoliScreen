# Changelog

Notable changes. Dates are release dates; the format follows [Keep a Changelog](https://keepachangelog.com).

## [1.0.1] — 2026-08-17

Polish on the container route, which is the one people install.

### Interface

- **Dark by default.** The theme is fixed at start, in `.streamlit/config.toml`. That file is read
  only when Streamlit starts with its directory as the working directory, which the container does
  and the Windows installer does not -- the reason the same build had an appearance switcher on one
  and not the other. `POLISCREEN_THEME=light` overrides it without editing anything.
- The 3D viewer fills its panel. It was reserving up to 300 px for controls that are one row tall,
  leaving a third of the panel empty.
- The results table and the summary panel no longer disagree. The summary read `ranking.csv` from
  disk, written with the weights of the run, while the table recomputed with the weights on screen;
  moving a slider gave one compound two effectiveness values and nothing said why.
- Default exhaustiveness back to 8, with the trade stated in the help text: 8 to explore, 24 for a
  run whose ranking will be reported.
- The IUPAC naming notice is one line. Explanations that long belong in Help.

### Distribution

- `PoliScreen-Docker.bat` is attached to the release, opens with a coloured banner naming the
  version, checks what it needs before it needs it, and looks for a newer image on every start --
  offline it uses the copy on disk and says so.
- `docs/REQUIREMENTS.md`: everything needed, including the obvious, and what is deliberately not.

### Repository

- The test suite runs in CI on every push and pull request, from the same environment file users
  install from.
- `CONTRIBUTING.md`, this changelog, and a Sponsor button.

### Fixed

- Fused complexes carried CONECT records naming receptor atoms, because a pose numbers its atoms
  from 1 and so does the receptor. Any viewer opening those files believed the bonds. Measured on
  20 complexes, correcting them changes nothing PLIP reports: a data fix, not a scoring change.

## [1.0.0] — 2026-08-17

First public release. The full cycle runs and has been used on real targets.

### Screening

- Design from a lead, filtered by real reaction feasibility, so what is docked can be made.
- Ligands from four sources that add up: reaction products, approved drugs from ChEMBL under
  property filters, enumerated peptides, and prepared files.
- Docking with AutoDock Vina 1.2.5; peptides of 5-20 residues routed to ADCP automatically.
- Interactions profiled with PLIP; cavities detected with fpocket; optional gnina re-scoring.
- Per-cavity objective score, normalised against the crystallographic ligand, with a confidence
  metric built from independent evidence and orthogonal to the score.
- ADMET and applicability domain through admelab, run out of process.

### Reproducibility

- Fixed seed, one thread per docking, and an exported Methods file with every parameter and version.
- The control's bond orders come from the PDB chemical component dictionary rather than being
  inferred from geometry, which gave different molecules on different machines.
- `poliscreen fingerprint` hashes each stage in pipeline order so two machines can be compared;
  line endings are normalised, since CRLF alone made identical files look different.
- Interface default exhaustiveness raised to 24, matching the CLI: at 8, one seed in five landed
  0.29 kcal/mol away from the rest, enough to move the top of the ranking.
- A published container image is the reproducible reference; cite its digest.

### Distribution

- Container image on GHCR, and a one-click `PoliScreen-Docker.bat` for Windows.
- One-click installers for Windows and Linux, with fpocket cross-compiled for Windows.
- conda for development.

### Known limits

- Rigid receptor, non-covalent docking only, Ki estimated from the docking score.
- Two different installs do not produce identical numbers: generating the ligand's 3D coordinates,
  converting to PDBQT and docking each do floating-point arithmetic in libraries compiled per
  platform. Inside one image they agree.
