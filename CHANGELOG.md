# Changelog

Notable changes. Dates are release dates; the format follows [Keep a Changelog](https://keepachangelog.com).

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
