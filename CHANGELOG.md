# Changelog

Notable changes. Dates are release dates; the format follows [Keep a Changelog](https://keepachangelog.com).

## [1.0.3] — 2026-08-17

### Fixed

- **The installers did not build for 1.0.1 or 1.0.2.** `installer/construct.yaml` names the wheel
  by file name and still asked for `poliscreen-1.0.0-py3-none-any.whl`, so both tag builds failed
  with a FileNotFoundError after the container image had already published. The releases went out
  with an image and no installers, and without the launcher attached.
- **The interface reported 1.0.0 from a 1.0.2 image.** `__init__.py` carried the version as a
  literal and was never bumped; it now reads the installed distribution metadata.
- Tests hold the version identical across pyproject, CITATION.cff, the installer recipe, the
  launcher banner and this changelog. Five hand-written copies with nothing checking them is what
  caused all of the above.

## [1.0.2] — 2026-08-17

### Fixed

- **The Windows launcher did not run.** Two causes, one message. The file was committed with
  LF-only line endings, which cmd.exe handles unreliably: lines split mid-command and each fragment
  came back as *"is not recognized as an internal or external command"*. And the banner used
  block-drawing characters, which need `chcp 65001`; changing the code page part-way through a
  batch file shifts the parser's byte offset, splitting lines the same way.

  The banner is seven-bit ASCII now, the file is CRLF, and `.gitattributes` keeps it that way
  through any clone. The 24-bit colour gradient stays -- it is plain ASCII and it was carrying most
  of the look. Tests fail on a non-ASCII byte, on a bare LF, and if the attributes go missing.

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
