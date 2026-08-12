# CLAUDE.md

Guidance for Claude Code when working in this repository. These instructions override default
behaviour; follow them exactly.

## What PoliScreen is

Reproducible virtual screening that closes the loop **design → synthesizability filter → docking →
interaction-quality scoring → ADMET**, with a per-cavity objective scoring function and an
orthogonal confidence metric. Author: Diego Cesar Anaya Guerrero. Being published (paper +
congress).

## Working notes

Longer-lived context — objectives, decisions already made and their cost, and a log of what
changed and why — is kept **outside this repository**, in a `poliscreen_notes` folder beside the
checkout. Read it before proposing work, and whenever a question feels like one that has already
been settled. If it is out of date, correct it there first, then act.

## Hard rules

1. **Never bend the code toward a result.** Never write code so that a chosen answer comes out.
   Correcting a fault, optimising, adding a method, or changing a parameter for a stated
   methodological reason is encouraged. Seeding a stochastic step is allowed: it fixes *which*
   structure is sampled, never *what score it earns*. The line: a change may decide what the
   method is, never what the answer is. If it would be embarrassing to describe in the Methods
   section, it is the wrong change.
2. **Never leak anything personal.** No absolute paths from a personal machine, no usernames, no
   machine names, no credentials — in code, comments, tests, docs or commit messages. The one
   exception is the author's name, for citation and authorship.
3. **Never hardcode.** No hardcoded paths or lists. Discover from disk; honour
   `POLISCREEN_PROJECTS`, `POLISCREEN_VINA`, `POLISCREEN_GNINA`, `POLISCREEN_OPSIN`,
   `POLISCREEN_ADME_PYTHON`, `POLISCREEN_ADME_ROOT`. A path that exists only on the machine the
   code was written on does not error — it silently reports that everything is absent.
4. **English, always** — code, comments, docstrings, identifiers, on-disk names, tests, docs and
   commit messages.
5. **Keep the workspace clean.** Announce every new directory: where, what for, and whether it can
   be deleted. Temporary work goes outside the repository.

## Deliberately Spanish, do not translate

The Spanish half of the `ui/i18n.py` catalogue; admelab's own keys where they are still Spanish;
the user-CSV header aliases that real input files use (`nombre`, `compuesto`); the `LEGACY_*` maps
and legacy on-disk names; the author's name.

## Environment and running

Conda env **`cribado`** (Python 3.11), editable install (`pip install -e .`).

```bash
conda activate cribado
python -m pytest tests/ -q      # the safety net — currently 119 tests, keep it green
poliscreen info                 # external-tool check: vina, obabel, obrms, plip, fpocket
poliscreen ui                   # Streamlit interface on 127.0.0.1:8501
poliscreen run --help           # full cycle: design, docking, interactions, ranking
```

External tools are **not** Python packages: `vina` (via `scripts/get_vina.*`), `obabel`/`obrms`
(OpenBabel), `plip`, `fpocket`, and the optional `adcp` and `gnina`. `poliscreen info` reports
which are present.

## Architecture

- **`core/`** — the engine, UI-agnostic. Everything the interface does is callable here.
  - `pipeline.py` orchestrates a run; `docking.py` (Vina), `adcp.py` (peptides), `pockets.py`
    (fpocket), `interactions.py` (PLIP), `screening.py` (scoring + `display_name`,
    `normalize_columns`), `validation.py` (redocking), `report.py` (methods section).
  - `receptor.py`, `ligands.py`, `peptides.py`, `reactions.py`, `reagents.py` prepare inputs.
  - `design.py` is the **admelab bridge** (analogue design + ADMET, run out-of-process);
    `_admelab_runner.py` is what it invokes; `naming.py` verifies IUPAC names with OPSIN.
  - `layout.py` (artifact file names), `session.py` (projects, paths, sessions), `viewer.py`
    (3D viewer + `assets/`).
- **`ui/streamlit_app.py`** — the interface; `ui/i18n.py` the EN/ES catalog; `ui/ayuda.py` the help.
- **`cli.py`** — every capability, scriptable. First-party surface for reproducibility.

## Conventions — do not break these

- **Back-compat is mandatory.** Older projects must keep opening: add to `LEGACY` /
  `layout.artifact()` / `screening.normalize_columns()` rather than renaming in place.
- **Scoring philosophy**: more contacts ≠ better. Bounded similarity to the control, catalytic
  residues weighted, confidence orthogonal to affinity. Do not add "bigger number wins" metrics.
- **Git**: do **not** add a `Co-Authored-By: Claude` trailer, and never list an AI as an author in
  `pyproject.toml` or a `CITATION`. (Institutional reason; the PI is not familiar with AI.)
- **Verify, don't assume.** Changes that only surface at runtime (Streamlit widget keys, packaging,
  installer layout, Windows paths) get a test — see `test_no_shadowing.py`, `test_packaging.py`,
  `test_windows_paths.py`, `test_installer_layout.py`, `test_reproducible_pdbqt.py`.

## Reproducibility

Every stage that touches coordinates must be deterministic, and each has a test. Three independent
sources of randomness were found in this pipeline, each of which alone changed the final ranking:

1. **OpenMM places added hydrogens at random positions** and relaxes them with a 50-step
   minimisation, from Python's global RNG. Seeded in `receptor.prepare`, with the minimisation
   pinned to the `Reference` platform — the multi-threaded ones sum in a different order each run.
2. **OpenBabel `-p` strips every hydrogen and re-adds them at random positions.** This undid the
   above, one step later, on the way to pdbqt. The flag now goes only to a receptor arriving
   without hydrogens. Note the pdbqt is **cached**: an old project keeps its old file.
3. **Vina** is run with a fixed seed and one CPU; it is not deterministic multi-threaded.

## Target platforms

**Windows** is the audience, by two routes: the one-click installer, and Docker (not yet built).
**Linux** natively and in the same container — it is what PoliScreen is developed on.
**macOS is out of scope**; the recipe still carries its selectors, leave them, do not extend them.

## admelab (ADMET / analogue design)

Invoked as a subprocess, never imported. What conflicts with the docking deps is **ADMET-AI**
(torch), not admelab: everything else in it is RDKit. So the installer **ships admelab without
ADMET-AI** (wheel built in the release workflow, `pip install --no-deps` in `post_install`), and
the ML predictions stay an opt-in second environment pointed at with `POLISCREEN_ADME_PYTHON` /
`POLISCREEN_ADME_ROOT`. `design._installed()` finds a copy sitting in PoliScreen's own environment.
The screening works without any of it. The **applicability domain** needs admelab **≥ 0.3**
(`domain` module); `design.AdmelabBridge.has_applicability()` gates it.

admelab renamed its verdict keys and vocabulary to English (`fischer_viability`: `good`,
`moderate`, `unfavorable (…)`, `difficult (…)`). Both vocabularies are accepted, so an older
admelab keeps working. Reading only the old key made every product read `unknown` **and** count as
synthesizable, which silently changed what got docked.

## Installer

`installer/construct.yaml` (conda constructor) builds a one-click installer per platform; the
`.github/workflows/installer.yml` matrix builds them and, on a tag, attaches them to the release.

- fpocket for Windows is cross-compiled in that workflow from `contrib/fpocket-windows/`
  (conda-forge has no win-64 build), and **Vina is vendored there too on Windows** (SHA256 verified
  in CI, Apache-2.0, licence shipped beside it): downloading it from `post_install` on the user's
  machine meant one refused request left an install with no docking engine, announced only by a
  warning that scrolls past. Linux still downloads it at install time.
- **OPSIN is not shipped**: it needs a JVM, about 120 MB per installer and 500 MB unpacked, to
  decide a label and never a number. `POLISCREEN_OPSIN` plus a Java on `PATH` re-enables it.
- **Actions artifact storage is 500 MB** and one build leaves about 1.5 GB, so `upload-artifact` is
  skipped on tag builds — the release gets them instead, and release assets do not count. Deleting
  artifacts does not free the quota immediately; GitHub recalculates every 6-12 hours.

Two hard-won fpocket facts: **never hand fpocket a path at all** — run it from the folder holding
the copy with a bare file name, because it takes the path apart itself and on Windows both halves
fail (backslashes → access violation; forward slashes → it mkdirs the drive letter and returns
silently, exit 0, having written nothing); and on Windows never rewrite `C:\...` to `/mnt/c/...`.

## Pending / roadmap

- **Docker route** for Windows, plus a setup script that installs what is needed and turns on
  Docker Desktop's WSL integration. Idempotent, and it must say what it is about to change.
- **Clean interrupt in the container**: `Ctrl+C` must stop everything and leave no orphans, the way
  the Windows launcher already does.
- Publish the **release tag** so the installers reach users.
- Confirm **Detect pockets** on Windows through the installed app.
- `m_mkdir_p` in `contrib/fpocket-windows/0001-*.patch` (PR #184) still aborts when an intermediate
  component cannot be created — the bare drive letter, whose result depends on whether the process
  has a current directory on that drive. `mkdir -p` never behaved that way: only the final mkdir
  should decide. PoliScreen no longer reaches that code, so this is for the PR, not for us.
- fpocket **PR #184** on `Discngine/fpocket`: merge → upstream release → conda-forge win-64, after
  which the `# [not win]` selector in `construct.yaml` can go.
- Bump the deprecated Node-20 GitHub Actions (`checkout@v4`, `upload-artifact@v4`,
  `setup-miniconda@v3`).
- Surface `fraction_in_domain` in the UI and the report (it already reaches the ranking table).
- Clean up the `tempfile.mkdtemp()` leak in `adcp.py` (diagnostic path, Linux-only).
