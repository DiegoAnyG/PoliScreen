# PoliScreen — Tutorial

The four-step cycle and the reasoning behind each control. Mirrors the app's **Help** menu.

**Receptors → Ligands → Run → Results.** The right-hand panel follows the active step.

What distinguishes the scoring: it measures *which* contacts are made and *with which residues*,
not affinity and contact count, and attaches a measure of how much to trust each result.

**Project folder.** Everything is written there — poses, complexes, PLIP XML, tables. Changing it
changes the analysis; anything already prepared in it is detected on entry. The path is a Linux
path; a pasted Windows path is translated and you are told which was used. In Docker, only the
mounted volume survives the container.

**Sessions.** *File › Save session* packs the analysis into a `.poliscreen` file. Restoring it
brings back tables, receptors and ligands, and lets you change the weighting **without re-running
the docking**.

**Downloads.** *File › Download results* builds one ZIP. The minimal reproducible set: Methods
(only exists once exported), ranking, interaction matrix, ligand table, redocking validation, and
the input receptors and ligands.

---

## 1. Receptors

Preparation removes waters, adds hydrogens, and **keeps the original residue numbering** — if that
changes, an active-site residue is identified with another's name and the pharmacophoric analysis
is wrong.

**The co-crystallised control is the most important part of the setup.** It defines the reference
everything is measured against, marks the real site, and enables redocking validation. Give its
**SMILES** when extracting: PDB stores no bond orders, and without that template heterocycles such
as a benzofuroxan *N*-oxide are read incorrectly.

**Cofactors.** Keep one that is part of the site (NADP in a reductase): the ligand competes or
cooperates with it and removing it changes the pocket. Ions and crystallisation molecules that do
not participate can go.

**Waters** are removed because rigid docking does not model them or the cost of displacing them.
The exception is conserved structural waters that mediate binding — a case-by-case decision, after
which PoliScreen detects water-mediated bridges automatically.

## 2. Ligands

Four sources, and they add up: one run can hold compounds from several, and the results table says
which is which.

- **Build by reaction** — a core plus a reagent library, filtered by real chemical feasibility, so
  what is docked is what can be synthesised.
- **Screen approved drugs** — compounds already approved as medicines, from ChEMBL, filtered by
  property (Lipinski, Veber, or any of MW/LogP/TPSA/HBD/HBA/RotB). Nothing is designed here, so
  there is no synthesis feasibility to judge. The library is cached in the project, so the run
  records which library it screened.
- **Generate peptides** — sequences enumerated under composition and property rules.
- **Upload ready ligands** — prepared structures, whose chemistry is read for ADMET and descriptors.

### Peptides

Filters: alphabet restricted by class, fixed prefix or suffix, repeats forbidden, consecutive
residues limited. For antimicrobials the discriminating properties are **net positive charge** (the
bacterial membrane is anionic) and moderate hydrophobicity.

- **Net charge** at pH 7.4 — most associated with antimicrobial activity.
- **GRAVY** — mean hydropathy.
- **Hydrophobic moment** — amphipathicity; whether, folded into a helix, hydrophobic residues sit
  on one face. This is what allows membrane insertion.
- **Boman index** — tendency to bind other proteins; above 2.5 kcal/mol is promiscuous.

Termini: *C-terminal amidation* removes the terminal negative charge (+1 net, usually more active);
*N-terminal acetylation* protects against aminopeptidases; *head-to-tail cyclisation* rigidifies —
fewer degrees of freedom, protease resistance, more reliable docking.

Peptides go to **ADCP**, not Vina, because Vina treats the ligand as a torsion tree and the
sampling stops covering the space (saFtsZ, 23 Å box, exhaustiveness 8, one thread):

| Residues | Rotatable bonds | Vina |
|---|---|---|
| 3 | 15 | ~98 s |
| 5 | 23 | > 2 min |
| 10 | 43 | does not finish |

Without ADCP they fall back to Vina with a warning. Peptide docking ranks candidates rather than
asserting a binding mode.

## 3. Run

**The box** is most reliable centred on the co-crystallised ligand, which marks the real site.
Otherwise the geometric centre, or a cofactor.

**Cavities.** `Cavity` is the pocket's real extent; `Box` is the search region assigned to it, with
a 14 Å minimum (below that a ligand does not fit), marked `*` when it applies. **Druggability**
estimates whether the pocket can bind a small molecule at all.

**Hybrid docking** screens the same compounds against several pockets of one receptor, each with
its own ranking and its own reference — site-selectivity information a single-site screen cannot
give.

**Parameters.** *Exhaustiveness* — how much the search explores. *Poses per ligand* — below 3 the
confidence metric loses resolution. *Energy range* — window for reporting alternative modes. *Seed
and one thread* — Vina is not deterministic with more than one thread per docking.

**gnina** re-scores poses Vina already found, with a network trained on crystallographic complexes.
The value is independence, not speed: Vina's score is empirical, gnina's is learned, and agreement
between two different assumptions is evidence neither gives alone. It cannot recover a pose the
sampling missed.

## 4. Results

**Effectiveness** — percentage against that site's reference ligand, which sits at 100 %. It is not
similarity to the control: each contact is weighted by bond type (salt bridge > H-bond > π >
hydrophobic) and residue role (catalytic, secondary, cavity, external). A compound beats the
control by making **more and better productive contacts**, not by copying it.

**Confidence** — orthogonal to effectiveness: not how good the compound is, but how much to trust
the number. A geometric mean, so one failing piece of evidence pulls it down: binding-mode
convergence, affinity↔interaction agreement, and Vina vs. network consensus when re-scoring ran. It
falls if the target fails its redocking. **High effectiveness with low confidence is an alarm.**

**Affinity.** *best_dock* — kcal/mol, more negative is better. *pKi* — sortable, but does **not**
enter the score: it derives from the docking score, so scoring it would count the same evidence
twice. *LE* and *LLE* reward binding per atom, guarding against size bias.

**Catalytic and secondary residues.** Catalytic are mandatory — not touching them penalises through
the catalytic gate — and are a property of the enzyme, determined experimentally. Secondary are
known anchors that add more than an ordinary contact without being required. Where the catalytic
site is unknown, PoliScreen suggests the residues the crystallographic ligand makes directional
interactions with: a starting point, not an answer. **This list moves the ranking more than any
weight.**

**Weighting** auto-normalises; all axes at 1.0 is a simple average. An axis with weight but no data
(toxicity without ADMET) is ignored and flagged, so the Methods section does not claim something
that never intervened.

**Percentage vs. percentile.** The percentage is measured against that target's control, so it
depends on how strong that control is and does not compare across targets. The percentile places
the compound within its own library and does.

## 5. Reproducibility

Fixed seed, one thread per docking, pinned environment, Vina pinned by SHA256. Two runs with the
same configuration give the same result.

`poliscreen fingerprint <project>` hashes each stage and lists tool versions — run it on two
machines and diff to find where they diverge.

*File › Export Methods* writes the parameters, box, weights, reference and exact versions for a
paper.

**Declared limitations.** Rigid docking (the receptor does not change conformation); non-covalent
only; Ki estimated from Vina's affinity, informative rather than measured; LD50 and toxicity from a
single model and tending optimistic — contrast with a second predictor before claiming low
toxicity.
