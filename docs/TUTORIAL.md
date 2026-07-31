# PoliScreen — Tutorial

A walkthrough of the four-step cycle and the scientific rationale behind each
control. The panels below mirror the app's **Help** menu.

> Images will be added per section. Placeholders point to `docs/img/`.

---

## 1. Getting started

**What PoliScreen does.** It chains the full screening cycle: prepare the target,
obtain compounds, dock them, and evaluate binding quality. The difference from an
ordinary docking panel is *how it scores* — instead of rewarding affinity and
contact count, it measures **which contacts** are made and **with which
residues**, and adds a measure of **how much to trust** each result.

**Order of work.**

1. **Receptors** — download or upload the structure, choose what to keep, and
   extract the co-crystallized ligand as the control.
2. **Ligands** — build the series by reaction, generate peptides, or upload ready
   structures.
3. **Run** — define where to search and launch docking.
4. **Results** — adjust the weighting and examine the ranking.

The right-hand panel always shows what corresponds to the active step.

![Workspace overview](img/overview.png)

**Save and restore.** *File › Save session* packs the analysis into a
`.poliscreen` file. Restoring it brings back tables, receptors and ligands, and
lets you **change the weighting without re-running the docking**. The light
session is a few megabytes; the full one adds poses and complexes.

**The project folder.** Everything is written there: poses, complexes, PLIP XML
and tables. Changing folder changes analysis, and anything already prepared in it
(receptors, controls, ligands) is auto-detected on entry. PoliScreen runs inside
Linux, so the path is a Linux path (`/home/user/...`); a pasted Windows path is
translated automatically and you are told which one was used.

**What to download.** *File › Download results* builds a single ZIP with what you
select, leaving no copies in the project folder. Most of it already lives in that
folder; downloading only makes sense to move the analysis elsewhere, attach it to
a manuscript, or archive it. The **recommended** items form the minimal
reproducible package: the **Methods** section (the only item that does not exist
until exported), the **ranking**, **interaction matrix** and **ligand table**,
the **redocking validation**, and the **input receptors and ligands**.

---

## 2. Receptors

**What preparation does.** Removes waters, adds hydrogens, and keeps the original
residue numbering. Preserving the numbering matters: if it changes, an active-site
residue can be identified with another's name and the whole pharmacophoric
analysis is wrong.

![Receptor preparation](img/receptor.png)

**Cofactors — when to keep them.** A cofactor that is part of the site (e.g. NADP
in a reductase) must be kept: the ligand competes or cooperates with it, and
removing it changes the pocket shape. Ions and crystallization molecules that do
not participate can be removed.

**The co-crystallized control.** The most important piece of the setup. It defines
the reference everything is measured against, marks the real binding site, and
enables protocol validation by redocking. When extracting it, provide its
**SMILES**: PDB does not store bond orders, and without that template heterocycles
such as the benzofuroxan *N*-oxide are read incorrectly.

**Why waters are removed.** Standard practice in rigid docking. AutoDock Vina does
not model explicit waters or the cost of displacing them, so leaving them produces
artifacts. The exception is **conserved structural waters** that mediate binding;
keeping them is a case-by-case decision, and PoliScreen then detects
water-mediated bridges automatically.

---

## 3. Ligands

**Three routes.** *Build by reaction* — start from a core and a reagent library,
filtered by real chemical feasibility, so what is docked is what can be
synthesized. *Generate peptides* — enumerate sequences under composition and
property rules. *Upload ready ligands* — prepared structures, whose chemistry is
read to compute ADMET and descriptors.

![Ligand design](img/ligands.png)

**Peptides — the rules.** The alphabet is restricted by class (hydrophobic,
cationic, aromatic…); a prefix or suffix can be fixed, repeats forbidden, or
consecutive residues limited. The physicochemical filters discriminate most in
antimicrobial peptides: **net positive charge** (the bacterial membrane is
anionic) and moderate hydrophobicity.

**Peptides — the descriptors.**

- **Net charge** at pH 7.4 — the trait most associated with antimicrobial activity.
- **GRAVY** — mean hydropathy; positive indicates overall hydrophobic character.
- **Hydrophobic moment** — amphipathicity: whether, folded into a helix,
  hydrophobic residues sit on one face and polar ones on the other. This is what
  enables membrane insertion.
- **Boman index** — tendency to bind other proteins; above 2.5 kcal/mol is
  considered promiscuous.

**Peptides — the termini.** *C-terminal amidation* removes the terminal negative
charge and adds +1 to the net charge, usually increasing antimicrobial activity.
*N-terminal acetylation* protects against aminopeptidases. *Head-to-tail
cyclization* rigidifies the peptide: fewer degrees of freedom, protease
resistance, and more reliable docking.

**Limit of peptide docking.** Measured on saFtsZ with AutoDock Vina (23 Å box,
exhaustiveness 8, one thread):

| Residues | Rotatable bonds | Time |
|---|---|---|
| 3 | 15 | ~98 s |
| 5 | 23 | > 2 min |
| 10 | 43 | does not finish |

This is a Vina limitation, not PoliScreen's: it treats the ligand as an
independent torsion tree, and with many bonds the sampling stops covering the
space. **For peptides, docking serves to rank candidates, not to propose a
binding mode.**

---

## 4. Run

**The search box.** Most reliable when centered on the co-crystallized ligand: it
marks the real site. The protein's geometric center or a cofactor point elsewhere.
The X, Y and Z axes are drawn in the viewer to show which direction each control
moves.

![Search box](img/run.png)

**Detected cavities.** `Cavity` is the pocket's real extent; `Box` is the search
region assigned to it, with a **14 Å minimum** because below that a ligand would
not fit. When that minimum applies it is marked with `*`, so cavities of different
volume can share the same box size. **Druggability** estimates whether the pocket
has suitable shape and chemistry to bind a small molecule.

**Hybrid docking.** Docks the same compounds into **several pockets** of the same
receptor, each with its own ranking — **site-selectivity** information that a
single-site screen does not give. Each site uses its own reference: the control if
present, a cofactor if inside the box, or the catalytic residues you designate.

**Docking parameters.** *Exhaustiveness* — how much the search explores; higher is
finer and slower. *Poses per ligand* — how many binding modes are kept; below 3
the confidence metric loses resolution. *Energy range* — window relative to the
best pose for reporting alternative modes. *Seed and one thread* — guarantee that
two identical runs give the same result; with more than one thread per docking,
Vina stops being deterministic.

**Second scoring with gnina.** Re-evaluates with a **neural network** the poses
Vina already generated, without repeating the search. Its value is not speed but
independence: Vina's score is empirical, gnina's is learned from crystallographic
complexes. Agreement between two methods with different assumptions is evidence
neither gives alone, and it feeds the confidence metric. *Limitation*: the network
can only judge poses Vina found; if sampling missed the correct pose, re-scoring
does not recover it.

---

## 5. Results

**Effectiveness.** Percentage relative to that site's reference ligand, which sits
at 100 %. Interaction quality does **not** measure similarity to the control: it
sums each contact weighted by bond type (salt bridge > H-bond > π > hydrophobic)
and residue role (catalytic, secondary, cavity, external). A compound beats the
control by making **more and better productive contacts**, not by copying it.

![Results ranking](img/results.png)

**Confidence — what it is and why it differs.** **Orthogonal to effectiveness**:
it measures not how good the compound is, but how much to trust the number. It is
a **geometric** mean — one failing piece of evidence pulls it down — of *conv*
(binding-mode convergence, interaction overlap between the best poses), *conc*
(affinity↔interaction agreement), and *consensus* (Vina vs. the neural network, if
re-scoring was enabled). It drops if the target fails its redocking. **High
effectiveness with low confidence is an alarm**: good score, disagreeing evidence.

**Affinity metrics.** *best_dock* — Vina energy in kcal/mol; more negative is
better. *pKi* — −log₁₀ of the estimated Ki, numeric and sortable; Ki does **not**
enter the score (it derives from the docking score, so scoring it would count the
same thing twice). *LE* and *LLE* reward binding **per atom**, guarding against
size bias.

**Catalytic and secondary residues.** *Catalytic* — mandatory: not touching them
penalizes via the catalytic gate. They are a property of the enzyme, determined
experimentally, not of the ligand. *Secondary* — known pocket anchors that add
more than an ordinary contact but are not required. If you do not know the
catalytic site, PoliScreen suggests the residues with which the crystallographic
ligand makes directional interactions — a starting point, not an answer: **this
list influences the ranking more than any weight**.

**Weighting.** Axis weights **auto-normalize**: they need not sum to 1, and all at
1.0 is a simple average. An axis with weight but **no data** (toxicity without
ADMET) is ignored and flagged, so the Methods section does not declare something
that did not intervene.

**Percentile vs. percentage.** The **percentage** is measured against that
target's control, so it depends on how strong that control is and is not
comparable across targets. The **percentile** places the compound within its own
library and does allow comparison.

---

## 6. Reproducibility

**What is fixed.** Constant seed and one thread per docking; pinned environment and
AutoDock Vina binary, verified by its SHA256 sum. Two runs with the same
configuration give the same result.

**Export the methods.** *File › Export Methods* generates a document with the
parameters, box, weights, reference used, and exact tool versions, ready for a
paper's Methods section.

**Declared limitations.** **Rigid** docking (the receptor does not change
conformation); **non-covalent** (covalent bonds with the target are not modeled);
**estimated Ki** from Vina's affinity (informative, not measured); **LD50 and
toxicity** from a single predictive model, which tend to be optimistic — contrast
with a second predictor before claiming low toxicity.
