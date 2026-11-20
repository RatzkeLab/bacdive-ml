# 08 Models

## Evaluation protocol (applies to every model)

- **Folds:** 5-fold CV, **grouped by genus** (`make_folds` in s12). Genera are shuffled with a
  fixed seed, then assigned greedily by size. The fold table is
  `data/tables/cv_folds_genus.csv`, and all models use the same folds. `--split family|strain`
  exist for comparison.
- **Metric: AUROC within each fold, averaged weighted by fold size** (`fold_auroc`). Folds with
  fewer than 3 positives or 3 negatives for a label are skipped.
  - **Do not pool out-of-fold predictions.** Under taxon-grouped folds, label prevalence differs
    between folds. A model trained on fewer positives predicts lower for the fold that holds more
    of them, and pooling turns that into anti-correlation. On random features, pooled AUROC was
    **0.43** for logistic regression, against 0.50 within folds. Pooled values are still written
    (`*_pooled` columns) for reference only.
- **Negative control:** random sparse features (density 0.3) on the real labels, taxonomy and
  folds (`--tag synthetic`). Median AUROC over 351 labels:

  | model | median AUROC |
  |---|---|
  | kNN | 0.504 |
  | logistic regression | 0.498 |
  | MLP | 0.498 |
  | taxonomy (real taxonomy) | 0.623 |

  A separate 12-seed × 12-label simulation gave 0.505 ± 0.006 (genus folds) and 0.508 ± 0.006
  (strain folds). Small labels are noisy, ±0.05 per label, so also read the subset of labels
  with ≥100 strains in each class. Files: `outputs/s12_synthetic_genus_*`, `outputs/s13_synthetic_genus_*`.

## Models

| model | script | idea |
|---|---|---|
| `tax` | s12 | label frequency among training strains in the nearest shared taxon (family → order → class → phylum). The test genus is never in training, so this is the genus-out taxonomy bar. |
| `knn` | s12 | 10 nearest training genomes by Jaccard on family presence, similarity-weighted. Relatedness at genome resolution, no gene → function learning. |
| `logreg` | s12 | per-label L2 logistic regression on binary families, liblinear. C chosen once by inner genus-grouped 3-fold CV on 40 labels inside fold 0's training part. |
| `mlp` | s13 | multi-task: sparse input → 512 → 256 → 351. Masked BCE with per-label pos_weight, feature dropout 0.2, hidden dropout 0.5, AdamW. Early stopping on 15% of training *genera*. `--final` trains on everything for transfer. |

**Reading the comparison:**
- logreg or MLP beating `knn` means the model learns which genes matter, not just who is related.
- `knn` beating `tax` means genome content resolves relatedness better than taxonomy.

## Transfer to lab genomes (s14)

- **Domain shift:** each lab genome's Jaccard to its nearest BacDive genome, compared with the same
  quantity for BacDive test strains under genus-out and strain-out CV.
- **Outputs:** predicted 351-label profiles (MLP and logreg) and the 256-d MLP embedding for every
  lab genome, with clone clusters attached (20261021 concordance) so clone-mates can be collapsed.

## Results: genus-held-out CV (id50_c80, 8,740 strains)

### s12 baselines

C = 0.003, chosen by inner CV; the grid was flat, 0.660–0.662 over 0.003–0.3. Features: 413,525
families in ≥10 genomes. Runtime 2.7 min. Within-fold AUROC:

| labels (min. strains per class) | n | tax | knn | logreg | logreg > tax |
|---|---|---|---|---|---|
| all ≥25 (threshold of the label matrix) | 351 | 0.612 | 0.655 | **0.680** | 87% |
| ≥100 | 126 | 0.642 | 0.693 | **0.743** | 99% |
| ≥300 | 58 | 0.680 | 0.733 | **0.793** | 98% |

- **Significance:** logreg beats knn on 77% of labels; Wilcoxon p = 9e-26 (logreg vs knn) and
  4e-44 (logreg vs tax).
- **Median by group, tax → logreg:**

  | group | tax | logreg |
  |---|---|---|
  | acid/gas | 0.66 | 0.75 |
  | growth/assimilation | 0.61 | 0.67 |
  | hydrolysis | 0.60 | 0.66 |
  | respiration | 0.58 | 0.63 |
  | reduction | 0.80 | 0.87 |

- **Best:** sulfate/thiosulfate/sulfite reduction, D-galacturonate acid 0.97, nitrate respiration
  0.90, L-arabinose acid 0.89, D-glucose acid 0.88, decanoate growth 0.88, nitrite/nitrate
  reduction 0.85.
- **Worst (≥100 per class):** cellulose 0.55, L-aspartate / L-glutamate growth 0.59–0.60, Tween 40,
  raffinose, hippurate ~0.61.
- **Biggest gains over taxonomy:** L-arabinose growth 0.62 → 0.83, malonate 0.59 → 0.76,
  tryptophan 0.63 → 0.80, decanoate 0.74 → 0.88.

**Interpretation:**
- Gene content beats taxonomy, and learning beats relatedness.
- Performance rises steeply with label size, so more or cleaner labels are a lever.
- The ~15% label noise ([05](05_label_noise_and_taxonomy_baseline.md)) caps attainable AUROC well
  below 1.

Files: `outputs/s12_id50_c80_genus_{auroc_per_label.csv,summary.json,oof.csv.gz}`,
`outputs/figures/s12_id50_c80_genus_logreg_vs_tax.png`.

### s15: do the models weight the right genes?

`s15_interpret_families.py` refits logreg on all strains for 12 labels with known genetics. It
names the top-coefficient families by MMseqs2 search of their representative protein against
Swiss-Prot 2026_04. Full tables: `outputs/s15_id50_c80_top_families.md`. The top families are
textbook pathway genes:

| label | top-weighted families |
|---|---|
| L-arabinose (growth and acid) | *araA*, *araB*, *araD*; the oxidative L-arabinose route (arabinose 1-dehydrogenase, *kdaD*); arabinose ABC transport (*ytfQRT*); α-L-arabinofuranosidase |
| nitrate reduction | *narG*, *narH*, *narV*, *narK*, *narX*/*narP*, *napA*, *napC*; molybdopterin synthesis (*mog*) |
| nitrite reduction | denitrification: *nirS*, heme-d1 synthesis *nirFJNGH*, *aniA* (NirK), *norB*, *nosZ* |
| D-xylose | *xylA*, *xylB*, *xylF*, *xylH*, *xylE*, *xylR*, *xynB* |
| L-rhamnose | *rhaA*, *rhaB*, *rhaD*, *rhaM*, *rhaT*, *fucO* |
| D-galacturonate | *exuT*, *uxaB* |
| sulfate reduction | *aprA*, *dsrA*, *dsrB*, *hdrA*. Only ~18 positives, so lineage markers also rank high. |

Decanoate and citrate pick up regulators and clade markers rather than a clean pathway.

### s13 multi-task MLP (genus-held-out)

- **Features:** 56,462 families present in ≥87 genomes (~1%). Rarer families are mostly
  genus-private and would make the first layer about 8× larger.
- **Training:** best epochs per fold were 6, 6, 6, 10, 5, so it overfits quickly. CV took 35.6 min.
  The final model trained on all strains for 6 epochs: `data/models/mlp_id50_c80_all.pt`.

| median AUROC | tax | knn | logreg | **mlp** |
|---|---|---|---|---|
| all 351 labels | 0.612 | 0.655 | 0.680 | 0.672 |
| labels ≥100 per class | 0.642 | 0.693 | 0.743 | 0.719 |

- **vs logreg:** the MLP beats logreg on only 38% of labels (Wilcoxon p = 3e-4).
- **Feature-matched check:** logreg on the same 56k families scores 0.683 (≥100: 0.747). The gap
  comes from the model, not the features, and rare families add nothing to logreg (413k → 0.680).
- **Where the MLP wins:** only respiration/oxidation (Biolog-type redox labels), 0.657 vs 0.631.
  Multi-task sharing helps those.
- **Untuned so far:** more dropout/weight decay, a smaller hidden layer, ensembling, or initializing
  from logreg are all untried.

### CV regime matters more than the model

`s16_summary_figure.py`; the s12 split runs are `--split strain|family`. The taxonomy baseline uses
species → genus → … levels, which only matter for the strain split (the genus split is unchanged).

| median AUROC (all 351 / ≥100 per class) | tax | knn | logreg |
|---|---|---|---|
| strain held out | 0.689 / 0.754 | 0.726 / 0.785 | **0.743 / 0.821** |
| genus held out | 0.612 / 0.642 | 0.655 / 0.693 | **0.680 / 0.743** |
| family held out | 0.583 / 0.616 | 0.630 / 0.665 | **0.653 / 0.719** |

Logreg beats tax on 88% / 87% / 78% of labels (strain / genus / family).

![by split](../outputs/figures/s16_id50_c80_auroc_by_split.png)

### s14 transfer to our 298 genomes

**Domain shift.** Jaccard of family presence to the nearest BacDive genome (10/25/50/75/90% quantiles):

| | 10% | 25% | median | 75% | 90% |
|---|---|---|---|---|---|
| **our lab genomes** | 0.66 | 0.73 | **0.85** | 0.90 | 0.92 |
| BacDive, strain held out | 0.48 | 0.58 | 0.71 | 0.83 | 0.91 |
| BacDive, genus held out (CV regime) | 0.29 | 0.35 | 0.42 | 0.51 | 0.61 |

**All 298 lab genomes are closer to BacDive than the median genus-held-out test strain.** Most of
our collection consists of species with a type strain in BacDive. So the realistic accuracy for our
strains is the strain-held-out row (~0.74 median, ~0.82 for well-populated labels). The genus-out
numbers are a conservative floor.

Nearest BacDive genera of lab genomes:

| genus | n | genus | n |
|---|---|---|---|
| Pseudomonas | 39 | Kluyvera | 15 |
| Klebsiella | 35 | Citrobacter | 14 |
| Brucella/Ochrobactrum | 22 | Chryseobacterium | 12 |
| Enterobacter | 19 | Rahnella | 11 |
| Bacillus | 18 | Microbacterium | 11 |

**IDs.**
- Internally a lab genome is `lab_<gbk stem>`, the GenBank file name in the panX folder
  `.../all_gbk_files_including_CeMbio43/input_GenBank` (s09).
- The stem is the **assembly name, not a 384 well**:
  - `P1_D6` = Or's sequencing plate 1, position D6, which is 384 well C24, strain 110;
  - Schulenburg genomes have NCBI stems `GCF_…_ASM…_genomic`;
  - CeMbio and others keep their isolate names (`BIGb0170`, `MYb…`).
- Since 2026-11-20 all s14 outputs are keyed by `strain_no`, the lab strain number from the
  20261113 rule (unique over the 298 genome wells).
- Chain: gbk stem → `20261021/03_genome_concordance/outputs/g01_genome_index.csv` (gbk →
  assembly_name, well) → `20261113/clean/strain_ids.csv` (well → strain_no).
- **Alignment check:** each lab genome's GC (from its gbk) against its nearest BacDive genome's
  NCBI GC gives r = 0.99, median |Δ| 0.4 percentage points; shuffled pairs give 8.9. So the rows
  belong to the right genomes.

**Nearest genus vs the master collection sheet (s17).**
- **Sheet:** `data_links/resources/Strain collection including Schulenberg's - collection.csv`
  (= `/home/rl/data/resources/karl/`).
  - It is keyed by stock (`Strain number`), so the join goes well → `stock_no` →
    `Strain number`.
  - Genus comes from its `phylogeny` column (a Greengenes-style 16S string).
  - Sanger-vs-genome 16S divergence comes from
    `external_reference_crosscheck/results/per_well_concordance.csv` (`sanger_vs_ngs_pct`).
- **Coverage:** 159 of 298 genomes have a sheet genus.

| verdict | n | meaning |
|---|---|---|
| agree | 101 | sheet genus = BacDive genus (or synonym: Ochrobactrum → Brucella, Raoultella → Klebsiella) |
| same 16S, name differs | 8 | the sheet's own Sanger 16S ≤1.5% from the genome's. Renames (Bacillus → *Priestia*, Ewingella → *Rahnella*, Trabulsiella → *Citrobacter*), 16S limits in Enterobacterales (Serratia → *K. planticola*), and one plain wrong string (164: "Bacillus", genome *Klebsiella*, 16S identical) |
| **different organism** | **44** | the sheet's Sanger 16S is **>1.5% (median 19%) from the genome's 16S**: different orders or phyla |
| no Sanger | 6 | genus differs, no sequence comparison |

- **Where the Sanger data comes from:** the sheet's `16S sequence (NCBI access number)` column.
  Despite the name it holds sequences (355 of 464 rows). They were loaded by
  `external_reference_crosscheck/scripts/01_assemble.py` (well → stock via
  `Full_384_strain_collection_no_seq_info.csv`) and aligned to the genome-derived 16S. Of the 154
  genome strains with both, **54 (35%) differ by >1.5%**.
- **The disagreeing Sanger reads are clean:** median 1,338 bp and no ambiguity codes, the same as
  the agreeing ones. Read quality cannot flag them.
- **Pattern:** a few clones account for many of the different-organism cases.
  - 32 of the 44 are clone-cluster members.
  - **C002** (*Citrobacter braakii*, the G02 invader) sits under 9 strain numbers that the sheet
    lists as 7 different genera.
  - **C004** (*Klebsiella ornithinolytica*) sits under 6 strain numbers listed as 5 genera.
  - **C005** (*K. planticola*) sits under 4.
  - The cases concentrate in strain numbers 51–150 (genus agreement 36–38%, against 65–83%
    elsewhere and 92% for MYb isolates).
- **40 of the 44 are concordance-`ok`:** the well's own reads match the genome. So genome and plate
  agree, and only the sheet differs.
- **Two explanations fit, and these data cannot tell them apart:**
  - **(a) The Sanger record is wrong** (sample or row mix-ups), and the clones are the same
    organism isolated repeatedly. The lab knows the Sanger sequences and the taxonomy derived from
    them have problems, which favours (a).
  - **(b) The record describes the original isolate,** and the stock was overgrown before it was
    genome-sequenced and plated. Earlier notes (external_reference_crosscheck: "original-isolate
    lineage") assumed this.
- **Either way,** concordance `ok` does not certify agreement with the sheet.
- **Phylum agreement:** 84%. The P1_D6 example is strain 110: the sheet says *Rhodococcus*, the
  genome has 55.9% GC and is nearly identical to *Klebsiella ornithinolytica*, J = 0.93, clone C004.
- **Correction:** an earlier version of this page blamed the join. The join is correct (it
  reproduces with the 20261113 ID table). For about a third of genome strains, the sheet's 16S and
  taxonomy do not describe what the genomes and plates contain.
- **Consequence:**
  - The profiles describe the organism actually in the wells and in the experiments, which is what
    the interaction model needs.
  - Never take a lab strain's taxonomy from the master sheet.
  - Treat clone-mates as one organism (`clone_cluster`).

Files: `outputs/s17_lab_vs_sheet_taxonomy.csv` (per genome, with `verdict`),
`outputs/s17_summary.json`.

**Outputs:**

| file | contents |
|---|---|
| `outputs/s14_lab_genome_ids.csv` | `strain_no`, well, isolate_name, genome_id, assembly_name, gbk, clone_cluster, clone_size |
| `outputs/s14_lab_profiles_logreg.csv` | 298 × 351 predicted probabilities, index `strain_no`. **Logreg is the better model.** |
| `outputs/s14_lab_profiles_mlp.csv` | same, from the MLP |
| `outputs/s14_lab_embeddings_mlp.csv` | 298 × 256 hidden layer, index `strain_no` |
| `outputs/s14_lab_nearest_bacdive.csv` | nearest BacDive genome per lab genome (strain_no, well, genome_id, clone cluster) |

![domain shift](../outputs/figures/s14_domain_shift.png)
