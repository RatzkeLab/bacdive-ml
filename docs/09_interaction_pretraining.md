# 09 Does BacDive pretraining help predict interactions?

Script: `scripts/s18_interaction_pretraining.py` (env karl_datasci). Question (2026-11-21): does
an MLP whose genome encoder was **pretrained on BacDive** (genome → 351 metabolic phenotypes,
s13 `--final`) predict pairwise interaction outcomes better than the **same MLP trained from
scratch**, especially for taxa never seen in training (genus- or family-held-out)?

## Data

- **Target:** `mean_log2_ratio_a_over_b` from the published modelling pairs.
  - 20260907: `analyses/20260923/01_model_sweep/outputs/g01_modeling_pairs.csv`
  - 20260630: `analyses/20260710_increased_tolerance/06_ml_analysis/genomic_ml/outputs_merged_ab/`
  - High-uncertainty and unstable pairs are already excluded there.
- **Wells:** wrong-genome wells are dropped using each experiment's own 20261113 verdict
  (`s01_wells_<exp>.csv`, `genome_status != wrong`). This is setting "D" of 20261113.

  | | pairs | strains | organisms | species | genera | families |
  |---|---|---|---|---|---|---|
  | 20260907 | 377 | 222 | 151 | 117 | 42 | 25 |
  | 20260630 | 1,021 | 61 | 46 | 38 | 19 | 13 |

- **Features:** the BacDive MLP's 56,462 protein families (MMseqs2 id50/c80, ≥87 BacDive genomes),
  restricted per experiment to families present in its strains (31,787 / 20,490). The restriction
  is exact: absent families are zero for every strain.
- **Taxonomy comes from genomes only; the master sheet is unreliable (s17):**

  | level | grouping |
  |---|---|
  | organism | clone cluster (`genome_unit`, ≥99.9% ANI) |
  | species | ≥95% ANI connected components (20261021 skani pairs, AF ≥50%). 139 clusters among 298 genomes; natural gap at 92–95% ANI; every cluster maps to one BacDive genus |
  | genus, family | LPSN genus/family of the nearest BacDive genome (s14) |

## CV

- **Folds:** 20261113 `group_folds` at each level. Groups are permuted and assigned round-robin
  into 5 folds. **Test = pairs with both strains in fold k; train = pairs with neither.**
- **Repeats:** 10, each with a different fold partition.
- **Caveat:** family folds are lopsided. Enterobacteriaceae is 70/222 strains (0907) and 31/61
  (0630), so one fold per repeat holds most of the test pairs and trains without that family.

## Variants

All variants use one architecture and one training recipe; only the encoder's initialization and
trainability differ.

- **Encoder** (= s13 body): sparse families → 512 → ReLU/Dropout(0.5) → 256 → ReLU/Dropout(0.5).
- **Pair head:** f([h_a, h_b]) = 512 → 64 → ReLU/Dropout(0.3) → 1, with
  ŷ(a,b) = f(a,b) − f(b,a), antisymmetric by construction.

| variant | encoder init | encoder trained? | strain vector |
|---|---|---|---|
| `naive_lr1e-3`, `naive_lr1e-4` | random (same distribution as the BacDive MLP's init) | yes (encoder lr 1e-3 / 1e-4) | 256-d |
| `pretrained_ft_lr1e-3`, `pretrained_ft_lr1e-4` | BacDive MLP | yes (fine-tuned) | 256-d |
| `pretrained_frozen` | BacDive MLP | no | 256-d embedding, z-scored |
| `pretrained_profiles_frozen` | BacDive MLP | no | 351 predicted phenotypes, z-scored |
| `random_frozen` | random | no | 256-d, z-scored. **Control:** does frozen-pretrained win through pretraining or just through fewer trainable weights? |
| `SHUFFLED_pretrained_frozen` | BacDive MLP | no | training targets permuted (noise floor) |

**Training:**
- Full-batch AdamW: head lr 1e-3, weight decay 1e-4, input dropout 0.2 on trainable encoders.
- **Early stopping:** on 20% of the *training* groups at the same level. Validation pairs are
  training pairs touching a validation group. Patience 40 steps, max 300.
- **Refit** on all training pairs for the chosen number of steps. The test fold is never used
  for any choice.
- **Both encoder learning rates** are reported for the naive and fine-tuned variants, so neither
  side gets a tuning advantage.

**Metrics:**
- `genomic_ml._metrics` per fold: Spearman ρ (primary), R² vs predicting 0, and sign accuracy on
  |y| > 1.
- Means over folds; error bar = SD of the per-repeat mean across the 10 repeats.
- Paired Δρ vs naive on identical folds, with a Wilcoxon test over repeats.

## Smoke-test lessons

- **Frozen embeddings need a floor on the z-score scale.** Dead ReLU units have SD ≈ 0 on the
  training strains, but they fire on test strains. Plain z-scoring gave R² ≈ −10⁷. Fixed: the
  SD is floored at 10% of the mean SD.
- **Build the sparse input tensor once.** Re-coalescing it every step cost ~1 s/step. It is now
  built once, and dropout is applied as a mask on its values: ~0.35 s/step for a trainable
  encoder under load.

- **Forked pool workers came up pinned to one physical core** (CPUs 0,32), although the parent
  reported 0–63. The first full run put all 20 workers on one core, ~9% CPU each, and was stopped.
  Fixed: a pool initializer resets the mask (`_init_worker`). Verified in a test pool: workers
  keep 0–63 through a training step.

## Status

Run on 2026-11-23: 346 fold jobs × 8 variants, 78.5 min on 48 workers.

```
/home/rl/mambaforge/envs/karl_datasci/bin/python -u scripts/s18_interaction_pretraining.py --reps 10 --processes 48
/home/rl/mambaforge/envs/karl_datasci/bin/python -u scripts/s19_pretraining_bootstrap.py --boots 2000 --processes 48
```

## Results

**Pretraining did not detectably help, at any level.** The hoped-for gain on unseen genera did
not appear. Family level cannot be settled with this collection (the intervals are very wide).

Spearman ρ on held-out taxa (mean over folds and 10 repeats). "Naive" and "fine-tuned" are each
the mean of their two encoder learning rates, so no learning rate was picked after seeing test
scores.

| held out | 0907 naive | 0907 fine-tuned | 0907 frozen profiles | 0630 naive | 0630 fine-tuned | 0630 frozen profiles |
|---|---|---|---|---|---|---|
| organism | 0.748 | 0.761 | 0.790 | 0.773 | 0.761 | 0.750 |
| species | 0.731 | 0.747 | 0.776 | 0.723 | 0.724 | 0.681 |
| genus | 0.656 | 0.673 | 0.686 | 0.650 | 0.595 | 0.572 |
| family | 0.536 | 0.541 | 0.611 | 0.365 | 0.385 | 0.429 |

Difference from naive, both experiments combined, with a 95% CI from resampling **taxa**
(`s19`, 2,000 draws):

| held out | taxa | fine-tuned − naive | frozen profiles − naive |
|---|---|---|---|
| organism | 151 | 0.000 (−0.044 to +0.031) | +0.009 (−0.048 to +0.064) |
| species | 117 | +0.008 (−0.046 to +0.049) | +0.001 (−0.059 to +0.067) |
| genus | 42 | −0.019 (−0.108 to +0.078) | −0.024 (−0.133 to +0.090) |
| family | 25 | +0.012 (−0.150 to +0.308) | +0.069 (−0.105 to +0.320) |

All eight intervals include zero; Holm-adjusted p = 1.0 for all eight.

![pretraining](../outputs/figures/s19_r10_pretraining.png)

**Reading it:**
- **Fine-tuning from BacDive ≈ training from scratch.** The combined difference is within ±0.02 ρ
  at every level.
- **Frozen phenotype profiles are the only arm with a hint.** At family level they are ahead in
  both experiments (+0.075 and +0.064), but the interval is about −0.1 to +0.3. At genus the
  experiments disagree (+0.030 in 0907, −0.078 in 0630).
- **The frozen 256-d embedding is the weakest pretrained arm** (combined −0.04 to −0.11, n.s.).
- **The naive MLP is already good:** ρ 0.75–0.77 with the organism held out, 0.72–0.73 species,
  0.65 genus, 0.37–0.54 family. Organism-level performance is similar to 20261113's XGBoost
  (0.78–0.81; different features and fold draws, so not a head-to-head).
- **Controls:** shuffled labels score ρ ≈ −0.07 to +0.12. A frozen random encoder is close to naive
  at organism and species level but worse at genus (−0.105, CI −0.195 to −0.025), so training the
  encoder does add something once the taxon gap is large.

**Limits of this test:**
1. **Power at family level is low.** There are 25 families (13 in 0630) and Enterobacteriaceae
   dominates. A real gain of +0.1 or more could be hidden. More CV repeats would not fix this;
   more distinct families would.
2. **"Unseen" means unseen by the interaction model only.** BacDive pretraining saw all these
   taxa. That favours pretraining, and it still showed no gain. Holding the same taxa out of
   pretraining has not been tested.
3. **The training recipe was not tuned for transfer.** Early stopping chose only ~20–50 full-batch
   steps (floor 20), so fine-tuning barely moves the encoder. A slower, more regularized recipe
   could behave differently.
4. **The pretrained encoder is the MLP**, which was slightly weaker than logistic regression on
   BacDive (0.672 vs 0.680). The logistic-regression profiles are tested in the follow-up below.

## Follow-up: frozen logistic-regression profiles (s20, 2026-11-23)

**Question:** does swapping the MLP's 351 predicted phenotypes for logistic regression's
(`s14_lab_profiles_logreg.csv`, fitted on BacDive genomes only) produce a gain on unseen taxa?

**Design:** same pairs, folds, seeds, head, early stopping and z-scoring as the s18 frozen-profile
arm. Probabilities are used as is, and no other transform was tried. **Decision rule, fixed before
running:** a gain requires the combined logreg − naive 95% CI to exclude 0 at genus or family level
after Holm correction over the 4 levels.

**Reproduction checks passed:** the re-run MLP-profile arm matches s18's 12,919 saved predictions
(max difference 8e-7, identical step counts), and all 60 s19 contrasts are reproduced.

**Result: no detectable gain; the decision rule is not met.**

Spearman ρ on held-out taxa:

| held out | 0907 naive | 0907 MLP profiles | 0907 logreg profiles | 0630 naive | 0630 MLP profiles | 0630 logreg profiles |
|---|---|---|---|---|---|---|
| organism | 0.748 | 0.790 | 0.776 | 0.773 | 0.750 | 0.738 |
| species | 0.731 | 0.776 | 0.762 | 0.723 | 0.681 | 0.701 |
| genus | 0.656 | 0.686 | 0.678 | 0.650 | 0.572 | 0.606 |
| family | 0.536 | 0.611 | 0.570 | 0.365 | 0.429 | 0.448 |

Logreg profiles − naive, both experiments combined (95% CI from resampling taxa, 2,000 draws):

| held out | taxa | logreg − naive | Holm p |
|---|---|---|---|
| organism | 151 | −0.004 (−0.055 to +0.055) | 1.0 |
| species | 117 | +0.004 (−0.057 to +0.062) | 1.0 |
| genus | 42 | −0.011 (−0.090 to +0.105) | 1.0 |
| family | 25 | +0.058 (−0.022 to +0.258) | 0.42 |

**Reading it:**
- Logreg profiles are indistinguishable from MLP profiles (combined difference −0.013 to +0.013).
- Family level is again ahead in both experiments (+0.033 and +0.083), as with the MLP profiles.
  It remains a hint only. The two profile arms share the same labels and genomes, so they are not
  independent evidence.
- **Conclusion:** stop pursuing BacDive pretraining with this interaction dataset. The family-level
  question needs more distinct families, not another model variant.

**Files:**
- `outputs/s18_r10_{predictions.csv.gz,fold_metrics.csv,summary.csv}`
- `outputs/s19_r10_{bootstrap.csv,rho_by_level.csv,summary.json}`
- `outputs/s20_r10_{predictions.csv.gz,rho_by_level.csv,bootstrap.csv}`
- `outputs/figures/s19_r10_pretraining.png`
