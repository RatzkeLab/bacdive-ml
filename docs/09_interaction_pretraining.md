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

## Results

*(pending: full run with 10 repeats)*
