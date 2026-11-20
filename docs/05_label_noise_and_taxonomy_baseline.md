# 05 Label noise & taxonomy baseline

## Label noise: `s07_label_reliability.py`

Method: pick a strain × metabolite × assay label measured by **≥2 different `@ref`s**. Collapse
each ref to one call and count how often the refs disagree. The headline s03 "conflict rate"
(0.7%) is misleading, because most labels have only one observation.

| comparison | replicated labels | disagreement |
|---|---|---|
| exact kind, all | 31,493 | **15.3%** |
| two DSMZ API kits | 23,732 | 14.4% |
| DSMZ kit vs literature | 6,774 | 17.9% |
| two literature refs | 987 | 19.4% |
| grouped kind | 53,032 | 19.2% |

- **Source bias:** when the literature and a DSMZ kit disagree, the literature says `+` 957
  times vs 338 for the kit. Species descriptions report positives more readily, or the kits are
  less sensitive.
- Per-label disagreement is in `outputs/s07_disagreement_per_label.csv`.
  - Some labels are near-perfect: D-fucose, L-sorbose, galactitol acid production.
  - Some are poor: starch acid 41%, citrate as carbon source 47%, D-mannose fermentation 65%.
- Implication: accuracy against these labels is capped at roughly 85%. Consider keeping
  source/assay in the label, or adding a source covariate, or using per-source heads.

## Taxonomy-only baseline: `s05_taxonomy_baseline.py`

A nearest-taxon frequency predictor, scored by AUROC per label (351 labels):

| regime | median | IQR |
|---|---|---|
| strain held out (genus-mates visible; fallback family → order → global) | **0.68** | 0.62–0.76 |
| whole genus held out (predict from rest of family → order → class) | **0.62** | 0.56–0.68 |

Median by group, strain-out / genus-out:

| group | strain-out | genus-out |
|---|---|---|
| acid/gas | 0.79 | 0.69 |
| reduction | 0.79 | 0.78 |
| hydrolysis | 0.69 | 0.62 |
| growth/assimilation | 0.67 | 0.61 |
| respiration | 0.64 | 0.60 |

**Positive controls** (same code, conserved traits), strain-out / genus-out:

| trait | strain-out | genus-out |
|---|---|---|
| Gram stain | 0.99 | 0.99 |
| oxygen tolerance | 0.98 | 0.96 |
| sporulation | 0.97 | 0.96 |
| motility | 0.88 | 0.86 |

The predictor works, so **metabolic labels really are only weakly phylogenetic**, partly
because of the ~15% noise. This is the bar a genome model has to clear, and there is room to
clear it. That contrasts with our interaction data, where KO ≈ 16S once leaks were fixed.

![taxonomy baseline](../outputs/figures/s05_taxonomy_baseline_auroc.png)
