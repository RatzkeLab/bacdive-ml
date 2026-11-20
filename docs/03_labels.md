# 03 Labels

Code: `s03_coverage_report.py`: `merge_unspecified_stereo`, `aggregate_labels`, `KIND_GROUP`.

## Definition

A label is **metabolite × kind of test** (or × kind *group*). Per strain:

| value | meaning |
|---|---|
| 1 | only `+` observations |
| 0 | only `-` observations |
| NaN | untested, only `+/-` (variable/weak), or both `+` and `-` |

**Missing ≠ negative.** The matrix is 7.5% filled, so models need a masked loss or per-label
fitting on labelled strains only.

## Metabolite key: normalized name, NOT ChEBI

- BacDive uses **ChEBI `0` as a placeholder**, which merges "Potassium 5-ketogluconate" (7.6k rows)
  with "Potassium 2-ketogluconate" (6.7k rows).
- Some ChEBI IDs are wrong: CHEBI:16411 holds both "acetic acid" and "indole-3-acetic acid".
- 19,113 rows have no ChEBI ID, or ChEBI 0.
- Names are consistent apart from case and whitespace, so the key is `lower(strip(name))`.

### D/L-unspecified names

"glucose" and "D-glucose" are separate entries from different sources. Rule: merge `x` → `d-x`
if `l-x` is never tested, and `x` → `l-x` if `d-x` is never tested. If both enantiomers are
tested (xylose, arabinose, fucose, alanine, ...), API 50CH tests both forms, so the unprefixed
name is ambiguous and stays a separate label. 27 merges, listed in
`outputs/s03_stereo_name_merges.csv` (e.g. glucose → d-glucose 8.5k rows, trehalose →
d-trehalose 16k, arginine → l-arginine 16k).

## Kinds of test

Rows by kind:

| kind | rows | kind | rows |
|---|---|---|---|
| builds acid from | 457k | reduction | 33k |
| hydrolysis | 85k | degradation | 19k |
| assimilation | 77k | energy source | 17k |
| carbon source | 76k | growth | 12k |
| fermentation | 74k | respiration | 12k |
| (none) | 8.6k | oxidation | 5.6k |

A long tail follows. Groups (`KIND_GROUP`):

| group | kinds | closest to resource competition? |
|---|---|---|
| growth/assimilation | carbon source, energy source, nitrogen source, growth, assimilation | **yes** |
| acid/gas (fermentation) | fermentation, builds acid from, builds gas from | partly (API 50CH acid ≠ growth) |
| respiration/oxidation | respiration, oxidation (Biolog-type redox) | partly |
| hydrolysis/degradation | hydrolysis, degradation | extracellular enzymes |
| reduction/e- acceptor | reduction, electron acceptor | nitrate etc. |
| other | everything else, incl. missing kind | |

Labels are produced at both levels. The grouped matrix is what s03 writes; per-kind stats are in
`outputs/s03_label_stats_mu_by_kind.csv`. Grouping trades purity for coverage: disagreement is
15% within an exact kind vs 19% within a group ([05](05_label_noise_and_taxonomy_baseline.md)).

## API kits are already in the utilization list

DSMZ files each kit's results in `metabolite utilization` / `enzymes` under a kit-specific
`@ref`:

| kit | MU ref | kit | MU ref |
|---|---|---|---|
| 50CHac | 68371 | rID32STR | 68381 |
| 20NE | 68369 | coryne | 68379 |
| 20E | 68368 (enzymes) | ID32E | 68374 |
| 20A | 68367 | ID32STA | 68375 |
| NH | 68377 | CAM | 68373 |

API zym is in enzymes under 68382. Coverage of kit strains is ~100% for each
(`outputs/s03_api_kits_folded_into_mu.csv`), so the kit code tables need no separate parsing.
These DSMZ kit refs supply **43% of all utilization rows on genome strains**. That's the most
standardized slice of the data.

## Usable-label threshold

A label is "usable" if it has ≥25 positive **and** ≥25 negative strains *with a genome*. That
gives 351 grouped labels: 173 growth/assimilation, 67 respiration/oxidation, 65 acid/gas,
30 hydrolysis, 11 other, 5 reduction. At ≥100 per class: 149 labels.
