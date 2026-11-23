# bacdive-ml wiki

Pretraining a genome → metabolic-capability model on BacDive. The model's predictions or
embedding will become strain features for the pairwise-interaction model, so it reaches
strains outside our ~300 genomes.

| page | what's in it |
|---|---|
| [01 BacDive API & download](01_bacdive_api_and_download.md) | how the API behaves, how the bulk download works, rate/size facts |
| [02 Data tables](02_data_tables.md) | every table in `data/tables/`, its columns and gotchas |
| [03 Labels](03_labels.md) | how a phenotype label is defined; ChEBI and naming traps; API kits |
| [04 Feasibility results](04_feasibility_results.md) | coverage numbers, overlap with our collection, figures |
| [05 Label noise & taxonomy baseline](05_label_noise_and_taxonomy_baseline.md) | how noisy the labels are; how well taxonomy alone predicts them |
| [06 Genomes (NCBI)](06_genomes_ncbi.md) | accession QC, RefSeq partners, pilot download, sizes |
| [07 Gene-family features (MMseqs2)](07_gene_family_features.md) | gene calling + protein-family clustering pipeline, transfer coverage |
| [08 Models](08_models.md) | evaluation protocol, baselines, MLP, CV regimes, gene-level sanity check, transfer to lab genomes |
| [09 Interaction pretraining test](09_interaction_pretraining.md) | naive vs BacDive-pretrained MLP for genome → interaction, organism to family held out: **no detectable gain** |
| [99 Decisions & TODO](99_decisions_and_todo.md) | decision log, open questions, deferred work (eggNOG) |

## Status (2026-11-20)

- **Done:**
  - BacDive download, parse and feasibility survey
  - genome QC + download (8,740 genomes)
  - Prodigal + MMseqs2 protein families (6.4M families; 413k used)
  - genome → phenotype models
  - predicted profiles and embeddings for our 298 genomes
- **Headline:**
  - Logistic regression on protein-family presence beats taxonomy on 87% of labels (genus held
    out). Median AUROC: taxonomy 0.61 → genome kNN 0.66 → logreg 0.68 (0.74 for labels with
    ≥100 per class).
  - With only the strain held out, the regime that matches our strains: 0.74 (0.82).
  - The top-weighted genes are the textbook pathways (ara, xyl, rha, nar/nap, nir/nor/nos).
  - The multi-task MLP is slightly worse than logreg so far.
- **Pretraining test (2026-11-23, [09](09_interaction_pretraining.md)):** BacDive pretraining did
  not detectably improve interaction prediction at any held-out level, with MLP or logreg
  profiles; family level is underpowered.
- **Not started:**
  - genome autoencoder
  - the eggNOG/KO route

## One-glance figures

![by split](../outputs/figures/s16_id50_c80_auroc_by_split.png)

## Folder layout

```
bacdive-ml/
  README.md             entry point + run order
  docs/                 this wiki
  scripts/              s01..sNN, run in order, env karl_datasci unless noted
  data/raw/             caches of downloaded data (safe to delete, scripts re-fetch)
  data/interim/         large intermediates (proteins, mmseqs DBs)
  data/tables/          parsed / derived tables (csv.gz, npz)
  outputs/              summaries (json/csv), logs, figures/
```
