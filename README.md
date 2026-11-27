# bacdive-ml

Pretrain a genome → metabolic-capability model on BacDive (11k strains with genome + phenotype
labels). Then reuse it as strain features for the pairwise-interaction model.

**Status (2026-11-27): shelved.** ~51 GB of regenerable data (genomes, proteins, MMseqs2 DBs,
family matrix, MLP weights, logs) was deleted; see [data/DELETED_DATA.md](data/DELETED_DATA.md)
for what is gone and how to rebuild it. Results, docs and small tables are kept.

**Documentation:** [docs/00_index.md](docs/00_index.md) is the mini-wiki home (data, labels,
results, decisions, TODO).

## Run order

Use env `karl_datasci` unless noted: `/home/rl/mambaforge/envs/karl_datasci/bin/python scripts/<s>.py`.

| script | does | runtime |
|---|---|---|
| `s01_download_bacdive.py` | sweep BacDive IDs via `/v2/fetch`, cache gz JSON (resumable) | 3 min |
| `s02_parse_tables.py` | flatten JSON → `data/tables/*.csv.gz` | 20 s |
| `s03_coverage_report.py` | labels, coverage, genome choice, label matrix, collection overlap, figures | 15 s |
| `s04_ncbi_genome_check.py` | NCBI QC metadata for every chosen genome | 2 min |
| `s05_taxonomy_baseline.py` | taxonomy-only AUROC per label + positive controls | 10 s |
| `s06_pilot_genome_download.py` | 12-genome pilot download (`datasets` CLI from env bw25113_panx) | 5 s |
| `s07_label_reliability.py` | label noise from independent references | 30 s |
| `s08_download_genomes.py` | select QC'd genomes (8,740) + NCBI dehydrated download | ~45 min |
| `s09`–`s14` via `run_s09_to_s14.sh` | Prodigal → MMseqs2 families → matrix → baselines → MLP → lab transfer | ~1.5 h |
| `s09_call_genes.py` | Prodigal on BacDive + our 298 genomes (60 processes) | 30 min |
| `s10_mmseqs_cluster.py` | MMseqs2 v16 (env eggnog_mapper) 50% id / 80% cov families | 35 min |
| `s11_family_matrix.py` | genome × family sparse matrix, transfer coverage | 2 min |
| `s12_baseline_models.py` | taxonomy / kNN / logreg per label; `--split genus\|strain\|family`, `--min-prev` | 3 min |
| `s13_multitask_mlp.py` | multi-task MLP, same folds; `--final` saves the transfer model | 40 min |
| `s14_transfer_lab_genomes.py` | domain shift + predicted profiles/embeddings for our genomes, keyed by `strain_no` | 2 min |
| `s15_interpret_families.py` | top-weighted families named via Swiss-Prot | 1 min |
| `s16_summary_figure.py` | AUROC by CV regime figure | s |
| `s17_sheet_taxonomy_check.py` | lab genomes' nearest BacDive genus vs the master sheet's 16S taxonomy | s |
| `s18_interaction_pretraining.py` | naive vs BacDive-pretrained MLP for genome → interaction, organism/species/genus/family held out | 80 min |
| `s19_pretraining_bootstrap.py` | taxon-level bootstrap CIs and figure for s18 | 10 s |
| `s20_logreg_profiles.py` | follow-up: frozen logreg BacDive profiles on s18's folds, with reproduction checks | 30 s |

## Headline numbers

| | |
|---|---|
| BacDive records | 102,187 |
| strains with genome + metabolic labels | 11,143 (~9k pass genome QC) |
| usable labels (≥25 pos and ≥25 neg) | 351 |
| label noise between independent sources | ~15% |
| taxonomy-only AUROC (median), strain out / genus out | 0.68 / 0.62 |
| genome (MMseqs2 families, logreg) AUROC, genus out: all / ≥100 per class | **0.68 / 0.74** (taxonomy 0.61 / 0.64) |
| same, strain out (the regime matching our strains) | **0.74 / 0.82** |
| BacDive pretraining → interaction model, Δρ vs naive (genus / family held out) | −0.02 / +0.01, CIs include 0: **no detectable gain** |
| frozen logreg BacDive profiles, Δρ vs naive (genus / family held out) | −0.01 / +0.06, CIs include 0: **no detectable gain** |
