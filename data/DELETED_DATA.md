# Deleted data (project shelved 2026-11-27)

bacdive-ml went on the backburner on 2026-11-27. ~51 GB of regenerable data was deleted to free
workstation space. Results, docs, scripts and small tables were kept.
Everything below can be rebuilt with the scripts in `scripts/` (envs as in the top-level README).

## Kept (small, or not exactly reproducible)

- `data/raw/ncbi_reports/` (8 MB): NCBI genome QC metadata from `s04`.
- `data/raw/ncbi_genomes/{accessions.txt, dehydrated.zip, md5sum.txt, missing_accessions.txt}`:
  the exact genome set (8,957 accessions) and NCBI's dehydrated package for it.
- `data/raw/ncbi_pilot/accessions.txt`, `data/raw/swissprot/reldate.txt`.
- `data/tables/*` except the family matrix (labels, strains, genome choice, ...).
- `data/tables/families_id50_c80/{genomes.txt, summary.json}`.
- `outputs/` except logs, `s15_id50_c80_tmp/` and `s13_id50_c80_genus_oof.csv.gz`. In particular
  the lab-strain features `outputs/s14_lab_*.csv` (referenced by pairs_to_communities_ml) are kept.

## Deleted, and how to get it back

Run from the repo root with `PY=/home/rl/mambaforge/envs/karl_datasci/bin/python`.

| deleted | size | regenerate with | time |
|---|---|---|---|
| `data/raw/bacdive_v2/` (BacDive API snapshot, 2,601 gz batches, Nov 2026) | 117 MB | `$PY scripts/s01_download_bacdive.py`, then `s02_parse_tables.py` to refresh `data/tables/`. Deleted at the user's request: a fresh download is a *newer* BacDive, so labels/counts will drift a little from the documented results (the kept `data/tables/*.csv.gz` are the Nov 2026 parse) | 3 min |
| `data/raw/ncbi_genomes/ncbi_dataset/` (8,957 BacDive genome FASTAs) | 12 GB | `datasets rehydrate --directory data/raw/ncbi_genomes --gzip` (CLI in env `bw25113_panx`; `dehydrated.zip` is still there). Or rerun `$PY scripts/s08_download_genomes.py` | ~45 min |
| `data/raw/ncbi_pilot/{pilot.zip, unzipped/}` | 125 MB | `$PY scripts/s06_pilot_genome_download.py` (only a 12-genome pilot, not needed downstream) | s |
| `data/raw/swissprot/uniprot_sprot.fasta.gz` | 90 MB | `wget -P data/raw/swissprot https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz`. The run used release 2026_04; a newer release only changes s15's family names | 1 min |
| `data/interim/proteins/` (Prodigal proteins, 9,254 genomes incl. 298 lab) | 7.2 GB | `$PY scripts/s09_call_genes.py`. Lab genomes are read from the panX/Prokka GBKs under `/home/rl/Documents/Workstation_Data/Or/KEIO_collection_panX_mapping/...` (path in s09) | 30 min |
| `data/interim/all_proteins.faa`, `data/interim/mmseqs_db/`, `data/interim/mmseqs_id50_c80/` | 32 GB | `$PY scripts/s10_mmseqs_cluster.py` (MMseqs2 v16 from env `eggnog_mapper`) | 35 min |
| `data/tables/families_id50_c80/{counts.npz, families.csv.gz}` (genome x family matrix) | 126 MB | `$PY scripts/s11_family_matrix.py --tag id50_c80` | 2 min |
| `data/models/mlp_id50_c80_all.pt` (BacDive multi-task MLP) | 112 MB | `$PY scripts/s13_multitask_mlp.py --tag id50_c80 --final` | 40 min |
| `outputs/s13_id50_c80_genus_oof.csv.gz` (MLP out-of-fold preds) | 38 MB | same s13 run | |
| `outputs/*.log`, `outputs/s15_id50_c80_tmp/`, stray `tmp.prodigal.stdin.*` | ~2 MB | logs come back with any rerun; headline numbers are in the `*_summary.json` files and `docs/` | |

The whole chain after the genome download is `bash scripts/run_s09_to_s14.sh` (~1.5 h), then
`s15` and `s18`-`s20` individually if needed. A rebuild should reproduce `summary.json` in
`data/tables/families_id50_c80/` (39,867,422 proteins, 6,371,151 families); if NCBI has
suppressed or replaced genomes in the meantime, the counts will differ slightly.
