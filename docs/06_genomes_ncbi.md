# 06 Genomes (NCBI)

## Accession check: `s04_ncbi_genome_check.py`

- **API:** NCBI Datasets v2 REST, `POST https://api.ncbi.nlm.nih.gov/datasets/v2/genome/dataset_report`
  with `{"accessions": [...], "page_size": 1000}`, sent in batches of 200 accessions.
- **Accessions:** BacDive gives versionless accessions (`GCA_003269275`). NCBI resolves them to
  the current version (378 are not `.1`).
- **Rate:** no API key is needed at our rate (stay under 5 requests/s). The run took about 2 min.
  Raw reports are cached in `data/raw/ncbi_reports/`.

Results for all 19,639 chosen genomes (`outputs/s04_ncbi_summary.json`):

| | n |
|---|---|
| resolve at NCBI | 19,598 (all `current`); 41 don't |
| current RefSeq (GCF) partner | 18,664 |
| CheckM available | 19,207 |
| CheckM ≥90% complete, ≤5% contamination | 14,138 |
| ANI taxonomy check OK / inconclusive / failed | 18,379 / 1,064 / 153 |
| total size | 88.6 Gb; median genome 4.16 Mb |

Among the 11,105 phenotyped strains that resolve:

| QC rule (CheckM NA kept, ANI-failed dropped) | strains | labels |
|---|---|---|
| CheckM ≥90 / ≤5 | 8,236 | 216k |
| CheckM ≥90 / ≤10 | 8,957 | 235k |
| CheckM ≥80 / ≤10 | 9,467 | 249k |
| no QC | 11,105 | 291k |

≤5% contamination is too strict for CheckM v1 on finished genomes; the complete Phaeobacter
genome scores 5.3%.

**Final rule: CheckM ≥90 / ≤10 AND CheckM must exist → 8,740 genomes**
(`data/tables/genome_set_hq.csv`). The 217 genomes without CheckM were first kept and
downloaded, but they include partial assemblies:
- bd16333 *Streptomyces calvus*, 5 kb (Prodigal refused it)
- a 68 kb *Leucobacter*
- a 721 kb *E. coli*
- a 902 kb *Kineosporia*

Another 10 are more than 1.6× their genus median size (possibly mixed), so all 217 are dropped.
They remain in the MMseqs2 clustering (harmless: extra members) but are excluded from features
and models.

## Pilot download: `s06_pilot_genome_download.py`

- **Tool:** NCBI `datasets` CLI v18.33.1. It is **not in karl_datasci**, so we call
  `/home/rl/mambaforge/envs/bw25113_panx/bin/datasets`.
- **Test:** 12 genomes (`--include genome,protein,gff3`) took **5 s**.
- **Size per genome:** about 1.23 MB fna, 0.36 MB gff and 0.81 MB faa (zipped), 2.4 MB total.
- **Extrapolated to all ~11.1k phenotyped genomes:** about 26 GB zipped / 94 GB unzipped;
  proteins alone about 9 GB zipped / 17 GB unzipped.
- **Bulk downloads:** NCBI recommends `--dehydrated` + `datasets rehydrate` above ~1,000 genomes.
- **Disk:** `/` has ~480 GB free and is 87% full. Keep fna gzipped and delete what isn't needed.

## Our own genomes

- **Location:** 298 Prokka GenBank files at
  `/home/rl/Documents/Workstation_Data/Or/KEIO_collection_panX_mapping/STRAIN_LIBRARY_SEQ_PROCCSSED_DATA/all_gbk_files_including_CeMbio43/input_GenBank/`.
  This is the same axis as the KO tables, and the files contain the contig sequence.
- **Same organism, several genomes:** 298 genomes = 181 organisms (clone clusters), and about 1 in
  5 wells holds a different organism than its genome. See the project memories on genome-well
  concordance and the 20261113 genomic-ML re-check before using these in interaction ML.
