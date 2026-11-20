# 07 Gene-family features (MMseqs2)

Lightweight, annotation-free gene content. Every genome becomes a presence/absence vector over
protein families defined by clustering **all** proteins together. "All" means the BacDive QC set
plus our 298 lab genomes, so both sides share one feature space. The eggNOG/KO route is
deferred ([99](99_decisions_and_todo.md)).

## Pipeline

| step | script | what | runtime |
|---|---|---|---|
| genome set | `s08_download_genomes.py` | phenotyped strains passing QC (current at NCBI, CheckM ≥90/≤10 and present, ANI not failed): **8,740** (8,957 downloaded before the CheckM-present rule); `datasets download --dehydrated` + `rehydrate --gzip` | ~45 min (NCBI-bound) |
| gene calling | `s09_call_genes.py` | Prodigal 2.6.3 `-p single` on every genome. Lab genomes: contigs from the Prokka `.gbk`. 60 processes. | ~1 min per 300 genomes |
| clustering | `s10_mmseqs_cluster.py` | MMseqs2 v16 cascaded `cluster`, `--min-seq-id 0.5 -c 0.8 --cov-mode 0` | see results |
| matrix | `s11_family_matrix.py` | genome × family copy-number matrix (sparse) + family table + transfer coverage | minutes |

## Design choices

- **Same gene caller everywhere.** NCBI ships PGAP proteins for 95% of the BacDive genomes, but
  our genomes would then be called differently (Prodigal, via Prokka). Prodigal is re-run on
  BacDive *and* lab genomes, so a family's presence never depends on which caller saw it.
  Prokka itself uses Prodigal, and the lab gene counts match what Prokka reported.
- **IDs:**
  - Genomes are `bd<bacdive_id>` or `lab_<gbk stem>`.
  - Proteins are `<genome_id>|<n>`, so cluster members map back to genomes with a string split.
- **50% identity / 80% bidirectional coverage** is roughly UniRef50-sized families. It is
  specific enough to keep paralog groups apart and loose enough to cross genera. The CLI flags
  allow 30% (more functional) or 90% (strain-level) runs side by side (`mmseqs_id30_c80/`, ...).
- **Feature filter (in s12):** keep families present in ≥10 BacDive genomes and ≤99% of them. The
  filter is computed on BacDive genomes only and is label-free, so lab genomes are projected into
  exactly the training feature space.
- **Kept for later:** `mmseqs_db/proteinsDB` (re-clustering without re-reading FASTA) and
  `reps.faa.gz` (representative sequences). For eggNOG, annotating representatives only and
  propagating to members is the cheap route.

## Files

```
data/interim/proteins/<genome_id>.faa.gz      Prodigal proteins (kept)
data/interim/mmseqs_db/proteinsDB*            mmseqs sequence DB
data/interim/mmseqs_id50_c80/clusters.tsv.gz  rep \t member
data/interim/mmseqs_id50_c80/reps.faa.gz
data/tables/families_id50_c80/counts.npz      genomes x families (CSR int32 copy number)
data/tables/families_id50_c80/genomes.txt     row order
data/tables/families_id50_c80/families.csv.gz fam, rep, n_members, n_genomes, n_bd, n_lab
data/tables/genomes_all.csv                   genome_id, source, n_proteins, n_aa, path
```

## Results (2026-11-20, id50_c80)

**Gene calling (s09):**
- 9,254 genomes: 8,956 BacDive + 298 lab. bd16333 (5 kb) was refused by Prodigal.
- Proteins per genome: BacDive median 4,010 (IQR 2,990–5,185), lab median 4,678. The two sides
  are comparable.
- 39.87M proteins, 13.3 GB FASTA.
- The NCBI rehydrate dir also holds `*_cds_from_genomic` / `*_rna_from_genomic` files for one
  accession. s09 now excludes them explicitly; an early glob would have taken the CDS file.

**Clustering (s10):** MMseqs2 cascaded clustering took **33 min** on 60 threads (plus 1 min createdb).

**Families (s11, `data/tables/families_id50_c80/summary.json`):**

| | n |
|---|---|
| families | 6,371,151 |
| singletons | 4,227,367 |
| in ≥2 genomes | 2,057,389 |
| in ≥10 BacDive genomes (logreg/kNN features) | 420,300 (413,525 after the 8,740-genome QC filter and ≤99% cap) |
| in ≥1% of BacDive genomes (MLP features, ≥87) | ~56k (56,462) |
| lab-only families | 46,658 |

**Transfer coverage.** This is the share of a lab genome's proteins that fall in families also seen
in BacDive genomes. Quantiles are over the 298 lab genomes:

| seen in ≥k BacDive genomes | 5% | 25% | median | 75% | 95% |
|---|---|---|---|---|---|
| k = 1 | 0.88 | 0.93 | **0.955** | 0.98 | 0.99 |
| k = 5 | 0.75 | 0.83 | 0.89 | 0.94 | 0.96 |
| k = 20 | 0.55 | 0.67 | 0.78 | 0.89 | 0.93 |

Our strains' gene content is well represented in the BacDive feature space.

**Disk:**

| path | size |
|---|---|
| `data/interim/proteins` | 7.2 GB |
| `mmseqs_db` | 16 GB |
| `mmseqs_id50_c80` | 3.5 GB |
| `all_proteins.faa` | 13 GB, deletable: rebuild from `proteins/` with s09 |
| `data/raw/ncbi_genomes` | 12 GB |
