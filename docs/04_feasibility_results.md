# 04 Feasibility results

From `s03_coverage_report.py` → `outputs/s03_coverage_summary.json`, `outputs/figures/s03_*.png`.

## Coverage

| | strains |
|---|---|
| BacDive records | 102,187 (22,126 type strains) |
| ≥1 resolved metabolite-utilization label | 30,312 |
| genome accession (valid GCA/GCF) | 19,639 |
| **genome + labels** | **11,143** (9,647 type strains) |
| genome + ≥10 / ≥20 / ≥50 usable labels | 8,360 / 5,340 / 2,053 |
| genome + labels + good genome (CheckM ≥90/≤10, ANI not failed) | 8,957 (235k labels) |

- Genomes per strain: 14,648 strains have one; the rest have 2–10. The chosen assemblies are
  4,502 complete, 373 chromosome, 7,105 scaffold and 7,659 contig.
- **Taxonomic breadth** (genome + labels): 9,793 species, 2,236 genera, 466 families, 36 phyla.
  - Phyla: Pseudomonadota 4,446; Actinomycetota 2,793; Bacillota 2,187; Bacteroidota 1,208.
  - **Median one strain per species**, 99th percentile three.
- Top genera: Streptomyces 463, Pseudomonas 245, Streptococcus 183, Mycobacterium 161, Paenibacillus 149.

![labels per strain](../outputs/figures/s03_labels_per_genome_strain.png)
![usable labels by kind](../outputs/figures/s03_usable_labels_by_kind.png)

## Overlap with our lab collection

`outputs/s03_collection_genus_coverage.csv`. Genus comes from the collection sheet's
`phylogeny` column (Greengenes names), matched to LPSN **or** BacDive genus.

> **Caveat (2026-11-20):** the sheet's taxonomy is unreliable for about a third of genome strains
> ([08](08_models.md), s17). The genome-based nearest-BacDive genera in s14 are the better guide;
> numbers below are a rough survey only.

- 213 of 257 genus-assigned collection strains sit in genera with ≥10 BacDive genome+label strains.
- Well covered: Pseudomonas 287, Paenibacillus 149, Corynebacterium 136, Staphylococcus 123,
  Acinetobacter 101, Bacillus 98, Microbacterium 80, Chryseobacterium 80.
- Thin: Rahnella 8, Ochrobactrum 10, Lysinibacillus 12, Gluconobacter 7, Erwinia 7, Ewingella 1,
  Wautersiella 0, Cryocola 0.
- Our strains in BacDive (`outputs/s03_lab_strain_name_hits.csv`): only MYb115 (DSM 119730,
  Pseudomonas) and MYb71 (DSM 119732, Brucella/Ochrobactrum). Neither has a genome or
  phenotype data there. A "JUB59" hit is an unrelated Bizionia strain.

![top genera](../outputs/figures/s03_top_genera.png)

## Other phenotype blocks (genome strains)

These are potential extra multi-task heads:

| block | genome strains |
|---|---|
| culture temperature | 17,061 |
| isolation category | 13,516 |
| oxygen tolerance | 10,447 |
| enzymes (62 usable labels) | 10,598 |
| Gram stain | 9,278 |
| motility | 8,703 |
| halophily | 6,915 |
| culture pH | 5,345 |
| spore formation | 4,178 |
| antibiotic resistance | 2,494 |
