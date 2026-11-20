# 02 Data tables

`scripts/s02_parse_tables.py` flattens the JSON without interpreting it. Each list-of-dicts
section becomes one row per item, with every key kept as a column, plus `bacdive_id`.
All interpretation happens later (s03+). Files are `data/tables/<name>.csv.gz`.

| table | rows | key columns | notes |
|---|---|---|---|
| `strains` | 102,187 | `bacdive_id`, `dsm_number`, `type_strain`, `ncbi_taxid`, `strain_designation`, `culture_collection_nos`, `domain..species` (+`*_bacdive`) | taxonomy columns = LPSN where available; `*_bacdive` = BacDive's own (older) names |
| `metabolite_utilization` | 883,830 | `@ref`, `Chebi-ID`, `metabolite`, `utilization activity` (+, -, +/-), `kind of utilization tested` | **the main label source**, see [03](03_labels.md) |
| `enzymes` | 579,801 | `value` (enzyme), `activity`, `ec` | API zym folded in here |
| `api_kits` | 1,560,882 | `kit`, `test` (kit code), `result` | redundant: already in MU/enzymes ([03](03_labels.md)) |
| `metabolite_production` | 39,238 | `metabolite`, `production` | |
| `metabolite_tests` | 30,563 | indole / VP / methyl red / citrate tests | |
| `genomes` | 28,880 | `INSDC accession` (no version), `assembly level`, `score`, `BV-BRC accession`, `IMG accession` | 1,332 rows have no valid GCA/GCF accession |
| `seq16s` | 44,899 | `accession`, `length` | |
| `traits` | 87,258 | `trait`, `value` | Gram stain, cell shape, motility, oxygen tolerance, spore formation, nutrition type, pigment |
| `culture_temp` / `culture_ph` / `halophily` / `tolerance` | | growth / optimum ranges | strings like `"25-30"`; still need parsing |
| `antibiotic_resistance` | 26,151 | | |
| `isolation_cat` | 93,118 | `Cat1..Cat3` | e.g. `#Host / #Mammals` |

Not parsed yet: `culture medium` (DSMZ medium names, MediaDive links, full compositions, growth
yes/no), `Literature`, `Isolation` free text, `multicellular morphology`, `fatty acid profile`.

## Derived tables (from s03 / s04)

| table | what |
|---|---|
| `genome_choice.csv.gz` | one genome per strain: best assembly level (complete > chromosome > scaffold > contig), then highest BacDive `score`. No accession is shared between strains. |
| `genome_choice_ncbi.csv.gz` | the same, plus NCBI QC: `ncbi_accession` (versioned), `paired_accession` (GCF), `paired_status`, `checkm_completeness`/`checkm_contamination`, `ani_taxonomy_check`, `genome_size`, ... |
| `label_matrix_mu_grouped_minority25_genome_strains.csv.gz` | 11,115 strains × 351 labels, values 1/0/NaN; columns named `<metabolite>\|<kind group>` |

## Gotchas

- **`@ref`** is the source document. DSMZ's in-house API-kit series use refs 68367–68382.
- **Taxonomy renames:** LPSN moved several of our families (Xanthomonadaceae → Lysobacteraceae,
  Leuconostocaceae → Lactobacillaceae, Planococcaceae → Caryophanaceae, Brucellaceae partly
  → Bartonellaceae) and genera (Ochrobactrum → Brucella). **Match on `x` OR `x_bacdive`.**
