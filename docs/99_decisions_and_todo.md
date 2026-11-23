# 99 Decisions & TODO

## Decision log

| date | decision | why |
|---|---|---|
| 2026-11-20 | BacDive API v2 with `predictions=0`, full ID sweep | v1 frozen; predicted phenotypes would make a genome model circular; no list endpoint |
| 2026-11-20 | label key = normalized metabolite name (+ D/L merge rule), not ChEBI | ChEBI 0 placeholder and mis-assigned IDs ([03](03_labels.md)) |
| 2026-11-20 | keep assay kinds separate / grouped, never pool all kinds | acid production ≠ growth; disagreement rises 15% → 19% when grouped |
| 2026-11-20 | usable label = ≥25 pos and ≥25 neg genome strains | 351 labels; threshold is a knob, not a finding |
| 2026-11-20 | genome QC: NCBI current, CheckM ≥90/≤10, ANI check not failed | ≤5% contamination too strict for CheckM v1 ([06](06_genomes_ncbi.md)) |
| 2026-11-20 | CheckM must exist (drops 217 → 8,740 genomes) | CheckM-less genomes included 5 kb–900 kb partial assemblies and oversized ones |
| 2026-11-20 | **gene-family features via MMseqs2 clustering first; eggNOG later** | lighter weight: no 50 GB DB, hours not days. User's call. |
| 2026-11-20 | re-call genes with Prodigal 2.6.3 on **all** genomes (BacDive + ours) | same gene caller on both sides of the transfer; NCBI PGAP proteins would differ systematically from ours |
| 2026-11-20 | evaluate with genus-grouped CV against a fold-matched taxonomy baseline | random strain splits overstate generalization (cf. interaction-ML leak) |
| 2026-11-20 | **metric = within-fold AUROC**, never pooled OOF | pooled AUROC is biased under taxon-grouped folds: 0.43 on random features ([08](08_models.md)) |
| 2026-11-20 | MMseqs2 50% id / 80% bidirectional coverage; features = families in ≥10 BacDive genomes (logreg) or ≥87 (MLP) | rare families add nothing (0.680 vs 0.683) |
| 2026-11-20 | logreg (C=0.003) is the current best model; use `s14_lab_profiles_logreg.csv` for transfer | beats MLP on 62% of labels |
| 2026-11-20 | report strain-out as the realistic regime for OUR strains | lab genomes sit closer to BacDive (median J 0.85) than even strain-out test strains (0.71) |
| 2026-11-20 | s14 outputs keyed by **`strain_no`** (20261113 rule); id map in `s14_lab_genome_ids.csv` | user's convention in the interaction project; 1:1 over the 298 genome wells. Fixed a bug on the way: clone clusters were looked up by gbk stem, losing 38 Schulenburg genomes (124 → 162 attached) |
| 2026-11-20 | taxonomy for lab strains comes from the genome, not the master sheet | s17: for 44 of 159 checkable strains the sheet's own Sanger 16S is a different organism from the genome (median 19% divergent), mostly invader clones ([08](08_models.md)) |
| 2026-11-23 | do not adopt BacDive pretraining as the default for interaction models yet | fine-tuned − naive is within ±0.02 ρ at every held-out level; all CIs include 0 ([09](09_interaction_pretraining.md)) |
| 2026-11-23 | judge pretraining with a taxon-level bootstrap, not SD across CV repeats | repeats reuse the same strains; family level has only 13–25 taxa |
| 2026-11-23 | stop pursuing BacDive pretraining on the current interaction data | frozen logreg profiles also miss the pre-set gain rule (family +0.058, CI −0.022 to +0.258, Holm p 0.42; [09](09_interaction_pretraining.md)) |

## TODO

### Deferred: eggNOG-mapper / KEGG KO route (heavier duty)
- [ ] Download the eggNOG 5 database (~50 GB with the diamond DB; `download_eggnog_data.py` in
  env `eggnog_mapper`, emapper v2.1.15). Disk is 87% full, so check space first.
- [ ] Annotate the same Prodigal proteins used for MMseqs2 (BacDive HQ set + our 298). Rough
  estimate: ~11k genomes × 1–2 min at 64 threads is days. Ways to speed it up:
  - annotate only MMseqs2 cluster representatives and propagate to members (much faster), or
  - use `--dmnd_iterate no`, or MMseqs2 mode (`-m mmseqs`).
- [ ] Build KO / KEGG-module / CAZy / PFAM presence tables. Rerun the s-models on them and compare
  with the MMseqs2 families under the **same folds**.
- [ ] Our 298 genomes must be re-annotated by the same run. The existing
  `final_genomic_tables/KEGG_ko_and_strains_table.csv` came from a different pipeline and is
  not guaranteed feature-compatible.

### Next
- [x] (done 2026-11-23, s20: no detectable gain) frozen **logistic-regression** profiles
  (`s14_lab_profiles_logreg.csv`) as the strain vector
- [ ] Remaining pretraining follow-ups, only if still of interest:
  - a slower, regularized fine-tuning recipe
  - if a gain appears: repeat with the held-out taxa removed from BacDive pretraining too
  - the family-level question needs more distinct families in the interaction data
- [x] (done 2026-11-23: no detectable gain, see [09](09_interaction_pretraining.md)) **Transfer to the interaction model:** `outputs/s14_lab_profiles_logreg.csv` (351-d predicted
  phenotypes) and/or `s14_lab_embeddings_mlp.csv` as strain features in the 20261113 leak-free
  setup (organism-grouped folds, wrong-genome wells dropped). Compare against KO and 16S on the
  same folds. Files are keyed by `strain_no`; group folds by `clone_cluster`
  (`s14_lab_genome_ids.csv`). Never take taxonomy from the master sheet (s17).
- [ ] Improve the MLP: stronger regularization, smaller hidden layer, ensembling, or a linear
  multi-task model (shared low-rank weights) that keeps logreg's strength.
- [ ] More labels per strain: the per-label AUROC climbs steeply with label size.
- [ ] Name families for other labels / build an interpretable "pathway completeness" feature set
  from the s15 hits.
- [ ] Cleanup option: `data/interim/all_proteins.faa` (13 GB) can be deleted. It is rebuilt from
  `proteins/` by s09 in about 1 min, and only s10 needs it.
- [ ] Mechanistic pair features from predicted carbon-source profiles: niche overlap
  (competition) and complementarity.
- [ ] Source/assay-aware labels: keep DSMZ kit refs as their own label family, or add a source covariate.
- [ ] Extra heads: enzymes (62 usable), temperature/pH/salt ranges (need parsing), oxygen, Gram.
- [ ] Parse `culture medium` (MediaDive compositions + growth yes/no) for growth-on-defined-media labels.
- [ ] Autoencoder: semi-supervised encoder (reconstruct gene families for all genomes plus phenotype
  heads where labelled). Unlabelled corpus: all 19.6k BacDive genomes, later GTDB species reps.

### Open questions
- Within-species variation is nearly absent (about one strain per species). Does that limit transfer to
  our collection, which has many strains per species?
- Identity/coverage thresholds for families (start: 50% id / 80% coverage). Try 30% (more
  functional) and 90% (more strain-specific)?
