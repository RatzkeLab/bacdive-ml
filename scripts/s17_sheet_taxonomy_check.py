"""Does the master collection sheet's taxonomy describe the genomes we predict for?

For each lab genome with a strain number (s14_lab_genome_ids.csv) this compares:
  - genome: nearest BacDive genome by protein-family Jaccard (s14), genus/phylum from BacDive
  - sheet:  `phylogeny` (Greengenes-style 16S string) of that strain's STOCK in the master
            sheet "Strain collection including Schulenberg's - collection.csv"
  - sheet Sanger 16S vs the genome's own 16S: `sanger_vs_ngs_pct`, % divergence, from
            external_reference_crosscheck/results/per_well_concordance.csv
  - genome_status from the 20261021 concordance: does the WELL's own 16S reads match the
            genome? ('ok' = yes, 'wrong' = the well holds a different organism)

Join: genome_id -> strain_no/well (s14 id table) -> stock_no (20261113 strain_ids.csv)
      -> master sheet 'Strain number' (the sheet is keyed by stock).

Verdict per genome:
  agree               sheet genus == BacDive genus (or a listed synonym)
  same_16S_name_differs  genus differs, but the sheet's own Sanger 16S is <=1.5% from the
                      genome's 16S -> same organism; rename, 16S-classifier limit or a wrong
                      phylogeny string
  different_organism  sheet's Sanger 16S >1.5% from the genome's 16S
  no_sanger           genus differs and no Sanger-vs-genome comparison exists
  no_sheet_genus      sheet has no genus for this stock

Result (2026-11-20): 44 different_organism, 40 of them with genome_status 'ok'. 'ok' only
means the well matches the genome. 32/44 are clone-cluster members: C002 Citrobacter braakii
under 9 strain numbers, C004 Klebsiella ornithinolytica under 6. Either the Sanger record is
wrong (the lab knows it has problems; repeated isolation of one clone) or the stock was
overgrown before genome sequencing; these data can't tell which. The disagreeing Sanger reads
are full-length and clean, so read quality does not flag them.

Outputs: outputs/s17_lab_vs_sheet_taxonomy.csv, outputs/s17_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
TAB = ROOT / "data" / "tables"
KARL = Path("/home/rl/scripts/karl")
SHEET = KARL / "data_links" / "resources" / "Strain collection including Schulenberg's - collection.csv"
STRAIN_IDS = (KARL / "sequencing_experiments/pairwise_interaction_experiments/analyses/"
              "20261113/clean/strain_ids.csv")
STATUS = (KARL / "sequencing_experiments/pairwise_interaction_experiments/analyses/"
          "20261021/03_genome_concordance/outputs/g02_strain_genome_status.csv")
XCHECK = (KARL / "sequencing_experiments/initial_database_creation/merge_consensus_sequences/"
          "external_reference_crosscheck/results/per_well_concordance.csv")
SYN = {"Ochrobactrum": "Brucella", "Wautersiella": "Empedobacter", "Raoultella": "Klebsiella"}
PHYLUM = {"Firmicutes": "Bacillota", "Proteobacteria": "Pseudomonadota",
          "Actinobacteria": "Actinomycetota", "Bacteroidetes": "Bacteroidota"}
SAME_16S_PCT = 1.5


def main() -> None:
    near = pd.read_csv(OUT / "s14_lab_nearest_bacdive.csv", dtype={"strain_no": str})
    hq = pd.read_csv(TAB / "genome_set_hq.csv", usecols=["bacdive_id", "phylum"])
    near = near.merge(hq.rename(columns={"bacdive_id": "nearest_bacdive_id",
                                         "phylum": "bacdive_phylum"}),
                      on="nearest_bacdive_id", how="left")
    sid = pd.read_csv(STRAIN_IDS, dtype=str)[["well", "stock_no"]]
    sheet = pd.read_csv(SHEET, dtype=str)
    sheet["stock_no"] = sheet["Strain number"].str.strip()
    ph = sheet.phylogeny.fillna("")
    sheet["sheet_genus"] = ph.str.extract(r"g__([^;]*)")[0].str.strip().replace("", np.nan)
    sheet["sheet_phylum"] = ph.str.extract(r"p__([^;]*)")[0].str.strip().replace(PHYLUM)
    st = pd.read_csv(STATUS)[["strain", "genome_status"]].rename(columns={"strain": "well"})
    xc = pd.read_csv(XCHECK)[["well", "sanger_vs_ngs_pct"]]

    t = (near.merge(sid, on="well", how="left", validate="1:1")
         .merge(sheet[["stock_no", "sheet_genus", "sheet_phylum", "phylogeny"]]
                .drop_duplicates("stock_no"), on="stock_no", how="left")
         .merge(st, on="well", how="left").merge(xc, on="well", how="left"))
    same = [isinstance(a, str) and (a == b or SYN.get(a) == b)
            for a, b in zip(t.sheet_genus, t.genus)]
    t["verdict"] = np.select(
        [t.sheet_genus.isna(), same, t.sanger_vs_ngs_pct <= SAME_16S_PCT,
         t.sanger_vs_ngs_pct > SAME_16S_PCT],
        ["no_sheet_genus", "agree", "same_16S_name_differs", "different_organism"], "no_sanger")
    t["phylum_agree"] = np.where(t.sheet_phylum.isna(), np.nan,
                                 (t.sheet_phylum == t.bacdive_phylum).astype(float))
    cols = ["strain_no", "well", "stock_no", "genome_id", "clone_cluster", "verdict",
            "sheet_genus", "genus", "species", "jaccard", "sheet_phylum", "bacdive_phylum",
            "phylum_agree", "sanger_vs_ngs_pct", "genome_status", "phylogeny"]
    t[cols].to_csv(OUT / "s17_lab_vs_sheet_taxonomy.csv", index=False)

    h = t[t.sheet_genus.notna()]
    num = pd.to_numeric(h.strain_no, errors="coerce")
    S = {
        "n_lab_genomes": len(t),
        "n_with_sheet_genus": len(h),
        "verdicts": t.verdict.value_counts().to_dict(),
        "genus_agreement": round(float(h.verdict.eq("agree").mean()), 3),
        "phylum_agreement": round(float(h.phylum_agree.mean()), 3),
        "genus_agreement_by_genome_status": h.assign(a=h.verdict.eq("agree"))
            .groupby(h.genome_status.fillna("NA")).a.agg(["size", "mean"]).round(3)
            .to_dict(orient="index"),
        "different_organism_by_genome_status": t[t.verdict == "different_organism"]
            .genome_status.fillna("NA").value_counts().to_dict(),
        "genus_agreement_by_strain_no_bin": h.assign(a=h.verdict.eq("agree"))
            .groupby(pd.cut(num, [0, 50, 100, 150, 200, 250, 300]), observed=True).a
            .agg(["size", "mean"]).round(2).rename(index=str).to_dict(orient="index"),
        "median_sanger_vs_genome_pct_if_different": float(
            t.loc[t.verdict == "different_organism", "sanger_vs_ngs_pct"].median()),
    }
    (OUT / "s17_summary.json").write_text(json.dumps(S, indent=1, default=str))
    print(json.dumps(S, indent=1, default=str))
    print(t[t.verdict.isin(["different_organism", "same_16S_name_differs"])]
          .sort_values(["verdict", "strain_no"])[
              ["strain_no", "well", "verdict", "sheet_genus", "species", "sanger_vs_ngs_pct",
               "genome_status"]].to_string(index=False))


if __name__ == "__main__":
    main()
