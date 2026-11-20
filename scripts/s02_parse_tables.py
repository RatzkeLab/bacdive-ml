"""Flatten cached BacDive v2 JSON batches into long-format tables.

Nothing is collapsed or interpreted here: every list-of-dicts section becomes
one row per item with all its keys kept as columns, plus bacdive_id. Semantics
(what counts as positive, which assay kinds to merge, ...) are decided
downstream in s03, so this step only needs re-running if the download changes.

Outputs (data/tables/*.csv.gz):
    strains                 one row per strain: taxonomy, DSM no., type strain, ...
    metabolite_utilization  ChEBI, metabolite, utilization activity, kind tested
    metabolite_production   ChEBI, metabolite, production
    metabolite_tests        indole / VP / MR ... tests (ChEBI-linked)
    enzymes                 enzyme, activity, EC
    api_kits                one row per (strain, ref, kit, test code, result)
    genomes                 genome assembly accessions + BacDive assembly score
    seq16s                  16S accessions
    traits                  simple categorical traits, long format (trait, value)
    culture_temp / culture_ph / halophily / antibiotic_resistance / isolation_cat
"""
from __future__ import annotations

import gzip
import json
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "bacdive_v2"
OUT = ROOT / "data" / "tables"

PM = "Physiology and metabolism"
SECTION_TABLES = {
    # (section, key) -> table name ; each item -> one row with all keys
    (PM, "metabolite utilization"): "metabolite_utilization",
    (PM, "metabolite production"): "metabolite_production",
    (PM, "metabolite tests"): "metabolite_tests",
    (PM, "enzymes"): "enzymes",
    (PM, "halophily"): "halophily",
    (PM, "antibiotic resistance"): "antibiotic_resistance",
    (PM, "tolerance"): "tolerance",
    ("Sequence information", "Genome sequences"): "genomes",
    ("Sequence information", "16S sequences"): "seq16s",
    ("Culture and growth conditions", "culture temp"): "culture_temp",
    ("Culture and growth conditions", "culture pH"): "culture_ph",
    ("Isolation, sampling and environmental information",
     "isolation source categories"): "isolation_cat",
}

# simple traits -> long table (trait, value)
TRAIT_FIELDS = [
    (PM, "oxygen tolerance", "oxygen tolerance"),
    (PM, "spore formation", "spore formation"),
    (PM, "nutrition type", "type"),
    ("Morphology", "cell morphology", "gram stain"),
    ("Morphology", "cell morphology", "cell shape"),
    ("Morphology", "cell morphology", "motility"),
    ("Morphology", "pigmentation", "color"),
]


def as_list(v):
    if isinstance(v, list):
        return [x for x in v if isinstance(x, dict)]
    if isinstance(v, dict):
        return [v]
    return []


def strain_row(bid: int, rec: dict) -> dict:
    g = rec.get("General", {}) or {}
    tx = rec.get("Name and taxonomic classification", {}) or {}
    lpsn = tx.get("LPSN", {}) if isinstance(tx.get("LPSN"), dict) else {}
    ncbi = g.get("NCBI tax id")
    ncbi = as_list(ncbi)
    lit = rec.get("Literature", {}) or {}
    ccn = lit.get("culture collection no.") if isinstance(lit, dict) else None
    row = {
        "bacdive_id": bid,
        "dsm_number": g.get("DSM-Number"),
        "description": g.get("description"),
        "keywords": ";".join(g["keywords"]) if isinstance(g.get("keywords"), list)
        else g.get("keywords"),
        "ncbi_taxid": ncbi[0].get("NCBI tax id") if ncbi else None,
        "ncbi_taxid_level": ncbi[0].get("Matching level") if ncbi else None,
        "type_strain": tx.get("type strain"),
        "strain_designation": tx.get("strain designation"),
        "culture_collection_nos": ccn,
    }
    for lvl in ["domain", "phylum", "class", "order", "family", "genus", "species"]:
        # LPSN is the curated, up-to-date taxonomy; fall back to BacDive's own
        row[lvl] = lpsn.get(lvl) or tx.get(lvl)
        row[f"{lvl}_bacdive"] = tx.get(lvl)
    return row


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    files = sorted(RAW.glob("batch_*.json.gz"))
    print(f"{len(files)} batch files")
    strains = []
    tables: dict[str, list] = defaultdict(list)
    api_rows = []
    trait_rows = []
    pm_keys = defaultdict(int)
    for i, fp in enumerate(files):
        with gzip.open(fp, "rt") as f:
            res = json.load(f).get("results") or {}
        for sid, rec in res.items():
            bid = int(sid)
            strains.append(strain_row(bid, rec))
            for (sec, key), tname in SECTION_TABLES.items():
                for it in as_list((rec.get(sec) or {}).get(key)):
                    tables[tname].append({"bacdive_id": bid, **it})
            pm = rec.get(PM) or {}
            for k, v in pm.items():
                pm_keys[k] += 1
                if k.startswith("API "):
                    for it in as_list(v):
                        ref = it.get("@ref")
                        for code, val in it.items():
                            if code == "@ref":
                                continue
                            api_rows.append((bid, ref, k, code, val))
            for sec, key, field in TRAIT_FIELDS:
                for it in as_list((rec.get(sec) or {}).get(key)):
                    if field in it and it[field] not in (None, ""):
                        trait_rows.append((bid, it.get("@ref"), f"{key}:{field}"
                                           if key != field else key, it[field]))
        if (i + 1) % 500 == 0:
            print(f"  parsed {i+1}/{len(files)} files, {len(strains)} strains")

    def write(df: pd.DataFrame, name: str) -> None:
        df.to_csv(OUT / f"{name}.csv.gz", index=False)
        print(f"{name:24s} {df.shape}")

    write(pd.DataFrame(strains).sort_values("bacdive_id"), "strains")
    for tname, rows in tables.items():
        write(pd.DataFrame(rows), tname)
    write(pd.DataFrame(api_rows, columns=["bacdive_id", "ref", "kit", "test", "result"]),
          "api_kits")
    write(pd.DataFrame(trait_rows, columns=["bacdive_id", "ref", "trait", "value"]),
          "traits")
    pd.Series(pm_keys).sort_values(ascending=False).rename("n_strains").to_csv(
        OUT / "_physiology_section_counts.csv")


if __name__ == "__main__":
    main()
