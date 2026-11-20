"""Pilot: download a handful of BacDive genomes from NCBI to test the route.

Picks N high-quality phenotyped strains (current assembly, CheckM >=90/<=10,
ANI check not failed, RefSeq partner present) preferring genera in our lab
collection, writes their RefSeq (GCF) accessions to a list and downloads
genome FASTA + protein FASTA + GFF via the NCBI `datasets` CLI. Reports bytes
per genome to extrapolate the full download.

The datasets CLI is not in karl_datasci; we call the binary from the
bw25113_panx env (v18.33.1).

Usage: python s06_pilot_genome_download.py [--n 12]
"""
from __future__ import annotations

import argparse
import json
import subprocess
import time
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "data" / "tables"
OUT = ROOT / "outputs"
DL = ROOT / "data" / "raw" / "ncbi_pilot"
DATASETS = "/home/rl/mambaforge/envs/bw25113_panx/bin/datasets"
PREFER = ["Pseudomonas", "Bacillus", "Microbacterium", "Enterobacter", "Lactococcus",
          "Paenibacillus", "Stenotrophomonas", "Rhodococcus", "Corynebacterium",
          "Acinetobacter", "Leuconostoc", "Serratia"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=12)
    args = ap.parse_args()
    DL.mkdir(parents=True, exist_ok=True)

    m = pd.read_csv(TAB / "genome_choice_ncbi.csv.gz", low_memory=False)
    st = pd.read_csv(TAB / "strains.csv.gz", low_memory=False,
                     usecols=["bacdive_id", "genus", "species", "type_strain"])
    mat = pd.read_csv(TAB / "label_matrix_mu_grouped_minority25_genome_strains.csv.gz",
                      index_col=0)
    nlab = mat.notna().sum(axis=1).rename("n_labels")
    m = m.merge(st, on="bacdive_id").merge(nlab, left_on="bacdive_id", right_index=True)
    hq = m[(m.assembly_status == "current") & (m.checkm_completeness >= 90)
           & (m.checkm_contamination <= 10) & (m.ani_taxonomy_check != "Failed")
           & (m.paired_status == "current")]
    pick = (hq[hq.genus.isin(PREFER)].sort_values("n_labels", ascending=False)
              .drop_duplicates("genus").head(args.n))
    pick.to_csv(OUT / "s06_pilot_genomes.csv", index=False)
    acc_file = DL / "accessions.txt"
    acc_file.write_text("\n".join(pick.paired_accession) + "\n")
    print(pick[["bacdive_id", "species", "paired_accession", "n_labels",
                "genome_size", "checkm_completeness", "checkm_contamination"]].to_string())

    zip_path = DL / "pilot.zip"
    t0 = time.time()
    subprocess.run([DATASETS, "download", "genome", "accession", "--inputfile", str(acc_file),
                    "--include", "genome,protein,gff3", "--filename", str(zip_path),
                    "--no-progressbar"], check=True)
    el = time.time() - t0
    with zipfile.ZipFile(zip_path) as z:
        info = z.infolist()
        z.extractall(DL / "unzipped")
    files = [(i.filename, i.file_size, i.compress_size) for i in info]
    df = pd.DataFrame(files, columns=["file", "bytes", "zipped_bytes"])
    df["type"] = df.file.str.extract(r"(genomic\.fna|protein\.faa|genomic\.gff)")[0]
    per = df.groupby("type")[["bytes", "zipped_bytes"]].sum() / len(pick)
    n_full = int((m.assembly_status == "current").sum())  # m = phenotyped (label-matrix) strains only
    summ = {
        "n_genomes": len(pick), "download_seconds": round(el, 1),
        "zip_bytes": zip_path.stat().st_size,
        "mean_bytes_per_genome_by_type": per.round(0).to_dict(orient="index"),
        "extrapolated_GB_for_all_phenotyped_genome_strains(~11.1k, uncompressed)":
            round(per.bytes.sum() * 11105 / 1e9, 1),
        "extrapolated_GB_proteins_only_11.1k_uncompressed":
            round(per.loc["protein.faa", "bytes"] * 11105 / 1e9, 1) if "protein.faa" in per.index else None,
        "n_phenotyped_strains_with_current_genome": n_full,
    }
    n_prot = {}
    for faa in (DL / "unzipped").rglob("protein.faa"):
        n_prot[faa.parent.name] = sum(1 for line in open(faa) if line.startswith(">"))
    summ["proteins_per_genome"] = n_prot
    (OUT / "s06_pilot_summary.json").write_text(json.dumps(summ, indent=1, default=str))
    print(json.dumps(summ, indent=1, default=str))


if __name__ == "__main__":
    main()
