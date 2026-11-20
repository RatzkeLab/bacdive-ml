"""Download genome FASTA for the QC-passing phenotyped BacDive strains.

Selection (docs/06): strain is in the s03 label matrix, its chosen accession is
current at NCBI, CheckM >=90 completeness and <=10 contamination (CheckM must
exist), NCBI ANI taxonomy check not 'Failed'.

Uses the NCBI datasets CLI (env bw25113_panx) with the dehydrated workflow NCBI
recommends for >1000 genomes, and rehydrates with --gzip so the FASTA stays
compressed on disk (~1.2 MB per genome). Downloads the versioned GCA accession
(BacDive's own); the sequence is identical to its GCF partner.

Outputs
  data/tables/genome_set_hq.csv            selected strains + accession + taxonomy
  data/raw/ncbi_genomes/ncbi_dataset/data/<acc>/*_genomic.fna.gz
"""
from __future__ import annotations

import subprocess
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "data" / "tables"
DL = ROOT / "data" / "raw" / "ncbi_genomes"
DATASETS = "/home/rl/mambaforge/envs/bw25113_panx/bin/datasets"


def select() -> pd.DataFrame:
    m = pd.read_csv(TAB / "genome_choice_ncbi.csv.gz", low_memory=False)
    mat_ids = pd.read_csv(TAB / "label_matrix_mu_grouped_minority25_genome_strains.csv.gz",
                          usecols=[0]).bacdive_id
    st = pd.read_csv(TAB / "strains.csv.gz", low_memory=False,
                     usecols=["bacdive_id", "type_strain", "phylum", "class", "order",
                              "family", "genus", "species"])
    m = m[m.bacdive_id.isin(mat_ids)]
    # CheckM must be PRESENT: the 217 CheckM-less genomes include partial assemblies
    # (5 kb Streptomyces, 68 kb Leucobacter, 721 kb E. coli) and 10 that are >1.6x their
    # genus median size; without CheckM contamination can't be judged. (They were
    # downloaded in the first run, before this rule; downstream steps ignore them.)
    ck_ok = (m.checkm_completeness >= 90) & (m.checkm_contamination <= 10)
    hq = m[(m.assembly_status == "current") & ck_ok & (m.ani_taxonomy_check != "Failed")]
    hq = hq.merge(st, on="bacdive_id")
    cols = ["bacdive_id", "ncbi_accession", "paired_accession", "assembly_level_ncbi",
            "genome_size", "checkm_completeness", "checkm_contamination",
            "ani_taxonomy_check", "type_strain", "phylum", "class", "order", "family",
            "genus", "species"]
    return hq[cols].sort_values("bacdive_id")


def main() -> None:
    DL.mkdir(parents=True, exist_ok=True)
    hq = select()
    hq.to_csv(TAB / "genome_set_hq.csv", index=False)
    print(f"selected {len(hq)} genomes; "
          f"{hq.ncbi_accession.nunique()} unique accessions")
    acc_file = DL / "accessions.txt"
    acc_file.write_text("\n".join(hq.ncbi_accession) + "\n")

    zip_path = DL / "dehydrated.zip"
    if not zip_path.exists():
        subprocess.run([DATASETS, "download", "genome", "accession", "--inputfile",
                        str(acc_file), "--include", "genome", "--dehydrated",
                        "--filename", str(zip_path), "--no-progressbar"], check=True)
    if not (DL / "ncbi_dataset" / "fetch.txt").exists():
        with zipfile.ZipFile(zip_path) as z:
            z.extractall(DL)
    subprocess.run([DATASETS, "rehydrate", "--directory", str(DL), "--gzip",
                    "--max-workers", "12", "--no-progressbar"], check=True)

    fna = sorted((DL / "ncbi_dataset" / "data").glob("GC*/*_genomic.fna.gz"))
    got = {p.parent.name for p in fna}
    missing = sorted(set(hq.ncbi_accession) - got)
    print(f"genome files on disk: {len(fna)}; missing: {len(missing)}")
    (DL / "missing_accessions.txt").write_text("\n".join(missing) + "\n")


if __name__ == "__main__":
    main()
