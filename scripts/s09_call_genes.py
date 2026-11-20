"""Call genes with Prodigal 2.6.3 on every genome: BacDive HQ set + our 298 lab genomes.

Same gene caller on both sides so protein families transfer cleanly (NCBI PGAP
proteins would differ systematically from Prodigal calls on our genomes).

Genome IDs
  bd<bacdive_id>    BacDive genome (data/tables/genome_set_hq.csv)
  lab_<gbk stem>    our genome; contigs extracted from the panX/Prokka GenBank files

Protein headers are rewritten to  >{genome_id}|{index}  (trailing '*' stripped),
so MMseqs2 cluster members map straight back to genomes.

Outputs
  data/interim/proteins/<genome_id>.faa.gz
  data/tables/genomes_all.csv               genome_id, source, n_proteins, n_aa, path
  data/interim/all_proteins.faa             concatenation for mmseqs createdb
"""
from __future__ import annotations

import gzip
import io
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "data" / "tables"
NCBI = ROOT / "data" / "raw" / "ncbi_genomes" / "ncbi_dataset" / "data"
PROT = ROOT / "data" / "interim" / "proteins"
ALL = ROOT / "data" / "interim" / "all_proteins.faa"
LAB_GBK = Path("/home/rl/Documents/Workstation_Data/Or/KEIO_collection_panX_mapping/"
               "STRAIN_LIBRARY_SEQ_PROCCSSED_DATA/all_gbk_files_including_CeMbio43/input_GenBank")
PRODIGAL = "/home/rl/mambaforge/bin/prodigal"
WORKERS = 60


def gbk_to_fasta(path: Path) -> bytes:
    from Bio import SeqIO
    out = io.StringIO()
    for i, rec in enumerate(SeqIO.parse(str(path), "genbank")):
        out.write(f">{rec.id or f'contig{i}'}\n{str(rec.seq)}\n")
    return out.getvalue().encode()


def run_one(genome_id: str, src: str) -> dict:
    out = PROT / f"{genome_id}.faa.gz"
    if out.exists() and out.stat().st_size > 0:
        n = n_aa = 0
        with gzip.open(out, "rt") as f:
            for line in f:
                if line.startswith(">"):
                    n += 1
                else:
                    n_aa += len(line.strip())
        return {"genome_id": genome_id, "n_proteins": n, "n_aa": n_aa, "cached": True}
    if src.endswith(".gbk"):
        fasta = gbk_to_fasta(Path(src))
    else:
        with gzip.open(src, "rb") as f:
            fasta = f.read()
    tmp_faa = PROT / f".{genome_id}.faa.tmp"
    r = subprocess.run([PRODIGAL, "-p", "single", "-q", "-a", str(tmp_faa), "-o", "/dev/null"],
                       input=fasta, capture_output=True)
    if r.returncode != 0:
        raise RuntimeError(f"{genome_id}: prodigal failed: {r.stderr.decode()[:500]}")
    n = n_aa = 0
    seq: list[str] = []
    tmp_gz = PROT / f".{genome_id}.faa.gz.tmp"
    with open(tmp_faa) as fi, gzip.open(tmp_gz, "wt", compresslevel=3) as fo:
        def flush():
            nonlocal n_aa
            s = "".join(seq).rstrip("*")
            if s:
                fo.write(f">{genome_id}|{n}\n{s}\n")
                n_aa += len(s)
        for line in fi:
            if line.startswith(">"):
                if seq:
                    flush()
                    seq.clear()
                n += 1
            else:
                seq.append(line.strip())
        if seq:
            flush()
    tmp_faa.unlink()
    tmp_gz.rename(out)
    return {"genome_id": genome_id, "n_proteins": n, "n_aa": n_aa, "cached": False}


def main() -> None:
    PROT.mkdir(parents=True, exist_ok=True)
    hq = pd.read_csv(TAB / "genome_set_hq.csv")
    jobs = []
    for bid, acc in zip(hq.bacdive_id, hq.ncbi_accession):
        # some accessions also ship *_cds_from_genomic / *_rna_from_genomic files
        f = [p for p in sorted((NCBI / acc).glob("*_genomic.fna.gz"))
             if "_cds_from_genomic" not in p.name and "_rna_from_genomic" not in p.name]
        if len(f) > 1:
            raise RuntimeError(f"{acc}: ambiguous genome files {f}")
        if f:
            jobs.append((f"bd{bid}", str(f[0]), "bacdive", bid, acc))
    n_missing = len(hq) - len(jobs)
    for g in sorted(LAB_GBK.glob("*.gbk")):
        jobs.append((f"lab_{g.stem}", str(g), "lab", None, g.stem))
    print(f"{len(jobs)} genomes ({n_missing} BacDive genomes missing on disk)", flush=True)

    meta = {j[0]: {"source": j[2], "bacdive_id": j[3], "accession_or_file": j[4], "path": j[1]}
            for j in jobs}
    rows = []
    with ProcessPoolExecutor(WORKERS) as ex:
        futs = {ex.submit(run_one, j[0], j[1]): j[0] for j in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            try:
                rows.append(fut.result())
            except Exception as e:  # keep going; report at the end
                print("FAILED", futs[fut], e, file=sys.stderr, flush=True)
            if i % 500 == 0:
                print(f"  {i}/{len(jobs)}", flush=True)
    df = pd.DataFrame(rows)
    df = df.assign(**{k: df.genome_id.map(lambda g, k=k: meta[g][k])
                      for k in ["source", "bacdive_id", "accession_or_file", "path"]})
    df = df.sort_values("genome_id")
    df.to_csv(TAB / "genomes_all.csv", index=False)
    print(df.groupby("source")[["n_proteins", "n_aa"]].describe().T.round(0).to_string())

    # concatenate for mmseqs (plain text; deleted after createdb by s10)
    with open(ALL, "wb") as fo:
        for g in df.genome_id:
            with gzip.open(PROT / f"{g}.faa.gz", "rb") as fi:
                fo.write(fi.read())
    print(f"wrote {ALL} ({ALL.stat().st_size/1e9:.1f} GB)")


if __name__ == "__main__":
    main()
