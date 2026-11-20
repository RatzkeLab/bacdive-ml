"""Cluster all Prodigal proteins (BacDive HQ + lab genomes) into protein families.

MMseqs2 v16 (env eggnog_mapper), cascaded `cluster`:
  --min-seq-id 0.5 -c 0.8 --cov-mode 0   (>=50% identity, >=80% coverage of BOTH
  sequences; roughly UniRef50-like families). Parameters are CLI args so other
  thresholds (0.3 functional / 0.9 strain-level) can be run side by side.

Outputs (data/interim/mmseqs_<tag>/):
  proteinsDB*        mmseqs sequence DB (kept: reusable for re-clustering / eggNOG on reps)
  clu*               cluster DB
  clusters.tsv.gz    representative \t member
  reps.faa.gz        representative sequences (for later annotation of families)

Usage: python s10_mmseqs_cluster.py [--id 0.5] [--cov 0.8] [--threads 60]
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INTERIM = ROOT / "data" / "interim"
ALL = INTERIM / "all_proteins.faa"
MMSEQS = "/home/rl/mambaforge/envs/eggnog_mapper/bin/mmseqs"


def run(args: list[str]) -> None:
    t0 = time.time()
    print("+", " ".join(args), flush=True)
    subprocess.run(args, check=True)
    print(f"  done in {(time.time() - t0) / 60:.1f} min", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--id", type=float, default=0.5)
    ap.add_argument("--cov", type=float, default=0.8)
    ap.add_argument("--threads", type=int, default=60)
    a = ap.parse_args()
    tag = f"id{int(a.id * 100)}_c{int(a.cov * 100)}"
    out = INTERIM / f"mmseqs_{tag}"
    out.mkdir(parents=True, exist_ok=True)
    db_dir = INTERIM / "mmseqs_db"
    db_dir.mkdir(exist_ok=True)
    db = db_dir / "proteinsDB"
    tmp = out / "tmp"

    if not db.exists():
        run([MMSEQS, "createdb", str(ALL), str(db), "--dbtype", "1"])
    clu = out / "clu"
    if not (out / "clu.index").exists():
        run([MMSEQS, "cluster", str(db), str(clu), str(tmp),
             "--min-seq-id", str(a.id), "-c", str(a.cov), "--cov-mode", "0",
             "--threads", str(a.threads), "--split-memory-limit", "150G",
             "--remove-tmp-files", "1"])
    tsv = out / "clusters.tsv"
    if not (out / "clusters.tsv.gz").exists():
        run([MMSEQS, "createtsv", str(db), str(db), str(clu), str(tsv),
             "--threads", str(a.threads)])
        run(["gzip", "-f", "-3", str(tsv)])
    if not (out / "reps.faa.gz").exists():
        run([MMSEQS, "createsubdb", str(clu), str(db), str(out / "repDB")])
        run([MMSEQS, "convert2fasta", str(out / "repDB"), str(out / "reps.faa")])
        run(["gzip", "-f", "-3", str(out / "reps.faa")])
    shutil.rmtree(tmp, ignore_errors=True)
    print("finished", out)


if __name__ == "__main__":
    main()
