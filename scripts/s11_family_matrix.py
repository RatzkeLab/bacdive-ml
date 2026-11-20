"""Genome x protein-family matrix from MMseqs2 clusters.

Reads data/interim/mmseqs_<tag>/clusters.tsv.gz (rep \t member; member ids are
'<genome_id>|<n>') and writes copy-number counts as a sparse matrix.

Outputs (data/tables/families_<tag>/):
  counts.npz          scipy.sparse CSR int32, genomes x families (copy number)
  genomes.txt         row order (genome_id)
  families.csv.gz     fam (col index), rep, n_members, n_genomes, n_bd, n_lab
  summary.json        family-size spectrum + transfer coverage of lab genomes
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
INTERIM = ROOT / "data" / "interim"
TAB = ROOT / "data" / "tables"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="id50_c80")
    a = ap.parse_args()
    src = INTERIM / f"mmseqs_{a.tag}" / "clusters.tsv.gz"
    out = TAB / f"families_{a.tag}"
    out.mkdir(parents=True, exist_ok=True)

    print("reading", src, flush=True)
    df = pd.read_csv(src, sep="\t", header=None, names=["rep", "member"],
                     dtype=str, engine="c")
    print(f"{len(df):,} proteins", flush=True)
    genome = df.member.str.split("|", n=1, regex=False).str[0]
    fam_idx, fam_rep = pd.factorize(df.rep, sort=False)
    gen_idx, gen_ids = pd.factorize(genome, sort=True)
    del df, genome
    counts = sp.coo_matrix((np.ones(len(fam_idx), dtype=np.int32), (gen_idx, fam_idx)),
                           shape=(len(gen_ids), len(fam_rep))).tocsr()
    counts.sum_duplicates()
    sp.save_npz(out / "counts.npz", counts)
    (out / "genomes.txt").write_text("\n".join(gen_ids) + "\n")

    is_lab = np.array([g.startswith("lab_") for g in gen_ids])
    pres = (counts > 0).astype(np.int8).tocsc()
    n_gen = np.asarray(pres.sum(axis=0)).ravel()
    n_lab = np.asarray(pres[is_lab].sum(axis=0)).ravel()
    n_bd = n_gen - n_lab
    n_mem = np.asarray(counts.sum(axis=0)).ravel()
    fams = pd.DataFrame({"fam": np.arange(len(fam_rep)), "rep": fam_rep, "n_members": n_mem,
                         "n_genomes": n_gen, "n_bd": n_bd, "n_lab": n_lab})
    fams.to_csv(out / "families.csv.gz", index=False)

    # transfer coverage: for each lab genome, share of its proteins / families that
    # fall in families seen in >= k BacDive genomes (only those can carry learned signal)
    lab_rows = counts[is_lab].tocsr()
    cov = {}
    for k in [1, 5, 20]:
        seen = n_bd >= k
        prot_in = np.asarray(lab_rows[:, seen].sum(axis=1)).ravel()
        prot_all = np.asarray(lab_rows.sum(axis=1)).ravel()
        cov[f"lab_protein_frac_in_fams_seen_in_>={k}_bacdive"] = \
            np.quantile(prot_in / prot_all, [0.05, 0.25, 0.5, 0.75, 0.95]).round(3).tolist()
    summ = {
        "n_proteins": int(n_mem.sum()), "n_genomes": len(gen_ids),
        "n_lab_genomes": int(is_lab.sum()), "n_families": len(fam_rep),
        "singleton_families": int((n_mem == 1).sum()),
        "families_in_>=2_genomes": int((n_gen >= 2).sum()),
        "families_in_>=10_bd_genomes": int((n_bd >= 10).sum()),
        "families_in_>=1pct_bd_genomes": int((n_bd >= 0.01 * (~is_lab).sum()).sum()),
        "families_lab_only": int(((n_lab > 0) & (n_bd == 0)).sum()),
        "transfer_coverage_quantiles_[5,25,50,75,95]": cov,
    }
    (out / "summary.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
