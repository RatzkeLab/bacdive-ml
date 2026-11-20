"""Sanity check: do the logistic-regression models weight the RIGHT genes?

For selected labels with known genetics (L-arabinose -> araBAD, nitrate reduction ->
nar/nap, sulfate reduction -> dsr/apr/sat, ...), fit the s12 logistic regression on
all QC'd BacDive strains, take the top positive / negative coefficient families, and
name each family's representative protein by MMseqs2 search against Swiss-Prot
(data/raw/swissprot/uniprot_sprot.fasta.gz, release in reldate.txt).

Also reports how often the family occurs in label-positive vs label-negative strains.

Usage: python s15_interpret_families.py [--tag id50_c80] [--top 12]
Outputs: outputs/s15_<tag>_top_families.csv, outputs/s15_<tag>_top_families.md
"""
from __future__ import annotations

import argparse
import gzip
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegression

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s12_baseline_models import load_presence  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "data" / "tables"
OUT = ROOT / "outputs"
SPROT = ROOT / "data" / "raw" / "swissprot" / "uniprot_sprot.fasta.gz"
MMSEQS = "/home/rl/mambaforge/envs/eggnog_mapper/bin/mmseqs"
LABELS = [
    "l-arabinose|growth/assimilation",
    "l-arabinose|acid/gas (fermentation)",
    "nitrate|reduction/e- acceptor",
    "nitrite|reduction/e- acceptor",
    "sulfate|reduction/e- acceptor",
    "decanoate|growth/assimilation",
    "d-xylose|acid/gas (fermentation)",
    "l-rhamnose|acid/gas (fermentation)",
    "d-galacturonic acid|acid/gas (fermentation)",
    "tryptophan|growth/assimilation",
    "urea|hydrolysis/degradation",
    "citrate|growth/assimilation",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="id50_c80")
    ap.add_argument("--top", type=int, default=12)
    a = ap.parse_args()
    P, genomes, is_bd, keep = load_presence(a.tag)
    X = P[is_bd]
    ids = np.array([int(g[2:]) for g in genomes[is_bd]])
    Y = pd.read_csv(TAB / "label_matrix_mu_grouped_minority25_genome_strains.csv.gz",
                    index_col=0).reindex(ids)
    fams = pd.read_csv(TAB / f"families_{a.tag}" / "families.csv.gz")
    C = json.loads((OUT / f"s12_{a.tag}_genus_summary.json").read_text())["C"]
    labels = [l for l in LABELS if l in Y.columns]

    def fit(col):
        y = Y[col].to_numpy()
        m = ~np.isnan(y)
        mdl = LogisticRegression(C=C, solver="liblinear", max_iter=1000).fit(X[m], y[m])
        coef = mdl.coef_.ravel()
        order = np.argsort(coef)
        sel = list(order[::-1][:a.top]) + list(order[:4])
        Xm = X[m]
        pos, neg = y[m] == 1, y[m] == 0
        rows = []
        for j in sel:
            col_j = Xm[:, j].toarray().ravel()
            rows.append({"label": col, "coef": coef[j], "fam": int(keep[j]),
                         "frac_in_pos": col_j[pos].mean(), "frac_in_neg": col_j[neg].mean()})
        return rows

    res = pd.DataFrame([r for rows in Parallel(n_jobs=len(labels))(
        delayed(fit)(c) for c in labels) for r in rows])
    res = res.merge(fams[["fam", "rep", "n_genomes"]], on="fam")

    # representative sequences -> Swiss-Prot
    work = OUT / f"s15_{a.tag}_tmp"
    work.mkdir(exist_ok=True)
    want = set(res.rep)
    q = work / "query.faa"
    with gzip.open(ROOT / "data" / "interim" / f"mmseqs_{a.tag}" / "reps.faa.gz", "rt") as fi, \
            open(q, "w") as fo:
        write = False
        for line in fi:
            if line.startswith(">"):
                write = line[1:].split()[0] in want
            if write:
                fo.write(line)
    hits = work / "hits.m8"
    subprocess.run([MMSEQS, "easy-search", str(q), str(SPROT), str(hits), str(work / "tmp"),
                    "--format-output", "query,target,pident,evalue,qcov,theader",
                    "-e", "1e-5", "--max-seqs", "50", "--threads", "16", "-v", "1"], check=True)
    h = pd.read_csv(hits, sep="\t", header=None,
                    names=["rep", "sp_id", "pident", "evalue", "qcov", "header"])
    h = h.sort_values(["rep", "evalue"]).drop_duplicates("rep")
    h["gene"] = h.header.str.extract(r"GN=(\S+)")[0]
    h["protein"] = h.header.str.extract(r"^\S+ (.*?) OS=")[0]
    res = res.merge(h[["rep", "gene", "protein", "pident", "evalue"]], on="rep", how="left")
    res.to_csv(OUT / f"s15_{a.tag}_top_families.csv", index=False)

    lines = ["# Top-weighted protein families per label (logreg on all QC strains)\n",
             f"C={C}; Swiss-Prot best hit (e<=1e-5); frac = share of +/- strains carrying it.\n"]
    for lab, d in res.groupby("label", sort=False):
        lines.append(f"\n## {lab}\n")
        lines.append("| coef | gene | protein | %id | frac in + | frac in - | n genomes |")
        lines.append("|---|---|---|---|---|---|---|")
        for _, r in d.iterrows():
            lines.append(f"| {r.coef:+.3f} | {r.gene if isinstance(r.gene, str) else ''} | "
                         f"{r.protein if isinstance(r.protein, str) else '(no hit)'} | "
                         f"{'' if pd.isna(r.pident) else f'{r.pident:.0f}'} | "
                         f"{r.frac_in_pos:.2f} | {r.frac_in_neg:.2f} | {r.n_genomes} |")
    (OUT / f"s15_{a.tag}_top_families.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
