"""Genome -> phenotype baselines on protein-family presence, taxon-grouped CV.

Per label (351 usable BacDive labels), all on the SAME folds:
  tax      frequency of the label among training strains in the nearest shared taxon
           (species -> genus -> family -> order -> class -> phylum -> global). With genus-grouped folds the
           test genus is never in training, so this is the genus-out taxonomy baseline.
  knn      k=10 nearest training genomes by Jaccard on family presence, similarity-
           weighted label mean. Relatedness only, at genome resolution.
  logreg   L2 logistic regression on binary family presence. C picked ONCE by inner
           genus-grouped CV inside the training part of fold 0 on 40 random labels.
If logreg > knn, the model learns gene -> function beyond relatedness.

Metric: AUROC of pooled out-of-fold predictions per label.

Usage: python s12_baseline_models.py [--tag id50_c80] [--split genus|family|strain]
Outputs: outputs/s12_<tag>_<split>_auroc_per_label.csv, ..._summary.json,
         ..._oof.csv.gz (long format), data/tables/cv_folds_<split>.csv, figure.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from joblib import Parallel, delayed
from scipy.stats import wilcoxon
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "data" / "tables"
OUT = ROOT / "outputs"
SEED = 20261120
N_FOLDS = 5
MIN_PREV = 10          # family must occur in >= this many BacDive genomes
K_NN = 10
C_GRID = [0.003, 0.01, 0.03, 0.1, 0.3]
# nearest shared taxon, finest first. Under genus/family-grouped folds the test strain's
# species+genus (+family) never occur in training, so those levels fall through.
TAX_LEVELS = ["species", "genus", "family", "order", "class", "phylum"]
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"


# ---------------------------------------------------------------- data
def load_presence(tag: str):
    """All genomes (BacDive + lab) x families, binary, restricted to families with
    BacDive prevalence in [MIN_PREV, 99%]. The column filter uses BacDive genomes
    only (label-free), so lab genomes get exactly the training feature space."""
    fdir = TAB / f"families_{tag}"
    counts = sp.load_npz(fdir / "counts.npz").tocsr()
    genomes = np.array((fdir / "genomes.txt").read_text().split())
    # keep lab genomes + BacDive genomes in the current QC set (genome_set_hq.csv);
    # genomes that were clustered but later failed QC are dropped here
    hq_ids = set("bd" + pd.read_csv(TAB / "genome_set_hq.csv").bacdive_id.astype(str))
    rows = np.array([g.startswith("lab_") or g in hq_ids for g in genomes])
    counts, genomes = counts[rows], genomes[rows]
    is_bd = np.char.startswith(genomes, "bd")
    P = (counts > 0).astype(np.float32).tocsr()
    prev = np.asarray(P[is_bd].sum(axis=0)).ravel()
    keep = (prev >= MIN_PREV) & (prev <= 0.99 * is_bd.sum())
    return P[:, keep].tocsr(), genomes, is_bd, np.where(keep)[0]


def load(tag: str):
    P, genomes, is_bd, keep = load_presence(tag)
    X = P[is_bd]
    ids = np.array([int(g[2:]) for g in genomes[is_bd]])
    labels = pd.read_csv(TAB / "label_matrix_mu_grouped_minority25_genome_strains.csv.gz",
                         index_col=0)
    hq = pd.read_csv(TAB / "genome_set_hq.csv").set_index("bacdive_id")
    order = [i for i in ids if i in labels.index and i in hq.index]
    pos = {b: k for k, b in enumerate(ids)}
    X = X[[pos[i] for i in order]]
    Y = labels.loc[order]
    tax = hq.loc[order, TAX_LEVELS].copy()
    return X, Y, tax, len(keep)


def make_folds(tax: pd.DataFrame, split: str) -> np.ndarray:
    rng = np.random.default_rng(SEED)
    if split == "strain":
        return rng.permutation(np.arange(len(tax)) % N_FOLDS)
    grp = tax[split].fillna(pd.Series(["_na_" + str(i) for i in tax.index], index=tax.index))
    sizes = grp.value_counts()
    keys = sizes.index.to_numpy()
    keys = keys[rng.permutation(len(keys))]
    keys = sorted(keys, key=lambda k: -sizes[k])  # greedy: biggest groups first
    load_ = np.zeros(N_FOLDS)
    assign = {}
    for k in keys:
        f = int(np.argmin(load_))
        assign[k] = f
        load_[f] += sizes[k]
    return grp.map(assign).to_numpy()


# ---------------------------------------------------------------- models
def tax_predict(tax: pd.DataFrame, y: np.ndarray, tr: np.ndarray, te: np.ndarray) -> np.ndarray:
    pred = np.full(len(te), np.nan)
    ttr = tax.iloc[tr]
    for lvl in TAX_LEVELS:
        rate = pd.Series(y[tr]).groupby(ttr[lvl].to_numpy()).mean()
        v = tax.iloc[te][lvl].map(rate).to_numpy(dtype=float)
        fill = np.isnan(pred) & ~np.isnan(v)
        pred[fill] = v[fill]
    pred[np.isnan(pred)] = y[tr].mean()
    return pred


def knn_predict(J: np.ndarray, y: np.ndarray, tr: np.ndarray, te: np.ndarray) -> np.ndarray:
    s = J[np.ix_(te, tr)]
    k = min(K_NN, len(tr))
    idx = np.argpartition(-s, k - 1, axis=1)[:, :k]
    w = np.take_along_axis(s, idx, axis=1) + 1e-6
    return (w * y[tr][idx]).sum(axis=1) / w.sum(axis=1)


def logreg_predict(X, y, tr, te, C) -> np.ndarray:
    if len(np.unique(y[tr])) < 2:
        return np.full(len(te), y[tr].mean())
    m = LogisticRegression(C=C, solver="liblinear", max_iter=500)
    m.fit(X[tr], y[tr])
    return m.predict_proba(X[te])[:, 1]


def run_label(j, X, Yv, tax, folds, J, C):
    mask = ~np.isnan(Yv[:, j])
    rows = np.where(mask)[0]
    y_full = np.zeros(len(Yv))
    y_full[rows] = Yv[rows, j]
    preds = {m: np.full(len(rows), np.nan) for m in ["tax", "knn", "logreg"]}
    for f in range(N_FOLDS):
        tr = rows[folds[rows] != f]
        te = rows[folds[rows] == f]
        if len(te) == 0:
            continue
        sel = folds[rows] == f
        preds["tax"][sel] = tax_predict(tax, y_full, tr, te)
        preds["knn"][sel] = knn_predict(J, y_full, tr, te)
        preds["logreg"][sel] = logreg_predict(X, y_full, tr, te, C)
    return j, rows, preds


def fold_auroc(y: np.ndarray, p: np.ndarray, f: np.ndarray, min_class: int = 3) -> float:
    """AUROC within each fold, averaged weighted by fold size.

    PRIMARY METRIC. Pooling out-of-fold predictions is biased under taxon-grouped
    folds: label prevalence differs between folds, a model trained on fewer
    positives predicts lower for the fold that holds more of them, and pooling
    turns that into anti-correlation. On pure-noise features pooled AUROC came out
    0.43 (s12 --tag synthetic). Within one fold the intercept is constant, so
    this artifact cancels."""
    aucs, ws = [], []
    for k in np.unique(f):
        m = f == k
        npos = int(y[m].sum())
        if npos >= min_class and m.sum() - npos >= min_class:
            aucs.append(roc_auc_score(y[m], p[m]))
            ws.append(m.sum())
    return float(np.average(aucs, weights=ws)) if aucs else np.nan


def jaccard(X) -> np.ndarray:
    Xb = X.astype(np.float32)
    inter = (Xb @ Xb.T).toarray()
    sz = np.asarray(Xb.sum(axis=1)).ravel()
    J = inter / (sz[:, None] + sz[None, :] - inter)
    np.fill_diagonal(J, -1)  # never your own neighbour
    return J.astype(np.float32)


def pick_C(X, Yv, tax, folds) -> tuple[float, dict]:
    """Inner genus-grouped 3-fold CV inside training part of outer fold 0."""
    rng = np.random.default_rng(SEED)
    dev = rng.choice(Yv.shape[1], size=min(40, Yv.shape[1]), replace=False)
    outer_tr = np.where(folds != 0)[0]
    g = tax.iloc[outer_tr]["genus"].fillna("_na").to_numpy()
    ug = rng.permutation(np.unique(g))
    inner = pd.Series(np.arange(len(ug)) % 3, index=ug).loc[g].to_numpy()

    def score(C, j):
        y = Yv[outer_tr, j]
        m = ~np.isnan(y)
        r, yy, ff = outer_tr[m], y[m], inner[m]
        p = np.full(len(r), np.nan)
        for f in range(3):
            tr, te = ff != f, ff == f
            if len(np.unique(yy[tr])) < 2 or te.sum() == 0:
                continue
            mdl = LogisticRegression(C=C, solver="liblinear", max_iter=500).fit(X[r[tr]], yy[tr])
            p[te] = mdl.predict_proba(X[r[te]])[:, 1]
        ok = ~np.isnan(p)
        return fold_auroc(yy[ok], p[ok], ff[ok])

    res = Parallel(n_jobs=50)(delayed(score)(C, j) for C in C_GRID for j in dev)
    res = np.array(res).reshape(len(C_GRID), len(dev))
    med = {C: float(np.nanmedian(r)) for C, r in zip(C_GRID, res)}
    return max(med, key=med.get), med


# ---------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="id50_c80")
    ap.add_argument("--split", default="genus", choices=["genus", "family", "strain"])
    ap.add_argument("--max-labels", type=int, default=None, help="quick test on first N labels")
    ap.add_argument("--min-prev", type=int, default=MIN_PREV,
                    help="family prevalence floor (BacDive genomes); 87 = the MLP's feature set")
    a = ap.parse_args()
    t0 = time.time()
    globals()["MIN_PREV"] = a.min_prev
    X, Y, tax, n_feat = load(a.tag)
    if a.max_labels:
        Y = Y.iloc[:, :a.max_labels]
    Yv = Y.to_numpy(dtype=float)
    folds = make_folds(tax, a.split)
    if not a.max_labels:
        pd.DataFrame({"bacdive_id": Y.index, "fold": folds}).to_csv(
            TAB / f"cv_folds_{a.split}.csv", index=False)
    print(f"X {X.shape} ({n_feat} families with prevalence >= {MIN_PREV}), "
          f"labels {Y.shape}, fold sizes {np.bincount(folds).tolist()}", flush=True)

    C, c_scores = pick_C(X, Yv, tax, folds)
    print(f"C grid (median inner AUROC): {c_scores} -> C={C}  "
          f"[{(time.time()-t0)/60:.1f} min]", flush=True)
    J = jaccard(X)
    print(f"jaccard done [{(time.time()-t0)/60:.1f} min]", flush=True)

    out = Parallel(n_jobs=50, verbose=0)(
        delayed(run_label)(j, X, Yv, tax, folds, J, C) for j in range(Yv.shape[1]))
    print(f"models done [{(time.time()-t0)/60:.1f} min]", flush=True)

    rows, long = [], []
    for j, r, preds in out:
        y = Yv[r, j]
        rec = {"label": Y.columns[j], "kind_group": Y.columns[j].split("|")[1],
               "n": len(r), "n_pos": int(y.sum())}
        for m, p in preds.items():
            rec[m] = fold_auroc(y, p, folds[r])            # primary
            rec[f"{m}_pooled"] = roc_auc_score(y, p)       # biased, for reference
        rows.append(rec)
        long.append(pd.DataFrame({"bacdive_id": Y.index[r], "label": Y.columns[j], "y": y,
                                  "fold": folds[r], **{m: p for m, p in preds.items()}}))
    res = pd.DataFrame(rows)
    n_all = len(res)
    res = res.dropna(subset=["tax", "knn", "logreg"])  # labels with no scorable fold
    tag = f"s12_{a.tag}_{a.split}" + (f"_first{a.max_labels}" if a.max_labels else "") \
        + (f"_minprev{a.min_prev}" if a.min_prev != 10 else "")
    res.to_csv(OUT / f"{tag}_auroc_per_label.csv", index=False)
    pd.concat(long).to_csv(OUT / f"{tag}_oof.csv.gz", index=False)

    summ = {"n_strains": int(len(Y)), "n_features": n_feat, "split": a.split, "C": C,
            "metric": "AUROC within fold, weighted mean over folds (pooled in *_pooled cols)",
            "n_labels_scored": len(res), "n_labels_unscorable": n_all - len(res),
            "median_auroc_pooled": res[["tax_pooled", "knn_pooled", "logreg_pooled"]]
                                      .median().round(3).to_dict(),
            "C_grid_inner_median_auroc": c_scores,
            "median_auroc": res[["tax", "knn", "logreg"]].median().round(3).to_dict(),
            "iqr_auroc": {m: res[m].quantile([.25, .75]).round(3).tolist()
                          for m in ["tax", "knn", "logreg"]},
            "median_by_kind_group": res.groupby("kind_group")[["tax", "knn", "logreg"]]
                                       .median().round(3).to_dict(orient="index"),
            "logreg_beats_tax_frac": round(float((res.logreg > res.tax).mean()), 3),
            "logreg_beats_knn_frac": round(float((res.logreg > res.knn).mean()), 3),
            "knn_beats_tax_frac": round(float((res.knn > res.tax).mean()), 3),
            "wilcoxon_p_logreg_vs_tax": float(wilcoxon(res.logreg, res.tax).pvalue),
            "wilcoxon_p_logreg_vs_knn": float(wilcoxon(res.logreg, res.knn).pvalue),
            "minutes": round((time.time() - t0) / 60, 1)}
    (OUT / f"{tag}_summary.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1))

    fig, ax = plt.subplots(figsize=(5.2, 5.0), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    for s in ["left", "bottom"]:
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    lo = min(res.tax.min(), res.logreg.min()) - 0.02
    ax.plot([lo, 1], [lo, 1], color=INK2, lw=1, ls="--")
    ax.scatter(res.tax, res.logreg, s=12, color=C1, alpha=0.7, edgecolor="none")
    ax.set_xlabel("taxonomy-only AUROC", color=INK2)
    ax.set_ylabel("logistic regression on protein families, AUROC", color=INK2)
    ax.set_title(f"{len(res)} BacDive labels, {a.split}-held-out CV", color=INK, fontsize=10.5,
                 loc="left")
    fig.tight_layout()
    fig.savefig(OUT / "figures" / f"{tag}_logreg_vs_tax.png", dpi=160)


if __name__ == "__main__":
    main()
