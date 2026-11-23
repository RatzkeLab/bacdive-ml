"""Does BacDive pretraining help predict pairwise interactions from genomes on unseen taxa?

Target: mean_log2_ratio_a_over_b of the published modelling pairs (20260907, 20260630), with
wrong-genome wells dropped using each experiment's own 20261113 verdict (setting "D" there).

Folds (as 20261113 lib.group_folds): groups permuted, round-robin into 5 folds; test = pairs
with BOTH strains in fold k, train = pairs with NEITHER. Groups at four levels:
  organism  clone cluster (genome_unit, >=99.9% ANI)
  species   >=95% ANI connected components among the 298 lab genomes (no sheet taxonomy)
  genus     LPSN genus of the nearest BacDive genome (protein-family Jaccard, s14)
  family    LPSN family of the nearest BacDive genome

One architecture for every variant:
  encoder  sparse protein-family presence (the BacDive MLP's 56,462 families, restricted to
           families present in this experiment's strains -- exact, absent columns are all zero)
           -> 512 -> ReLU/Dropout -> 256 -> ReLU/Dropout          (= s13 MTL body)
  head     f([h_a, h_b]) = Linear(512,64) -> ReLU -> Dropout(0.3) -> Linear(64,1)
           y_hat(a,b) = f(a,b) - f(b,a)                             (antisymmetric)
Variants differ ONLY in encoder init / trainability:
  naive_lr*                  random init, trained end-to-end (encoder lr 1e-3 or 1e-4)
  pretrained_ft_lr*          BacDive MLP weights (s13 --final), trained end-to-end
  pretrained_frozen          BacDive encoder frozen, head only
  pretrained_profiles_frozen frozen 351 predicted BacDive phenotypes as the strain vector
  random_frozen              random encoder frozen, head only (control: is "frozen pretrained"
                             better because of pretraining, or because fewer weights train?)
  SHUFFLED_pretrained_frozen training targets permuted (noise floor)
Frozen representations are z-scored on the training strains of each fold.

Training: full-batch AdamW (head lr 1e-3, wd 1e-4, input dropout 0.2 for trainable encoders).
Early stopping on 20% of the TRAINING groups (same level; val pairs = training pairs touching
a val group), then refit on all training pairs for the best step count. The test fold is
never used for any choice.

Metrics per fold (= genomic_ml._metrics): Spearman rho, R^2 vs predicting 0 (no winner),
sign accuracy on |y| > 1. Summaries: mean over folds; error bars = SD across CV repeats;
paired differences vs naive on identical folds.

Usage: python s18_interaction_pretraining.py [--reps 10] [--exps 20260907 20260630]
                                             [--levels organism species genus family]
                                             [--processes 56] [--quick]
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import json
import sys
import time
import zlib
from multiprocessing import get_context
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
import torch.nn as nn
from scipy.sparse.csgraph import connected_components
from scipy.stats import pearsonr, spearmanr, wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parent))
import s12_baseline_models as s12  # noqa: E402
from s12_baseline_models import load_presence  # noqa: E402
from s13_multitask_mlp import MTL, P_HID, to_torch  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
MODELS = ROOT / "data" / "models"
TAG = "id50_c80"
PW = Path("/home/rl/scripts/karl/sequencing_experiments/pairwise_interaction_experiments/analyses")
EXPS = {
    "20260907": PW / "20260923" / "01_model_sweep" / "outputs",
    "20260630": PW / "20260710_increased_tolerance" / "06_ml_analysis" / "genomic_ml" / "outputs_merged_ab",
}
WELLS = PW / "20261113" / "messy" / "outputs" / "s01_wells_{exp}.csv"
ANI = PW / "20261021" / "03_genome_concordance" / "outputs" / "g01_ani_pairs.csv"
TARGET = "mean_log2_ratio_a_over_b"
LEVELS = ["organism", "species", "genus", "family"]
N_SPLITS = 5

VARIANTS = {
    "naive_lr1e-3": dict(init="random", train_enc=True, enc_lr=1e-3),
    "naive_lr1e-4": dict(init="random", train_enc=True, enc_lr=1e-4),
    "pretrained_ft_lr1e-3": dict(init="pretrained", train_enc=True, enc_lr=1e-3),
    "pretrained_ft_lr1e-4": dict(init="pretrained", train_enc=True, enc_lr=1e-4),
    "pretrained_frozen": dict(init="pretrained", train_enc=False, rep="embed"),
    "pretrained_profiles_frozen": dict(init="pretrained", train_enc=False, rep="profiles"),
    "random_frozen": dict(init="random", train_enc=False, rep="embed"),
    "SHUFFLED_pretrained_frozen": dict(init="pretrained", train_enc=False, rep="embed",
                                       shuffle=True),
}
HEAD_HID, HEAD_DROP, HEAD_LR, WD, P_IN = 64, 0.3, 1e-3, 1e-4, 0.2
# best inner-validation steps in the smoke test were ~20-45 for every variant
MAX_STEPS, EVAL_EVERY, PATIENCE, MIN_STEPS, FALLBACK_STEPS, VAL_FRAC = 300, 5, 40, 20, 60, 0.2

_D: dict = {}   # filled in the parent before the pool forks


# ------------------------------------------------------------------ data
def species_clusters(wells: list[str]) -> dict:
    ani = pd.read_csv(ANI)
    e = ani[(ani.ani >= 95) & (ani[["af_a", "af_b"]].min(axis=1) >= 50)]
    ix = {w: i for i, w in enumerate(wells)}
    e = e[e.well_a.isin(ix) & e.well_b.isin(ix)]
    A = sp.coo_matrix((np.ones(len(e)), ([ix[w] for w in e.well_a], [ix[w] for w in e.well_b])),
                      shape=(len(wells), len(wells)))
    _, lab = connected_components(A, directed=False)
    return {w: f"S{l}" for w, l in zip(wells, lab)}


def prepare(exps: list[str]) -> None:
    s12.MIN_PREV = 87
    P, genomes, is_bd, _ = load_presence(TAG)
    ck = torch.load(MODELS / f"mlp_{TAG}_all.pt")
    assert P.shape[1] == ck["n_in"], (P.shape, ck["n_in"])
    Xlab = P[~is_bd].tocsr()
    row = {g: i for i, g in enumerate(genomes[~is_bd])}
    near = pd.read_csv(OUT / "s14_lab_nearest_bacdive.csv").set_index("well")
    tax = near[["genome_id", "genus", "family", "jaccard"]].copy()
    tax["species"] = pd.Series(species_clusters(sorted(tax.index)))

    # pretrained frozen representations (eval mode) for every lab genome
    mtl = MTL(ck["n_in"], len(ck["labels"]))
    mtl.load_state_dict(ck["state_dict"])
    mtl.eval()
    with torch.no_grad():
        xt = to_torch(Xlab)
        emb_all = mtl.embed(xt).numpy()
        prof_all = torch.sigmoid(mtl.head(torch.from_numpy(emb_all))).numpy()
    W = ck["state_dict"]

    for exp in exps:
        pairs = pd.read_csv(EXPS[exp] / "g01_modeling_pairs.csv")
        wells = pd.read_csv(str(WELLS).format(exp=exp)).set_index("well")
        keep = set(wells.index[wells.genome_status != "wrong"]) & set(tax.index)
        pairs = pairs[pairs.strain_a.isin(keep) & pairs.strain_b.isin(keep)].reset_index(drop=True)
        strains = sorted(set(pairs.strain_a) | set(pairs.strain_b))
        sidx = {s: i for i, s in enumerate(strains)}
        rows = [row[tax.loc[s, "genome_id"]] for s in strains]
        Xs = Xlab[rows]
        used = np.unique(Xs.indices)
        groups = {
            "organism": wells.loc[strains, "genome_unit"].astype(str).to_dict(),
            "species": tax.loc[strains, "species"].to_dict(),
            "genus": tax.loc[strains, "genus"].to_dict(),
            "family": tax.loc[strains, "family"].to_dict(),
        }
        _D[exp] = dict(
            pairs=pairs, strains=strains,
            ia=pairs.strain_a.map(sidx).to_numpy(), ib=pairs.strain_b.map(sidx).to_numpy(),
            y=pairs[TARGET].to_numpy(dtype=np.float32),
            Xs=Xs[:, used].tocsr(), groups=groups,
            emb=emb_all[rows], prof=prof_all[rows],
            w1=W["w1"][used].clone(), b1=W["b1"].clone(),
            w2=W["body.2.weight"].clone(), b2=W["body.2.bias"].clone(),
        )
        print(f"{exp}: {len(pairs)} pairs, {len(strains)} strains, {len(used)} families used; "
              + ", ".join(f"{lv} {len(set(g.values()))}" for lv, g in groups.items()), flush=True)


def group_folds(pairs, group_of, seed, n_splits=N_SPLITS):
    """20261113 lib.group_folds: (k, train_idx, test_idx), test both in k, train neither."""
    rng = np.random.default_rng(seed)
    strains = sorted(set(pairs["strain_a"]) | set(pairs["strain_b"]))
    groups = sorted({group_of[s] for s in strains})
    fold_of_group = dict(zip(groups, rng.permutation(len(groups)) % n_splits))
    fa = pairs["strain_a"].map(group_of).map(fold_of_group).values
    fb = pairs["strain_b"].map(group_of).map(fold_of_group).values
    for k in range(n_splits):
        test = np.where((fa == k) & (fb == k))[0]
        train = np.where((fa != k) & (fb != k))[0]
        if len(test) >= 5:
            yield k, train, test


# ------------------------------------------------------------------ model
class Encoder(nn.Module):
    """s13 MTL body on a column subset; random init has the same distribution as MTL's
    (kaiming_uniform on an (n_in, 512) tensor uses fan_in = 512 regardless of n_in)."""

    def __init__(self, n_in: int):
        super().__init__()
        self.w1 = nn.Parameter(torch.empty(n_in, 512))
        nn.init.kaiming_uniform_(self.w1, a=5 ** 0.5)
        self.b1 = nn.Parameter(torch.zeros(512))
        self.body = nn.Sequential(nn.ReLU(), nn.Dropout(P_HID), nn.Linear(512, 256), nn.ReLU(),
                                  nn.Dropout(P_HID))

    def forward(self, xs):
        return self.body(torch.sparse.mm(xs, self.w1) + self.b1)


class Head(nn.Module):
    def __init__(self, d: int):
        super().__init__()
        self.f = nn.Sequential(nn.Linear(2 * d, HEAD_HID), nn.ReLU(), nn.Dropout(HEAD_DROP),
                               nn.Linear(HEAD_HID, 1))

    def forward(self, ha, hb):
        return (self.f(torch.cat([ha, hb], 1)) - self.f(torch.cat([hb, ha], 1))).squeeze(-1)


class SparseRows:
    """Coalesced torch COO for a fixed set of rows, built once; input dropout is a mask on
    the stored values (no re-sorting per step -- rebuilding cost ~1 s/step in the smoke test)."""

    def __init__(self, X):
        t = torch.sparse_coo_tensor(*_coo_parts(X)).coalesce()
        self.idx, self.val, self.shape = t.indices(), t.values(), X.shape

    def get(self, drop=0.0, gen=None):
        val = self.val
        if drop > 0:
            val = val * ((torch.rand(val.shape, generator=gen) >= drop).float() / (1 - drop))
        return torch.sparse_coo_tensor(self.idx, val, self.shape, is_coalesced=True)


def _coo_parts(X):
    c = X.tocoo()
    return (torch.from_numpy(np.vstack([c.row, c.col]).astype(np.int64)),
            torch.from_numpy(c.data.astype(np.float32)), X.shape)


def make_encoder(d, init, seed):
    torch.manual_seed(seed)
    enc = Encoder(d["Xs"].shape[1])
    if init == "pretrained":
        with torch.no_grad():
            enc.w1.copy_(d["w1"])
            enc.b1.copy_(d["b1"])
            enc.body[2].weight.copy_(d["w2"])
            enc.body[2].bias.copy_(d["b2"])
    return enc


def frozen_rep(d, v, seed, train_strains):
    if v["rep"] == "profiles":
        H = d["prof"].copy()
    elif v["init"] == "pretrained":
        H = d["emb"].copy()
    else:
        enc = make_encoder(d, "random", seed)
        enc.eval()
        with torch.no_grad():
            H = enc(SparseRows(d["Xs"]).get()).numpy()
    mu = H[train_strains].mean(0)
    sd = H[train_strains].std(0)
    # Floor the scale: many ReLU units are dead (sd ~ 0) on the training strains but fire on a
    # test strain; dividing by ~1e-6 blew R^2 to -1e7 in the smoke test (same failure as the
    # 16S PCA whitening in 20261113).
    sd = np.maximum(sd, 0.1 * sd.mean() + 1e-6)
    return torch.from_numpy(((H - mu) / sd).astype(np.float32))


def run_variant(d, v, tr, val, te, seed, n_steps=None):
    """Train on pairs tr. If n_steps is None: early-stop on val, return best step.
    Else: train exactly n_steps, return predictions on te."""
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    ia, ib, y = d["ia"], d["ib"], d["y"].copy()
    if v.get("shuffle"):
        y[tr] = y[rng.permutation(tr)]
    tr_strains = np.unique(np.r_[ia[tr], ib[tr]])
    if v["train_enc"]:
        enc = make_encoder(d, v["init"], seed)
        H_fix = None
        dim = 256
    else:
        enc = None
        H_fix = frozen_rep(d, v, seed, tr_strains)
        dim = H_fix.shape[1]
    torch.manual_seed(seed + 1)
    head = Head(dim)
    groups = [{"params": head.parameters(), "lr": HEAD_LR}]
    if enc is not None:
        groups.append({"params": enc.parameters(), "lr": v["enc_lr"]})
    opt = torch.optim.AdamW(groups, weight_decay=WD)
    loc = np.full(len(d["strains"]), -1)
    loc[tr_strains] = np.arange(len(tr_strains))
    la, lb = torch.from_numpy(loc[ia[tr]]), torch.from_numpy(loc[ib[tr]])
    yt = torch.from_numpy(y[tr])
    gen = torch.Generator().manual_seed(seed)
    if enc is not None:
        Xtr_t = SparseRows(d["Xs"][tr_strains])
        Xall_t = SparseRows(d["Xs"]).get()

    def all_embeddings():
        if enc is None:
            return H_fix
        enc.eval()
        with torch.no_grad():
            return enc(Xall_t)

    def predict(idx):
        head.eval()
        H = all_embeddings()
        with torch.no_grad():
            return head(H[ia[idx]], H[ib[idx]]).numpy()

    best, best_step, bad = np.inf, None, 0
    steps = n_steps if n_steps is not None else MAX_STEPS
    for step in range(1, steps + 1):
        head.train()
        if enc is not None:
            enc.train()
            Ht = enc(Xtr_t.get(P_IN, gen))
        else:
            Ht = H_fix[torch.from_numpy(tr_strains)]
        loss = ((head(Ht[la], Ht[lb]) - yt) ** 2).mean()
        opt.zero_grad()
        loss.backward()
        opt.step()
        if n_steps is None and step % EVAL_EVERY == 0:
            vl = float(np.mean((predict(val) - d["y"][val]) ** 2))
            if vl < best - 1e-6:
                best, best_step, bad = vl, step, 0
            else:
                bad += EVAL_EVERY
                if bad >= PATIENCE and step >= MIN_STEPS:
                    break
    if n_steps is None:
        return max(best_step or FALLBACK_STEPS, MIN_STEPS)
    return predict(te)


def inner_split(d, level, tr, seed):
    """Validation = training pairs touching 20% of the training groups; inner-train = the
    training pairs touching none of them."""
    g = d["groups"][level]
    pr = d["pairs"].iloc[tr]
    ga, gb = pr.strain_a.map(g).to_numpy(), pr.strain_b.map(g).to_numpy()
    ug = np.unique(np.r_[ga, gb])
    for attempt in range(10):
        rng = np.random.default_rng(seed * 100 + attempt)
        vg = set(rng.choice(ug, size=max(1, int(round(VAL_FRAC * len(ug)))), replace=False))
        in_v = np.array([a in vg or b in vg for a, b in zip(ga, gb)])
        if in_v.sum() >= 5 and (~in_v).sum() >= 10:
            return tr[~in_v], tr[in_v]
    return None, None


def _init_worker():
    """Forked workers came up pinned to ONE physical core (CPUs 0,32) although the parent
    reports 0-63 -- 20 workers shared one core (first full run, 2026-11-21). Reset the mask."""
    os.sched_setaffinity(0, range(os.cpu_count()))
    torch.set_num_threads(1)


def run_job(job):
    torch.set_num_threads(1)
    d = _D[job["exp"]]
    tr, te = np.asarray(job["train"]), np.asarray(job["test"])
    seed = job["seed"]
    itr, val = inner_split(d, job["level"], tr, seed)
    out = []
    for name, v in VARIANTS.items():
        t0 = time.time()
        if itr is not None:
            n = run_variant(d, v, itr, val, None, seed)
        else:
            n = FALLBACK_STEPS
        yhat = run_variant(d, v, tr, None, te, seed, n_steps=n)
        out.append(pd.DataFrame({"variant": name, "pair_idx": te, "y_true": d["y"][te],
                                 "y_pred": yhat, "steps": n, "n_train": len(tr),
                                 "secs": round(time.time() - t0, 1)}))
    res = pd.concat(out, ignore_index=True)
    for k in ["exp", "level", "rep", "fold"]:
        res[k] = job[k]
    return res


# ------------------------------------------------------------------ metrics
def metrics(y, yhat, sign_threshold=1.0):
    """= genomic_ml._metrics."""
    y, yhat = np.asarray(y, float), np.asarray(yhat, float)
    dec = np.abs(y) > sign_threshold
    out = {"n": len(y), "r2": float(1 - np.sum((y - yhat) ** 2) / np.sum(y ** 2))}
    if len(y) < 3 or np.std(yhat) == 0:
        return {**out, "spearman_rho": np.nan, "pearson_r": np.nan, "sign_accuracy": np.nan}
    return {**out, "spearman_rho": float(spearmanr(y, yhat)[0]),
            "pearson_r": float(pearsonr(y, yhat)[0]),
            "sign_accuracy": float(np.mean(np.sign(y[dec]) == np.sign(yhat[dec]))) if dec.sum() else np.nan}


def summarise(pred: pd.DataFrame, tag: str) -> pd.DataFrame:
    by = ["exp", "level", "variant", "rep", "fold"]
    fm = pd.DataFrame([{**dict(zip(by, k)), **metrics(g.y_true, g.y_pred)}
                       for k, g in pred.groupby(by, sort=False)])
    fm.to_csv(OUT / f"s18_{tag}_fold_metrics.csv", index=False)
    num = ["spearman_rho", "r2", "sign_accuracy"]
    per_rep = fm.groupby(["exp", "level", "variant", "rep"])[num].mean().reset_index()
    rows = []
    for (exp, level, var), g in per_rep.groupby(["exp", "level", "variant"]):
        f = fm[(fm.exp == exp) & (fm.level == level) & (fm.variant == var)]
        r = {"exp": exp, "level": level, "variant": var, "n_reps": g.rep.nunique(),
             "n_folds": len(f), "median_test_pairs_per_fold": float(f.n.median())}
        for c in num:
            r[f"{c}"] = g[c].mean()
            r[f"{c}_rep_sd"] = g[c].std()
        r["r2_median_fold"] = f.r2.median()
        # paired vs naive (same encoder lr) and vs best naive, per repeat
        for ref in ["naive_lr1e-3", "naive_lr1e-4"]:
            gr = per_rep[(per_rep.exp == exp) & (per_rep.level == level) & (per_rep.variant == ref)]
            m = g.merge(gr, on="rep", suffixes=("", "_ref"))
            dlt = m.spearman_rho - m.spearman_rho_ref
            r[f"d_rho_vs_{ref}"] = dlt.mean()
            r[f"d_rho_vs_{ref}_sd"] = dlt.std()
            r[f"p_rho_vs_{ref}"] = (wilcoxon(dlt).pvalue if len(dlt) >= 5 and dlt.abs().sum() > 0
                                    else np.nan)
        rows.append(r)
    s = pd.DataFrame(rows)
    order = {v: i for i, v in enumerate(VARIANTS)}
    s["_o"] = s.variant.map(order)
    s["_l"] = s.level.map({l: i for i, l in enumerate(LEVELS)})
    s = s.sort_values(["exp", "_l", "_o"]).drop(columns=["_o", "_l"])
    s.to_csv(OUT / f"s18_{tag}_summary.csv", index=False)
    return s


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--exps", nargs="+", default=list(EXPS))
    ap.add_argument("--levels", nargs="+", default=LEVELS)
    ap.add_argument("--processes", type=int, default=56)
    ap.add_argument("--quick", action="store_true", help="smoke test: 1 rep, 2 folds/level")
    a = ap.parse_args()
    global MAX_STEPS
    if a.quick:
        MAX_STEPS = 60
    t0 = time.time()
    prepare(a.exps)
    jobs = []
    for exp in a.exps:
        d = _D[exp]
        for level in a.levels:
            for rep in range(1 if a.quick else a.reps):
                for i, (k, tr, te) in enumerate(group_folds(d["pairs"], d["groups"][level], seed=rep)):
                    if a.quick and i >= 2:
                        break
                    jobs.append(dict(exp=exp, level=level, rep=rep, fold=k, train=tr, test=te,
                                     seed=zlib.crc32(f"{exp}|{level}|{rep}|{k}".encode()) % 2**31))
    print(f"{len(jobs)} jobs x {len(VARIANTS)} variants [{time.time()-t0:.0f}s]", flush=True)
    # big (slow) jobs first for better packing
    jobs.sort(key=lambda j: -len(j["train"]))
    parts = []
    with get_context("fork").Pool(a.processes, initializer=_init_worker) as pool:
        for i, part in enumerate(pool.imap_unordered(run_job, jobs, chunksize=1), 1):
            parts.append(part)
            if i % 20 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)} jobs [{(time.time()-t0)/60:.1f} min]", flush=True)
    pred = pd.concat(parts, ignore_index=True)
    tag = "quick" if a.quick else f"r{a.reps}"
    pred.to_csv(OUT / f"s18_{tag}_predictions.csv.gz", index=False)
    s = summarise(pred, tag)
    cols = ["exp", "level", "variant", "n_folds", "spearman_rho", "spearman_rho_rep_sd", "r2",
            "r2_median_fold", "sign_accuracy", "d_rho_vs_naive_lr1e-3", "p_rho_vs_naive_lr1e-3"]
    with pd.option_context("display.width", 250, "display.max_rows", 200):
        print(s[cols].round(3).to_string(index=False))
    print(f"done [{(time.time()-t0)/60:.1f} min]")


if __name__ == "__main__":
    main()
