"""Honest uncertainty for s18: does pretraining beat the naive MLP?

s18's error bars (SD across CV repeats) and Wilcoxon tests treat repeats as independent, but
all repeats reuse the same strains, and at family level there are only 13-25 taxa. Here the
resampling unit is the TAXON at the level being held out (dyadic cluster bootstrap):

  draw taxa with replacement -> a pair counts c[taxon_a] * c[taxon_b] times -> recompute the
  s18 metric (Spearman rho within each fold, mean over folds, mean over repeats) for every
  variant on that same resample -> paired differences.

No model is refit: this uses s18's saved out-of-fold predictions, so it captures "which taxa
happened to be in the collection", not training randomness.

Contrasts (positive = pretraining better). "naive" and "pretrained_ft" are each the mean of
their two encoder learning rates, so no learning rate is picked after seeing test scores:
  ft_vs_naive          pretrained fine-tuned - naive            (the main question)
  profiles_vs_naive    frozen 351 BacDive phenotypes - naive
  profiles_vs_bestnaive  same, against whichever naive lr scored higher in that cell
                         (deliberately favours naive)
  frozen_vs_naive      frozen 256-d BacDive embedding - naive
  random_vs_naive      frozen RANDOM encoder - naive (control)

exp = "both" averages the two experiments' contrasts under one joint draw of taxa (they share
strains, so they are not independent). Holm correction is applied to the 8 primary tests
(ft_vs_naive and profiles_vs_naive x 4 levels, exp = both).

Usage: python s19_pretraining_bootstrap.py [--tag r10] [--boots 2000] [--processes 48]
Outputs: outputs/s19_<tag>_bootstrap.csv, s19_<tag>_rho_by_level.csv, s19_<tag>_summary.json,
         outputs/figures/s19_<tag>_pretraining.png
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import json
import sys
import warnings
from multiprocessing import get_context
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
import s18_interaction_pretraining as s18  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
LEVELS = s18.LEVELS
EXPS = ["20260907", "20260630"]
VARS = ["naive_lr1e-3", "naive_lr1e-4", "pretrained_ft_lr1e-3", "pretrained_ft_lr1e-4",
        "pretrained_profiles_frozen", "pretrained_frozen", "random_frozen"]
SEED = 20261123
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
SERIES = [("naive", "naive MLP", "#2a78d6"),
          ("pretrained_ft", "BacDive-pretrained, fine-tuned", "#eb6834"),
          ("profiles", "frozen BacDive phenotype profiles", "#1baf7a")]

_C: dict = {}
warnings.filterwarnings("ignore", message="Mean of empty slice")


def build_cells(pred: pd.DataFrame) -> None:
    for level in LEVELS:
        taxa = sorted({t for e in EXPS for t in s18._D[e]["groups"][level].values()})
        code = {t: i for i, t in enumerate(taxa)}
        cell = {"n_taxa": len(taxa), "exps": {}, "exp_taxa": {}}
        for exp in EXPS:
            d = s18._D[exp]
            g = d["groups"][level]
            cell["exp_taxa"][exp] = np.unique([code[x] for x in g.values()])
            ca = d["pairs"].strain_a.map(g).map(code).to_numpy()
            cb = d["pairs"].strain_b.map(g).map(code).to_numpy()
            p = pred[(pred.exp == exp) & (pred.level == level) & pred.variant.isin(VARS)]
            w = p.pivot_table(index=["rep", "fold", "pair_idx", "y_true"], columns="variant",
                              values="y_pred").reset_index()
            folds = []
            for (rep, fold), f in w.groupby(["rep", "fold"], sort=True):
                ix = f.pair_idx.to_numpy()
                folds.append(dict(rep=rep, a=ca[ix], b=cb[ix], y=f.y_true.to_numpy(),
                                  P=f[VARS].to_numpy()))
            cell["exps"][exp] = folds
        _C[level] = cell


def spearman_columns(y, P):
    ry = rankdata(y)
    rp = rankdata(P, axis=0)
    ry = ry - ry.mean()
    rp = rp - rp.mean(0)
    den = np.sqrt((ry ** 2).sum() * (rp ** 2).sum(0))
    with np.errstate(invalid="ignore", divide="ignore"):
        return (ry @ rp) / den


def score(folds, c=None):
    """s18 metric per variant: rho per fold, mean over folds per repeat, mean over repeats."""
    by_rep: dict = {}
    for f in folds:
        if c is None:
            y, P = f["y"], f["P"]
        else:
            k = c[f["a"]] * c[f["b"]]
            if k.sum() < 5:
                continue
            ix = np.repeat(np.arange(len(k)), k)
            y, P = f["y"][ix], f["P"][ix]
        if np.std(y) == 0:
            continue
        by_rep.setdefault(f["rep"], []).append(spearman_columns(y, P))
    if not by_rep:
        return np.full(len(VARS), np.nan)
    return np.nanmean([np.nanmean(v, axis=0) for v in by_rep.values()], axis=0)


def contrasts(r, best_naive):
    i = {v: k for k, v in enumerate(VARS)}
    naive = (r[i["naive_lr1e-3"]] + r[i["naive_lr1e-4"]]) / 2
    ft = (r[i["pretrained_ft_lr1e-3"]] + r[i["pretrained_ft_lr1e-4"]]) / 2
    return {"ft_vs_naive": ft - naive,
            "profiles_vs_naive": r[i["pretrained_profiles_frozen"]] - naive,
            "profiles_vs_bestnaive": r[i["pretrained_profiles_frozen"]] - r[i[best_naive]],
            "frozen_vs_naive": r[i["pretrained_frozen"]] - naive,
            "random_vs_naive": r[i["random_frozen"]] - naive}


def boot_one(job):
    level, seed, best = job
    rng = np.random.default_rng(seed)
    cell = _C[level]
    # combined estimate: one joint draw over the union of taxa (experiments share strains)
    c = rng.multinomial(cell["n_taxa"], np.ones(cell["n_taxa"]) / cell["n_taxa"])
    joint = {exp: contrasts(score(cell["exps"][exp], c), best[exp]) for exp in EXPS}
    out = {"both": {k: np.nanmean([joint[e][k] for e in EXPS]) for k in joint[EXPS[0]]}}
    # per-experiment estimate: resample only that experiment's own taxa
    for exp in EXPS:
        codes = cell["exp_taxa"][exp]
        ce = np.zeros(cell["n_taxa"], dtype=int)
        ce[codes] = rng.multinomial(len(codes), np.ones(len(codes)) / len(codes))
        out[exp] = contrasts(score(cell["exps"][exp], ce), best[exp])
    return level, out


def holm(p: pd.Series) -> pd.Series:
    order = p.sort_values().index
    m = len(p)
    adj, running = {}, 0.0
    for j, ix in enumerate(order):
        running = max(running, min(1.0, (m - j) * p[ix]))
        adj[ix] = running
    return pd.Series(adj).reindex(p.index)


def figure(table: pd.DataFrame, boot: pd.DataFrame, tag: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 6.8), facecolor=SURFACE, sharex=True,
                             gridspec_kw={"height_ratios": [1.35, 1]})
    x = np.arange(len(LEVELS))
    for col, exp in enumerate(EXPS):
        ax, bx = axes[0, col], axes[1, col]
        for a in (ax, bx):
            a.set_facecolor(SURFACE)
            for s in ["top", "right"]:
                a.spines[s].set_visible(False)
            for s in ["left", "bottom"]:
                a.spines[s].set_color(GRID)
            a.tick_params(colors=INK2, labelsize=9)
            a.grid(axis="y", color=GRID, lw=0.6)
            a.set_axisbelow(True)
        for key, name, c in SERIES:
            y = [table[(table.exp == exp) & (table.level == l)][key].iloc[0] for l in LEVELS]
            ax.plot(x, y, color=c, lw=2, marker="o", ms=8, mec=SURFACE, mew=2, label=name)
        ax.set_ylim(0.25, 0.85)
        ax.set_title(f"{exp}: prediction on held-out taxa", color=INK, fontsize=10.5, loc="left")
        bx.axhline(0, color=INK2, lw=1)
        for off, (con, c, lab) in zip([-0.09, 0.09], [("ft_vs_naive", "#eb6834", "fine-tuned"),
                                                       ("profiles_vs_naive", "#1baf7a", "profiles")]):
            b = boot[(boot.exp == exp) & (boot.contrast == con)].set_index("level").loc[LEVELS]
            bx.errorbar(x + off, b.delta, yerr=[b.delta - b.ci_lo, b.ci_hi - b.delta], fmt="o",
                        color=c, ms=7, mec=SURFACE, mew=1.5, elinewidth=2, capsize=0, label=lab)
        bx.set_xticks(x, [f"{l}\nheld out" for l in LEVELS])
        bx.set_title("difference from naive (95% CI, resampling taxa)", color=INK, fontsize=9.5,
                     loc="left")
    axes[0, 0].set_ylabel("Spearman ρ (mean over folds)", color=INK2)
    axes[1, 0].set_ylabel("Δρ vs naive", color=INK2)
    axes[0, 0].legend(frameon=False, fontsize=8, labelcolor=INK2, loc="lower left")
    for bx in axes[1]:
        bx.set_ylim(-0.3, 0.55)  # same scale, so the experiments' intervals are comparable
    fig.tight_layout()
    fig.savefig(OUT / "figures" / f"s19_{tag}_pretraining.png", dpi=160)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="r10")
    ap.add_argument("--boots", type=int, default=2000)
    ap.add_argument("--processes", type=int, default=48)
    a = ap.parse_args()
    pred = pd.read_csv(OUT / f"s18_{a.tag}_predictions.csv.gz", dtype={"exp": str})
    s18.prepare(EXPS)
    build_cells(pred)

    # point estimates, checked against the s18 summary
    summ = pd.read_csv(OUT / f"s18_{a.tag}_summary.csv")
    rows, point, best = [], {}, {}
    for level in LEVELS:
        best[level], point[level] = {}, {}
        for exp in EXPS:
            r = score(_C[level]["exps"][exp])
            for v, val in zip(VARS, r):
                ref = summ[(summ.exp.astype(str) == exp) & (summ.level == level)
                           & (summ.variant == v)].spearman_rho.iloc[0]
                assert abs(val - ref) < 1e-6, (exp, level, v, val, ref)
            i = {v: k for k, v in enumerate(VARS)}
            best[level][exp] = ("naive_lr1e-3" if r[i["naive_lr1e-3"]] >= r[i["naive_lr1e-4"]]
                                else "naive_lr1e-4")
            point[level][exp] = contrasts(r, best[level][exp])
            rows.append({"exp": exp, "level": level,
                         "naive": (r[i["naive_lr1e-3"]] + r[i["naive_lr1e-4"]]) / 2,
                         "pretrained_ft": (r[i["pretrained_ft_lr1e-3"]] + r[i["pretrained_ft_lr1e-4"]]) / 2,
                         "profiles": r[i["pretrained_profiles_frozen"]],
                         "frozen": r[i["pretrained_frozen"]], "random": r[i["random_frozen"]]})
        point[level]["both"] = {k: np.mean([point[level][e][k] for e in EXPS])
                                for k in point[level][EXPS[0]]}
    table = pd.DataFrame(rows)
    table.to_csv(OUT / f"s19_{a.tag}_rho_by_level.csv", index=False)

    seeds = np.random.SeedSequence(SEED).generate_state(a.boots * len(LEVELS))
    jobs = [(level, int(seeds[i * len(LEVELS) + j]), best[level])
            for i in range(a.boots) for j, level in enumerate(LEVELS)]
    with get_context("fork").Pool(a.processes, initializer=s18._init_worker) as pool:
        res = pool.map(boot_one, jobs, chunksize=20)

    out = []
    for level in LEVELS:
        draws = [o for l, o in res if l == level]
        for exp in [*EXPS, "both"]:
            for con, delta in point[level][exp].items():
                v = np.array([d[exp][con] for d in draws], float)
                v = v[np.isfinite(v)]
                p = min(1.0, 2 * min((np.sum(v <= 0) + 1) / (len(v) + 1),
                                     (np.sum(v >= 0) + 1) / (len(v) + 1)))
                out.append({"exp": exp, "level": level, "contrast": con, "delta": delta,
                            "ci_lo": np.percentile(v, 2.5), "ci_hi": np.percentile(v, 97.5),
                            "p_boot": p, "n_boot": len(v),
                            "n_taxa": (_C[level]["n_taxa"] if exp == "both"
                                       else len(_C[level]["exp_taxa"][exp]))})
    boot = pd.DataFrame(out)
    prim = (boot.exp == "both") & boot.contrast.isin(["ft_vs_naive", "profiles_vs_naive"])
    boot["p_holm_primary"] = np.nan
    boot.loc[prim, "p_holm_primary"] = holm(boot.loc[prim, "p_boot"])
    boot.to_csv(OUT / f"s19_{a.tag}_bootstrap.csv", index=False)
    figure(table, boot, a.tag)

    steps = pred.groupby("variant").steps.median().to_dict()
    (OUT / f"s19_{a.tag}_summary.json").write_text(json.dumps(
        {"n_boot": a.boots, "median_training_steps_by_variant": steps,
         "n_taxa_by_level": {l: _C[l]["n_taxa"] for l in LEVELS}}, indent=1))
    with pd.option_context("display.width", 220, "display.max_rows", 200):
        print(table.round(3).to_string(index=False))
        print(boot.round(3).to_string(index=False))
    print("median training steps:", steps)


if __name__ == "__main__":
    main()
