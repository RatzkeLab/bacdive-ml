"""Follow-up to s18/s19: do frozen LOGISTIC-REGRESSION BacDive profiles help on unseen taxa?

s18's frozen-profile arm used the multi-task MLP's 351 predicted phenotypes. Logistic regression
was slightly better on BacDive (s12/s13), so this swaps in its predictions
(outputs/s14_lab_profiles_logreg.csv: fitted on BacDive genomes only, predicted for lab genomes)
and changes nothing else:
  same pairs, wrong-genome filter, folds, seeds, head, early stopping and z-scoring as
  s18 `pretrained_profiles_frozen`; probabilities are used as is (no transform was tried).

Decision rule, fixed before running: logreg profiles count as a gain only if the combined
(both experiments) logreg - naive 95% CI excludes 0 at genus or family level after Holm
correction over the 4 levels.

Reproduction checks (the script stops if either fails):
  1. re-running the MLP-profile arm here reproduces s18's saved predictions;
  2. the s19 contrasts recomputed here, with the same seeds, reproduce s19's bootstrap table.

Usage: python s20_logreg_profiles.py [--reps 10] [--boots 2000] [--processes 32]
Outputs: outputs/s20_r10_{predictions.csv.gz,rho_by_level.csv,bootstrap.csv}
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import sys
import time
import zlib
from multiprocessing import get_context
from pathlib import Path

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import s18_interaction_pretraining as s18  # noqa: E402
import s19_pretraining_bootstrap as s19  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
NEW = "logreg_profiles_frozen"
CHECK = "pretrained_profiles_frozen"
KEYS = ["exp", "level", "rep", "fold", "pair_idx"]
EXPS = s19.EXPS
LEVELS = s18.LEVELS

_base_contrasts = s19.contrasts


def contrasts(r, best_naive):
    out = _base_contrasts(r, best_naive)
    i = {v: k for k, v in enumerate(s19.VARS)}
    naive = (r[i["naive_lr1e-3"]] + r[i["naive_lr1e-4"]]) / 2
    out["logreg_vs_naive"] = r[i[NEW]] - naive
    out["logreg_vs_bestnaive"] = r[i[NEW]] - r[i[best_naive]]
    out["logreg_vs_mlp_profiles"] = r[i[NEW]] - r[i[CHECK]]
    return out


def attach_logreg() -> None:
    ids = pd.read_csv(OUT / "s14_lab_genome_ids.csv").set_index("well")
    L = pd.read_csv(OUT / "s14_lab_profiles_logreg.csv", index_col=0)
    M = pd.read_csv(OUT / "s14_lab_profiles_mlp.csv", index_col=0)
    assert list(L.columns) == list(M.columns)
    for exp, d in s18._D.items():
        sn = ids.loc[d["strains"], "strain_no"].to_numpy()
        # join check: the same well -> strain_no map must recover s18's own MLP profiles
        err = np.abs(M.loc[sn].to_numpy(np.float32) - d["prof"]).max()
        assert err < 1e-3, (exp, err)
        d["prof_logreg"] = L.loc[sn].to_numpy(np.float32)
        print(f"{exp}: logreg profiles {d['prof_logreg'].shape}, join check max err {err:.1e}",
              flush=True)


def run_job(job):
    torch.set_num_threads(1)
    d = s18._D[job["exp"]]
    tr, te = np.asarray(job["train"]), np.asarray(job["test"])
    seed = job["seed"]
    itr, val = s18.inner_split(d, job["level"], tr, seed)
    v = s18.VARIANTS[CHECK]
    out = []
    for name, dd in ((CHECK, d), (NEW, {**d, "prof": d["prof_logreg"]})):
        n = (s18.run_variant(dd, v, itr, val, None, seed) if itr is not None
             else s18.FALLBACK_STEPS)
        yhat = s18.run_variant(dd, v, tr, None, te, seed, n_steps=n)
        out.append(pd.DataFrame({"variant": name, "pair_idx": te, "y_pred": yhat, "steps": n}))
    res = pd.concat(out, ignore_index=True)
    for k in ["exp", "level", "rep", "fold"]:
        res[k] = job[k]
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--boots", type=int, default=2000)
    ap.add_argument("--processes", type=int, default=32)
    a = ap.parse_args()
    t0 = time.time()
    tag = f"r{a.reps}"
    s18.prepare(EXPS)
    attach_logreg()

    jobs = []
    for exp in EXPS:
        d = s18._D[exp]
        for level in LEVELS:
            for rep in range(a.reps):
                for k, tr, te in s18.group_folds(d["pairs"], d["groups"][level], seed=rep):
                    jobs.append(dict(exp=exp, level=level, rep=rep, fold=k, train=tr, test=te,
                                     seed=zlib.crc32(f"{exp}|{level}|{rep}|{k}".encode()) % 2**31))
    print(f"{len(jobs)} fold jobs [{time.time()-t0:.0f}s]", flush=True)
    parts = []
    with get_context("fork").Pool(a.processes, initializer=s18._init_worker) as pool:
        for i, part in enumerate(pool.imap_unordered(run_job, jobs, chunksize=1), 1):
            parts.append(part)
            if i % 40 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)} [{(time.time()-t0)/60:.1f} min]", flush=True)
    new = pd.concat(parts, ignore_index=True)
    new["exp"] = new["exp"].astype(str)

    # check 1: the MLP-profile arm reproduces s18
    old = pd.read_csv(OUT / f"s18_{tag}_predictions.csv.gz", dtype={"exp": str})
    ref = old[old.variant == CHECK]
    chk = new[new.variant == CHECK].merge(ref[KEYS + ["y_pred", "steps"]], on=KEYS,
                                          suffixes=("", "_s18"), validate="1:1")
    assert len(chk) == len(ref), (len(chk), len(ref))
    max_diff = float((chk.y_pred - chk.y_pred_s18).abs().max())
    print(f"check 1: {len(chk)} MLP-profile predictions vs s18, max |diff| {max_diff:.2e}, "
          f"step mismatches {(chk.steps != chk.steps_s18).sum()}", flush=True)
    assert max_diff < 1e-3 and (chk.steps == chk.steps_s18).all()

    lg = new[new.variant == NEW].merge(ref[KEYS + ["y_true"]], on=KEYS, validate="1:1")
    assert len(lg) == len(ref)
    cols = KEYS + ["variant", "y_true", "y_pred", "steps"]
    lg[cols].to_csv(OUT / f"s20_{tag}_predictions.csv.gz", index=False)
    pred = pd.concat([old[cols], lg[cols]], ignore_index=True)

    # bootstrap: s19's machinery and seeds, with the new arm added
    s19.VARS = [*s19.VARS, NEW]
    s19.contrasts = contrasts
    s19.build_cells(pred)
    C = s19._C
    i = {v: k for k, v in enumerate(s19.VARS)}
    rows, point, best = [], {}, {}
    for level in LEVELS:
        best[level], point[level] = {}, {}
        for exp in EXPS:
            r = s19.score(C[level]["exps"][exp])
            best[level][exp] = ("naive_lr1e-3" if r[i["naive_lr1e-3"]] >= r[i["naive_lr1e-4"]]
                                else "naive_lr1e-4")
            point[level][exp] = contrasts(r, best[level][exp])
            rows.append({"exp": exp, "level": level,
                         "naive": (r[i["naive_lr1e-3"]] + r[i["naive_lr1e-4"]]) / 2,
                         "mlp_profiles": r[i[CHECK]], "logreg_profiles": r[i[NEW]]})
        point[level]["both"] = {k: np.mean([point[level][e][k] for e in EXPS])
                                for k in point[level][EXPS[0]]}
    table = pd.DataFrame(rows)
    table.to_csv(OUT / f"s20_{tag}_rho_by_level.csv", index=False)

    seeds = np.random.SeedSequence(s19.SEED).generate_state(a.boots * len(LEVELS))
    bjobs = [(level, int(seeds[n * len(LEVELS) + j]), best[level])
             for n in range(a.boots) for j, level in enumerate(LEVELS)]
    with get_context("fork").Pool(a.processes, initializer=s18._init_worker) as pool:
        res = pool.map(s19.boot_one, bjobs, chunksize=20)

    out = []
    for level in LEVELS:
        draws = [o for l, o in res if l == level]
        for exp in [*EXPS, "both"]:
            for con, delta in point[level][exp].items():
                v = np.array([x[exp][con] for x in draws], float)
                v = v[np.isfinite(v)]
                p = min(1.0, 2 * min((np.sum(v <= 0) + 1) / (len(v) + 1),
                                     (np.sum(v >= 0) + 1) / (len(v) + 1)))
                out.append({"exp": exp, "level": level, "contrast": con, "delta": delta,
                            "ci_lo": np.percentile(v, 2.5), "ci_hi": np.percentile(v, 97.5),
                            "p_boot": p, "n_boot": len(v)})
    boot = pd.DataFrame(out)

    # check 2: contrasts shared with s19 reproduce its table
    s = pd.read_csv(OUT / f"s19_{tag}_bootstrap.csv", dtype={"exp": str})
    m = boot.merge(s, on=["exp", "level", "contrast"], suffixes=("", "_s19"), validate="1:1")
    assert len(m) == len(s), (len(m), len(s))
    worst = max(float((m[c] - m[f"{c}_s19"]).abs().max()) for c in ["delta", "ci_lo", "ci_hi"])
    print(f"check 2: {len(m)} s19 contrasts recomputed, max |diff| {worst:.2e}", flush=True)
    assert worst < 1e-6

    prim = (boot.exp == "both") & (boot.contrast == "logreg_vs_naive")
    boot["p_holm_logreg"] = np.nan
    boot.loc[prim, "p_holm_logreg"] = s19.holm(boot.loc[prim, "p_boot"])
    boot.to_csv(OUT / f"s20_{tag}_bootstrap.csv", index=False)

    show = boot[boot.contrast.str.startswith("logreg")]
    with pd.option_context("display.width", 200, "display.max_rows", 100):
        print(table.round(3).to_string(index=False))
        print(show.round(3).to_string(index=False))
    print(f"done [{(time.time()-t0)/60:.1f} min]")


if __name__ == "__main__":
    main()
