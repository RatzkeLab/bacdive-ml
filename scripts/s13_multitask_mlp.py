"""Multi-task MLP: protein-family presence -> all 351 BacDive labels at once.

This is the candidate *pretrained* model: its hidden layer is the genome
embedding we would transfer to the interaction model.

  input (binary families, prevalence >= 10)  --feature dropout-->
  Linear(F, 512) -> ReLU -> Dropout -> Linear(512, 256) -> ReLU -> Dropout -> Linear(256, L)
  loss = BCE over observed labels only (NaN = masked), pos_weight per label

Evaluation: SAME folds and SAME metric as s12 (within-fold AUROC, imported).
Early stopping on a genus-held-out validation split (15% of training genera)
inside each training fold; the test fold is never looked at for model selection.

Usage: python s13_multitask_mlp.py [--tag id50_c80] [--split genus] [--final]
  --final  additionally trains on ALL strains (epochs = median best epoch of the
           CV folds) and saves the model for transfer (s14).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
import torch.nn as nn
from scipy.stats import wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parent))
import s12_baseline_models as s12  # noqa: E402
from s12_baseline_models import fold_auroc, load, make_folds  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
MODELS = ROOT / "data" / "models"
SEED = 20261120
H1, H2 = 512, 256
P_IN, P_HID = 0.2, 0.5
LR, WD = 1e-3, 1e-5
BATCH = 128
MAX_EPOCHS, PATIENCE = 120, 12


class MTL(nn.Module):
    def __init__(self, n_in: int, n_out: int):
        super().__init__()
        self.w1 = nn.Parameter(torch.empty(n_in, H1))
        nn.init.kaiming_uniform_(self.w1, a=5 ** 0.5)
        self.b1 = nn.Parameter(torch.zeros(H1))
        self.body = nn.Sequential(nn.ReLU(), nn.Dropout(P_HID), nn.Linear(H1, H2), nn.ReLU(),
                                  nn.Dropout(P_HID))
        self.head = nn.Linear(H2, n_out)

    def embed(self, x_sparse: torch.Tensor) -> torch.Tensor:
        return self.body(torch.sparse.mm(x_sparse, self.w1) + self.b1)

    def forward(self, x_sparse: torch.Tensor) -> torch.Tensor:
        return self.head(self.embed(x_sparse))


def to_torch(X: sp.csr_matrix, drop: float = 0.0, rng=None) -> torch.Tensor:
    X = X.tocoo()
    data = X.data.astype(np.float32)
    if drop > 0:
        keep = rng.random(len(data)) >= drop
        data = data * keep / (1 - drop)
    idx = torch.from_numpy(np.vstack([X.row, X.col]).astype(np.int64))
    return torch.sparse_coo_tensor(idx, torch.from_numpy(data), X.shape).coalesce()


def masked_bce(logits, y, pos_weight):
    m = ~torch.isnan(y)
    yy = torch.nan_to_num(y)
    loss = nn.functional.binary_cross_entropy_with_logits(
        logits, yy, pos_weight=pos_weight, reduction="none")
    return (loss * m).sum() / m.sum().clamp(min=1)


def predict(model, X) -> np.ndarray:
    model.eval()
    out = []
    with torch.no_grad():
        for i in range(0, X.shape[0], 1024):
            out.append(torch.sigmoid(model(to_torch(X[i:i + 1024]))).numpy())
    return np.vstack(out)


def train(X, Yv, tr, va, n_epochs: int | None = None, seed: int = SEED):
    """Train on rows tr; early-stop on rows va (if n_epochs is None)."""
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    model = MTL(X.shape[1], Yv.shape[1])
    opt = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WD)
    ytr = Yv[tr]
    npos = np.nansum(ytr, axis=0)
    nneg = np.sum(~np.isnan(ytr), axis=0) - npos
    pw = torch.tensor(np.clip(nneg / np.maximum(npos, 1), 0.1, 10), dtype=torch.float32)
    best, best_ep, best_state, bad = np.inf, 0, None, 0
    Xva = to_torch(X[va]) if va is not None else None
    yva = torch.tensor(Yv[va], dtype=torch.float32) if va is not None else None
    epochs = n_epochs or MAX_EPOCHS
    for ep in range(1, epochs + 1):
        model.train()
        perm = rng.permutation(tr)
        for i in range(0, len(perm), BATCH):
            b = perm[i:i + BATCH]
            xb = to_torch(X[b], P_IN, rng)
            yb = torch.tensor(Yv[b], dtype=torch.float32)
            loss = masked_bce(model(xb), yb, pw)
            opt.zero_grad()
            loss.backward()
            opt.step()
        if va is None:
            continue
        model.eval()
        with torch.no_grad():
            vl = masked_bce(model(Xva), yva, pw).item()
        if vl < best - 1e-4:
            best, best_ep, bad = vl, ep, 0
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
        else:
            bad += 1
            if bad >= PATIENCE:
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    return model, (best_ep if va is not None else epochs)


def genus_val_split(tax, tr, frac=0.15, seed=SEED):
    rng = np.random.default_rng(seed)
    g = tax.iloc[tr]["genus"].fillna("_na").to_numpy()
    ug = np.unique(g)
    val_g = set(rng.choice(ug, size=max(1, int(frac * len(ug))), replace=False))
    is_va = np.array([x in val_g for x in g])
    return tr[~is_va], tr[is_va]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="id50_c80")
    ap.add_argument("--split", default="genus", choices=["genus", "family", "strain"])
    ap.add_argument("--threads", type=int, default=48)
    ap.add_argument("--final", action="store_true")
    ap.add_argument("--max-epochs", type=int, default=None, help="override (testing)")
    ap.add_argument("--min-prev", type=int, default=87,
                    help="family must occur in >= this many BacDive genomes (87 ~ 1%% of "
                         "8,740). Rarer families are mostly genus-private and carry no "
                         "signal under genus-held-out CV, and they would make the first "
                         "layer ~8x larger.")
    a = ap.parse_args()
    torch.set_num_threads(a.threads)
    global MAX_EPOCHS
    if a.max_epochs:
        MAX_EPOCHS = a.max_epochs
    t0 = time.time()
    s12.MIN_PREV = a.min_prev  # load()/load_presence() read the module global
    X, Y, tax, n_feat = load(a.tag)
    Yv = Y.to_numpy(dtype=float)
    folds = make_folds(tax, a.split)
    print(f"X {X.shape}, labels {Y.shape}", flush=True)

    oof = np.full(Yv.shape, np.nan)
    best_eps = []
    for f in range(5):
        tr_all = np.where(folds != f)[0]
        te = np.where(folds == f)[0]
        tr, va = genus_val_split(tax, tr_all, seed=SEED + f)
        model, ep = train(X, Yv, tr, va, seed=SEED + f)
        best_eps.append(ep)
        oof[te] = predict(model, X[te])
        print(f"fold {f}: best epoch {ep}, n_train {len(tr)}, n_val {len(va)} "
              f"[{(time.time()-t0)/60:.1f} min]", flush=True)

    rows = []
    for j, col in enumerate(Y.columns):
        m = ~np.isnan(Yv[:, j])
        rows.append({"label": col, "mlp": fold_auroc(Yv[m, j], oof[m, j], folds[m])})
    res = pd.DataFrame(rows)
    tag = f"s13_{a.tag}_{a.split}"
    res.to_csv(OUT / f"{tag}_auroc_per_label.csv", index=False)
    long = pd.DataFrame(oof, index=Y.index, columns=Y.columns).stack().rename("mlp")
    long.reset_index().rename(columns={"level_1": "label"}).to_csv(
        OUT / f"{tag}_oof.csv.gz", index=False)

    base = pd.read_csv(OUT / f"s12_{a.tag}_{a.split}_auroc_per_label.csv")
    cmp_ = base.merge(res, on="label").dropna(subset=["tax", "logreg", "mlp"])
    big = cmp_[cmp_.n - cmp_.n_pos >= 100]
    big = big[big.n_pos >= 100]
    summ = {
        "best_epochs": best_eps, "n_features": n_feat, "min_prev": a.min_prev,
        "median_auroc": cmp_[["tax", "knn", "logreg", "mlp"]].median().round(3).to_dict(),
        "median_auroc_labels_ge100_each_class": big[["tax", "knn", "logreg", "mlp"]]
                                                   .median().round(3).to_dict(),
        "n_labels_ge100_each_class": len(big),
        "mlp_beats_logreg_frac": round(float((cmp_.mlp > cmp_.logreg).mean()), 3),
        "mlp_beats_tax_frac": round(float((cmp_.mlp > cmp_.tax).mean()), 3),
        "wilcoxon_p_mlp_vs_logreg": float(wilcoxon(cmp_.mlp, cmp_.logreg).pvalue),
        "median_by_kind_group": cmp_.groupby("kind_group")[["tax", "knn", "logreg", "mlp"]]
                                    .median().round(3).to_dict(orient="index"),
        "minutes": round((time.time() - t0) / 60, 1),
    }
    (OUT / f"{tag}_summary.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1), flush=True)

    if a.final:
        n_ep = int(np.median(best_eps))
        model, _ = train(X, Yv, np.arange(len(Yv)), None, n_epochs=n_ep, seed=SEED)
        MODELS.mkdir(parents=True, exist_ok=True)
        torch.save({"state_dict": model.state_dict(), "labels": list(Y.columns),
                    "n_in": X.shape[1], "tag": a.tag, "epochs": n_ep,
                    "min_prev": a.min_prev},
                   MODELS / f"mlp_{a.tag}_all.pt")
        print(f"saved final model ({n_ep} epochs)", flush=True)


if __name__ == "__main__":
    main()
