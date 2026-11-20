"""Summary figure: median AUROC per model across CV regimes (strain / genus / family out).

Reads outputs/s12_<tag>_<split>_auroc_per_label.csv (+ s13 MLP for genus) and the s14
domain-shift summary, writes outputs/figures/s16_<tag>_auroc_by_split.png and
outputs/s16_<tag>_auroc_by_split.csv.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
COLORS = {"tax": "#2a78d6", "knn": "#eb6834", "logreg": "#1baf7a"}
NAMES = {"tax": "taxonomy only", "knn": "genome kNN (relatedness)",
         "logreg": "logistic regression on families"}
SPLITS = [("strain", "strain\nheld out"), ("genus", "genus\nheld out"),
          ("family", "family\nheld out")]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="id50_c80")
    a = ap.parse_args()
    rows = []
    for split, _ in SPLITS:
        r = pd.read_csv(OUT / f"s12_{a.tag}_{split}_auroc_per_label.csv")
        big = r[(r.n_pos >= 100) & (r.n - r.n_pos >= 100)]
        for m in COLORS:
            rows.append({"split": split, "model": m, "labels": "all 351", "median": r[m].median()})
            rows.append({"split": split, "model": m, "labels": ">=100 per class",
                         "median": big[m].median()})
    mlp = pd.read_csv(OUT / f"s13_{a.tag}_genus_auroc_per_label.csv").merge(
        pd.read_csv(OUT / f"s12_{a.tag}_genus_auroc_per_label.csv")[["label", "n", "n_pos"]])
    big = mlp[(mlp.n_pos >= 100) & (mlp.n - mlp.n_pos >= 100)]
    rows += [{"split": "genus", "model": "mlp", "labels": "all 351", "median": mlp.mlp.median()},
             {"split": "genus", "model": "mlp", "labels": ">=100 per class",
              "median": big.mlp.median()}]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / f"s16_{a.tag}_auroc_by_split.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2), facecolor=SURFACE, sharey=True)
    for ax, lab in zip(axes, ["all 351", ">=100 per class"]):
        ax.set_facecolor(SURFACE)
        for s in ["top", "right"]:
            ax.spines[s].set_visible(False)
        for s in ["left", "bottom"]:
            ax.spines[s].set_color(GRID)
        ax.tick_params(colors=INK2, labelsize=9)
        ax.grid(axis="y", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        x = range(len(SPLITS))
        for m, c in COLORS.items():
            y = [df[(df.split == s) & (df.model == m) & (df.labels == lab)]["median"].iloc[0]
                 for s, _ in SPLITS]
            ax.plot(x, y, color=c, lw=2, marker="o", ms=8, mec=SURFACE, mew=2, label=NAMES[m])
            ax.text(len(SPLITS) - 1 + 0.08, y[-1], f"{y[-1]:.2f}", color=INK2, fontsize=8,
                    va="center")
            ax.text(-0.08, y[0], f"{y[0]:.2f}", color=INK2, fontsize=8, va="center", ha="right")
        ax.set_xticks(list(x), [t for _, t in SPLITS])
        ax.set_xlim(-0.45, len(SPLITS) - 0.55)
        ax.set_title(f"labels: {lab}", color=INK, fontsize=10, loc="left")
    axes[0].set_ylabel("median AUROC across labels (within-fold)", color=INK2)
    axes[0].legend(frameon=False, fontsize=8, labelcolor=INK2, loc="lower left")
    fig.suptitle("Genome -> BacDive phenotype: harder splits cost accuracy, gene content "
                 "stays ahead of taxonomy", color=INK, fontsize=10.5, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(OUT / "figures" / f"s16_{a.tag}_auroc_by_split.png", dpi=160)
    print(df.pivot_table(index=["labels", "model"], columns="split", values="median")
          .round(3).to_string())


if __name__ == "__main__":
    main()
