"""How much of each BacDive label is predictable from taxonomy alone?

This is the bar any genome->phenotype model has to clear (our interaction ML
found KO ~ 16S once leaks were removed, so this matters).  For every usable
label (s03 label matrix: genome strains, >=25 pos and >=25 neg) score a
nearest-taxon frequency predictor with AUROC, in two regimes:

  loo_strain   predict strain i from the other strains of its genus
               (fallback family -> order -> global), i.e. a random split.
  loo_genus    hide i's whole genus, predict from the rest of its family
               (fallback order -> class -> global), i.e. a genus-held-out split,
               the relevant regime for strains outside the training set.

No genomes involved; uses LPSN taxonomy from the strain table.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "data" / "tables"
OUT = ROOT / "outputs"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
C1, C2 = "#2a78d6", "#eb6834"


def loo_rate(df: pd.DataFrame, levels: list[str], exclude: str | None) -> np.ndarray:
    """Hierarchical leave-out frequency of label==1.
    exclude=None   -> leave out the strain itself at every level
    exclude='genus'-> leave out the strain's whole genus at every level
    Uses the first level with >=1 remaining labelled strain, else global."""
    y = df.label.values
    pred = np.full(len(df), np.nan)
    if exclude:
        g_sum = df.groupby(exclude).label.transform("sum").values
        g_n = df.groupby(exclude).label.transform("size").values
    else:
        g_sum, g_n = y, np.ones(len(df))
    for lvl in levels:
        s = df.groupby(lvl).label.transform("sum").values - g_sum
        n = df.groupby(lvl).label.transform("size").values - g_n
        ok = np.isnan(pred) & (n > 0) & df[lvl].notna().values
        pred[ok] = s[ok] / n[ok]
    tot_s, tot_n = y.sum() - g_sum, len(y) - g_n
    rest = np.isnan(pred)
    pred[rest] = tot_s[rest] / tot_n[rest]
    return pred


def positive_controls(st: pd.DataFrame) -> dict:
    """Same predictor on traits known to be phylogenetically conserved. If these
    score ~0.95+, low AUROCs on metabolic labels are real, not a predictor bug."""
    t = pd.read_csv(TAB / "traits.csv.gz", low_memory=False)
    gid = set(pd.read_csv(TAB / "genome_choice.csv.gz").bacdive_id)
    out = {}
    for trait, pos, neg in [
            ("cell morphology:gram stain", ["positive"], ["negative"]),
            ("oxygen tolerance", ["anaerobe", "obligate anaerobe"], ["aerobe", "obligate aerobe"]),
            ("spore formation", ["yes"], ["no"]),
            ("cell morphology:motility", ["yes"], ["no"])]:
        x = t[(t.trait == trait) & t.bacdive_id.isin(gid) & t.value.isin(pos + neg)]
        y = x.groupby("bacdive_id").value.agg(
            lambda s: 1.0 if s.isin(pos).all() else (0.0 if s.isin(neg).all() else np.nan)).dropna()
        df = st.loc[y.index].assign(label=y.values)
        df = df[df.genus.notna()]
        out[trait] = {
            "n": len(df),
            "auroc_loo_strain": round(roc_auc_score(df.label, loo_rate(df, ["genus", "family", "order"], None)), 3),
            "auroc_loo_genus": round(roc_auc_score(df.label, loo_rate(df, ["family", "order", "class"], "genus")), 3),
        }
    return out


def main() -> None:
    mat = pd.read_csv(TAB / "label_matrix_mu_grouped_minority25_genome_strains.csv.gz",
                      index_col=0)
    st = pd.read_csv(TAB / "strains.csv.gz", low_memory=False,
                     usecols=["bacdive_id", "genus", "family", "order", "class", "species"]
                     ).set_index("bacdive_id")
    rows = []
    for col in mat.columns:
        y = mat[col].dropna()
        df = st.loc[y.index].assign(label=y.values)
        df = df[df.genus.notna()]
        if df.label.nunique() < 2:
            continue
        p_strain = loo_rate(df, ["genus", "family", "order"], exclude=None)
        p_genus = loo_rate(df, ["family", "order", "class"], exclude="genus")
        n_gen = df.groupby("genus").label.transform("size")
        rows.append({
            "label": col, "n": len(df), "pos_frac": df.label.mean(),
            "auroc_loo_strain": roc_auc_score(df.label, p_strain),
            "auroc_loo_genus": roc_auc_score(df.label, p_genus),
            "frac_with_genus_mate": float((n_gen > 1).mean()),
            "n_genera": df.genus.nunique(),
        })
    res = pd.DataFrame(rows).sort_values("auroc_loo_genus")
    res.to_csv(OUT / "s05_taxonomy_baseline.csv", index=False)
    res["kind_group"] = res.label.str.split("|").str[1]
    summ = {
        "n_labels": len(res),
        "auroc_loo_strain_quantiles": res.auroc_loo_strain.quantile([.1, .25, .5, .75, .9]).round(3).to_dict(),
        "auroc_loo_genus_quantiles": res.auroc_loo_genus.quantile([.1, .25, .5, .75, .9]).round(3).to_dict(),
        "median_by_kind_group": res.groupby("kind_group")[["auroc_loo_strain", "auroc_loo_genus"]]
                                   .median().round(3).to_dict(orient="index"),
        "labels_genus_out_auroc_lt_0.75": int((res.auroc_loo_genus < 0.75).sum()),
        "positive_controls_conserved_traits": positive_controls(st),
    }
    (OUT / "s05_taxonomy_baseline_summary.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1))

    fig, ax = plt.subplots(figsize=(6.4, 3.8), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    for s in ["left", "bottom"]:
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for c, colr, lab in [("auroc_loo_strain", C1, "strain held out (genus-mates visible)"),
                         ("auroc_loo_genus", C2, "whole genus held out")]:
        v = np.sort(res[c].values)
        ax.plot(v, np.arange(1, len(v) + 1) / len(v), color=colr, lw=2, label=lab)
    ax.set_xlabel("AUROC of taxonomy-only predictor", color=INK2)
    ax.set_ylabel("fraction of labels", color=INK2)
    ax.set_title(f"Taxonomy-only predictor per metabolic label (n={len(res)} labels)",
                 color=INK, fontsize=10.5, loc="left")
    ax.legend(frameon=False, fontsize=8, labelcolor=INK2, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "figures" / "s05_taxonomy_baseline_auroc.png", dpi=160)


if __name__ == "__main__":
    main()
