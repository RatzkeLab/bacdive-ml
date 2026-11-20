"""Apply the BacDive-trained models to our 298 lab genomes.

1. Domain shift: Jaccard (protein-family presence) of each lab genome to its nearest
   BacDive genome, vs. the nearest-neighbour similarity that BacDive test strains had
   in CV: genus held out (= what s12/s13 scores reflect) and strain held out. If lab
   genomes sit closer to BacDive than genus-held-out test strains did, the CV numbers
   are a conservative guide for our strains.
2. Predicted metabolic profiles (351 labels) + 256-d embedding from the final MLP
   (s13 --final), and per-label logistic-regression profiles trained on all BacDive.

IDs. Internally a lab genome is 'lab_<gbk stem>', the GenBank file name in the panX
input_GenBank folder (s09). The stem is the ASSEMBLY name, not a 384 well: 'P1_D6' is Or's
sequencing plate 1 / position D6 and sits in 384 well C24 (strain 110); Schulenburg genomes
carry NCBI stems 'GCF_..._ASM..._genomic'. Outputs are keyed by the lab strain number
(`strain_no`, 20261113 rule, unique across the 298 genome wells). The full map is in
s14_lab_genome_ids.csv: strain_no, well, genome_id, assembly_name, clone_cluster.
Chain: gbk stem -> g01_genome_index (20261021: gbk -> assembly_name, well)
       -> 20261113/clean/strain_ids.csv (well -> strain_no).

Outputs: outputs/s14_lab_genome_ids.csv, s14_lab_nearest_bacdive.csv, s14_lab_profiles_mlp.csv,
         s14_lab_profiles_logreg.csv, s14_lab_embeddings_mlp.csv, s14_summary.json, figure.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegression

sys.path.insert(0, str(Path(__file__).resolve().parent))
import s12_baseline_models as s12  # noqa: E402
from s12_baseline_models import load_presence  # noqa: E402
from s13_multitask_mlp import MTL, to_torch  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "data" / "tables"
OUT = ROOT / "outputs"
MODELS = ROOT / "data" / "models"
ANALYSES = Path("/home/rl/scripts/karl/sequencing_experiments/pairwise_interaction_experiments/"
                "analyses")
CONC = ANALYSES / "20261021" / "03_genome_concordance" / "outputs"
STRAIN_IDS = ANALYSES / "20261113" / "clean" / "strain_ids.csv"


def lab_id_table(lab_ids: np.ndarray) -> pd.DataFrame:
    """genome_id ('lab_<gbk stem>') -> strain_no, well, assembly_name, clone_cluster."""
    idx = pd.read_csv(CONC / "g01_genome_index.csv", dtype=str)
    cl = pd.read_csv(CONC / "g02_genome_clusters.csv", dtype=str)
    sid = pd.read_csv(STRAIN_IDS, dtype=str)
    t = pd.DataFrame({"genome_id": lab_ids, "gbk": [g[4:] for g in lab_ids]})
    t = t.merge(idx[["gbk", "assembly_name", "well"]], on="gbk", how="left", validate="1:1")
    t = t.merge(sid[["well", "strain_no", "isolate_name"]], on="well", how="left", validate="1:1")
    # clusters are keyed by assembly_name; Schulenburg gbk stems carry an extra
    # '_ASM..._genomic' suffix, so never look them up by the gbk stem
    t = t.merge(cl[["assembly_name", "clone_cluster", "clone_size"]].drop_duplicates("assembly_name"),
                on="assembly_name", how="left", validate="1:1")
    assert t.strain_no.notna().all() and t.strain_no.is_unique, "strain_no must map 1:1"
    return t[["strain_no", "well", "isolate_name", "genome_id", "assembly_name", "gbk",
              "clone_cluster", "clone_size"]]
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"


def jac(A, B) -> np.ndarray:
    inter = (A @ B.T).toarray()
    sa = np.asarray(A.sum(axis=1)).ravel()
    sb = np.asarray(B.sum(axis=1)).ravel()
    return inter / (sa[:, None] + sb[None, :] - inter)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="id50_c80")
    a = ap.parse_args()
    P, genomes, is_bd, _ = load_presence(a.tag)
    labels = pd.read_csv(TAB / "label_matrix_mu_grouped_minority25_genome_strains.csv.gz",
                         index_col=0)
    hq = pd.read_csv(TAB / "genome_set_hq.csv").set_index("bacdive_id")
    bd_ids = np.array([int(g[2:]) for g in genomes[is_bd]])
    Xbd, Xlab = P[is_bd], P[~is_bd]
    lab_ids = genomes[~is_bd]
    S = {"n_lab": int(len(lab_ids)), "n_bacdive": int(len(bd_ids))}

    # ---- ids: strain number, well, clone cluster for lab genomes
    ids = lab_id_table(lab_ids)
    ids.to_csv(OUT / "s14_lab_genome_ids.csv", index=False)
    strain_no = pd.Index(ids.strain_no, name="strain_no")
    S["n_lab_in_clone_cluster"] = int(ids.clone_cluster.notna().sum())

    # ---- 1. domain shift
    J_lab = jac(Xlab, Xbd)
    nn_i = J_lab.argmax(axis=1)
    near = ids[["strain_no", "well", "genome_id", "clone_cluster"]].assign(
        nearest_bacdive_id=bd_ids[nn_i], jaccard=J_lab.max(axis=1))
    near = near.join(hq[["species", "genus", "family"]], on="nearest_bacdive_id")
    near.to_csv(OUT / "s14_lab_nearest_bacdive.csv", index=False)

    J_bd = jac(Xbd, Xbd)
    np.fill_diagonal(J_bd, -1)
    gen = hq.reindex(bd_ids).genus.fillna("_na").to_numpy()
    nn_strain = J_bd.max(axis=1)
    J_bd[gen[:, None] == gen[None, :]] = -1
    nn_genus = J_bd.max(axis=1)
    del J_bd
    q = [0.1, 0.25, 0.5, 0.75, 0.9]
    S["nearest_neighbour_jaccard_quantiles"] = {
        "lab_to_bacdive": np.quantile(near.jaccard, q).round(3).tolist(),
        "bacdive_strain_held_out": np.quantile(nn_strain, q).round(3).tolist(),
        "bacdive_genus_held_out": np.quantile(nn_genus, q).round(3).tolist(),
    }
    S["lab_frac_closer_than_median_genus_held_out"] = round(
        float((near.jaccard > np.median(nn_genus)).mean()), 3)
    S["lab_nearest_genus_counts"] = near.genus.value_counts().head(25).to_dict()

    fig, ax = plt.subplots(figsize=(6.4, 3.8), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    for s in ["left", "bottom"]:
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for v, c, lab in [(nn_genus, C2, "BacDive test strain, genus held out (CV regime)"),
                      (nn_strain, C1, "BacDive test strain, only itself held out"),
                      (near.jaccard.to_numpy(), C3, "our lab genomes")]:
        v = np.sort(v)
        ax.plot(v, np.arange(1, len(v) + 1) / len(v), color=c, lw=2, label=lab)
    ax.set_xlabel("Jaccard similarity to nearest BacDive training genome", color=INK2)
    ax.set_ylabel("fraction of genomes", color=INK2)
    ax.set_title("How far are genomes from the BacDive training set?", color=INK,
                 fontsize=10.5, loc="left")
    ax.legend(frameon=False, fontsize=8, labelcolor=INK2, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "figures" / "s14_domain_shift.png", dpi=160)

    # ---- 2a. MLP profiles + embeddings
    mp = MODELS / f"mlp_{a.tag}_all.pt"
    if mp.exists():
        ck = torch.load(mp)
        # the MLP was trained on a stricter prevalence filter than logreg/kNN
        s12.MIN_PREV = ck.get("min_prev", s12.MIN_PREV)
        Pm, gm, bdm, _ = load_presence(a.tag)
        s12.MIN_PREV = 10
        assert Pm.shape[1] == ck["n_in"], (Pm.shape, ck["n_in"])
        assert list(gm[~bdm]) == list(lab_ids)
        Xlab_mlp = Pm[~bdm]
        model = MTL(ck["n_in"], len(ck["labels"]))
        model.load_state_dict(ck["state_dict"])
        model.eval()
        with torch.no_grad():
            xt = to_torch(Xlab_mlp)
            prof = torch.sigmoid(model(xt)).numpy()
            emb = model.embed(xt).numpy()
        pd.DataFrame(prof, index=strain_no, columns=ck["labels"]) \
            .to_csv(OUT / "s14_lab_profiles_mlp.csv")
        pd.DataFrame(emb, index=strain_no).add_prefix("emb") \
            .to_csv(OUT / "s14_lab_embeddings_mlp.csv")
        S["mlp_model"] = str(mp)
    else:
        print("no final MLP found; run s13 --final first")

    # ---- 2b. logistic-regression profiles (C from s12)
    C = json.loads((OUT / f"s12_{a.tag}_genus_summary.json").read_text())["C"]
    Y = labels.reindex(bd_ids)

    def fit_pred(col):
        y = Y[col].to_numpy()
        m = ~np.isnan(y)
        mdl = LogisticRegression(C=C, solver="liblinear", max_iter=500).fit(Xbd[m], y[m])
        return mdl.predict_proba(Xlab)[:, 1]

    cols = list(labels.columns)
    preds = Parallel(n_jobs=50)(delayed(fit_pred)(c) for c in cols)
    pd.DataFrame(np.column_stack(preds), index=strain_no, columns=cols) \
        .to_csv(OUT / "s14_lab_profiles_logreg.csv")
    (OUT / "s14_summary.json").write_text(json.dumps(S, indent=1, default=str))
    print(json.dumps(S, indent=1, default=str))


if __name__ == "__main__":
    main()
