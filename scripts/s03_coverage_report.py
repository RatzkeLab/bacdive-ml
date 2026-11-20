"""Coverage / feasibility report: how many usable genome->phenotype examples exist?

Answers, for BacDive v2 (predictions excluded):
  1. strains with metabolite-utilization observations
  2. strains with a genome assembly accession (GCA_/GCF_)
  3. strains with both
  4. +/- counts per (metabolite, kind-of-test) label, restricted to genome strains
  5. strains tested for >=10/20/50/100 labels
plus: conflicts between replicate observations, phylogenetic redundancy
(species / genus counts), coverage of the genera in our lab collection, whether
any lab strains (MYb/BIGb/JUb/CEent) are in BacDive, whether API-kit results are
already folded into the metabolite-utilization list, and a strain x label matrix.

Label convention (per strain, label = metabolite x kind):
    1  only '+' observations     0  only '-' observations
    NA untested, only '+/-' (variable/weak), or conflicting '+' and '-'
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "data" / "tables"
OUT = ROOT / "outputs"
FIG = OUT / "figures"
COLLECTION = Path("/home/rl/scripts/karl/data_links/resources/"
                  "Strain collection including Schulenberg's - collection.csv")

KIND_GROUP = {
    "carbon source": "growth/assimilation", "energy source": "growth/assimilation",
    "nitrogen source": "growth/assimilation", "growth": "growth/assimilation",
    "assimilation": "growth/assimilation",
    "fermentation": "acid/gas (fermentation)", "builds acid from": "acid/gas (fermentation)",
    "builds gas from": "acid/gas (fermentation)",
    "hydrolysis": "hydrolysis/degradation", "degradation": "hydrolysis/degradation",
    "respiration": "respiration/oxidation", "oxidation": "respiration/oxidation",
    "reduction": "reduction/e- acceptor", "electron acceptor": "reduction/e- acceptor",
}
LEVEL_RANK = {"complete": 0, "chromosome": 1, "scaffold": 2, "contig": 3}
MINORITY_THRESHOLDS = [10, 25, 50, 100]

# reference palette (dataviz skill, light mode)
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
C1, C2 = "#2a78d6", "#eb6834"


def style_ax(ax):
    ax.set_facecolor(SURFACE)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    for s in ["left", "bottom"]:
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(axis="both", color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def load(name, **kw):
    return pd.read_csv(TAB / f"{name}.csv.gz", low_memory=False, **kw)


def aggregate_labels(obs: pd.DataFrame, keys: list[str], val: str) -> pd.DataFrame:
    """obs has bacdive_id, keys..., val in {'+','-','+/-',...}; one row per
    strain x key with n_pos/n_neg/n_var and the resolved label."""
    o = obs.assign(pos=obs[val].eq("+"), neg=obs[val].eq("-"),
                   var=~obs[val].isin(["+", "-"]))
    g = o.groupby(["bacdive_id", *keys], dropna=False)[["pos", "neg", "var"]].sum()
    g.columns = ["n_pos", "n_neg", "n_var"]
    g = g.reset_index()
    g["label"] = np.where((g.n_pos > 0) & (g.n_neg == 0), 1.0,
                          np.where((g.n_neg > 0) & (g.n_pos == 0), 0.0, np.nan))
    g["conflict"] = (g.n_pos > 0) & (g.n_neg > 0)
    return g


def merge_unspecified_stereo(mu: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """'glucose' and 'D-glucose' are the same test when the other enantiomer is
    never tested: merge 'x' into 'd-x' if 'l-x' is absent (and vice versa).
    When both D- and L- forms are tested (xylose, arabinose, alanine, ...) the
    unprefixed name is ambiguous and is kept as its own label."""
    keys = set(mu.metabolite_key)
    mapping = {}
    for k in keys:
        d, l = "d-" + k, "l-" + k
        if d in keys and l not in keys:
            mapping[k] = d
        elif l in keys and d not in keys:
            mapping[k] = l
    rep = pd.DataFrame([{"unspecified_name": k, "merged_into": v,
                         "n_rows_moved": int((mu.metabolite_key == k).sum())}
                        for k, v in sorted(mapping.items())])
    mu = mu.assign(metabolite_key=mu.metabolite_key.replace(mapping))
    return mu, rep


def label_stats(lab: pd.DataFrame, keys: list[str], strains: pd.DataFrame) -> pd.DataFrame:
    x = lab.merge(strains[["bacdive_id", "genus", "family", "species"]], on="bacdive_id")
    rows = []
    for k, d in x.groupby(keys, dropna=False):
        dd = d[d.label.notna()]
        pos = dd[dd.label == 1]
        neg = dd[dd.label == 0]
        rows.append({
            **dict(zip(keys, k if isinstance(k, tuple) else (k,))),
            "n_strains_tested": len(d), "n_pos": len(pos), "n_neg": len(neg),
            "n_conflict": int(d.conflict.sum()),
            "n_variable_only": int(((d.n_var > 0) & (d.n_pos == 0) & (d.n_neg == 0)).sum()),
            "pos_frac": len(pos) / max(1, len(dd)),
            "minority": min(len(pos), len(neg)),
            "n_species": dd.species.nunique(), "n_genera": dd.genus.nunique(),
            "n_families": dd.family.nunique(),
            # label variation that is NOT just between genera: genera with both classes
            "n_genera_with_both": int((dd.groupby("genus").label.nunique() > 1).sum()),
        })
    return pd.DataFrame(rows).sort_values("minority", ascending=False)


def main() -> None:
    OUT.mkdir(exist_ok=True)
    FIG.mkdir(exist_ok=True)
    S = {}  # summary numbers

    strains = load("strains")
    mu = load("metabolite_utilization")
    gen = load("genomes")
    S["n_strains"] = len(strains)
    S["n_type_strains"] = int(strains.type_strain.eq("yes").sum())

    # ---------------- genomes ----------------
    acc_col = "INSDC accession"
    gen["acc"] = gen[acc_col].astype(str).str.strip()
    gen["acc_ok"] = gen.acc.str.match(r"^GC[AF]_\d{9}(\.\d+)?$")
    S["genome_rows"] = len(gen)
    S["genome_rows_with_valid_acc"] = int(gen.acc_ok.sum())
    S["genome_columns"] = list(gen.columns)
    g = gen[gen.acc_ok].copy()
    g["acc_base"] = g.acc.str.replace(r"\.\d+$", "", regex=True)
    g["level_rank"] = g["assembly level"].map(LEVEL_RANK).fillna(9)
    g["score"] = pd.to_numeric(g.get("score"), errors="coerce")
    S["genomes_per_strain_dist"] = g.groupby("bacdive_id").acc_base.nunique() \
        .value_counts().sort_index().head(10).to_dict()
    acc_strains = g.groupby("acc_base").bacdive_id.nunique()
    S["accessions_shared_by_multiple_strains"] = int((acc_strains > 1).sum())
    best = (g.sort_values(["bacdive_id", "level_rank", "score"],
                          ascending=[True, True, False])
             .drop_duplicates("bacdive_id"))
    best = best.merge(acc_strains.rename("n_strains_sharing_acc"), left_on="acc_base",
                      right_index=True)
    best.to_csv(TAB / "genome_choice.csv.gz", index=False)
    S["assembly_level_of_chosen"] = best["assembly level"].value_counts().to_dict()
    S["score_quantiles_of_chosen"] = best.score.quantile([.05, .25, .5, .75, .95]).round(1).to_dict()
    genome_ids = set(best.bacdive_id)

    # ---------------- metabolite utilization labels ----------------
    mu = mu.rename(columns={"Chebi-ID": "chebi", "utilization activity": "activity",
                            "kind of utilization tested": "kind"})
    S["mu_rows"] = len(mu)
    S["mu_activity_values"] = mu.activity.value_counts(dropna=False).head(10).to_dict()
    S["mu_kind_values"] = mu.kind.value_counts(dropna=False).to_dict()
    # Key labels by normalized NAME, not ChEBI: BacDive uses ChEBI 0 as a placeholder
    # (merges 2- and 5-ketogluconate) and has a few wrong ChEBI ids (acetic acid vs
    # indole-3-acetic acid), while names are consistent apart from case.
    S["mu_rows_without_chebi_or_chebi0"] = int((mu.chebi.isna() | mu.chebi.eq(0)).sum())
    mu["metabolite_key"] = mu.metabolite.astype(str).str.strip().str.lower() \
        .str.replace(r"\s+", " ", regex=True)
    mu, stereo_merges = merge_unspecified_stereo(mu)
    stereo_merges.to_csv(OUT / "s03_stereo_name_merges.csv", index=False)
    S["n_stereo_unspecified_names_merged"] = len(stereo_merges)
    names = (mu.groupby("metabolite_key").metabolite
               .agg(lambda s: s.value_counts().index[0]))
    mu["kind_group"] = mu.kind.map(KIND_GROUP).fillna("other")
    mu = mu[mu.activity.notna()]

    lab = aggregate_labels(mu, ["metabolite_key", "kind"], "activity")
    lab_grp = aggregate_labels(mu, ["metabolite_key", "kind_group"], "activity")
    lab["metabolite"] = lab.metabolite_key.map(names)
    lab_grp["metabolite"] = lab_grp.metabolite_key.map(names)

    mu_ids = set(lab.loc[lab.label.notna(), "bacdive_id"])
    both = mu_ids & genome_ids
    S["strains_with_mu_resolved_label"] = len(mu_ids)
    S["strains_with_genome"] = len(genome_ids)
    S["strains_with_both"] = len(both)
    S["strains_with_both_type_strain"] = int(
        strains[strains.bacdive_id.isin(both)].type_strain.eq("yes").sum())
    S["label_conflict_rate_all"] = round(float(lab.conflict.mean()), 4)
    S["label_conflict_rate_grouped"] = round(float(lab_grp.conflict.mean()), 4)

    # per-label stats on genome strains
    lab_g = lab[lab.bacdive_id.isin(genome_ids)]
    lab_grp_g = lab_grp[lab_grp.bacdive_id.isin(genome_ids)]
    st = label_stats(lab_g, ["metabolite_key", "kind"], strains)
    st.insert(1, "metabolite", st.metabolite_key.map(names))
    st.to_csv(OUT / "s03_label_stats_mu_by_kind.csv", index=False)
    stg = label_stats(lab_grp_g, ["metabolite_key", "kind_group"], strains)
    stg.insert(1, "metabolite", stg.metabolite_key.map(names))
    stg.to_csv(OUT / "s03_label_stats_mu_by_kindgroup.csv", index=False)
    S["n_labels_by_kind_minority_ge"] = {t: int((st.minority >= t).sum())
                                        for t in MINORITY_THRESHOLDS}
    S["n_labels_by_kindgroup_minority_ge"] = {t: int((stg.minority >= t).sum())
                                             for t in MINORITY_THRESHOLDS}
    S["usable_labels_per_kindgroup_minority_ge25"] = (
        stg[stg.minority >= 25].kind_group.value_counts().to_dict())

    # labels per strain (genome strains, resolved labels only)
    per_strain = lab_g[lab_g.label.notna()].groupby("bacdive_id").size()
    per_strain_grp = lab_grp_g[lab_grp_g.label.notna()].groupby("bacdive_id").size()
    S["genome_strains_with_ge_k_labels"] = {k: int((per_strain >= k).sum())
                                            for k in [1, 10, 20, 50, 100]}
    usable_keys = set(map(tuple, stg.loc[stg.minority >= 25,
                                         ["metabolite_key", "kind_group"]].values))
    lab_use = lab_grp_g[[k in usable_keys for k in
                         zip(lab_grp_g.metabolite_key, lab_grp_g.kind_group)]]
    per_strain_use = lab_use[lab_use.label.notna()].groupby("bacdive_id").size()
    S["genome_strains_with_ge_k_usable_grouped_labels(minority>=25)"] = {
        k: int((per_strain_use >= k).sum()) for k in [1, 5, 10, 20, 50]}

    # label matrix: genome strains x usable grouped labels
    lab_use = lab_use.assign(col=lab_use.metabolite_key + "|" + lab_use.kind_group)
    mat = lab_use.pivot_table(index="bacdive_id", columns="col", values="label",
                              aggfunc="first")
    mat.to_csv(TAB / "label_matrix_mu_grouped_minority25_genome_strains.csv.gz")
    S["label_matrix_shape"] = list(mat.shape)
    S["label_matrix_fill"] = round(float(mat.notna().values.mean()), 4)

    # ---------------- phylogenetic redundancy ----------------
    sb = strains[strains.bacdive_id.isin(both)]
    S["both_n_species"] = int(sb.species.nunique())
    S["both_n_genera"] = int(sb.genus.nunique())
    S["both_n_families"] = int(sb.family.nunique())
    S["both_n_phyla"] = int(sb.phylum.nunique())
    S["both_top_phyla"] = sb.phylum.value_counts().head(8).to_dict()
    S["both_top_genera"] = sb.genus.value_counts().head(25).to_dict()
    S["both_strains_per_species_quantiles"] = sb.species.value_counts() \
        .quantile([.5, .9, .99]).to_dict()

    # ---------------- other phenotype blocks on genome strains ----------------
    other = {}
    for name in ["enzymes", "api_kits", "metabolite_production", "metabolite_tests",
                 "traits", "culture_temp", "culture_ph", "halophily",
                 "antibiotic_resistance", "isolation_cat"]:
        try:
            t = load(name)
        except FileNotFoundError:
            continue
        other[name] = {"rows": len(t), "strains": int(t.bacdive_id.nunique()),
                       "genome_strains": int(t.bacdive_id.isin(genome_ids).sum() and
                                             t[t.bacdive_id.isin(genome_ids)].bacdive_id.nunique())}
    S["other_blocks"] = other

    enz = load("enzymes")
    enz = enz[enz.activity.notna()]
    enz_lab = aggregate_labels(enz, ["value"], "activity")
    enz_st = label_stats(enz_lab[enz_lab.bacdive_id.isin(genome_ids)], ["value"], strains)
    enz_st.to_csv(OUT / "s03_label_stats_enzymes.csv", index=False)
    S["n_enzyme_labels_minority_ge25"] = int((enz_st.minority >= 25).sum())

    api = load("api_kits")
    S["api_kit_strains"] = api.groupby("kit").bacdive_id.nunique().sort_values(
        ascending=False).to_dict()
    S["api_kit_genome_strains"] = api[api.bacdive_id.isin(genome_ids)].groupby(
        "kit").bacdive_id.nunique().sort_values(ascending=False).to_dict()
    # Are API-kit results already folded into metabolite utilization / enzymes?
    # DSMZ files each kit's results in MU under a kit-specific @ref (683xx), so for
    # each kit report the MU ref found in the most kit strains and its coverage.
    folded = {}
    enz_all = load("enzymes")
    for kit, d in api.groupby("kit"):
        ids = set(d.bacdive_id)
        sub = mu[mu.bacdive_id.isin(ids)]
        cov_mu = sub.groupby("@ref").bacdive_id.nunique().sort_values(ascending=False)
        esub = enz_all[enz_all.bacdive_id.isin(ids)]
        cov_enz = esub.groupby("@ref").bacdive_id.nunique().sort_values(ascending=False)
        folded[kit] = {
            "n_strains": len(ids),
            "best_mu_ref": int(cov_mu.index[0]) if len(cov_mu) else None,
            "best_mu_ref_strain_frac": round(cov_mu.iloc[0] / len(ids), 2) if len(cov_mu) else 0,
            "best_enzyme_ref": int(cov_enz.index[0]) if len(cov_enz) else None,
            "best_enzyme_ref_strain_frac": round(cov_enz.iloc[0] / len(ids), 2) if len(cov_enz) else 0,
        }
    pd.DataFrame(folded).T.to_csv(OUT / "s03_api_kits_folded_into_mu.csv")
    ref_share = mu[mu.bacdive_id.isin(genome_ids)]["@ref"].value_counts(normalize=True)
    S["mu_rows_genome_strains_top10_ref_share"] = round(float(ref_share.head(10).sum()), 3)
    S["mu_top10_refs"] = ref_share.head(10).round(3).to_dict()

    traits = load("traits")
    tg = traits[traits.bacdive_id.isin(genome_ids)]
    S["traits_genome_strains"] = {t: int(d.bacdive_id.nunique()) for t, d in tg.groupby("trait")}
    S["gram_stain_values"] = tg[tg.trait == "cell morphology:gram stain"].value \
        .value_counts().head(5).to_dict()
    S["oxygen_values"] = tg[tg.trait == "oxygen tolerance"].value.value_counts().head(8).to_dict()

    # ---------------- relevance to our collection ----------------
    col = pd.read_csv(COLLECTION)
    ph = col["phylogeny"].fillna("")
    col["genus"] = ph.str.extract(r"g__([^;]*)")[0].str.strip()
    col["family"] = ph.str.extract(r"f__([^;]*)")[0].str.strip().str.strip("[]")
    cg = col.genus.replace("", np.nan).dropna().value_counts()
    rows = []
    n_lab_use = per_strain_use.reindex(strains.bacdive_id).fillna(0).values
    st2 = strains.assign(has_genome=strains.bacdive_id.isin(genome_ids),
                         has_both=strains.bacdive_id.isin(both),
                         n_usable=n_lab_use)
    for gname, n in cg.items():
        m = st2[(st2.genus == gname) | (st2.genus_bacdive == gname)]
        rows.append({"genus": gname, "n_in_collection": n, "bacdive_strains": len(m),
                     "with_genome": int(m.has_genome.sum()), "with_genome_and_mu": int(m.has_both.sum()),
                     "with_genome_and_ge10_usable_labels": int((m.has_genome & (m.n_usable >= 10)).sum()),
                     "species_with_genome_and_mu": int(m[m.has_both].species.nunique())})
    cov = pd.DataFrame(rows)
    cov.to_csv(OUT / "s03_collection_genus_coverage.csv", index=False)
    S["collection_strains_in_genera_with_>=10_both"] = int(
        cov.loc[cov.with_genome_and_mu >= 10, "n_in_collection"].sum())
    S["collection_strains_with_genus"] = int(cg.sum())
    # match family on LPSN OR BacDive's own family: LPSN renamed several of ours
    # (Xanthomonadaceae->Lysobacteraceae, Leuconostocaceae->Lactobacillaceae, ...)
    sbf = strains[strains.bacdive_id.isin(both)]
    cf = col.family.replace("", np.nan).dropna().value_counts()
    S["collection_families_cov"] = {
        f: int(((sbf.family == f) | (sbf.family_bacdive == f)).sum()) for f in cf.index}

    # lab strain names in BacDive?
    pat = re.compile(r"\b(MYb|BIGb|JUb|CEent|MSPm)\s?-?\d+", re.I)
    hay = (strains.strain_designation.fillna("") + " | " +
           strains.culture_collection_nos.fillna("") + " | " + strains.description.fillna(""))
    hit = strains[hay.str.contains(pat)].copy()
    hit["match"] = hay[hit.index].str.findall(pat.pattern + r"\w*").map(lambda x: x)
    hit = hit.assign(match=hay[hit.index].map(lambda s: ",".join(
        sorted(set(m.group(0) for m in pat.finditer(s))))))
    hit["has_genome"] = hit.bacdive_id.isin(genome_ids)
    hit["has_mu"] = hit.bacdive_id.isin(mu_ids)
    hit[["bacdive_id", "match", "species", "strain_designation", "culture_collection_nos",
         "has_genome", "has_mu"]].to_csv(OUT / "s03_lab_strain_name_hits.csv", index=False)
    S["lab_strain_name_hits"] = len(hit)
    S["lab_strain_name_hits_with_mu"] = int(hit.has_mu.sum())

    # ---------------- figures ----------------
    # 1. how many genome strains have >= k labels
    fig, ax = plt.subplots(figsize=(6.4, 3.8), facecolor=SURFACE)
    style_ax(ax)
    ks = np.arange(1, 201)
    ax.plot(ks, [(per_strain >= k).sum() for k in ks], color=C1, lw=2,
            label="any metabolite x kind label")
    ax.plot(ks, [(per_strain_use >= k).sum() for k in ks], color=C2, lw=2,
            label="usable grouped labels (minority class >= 25)")
    ax.set_yscale("log")
    ax.set_xlabel("k  (labels per strain)", color=INK2)
    ax.set_ylabel("genome strains with >= k labels", color=INK2)
    ax.set_title("Phenotype depth of BacDive strains that have a genome", color=INK,
                 fontsize=11, loc="left")
    ax.legend(frameon=False, fontsize=8, labelcolor=INK2)
    fig.tight_layout()
    fig.savefig(FIG / "s03_labels_per_genome_strain.png", dpi=160)
    plt.close(fig)

    # 2. usable labels per kind group at minority >= 25
    vc = stg[stg.minority >= 25].kind_group.value_counts().sort_values()
    fig, ax = plt.subplots(figsize=(6.4, 3.0), facecolor=SURFACE)
    style_ax(ax)
    ax.grid(axis="y", visible=False)
    ax.barh(vc.index, vc.values, color=C1, height=0.6)
    for y, v in enumerate(vc.values):
        ax.text(v, y, f" {v}", va="center", color=INK2, fontsize=8)
    ax.set_xlabel("metabolite labels with >= 25 positive AND >= 25 negative genome strains",
                  color=INK2, fontsize=8)
    ax.set_title("Usable labels by assay type", color=INK, fontsize=11, loc="left")
    fig.tight_layout()
    fig.savefig(FIG / "s03_usable_labels_by_kind.png", dpi=160)
    plt.close(fig)

    # 3. genome+phenotype strains by genus, collection genera highlighted
    top = sb.genus.value_counts().head(30).sort_values()
    in_col = [g in set(cg.index) for g in top.index]
    fig, ax = plt.subplots(figsize=(6.4, 6.4), facecolor=SURFACE)
    style_ax(ax)
    ax.grid(axis="y", visible=False)
    ax.barh(top.index, top.values, color=[C2 if c else C1 for c in in_col], height=0.65)
    ax.set_xlabel("strains with genome + metabolite-utilization data", color=INK2)
    ax.set_title("Top genera (orange = genus present in lab collection)", color=INK,
                 fontsize=10, loc="left")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=C2, label="in lab collection"),
                       Patch(color=C1, label="not in collection")],
              frameon=False, fontsize=8, labelcolor=INK2, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "s03_top_genera.png", dpi=160)
    plt.close(fig)

    (OUT / "s03_coverage_summary.json").write_text(json.dumps(S, indent=1, default=str))
    print(json.dumps(S, indent=1, default=str))


if __name__ == "__main__":
    main()
