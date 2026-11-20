"""Label reliability: how often do independent references disagree?

A strain x metabolite x assay-kind label measured by >=2 different BacDive
references (@ref) gives a replicate pair. The disagreement rate among those is
an empirical label-noise estimate (a ceiling for any model), which the overall
"conflict rate" in s03 hides because most labels have a single observation.

Sources: @ref 68360-68399 are DSMZ's in-house API-kit series (see
s03_api_kits_folded_into_mu.csv); everything else is literature curation.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s03_coverage_report import KIND_GROUP, merge_unspecified_stereo  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "data" / "tables"
OUT = ROOT / "outputs"


def per_ref_calls(mu: pd.DataFrame, kind_col: str) -> pd.DataFrame:
    """Collapse to one call per (strain, metabolite, kind, ref); drop refs that
    contradict themselves."""
    r = (mu.groupby(["bacdive_id", "metabolite_key", kind_col, "@ref"])
           .agg(act=("activity", lambda s: "+" if (s == "+").all()
                     else ("-" if (s == "-").all() else "mix")),
                src=("src", "first")).reset_index())
    return r[r.act != "mix"]


def disagreement(r: pd.DataFrame, kind_col: str) -> pd.DataFrame:
    g = r.groupby(["bacdive_id", "metabolite_key", kind_col]).agg(
        n_ref=("@ref", "nunique"), n_pos=("act", lambda s: (s == "+").sum()),
        sources=("src", lambda s: "+".join(sorted(set(s)))))
    g = g[g.n_ref >= 2].copy()
    g["disagree"] = (g.n_pos > 0) & (g.n_pos < g.n_ref)
    return g


def main() -> None:
    mu = pd.read_csv(TAB / "metabolite_utilization.csv.gz", low_memory=False).rename(
        columns={"utilization activity": "activity", "kind of utilization tested": "kind"})
    mu = mu[mu.activity.isin(["+", "-"])]
    mu["metabolite_key"] = (mu.metabolite.astype(str).str.strip().str.lower()
                              .str.replace(r"\s+", " ", regex=True))
    mu, _ = merge_unspecified_stereo(mu)
    mu["kind_group"] = mu.kind.map(KIND_GROUP).fillna("other")
    mu["src"] = np.where(mu["@ref"].between(68360, 68399), "DSMZ_API_kit", "literature")

    S = {}
    r_exact = per_ref_calls(mu, "kind")
    d_exact = disagreement(r_exact, "kind")
    S["exact_kind_n_replicated"] = len(d_exact)
    S["exact_kind_disagreement"] = round(float(d_exact.disagree.mean()), 3)
    S["exact_kind_by_kind"] = (d_exact.groupby(level="kind").disagree.agg(["size", "mean"])
                               .query("size >= 200").round(3).to_dict(orient="index"))
    S["exact_kind_by_source_pair"] = (d_exact.groupby("sources").disagree
                                      .agg(["size", "mean"]).round(3).to_dict(orient="index"))
    # direction of disagreement between DSMZ kits and literature
    dis_idx = d_exact[d_exact.disagree & (d_exact.sources == "DSMZ_API_kit+literature")].index
    x = r_exact.set_index(["bacdive_id", "metabolite_key", "kind"]).loc[dis_idx]
    S["kit_vs_literature_disagreements_calls"] = {
        f"{s}:{a}": int(n) for (s, a), n in x.groupby(["src", "act"]).size().items()}

    r_grp = per_ref_calls(mu, "kind_group")
    d_grp = disagreement(r_grp, "kind_group")
    S["kind_group_n_replicated"] = len(d_grp)
    S["kind_group_disagreement"] = round(float(d_grp.disagree.mean()), 3)
    S["kind_group_by_group"] = (d_grp.groupby(level="kind_group").disagree
                                .agg(["size", "mean"]).round(3).to_dict(orient="index"))

    per_label = (d_exact.reset_index().groupby(["metabolite_key", "kind"])
                 .disagree.agg(n_replicated="size", disagreement="mean")
                 .query("n_replicated >= 30").sort_values("disagreement"))
    per_label.to_csv(OUT / "s07_disagreement_per_label.csv")
    (OUT / "s07_label_reliability_summary.json").write_text(json.dumps(S, indent=1))
    print(json.dumps(S, indent=1))


if __name__ == "__main__":
    main()
