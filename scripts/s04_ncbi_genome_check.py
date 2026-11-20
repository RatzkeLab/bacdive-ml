"""Check every chosen BacDive genome accession against NCBI Datasets (metadata only).

For each strain's chosen genome (data/tables/genome_choice.csv.gz from s03) fetch
the NCBI assembly report: does the accession still resolve, current version,
assembly level, size, CheckM completeness/contamination, NCBI's ANI taxonomy
check, and whether a RefSeq (GCF) partner exists -- that last one means PGAP
protein annotations can be downloaded directly instead of running gene calling.

Raw reports are cached as JSONL (data/raw/ncbi_reports/*.jsonl.gz, one file per
request batch) so re-runs are incremental.

NCBI rate limit without an API key: 5 requests/s; we stay well under it.
"""
from __future__ import annotations

import gzip
import json
import time
from pathlib import Path

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "data" / "tables"
RAW = ROOT / "data" / "raw" / "ncbi_reports"
OUT = ROOT / "outputs"
URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/dataset_report"
BATCH = 200


def session() -> requests.Session:
    s = requests.Session()
    s.mount("https://", HTTPAdapter(max_retries=Retry(
        total=6, backoff_factor=2, status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["POST", "GET"])))
    s.headers["User-Agent"] = "bacdive-coverage-survey/0.1 (python-requests)"
    return s


def fetch(s: requests.Session, accs: list[str]) -> list[dict]:
    reports, token = [], None
    while True:
        body = {"accessions": accs, "page_size": 1000}
        if token:
            body["page_token"] = token
        r = s.post(URL, json=body, timeout=180)
        r.raise_for_status()
        d = r.json()
        reports += d.get("reports", [])
        token = d.get("next_page_token")
        if not token:
            return reports
        time.sleep(0.3)


def flat(r: dict) -> dict:
    ai, st = r.get("assembly_info", {}), r.get("assembly_stats", {})
    ck, ani = r.get("checkm_info", {}), r.get("average_nucleotide_identity", {})
    ann = r.get("annotation_info", {})
    return {
        "ncbi_accession": r.get("accession"),
        "acc_base": (r.get("accession") or "").split(".")[0],
        "current_accession": r.get("current_accession"),
        "paired_accession": r.get("paired_accession"),
        "assembly_level_ncbi": ai.get("assembly_level"),
        "assembly_status": ai.get("assembly_status"),
        "refseq_category": ai.get("refseq_category"),
        "paired_status": (ai.get("paired_assembly") or {}).get("status"),
        "release_date": ai.get("release_date"),
        "genome_size": pd.to_numeric(st.get("total_sequence_length"), errors="coerce"),
        "n_contigs": st.get("number_of_contigs"),
        "contig_n50": st.get("contig_n50"),
        "gc_percent": st.get("gc_percent"),
        "checkm_completeness": ck.get("completeness"),
        "checkm_contamination": ck.get("contamination"),
        "ani_taxonomy_check": ani.get("taxonomy_check_status"),
        "ani_match_status": ani.get("match_status"),
        "ani_category": ani.get("category"),
        "n_protein_coding": ((ann.get("stats") or {}).get("gene_counts") or {}).get("protein_coding"),
        "annotation_provider": ann.get("provider"),
        "ncbi_organism": (r.get("organism") or {}).get("organism_name"),
        "ncbi_taxid": (r.get("organism") or {}).get("tax_id"),
    }


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    ch = pd.read_csv(TAB / "genome_choice.csv.gz")
    accs = sorted(ch.acc_base.dropna().unique())
    print(f"{len(accs)} unique chosen accessions")
    s = session()
    for i in range(0, len(accs), BATCH):
        fp = RAW / f"reports_{i:06d}.jsonl.gz"
        if fp.exists():
            continue
        reps = fetch(s, accs[i:i + BATCH])
        with gzip.open(fp.with_suffix(".tmp"), "wt") as f:
            for r in reps:
                f.write(json.dumps(r) + "\n")
        fp.with_suffix(".tmp").rename(fp)
        if (i // BATCH) % 10 == 0:
            print(f"  {i + BATCH}/{len(accs)}", flush=True)
        time.sleep(0.4)

    rows = []
    for fp in sorted(RAW.glob("reports_*.jsonl.gz")):
        with gzip.open(fp, "rt") as f:
            rows += [flat(json.loads(line)) for line in f]
    rep = pd.DataFrame(rows).drop_duplicates("acc_base")
    m = ch.merge(rep, on="acc_base", how="left")
    m.to_csv(TAB / "genome_choice_ncbi.csv.gz", index=False)

    found = m.ncbi_accession.notna()
    summ = {
        "chosen_accessions": len(accs),
        "resolved_at_ncbi": int(found.sum()),
        "not_resolved": int((~found).sum()),
        "assembly_status": m.assembly_status.value_counts(dropna=False).to_dict(),
        "has_refseq_partner_current": int(m.paired_accession.notna().sum()
                                          and (m.paired_status == "current").sum()),
        "checkm_available": int(m.checkm_completeness.notna().sum()),
        "checkm_ge90_le5": int(((m.checkm_completeness >= 90) &
                                (m.checkm_contamination <= 5)).sum()),
        "ani_taxonomy_check": m.ani_taxonomy_check.value_counts(dropna=False).to_dict(),
        "genome_size_total_gb": round(float(m.genome_size.sum()) / 1e9, 2),
        "genome_size_median_mb": round(float(m.genome_size.median()) / 1e6, 2),
    }
    (OUT / "s04_ncbi_summary.json").write_text(json.dumps(summ, indent=1, default=str))
    print(json.dumps(summ, indent=1, default=str))


if __name__ == "__main__":
    main()
