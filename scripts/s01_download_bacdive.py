"""Bulk-download every BacDive strain record via API v2.

The API has no "list all IDs" endpoint, so we sweep the BacDive-ID space in
batches of 100 (the per-call limit of /v2/fetch/{id;id;...}). Missing IDs are
simply absent from the response, so empty batches are cheap (~50 bytes).

Each batch is cached as data/raw/bacdive_v2/batch_<start>.json.gz -> the run is
resumable: re-running skips batches already on disk.

Predictions: API v2 excludes AI/genome-based predictions by default; we pass
predictions=0 explicitly anyway so the phenotype labels are not genome-derived
(would be circular for a genome->phenotype model).

Usage:
    python s01_download_bacdive.py [--max-id 260000] [--workers 4]
"""
from __future__ import annotations

import argparse
import gzip
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "bacdive_v2"
BASE = "https://api.bacdive.dsmz.de/v2/fetch/"
BATCH = 100


def make_session() -> requests.Session:
    s = requests.Session()
    retry = Retry(total=6, backoff_factor=2,
                  status_forcelist=[429, 500, 502, 503, 504],
                  allowed_methods=["GET"])
    s.mount("https://", HTTPAdapter(max_retries=retry, pool_maxsize=16))
    s.headers["User-Agent"] = "bacdive-coverage-survey/0.1 (python-requests)"
    return s


def fetch_batch(session: requests.Session, start: int) -> tuple[int, int, int]:
    out = RAW / f"batch_{start:07d}.json.gz"
    if out.exists():
        with gzip.open(out, "rt") as f:
            return start, len(json.load(f).get("results") or {}), 0
    ids = ";".join(str(i) for i in range(start, start + BATCH))
    r = session.get(BASE + ids, params={"predictions": 0}, timeout=180)
    r.raise_for_status()
    data = r.json()
    res = data.get("results") or {}
    # sanity: the batch endpoint must not paginate for <=100 ids
    if data.get("next"):
        raise RuntimeError(f"unexpected pagination for batch {start}: {data.get('next')}")
    tmp = out.with_suffix(".tmp")
    with gzip.open(tmp, "wt") as f:
        json.dump(data, f)
    tmp.rename(out)
    return start, len(res), len(r.content)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-id", type=int, default=260_000)
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)

    starts = list(range(1, args.max_id + 1, BATCH))
    session = make_session()
    t0 = time.time()
    n_rec = n_bytes = done = 0
    max_seen = 0
    with ThreadPoolExecutor(args.workers) as ex:
        futs = {ex.submit(fetch_batch, session, s): s for s in starts}
        for fut in as_completed(futs):
            start, n, b = fut.result()
            done += 1
            n_rec += n
            n_bytes += b
            if n:
                max_seen = max(max_seen, start)
            if done % 100 == 0 or done == len(starts):
                el = time.time() - t0
                print(f"[{done}/{len(starts)}] records={n_rec} "
                      f"downloaded={n_bytes/1e6:.0f}MB elapsed={el/60:.1f}min "
                      f"highest non-empty batch start={max_seen}", flush=True)
    summary = {"n_records": n_rec, "highest_nonempty_batch_start": max_seen,
               "max_id_swept": args.max_id,
               "finished": time.strftime("%Y-%m-%d %H:%M:%S")}
    (RAW / "_download_summary.json").write_text(json.dumps(summary, indent=1))
    print(summary)
    if max_seen + BATCH > args.max_id - 20_000:
        print("WARNING: records found near the end of the swept range; "
              "re-run with a larger --max-id")


if __name__ == "__main__":
    main()
