# 01 BacDive API & download

Script: `scripts/s01_download_bacdive.py` → `data/raw/bacdive_v2/batch_<start>.json.gz`

## API facts (checked 2026-11-20)

- Base: `https://api.bacdive.dsmz.de/v2/`. **No registration since Feb 2026.**
- v1 (`/fetch/...`) is frozen at April 2025. Use **v2**: the genome-sequence section has a new
  structure, and data added since Dec 2025 is v2-only.
- Endpoints:
  - `/v2/fetch/{id;id;...}`: full records, **≤100 IDs per call**, IDs separated by semicolons
  - `/v2/culturecollectionno/{DSM 26640}`
  - `/v2/taxon/{genus}/{species?}/{subsp?}`
  - `/v2/sequence_16s/{acc}`
  - `/v2/sequence_genome/{acc}`

  All except fetch return BacDive IDs, paginated at 100.
- **No "list all IDs" endpoint.** We sweep IDs 1..260000 in blocks of 100. Missing IDs are just
  absent from the response; an empty block is a 52-byte reply in about 0.15 s.
- ID space: two dense blocks, 1–~24.9k and ~73k–179.2k. Nothing above 179,200 (checked to 260k).
- `predictions=0` is passed explicitly. v2 already excludes AI/genome-based predictions by
  default; the `Genome-based predictions` section came back empty everywhere.
- Each record's `General.doi` carries the release: `10.13145/bacdive<ID>.20260601.11`.

## Download run

- 2,600 requests with 4 threads took 3.2 min.
- 102,187 records, 745 MB raw JSON, 117 MB gzipped.
- Resumable: existing batch files are skipped. Delete a batch file to re-fetch it.

## Record structure (top-level sections)

`General`, `Name and taxonomic classification` (BacDive's own taxonomy plus an `LPSN` sub-dict
with current names), `Morphology`, `Culture and growth conditions`, `Physiology and metabolism`,
`Isolation, sampling and environmental information`, `Interaction and safety`,
`Sequence information`, `Genome-based predictions`, `Literature`, `Reference`.

**A field is a dict when it has one item and a list when it has several.** Always normalize
with `as_list()` (see `s02`).

## Alternatives considered

- **SPARQL** (`sparql.dsmz.de`): can enumerate IDs, but may lag the API. The sweep is cheaper
  and authoritative.
- **Asking DSMZ for a bulk export:** unnecessary; the full sweep takes 3 minutes.
