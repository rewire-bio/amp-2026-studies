# EVEE source

Written by `lookup_scores.py` when it built `data/evee/evee_slim.parquet` (PROTOCOL.md section 10).

- Record: https://zenodo.org/records/19701997 (EVEE: Interpretable variant effect prediction from genomic foundation model embeddings), licence cc-by-4.0
- Version: not stated; revision 4; published 2026-04-23; last updated 2026-04-23T01:51:33.674271+00:00
- Downloaded: 2026-11-02, columns chrom, pos, ref, alt, pathogenicity, label, significance, stars, consequence, variation_id only, read by byte range
- Rows: 4252870 (825720, 826113, 871493, 849270, 880274 per shard), all keys distinct
- `evee_slim.parquet` SHA-256: f4e9d04415e1b9d79da099dd9c8cbcd0c46f82f920307a8e88bea626ea2d58fc

Only the ten columns were read, so the checksums below are Zenodo's for the full shards, not checked here.

| File | Size (bytes) | Checksum (Zenodo) |
|---|---|---|
| clean_shard_0.parquet | 6805237030 | md5:6b25a426954b7590c04e93b42bd0624f |
| clean_shard_1.parquet | 6758853203 | md5:5f9fe5fb06689b0b7b6f412d86975e05 |
| clean_shard_2.parquet | 7157309607 | md5:ffbf39a9d047897932bbd4c7d43c05a2 |
| clean_shard_3.parquet | 6985087155 | md5:4714443ef1ac68edf7cc92fb49587df7 |
| clean_shard_4.parquet | 7303028172 | md5:fa19cce9b5d64235eedce017efb62074 |
