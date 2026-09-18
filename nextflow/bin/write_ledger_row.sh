#!/usr/bin/env bash
# Appends one row to the file-based progress ledger (see ../../design/latest_design-doc.md section 5).
#
# Each Nextflow task writes its own small CSV file instead of all tasks writing into one
# shared DuckDB file, because many FETCH_LOGAN/SOURMASH_SKETCH/etc. tasks run concurrently
# (one per accession) and DuckDB does not support concurrent writers to a single file.
# The ledger is queried later by unioning these files, e.g.:
#   duckdb -c "SELECT * FROM read_csv_auto('results/ledger/*.csv', union_by_name=true)"
set -euo pipefail

outfile=${1:?"Usage: $0 OUTFILE ACCESSION SEQ_TYPE STAGE STATUS [SHARD_ID]"}
accession=${2:?}
seq_type=${3:?}
stage=${4:?}
status=${5:?}
shard_id=${6:-NA}

updated_at="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

{
    echo "accession,seq_type,stage,status,shard_id,updated_at"
    echo "${accession},${seq_type},${stage},${status},${shard_id},${updated_at}"
} > "$outfile"
