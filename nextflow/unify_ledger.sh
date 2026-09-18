#!/usr/bin/env bash
# Unifies the per-accession/per-stage ledger CSVs written by WRITE_LEDGER_ROW and the
# ANALYZE stages (see modules/ledger.nf, bin/write_ledger_row.sh) into one table.
#
# This is a standalone convenience script, not part of the Nextflow workflow: run it
# by hand after a run finishes (or at any point while it's in progress) to get one
# queryable table instead of many one-row-per-file ledger CSVs.
#
# Usage:
#   ./unify_ledger.sh [OUTDIR] [OUT_DB]
#
#   OUTDIR  Nextflow run's output directory, i.e. params.outdir (default: results).
#           Must contain an OUTDIR/ledger/*.csv directory, e.g. results_500 from a
#           benchmark run.
#   OUT_DB  Path to write the DuckDB database file to (default: OUTDIR/ledger.duckdb).
#           A CSV export is also written alongside it, replacing a .duckdb suffix
#           with _combined.csv.
#
# Example queries against the result:
#   duckdb results/ledger.duckdb -c "SELECT stage, status, count(*) FROM ledger GROUP BY 1, 2 ORDER BY 1, 2;"
#   duckdb results/ledger.duckdb -c "SELECT * FROM ledger WHERE status NOT IN ('DONE', 'SKIPPED_EMPTY', 'SKIPPED_NO_HITS');"
set -euo pipefail

outdir="${1:-results}"
out_db="${2:-${outdir%/}/ledger.duckdb}"
ledger_glob="${outdir%/}/ledger/*.csv"
out_csv="${out_db%.duckdb}_combined.csv"

if ! compgen -G "$ledger_glob" > /dev/null; then
    echo "No ledger CSVs found matching: ${ledger_glob}" >&2
    exit 1
fi

rm -f "$out_db"

duckdb "$out_db" -c "
CREATE TABLE ledger AS
SELECT *
FROM read_csv_auto('${ledger_glob}', union_by_name = true);

COPY ledger TO '${out_csv}' (HEADER, DELIMITER ',');
"

n_rows="$(duckdb "$out_db" -csv -noheader -c 'SELECT count(*) FROM ledger;')"
echo "Unified ${n_rows} ledger rows from ${ledger_glob}"
echo "  DuckDB table: ${out_db} (table 'ledger')"
echo "  CSV export:   ${out_csv}"
