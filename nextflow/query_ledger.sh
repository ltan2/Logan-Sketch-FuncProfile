#!/usr/bin/env bash
# Accessions whose sketch finished but whose KO profile never did -- i.e. work the next
# run still owes. See unify_ledger.sh for building one queryable table out of the ledger.

duckdb -c "
WITH ledger AS (
    SELECT *
    FROM read_csv_auto('results/ledger/*.csv', union_by_name = true)
)
SELECT accession, seq_type
FROM ledger
WHERE stage = 'sketch' AND status = 'DONE'
EXCEPT
SELECT accession, seq_type
FROM ledger
WHERE stage = 'funprofiler';
"
