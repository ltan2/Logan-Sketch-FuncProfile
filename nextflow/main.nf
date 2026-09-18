#!/usr/bin/env nextflow

// Usage (from the repository root):
//   nextflow run nextflow --accessions path/to/accessions.txt [-profile test] [-resume]

nextflow.enable.dsl = 2

include { FETCH_LOGAN }      from './modules/fetch.nf'
include { WRITE_LEDGER_ROW } from './modules/ledger.nf'
include { ANALYZE }          from './subworkflows/analyze.nf'

workflow {

    // One accession per line; blank lines and '#' comments are skipped, matching
    // the original run_logan_subworkflow.sh convention.
    accession_ch = Channel.fromPath(params.accessions, checkIfExists: true)
        .splitText() { it.trim() }
        .filter { it && !it.startsWith('#') }

    fetched  = FETCH_LOGAN(accession_ch)
    /*
    *    fetched.unitigs ─┐
    *                     ├─> fetch_ch
    *    fetched.contigs ─┘
    */
    fetch_ch = fetched.unitigs.mix(fetched.contigs)   // tuple(accession, seq_type, zst_path, status)

    // Turn the manual "unitig/contig can sometimes be empty, no need to
    // conduct downstream analysis" rule into an explicit channel filter (design doc section 3.3).
    fetch_ch.branch {
        ok:    it[3] == 'OK'
        empty: it[3] == 'SKIPPED_EMPTY'
    }.set { gated }

    // Log skips to the ledger instead of silently dropping the accession.
    WRITE_LEDGER_ROW(
        gated.empty.map { row -> tuple(row[0], row[1], 'fetch', row[3]) }
    )

    // for files that are found, run analysis
    ANALYZE(
        gated.ok.map { row -> tuple(row[0], row[1], row[2]) }
    )
}
