// Generic progress-ledger writer, used for events that don't otherwise produce a
// ledger row of their own (e.g. logging a SKIPPED_EMPTY fetch). See bin/write_ledger_row.sh
// and ../../design/latest_design-doc.md section 5.

process WRITE_LEDGER_ROW {
    tag "${accession}.${seq_type}.${stage}"
    publishDir "${params.outdir}/ledger", mode: 'copy'

    input:
    tuple val(accession), val(seq_type), val(stage), val(status)

    output:
    path("${accession}.${seq_type}.${stage}.ledger.csv")

    script:
    """
    write_ledger_row.sh "${accession}.${seq_type}.${stage}.ledger.csv" "${accession}" "${seq_type}" "${stage}" "${status}"
    """
}
