// Build a sourmash sketch from the FASTA produced by the shared decompression process.

process SOURMASH_SKETCH {
    tag "${accession}.${seq_type}"
    // NOTE: path must be a closure here, not a plain interpolated string -- publishDir's
    // path is evaluated once at process-definition time unless wrapped in `{ ... }` to
    // defer it until a task actually runs and `seq_type` (an input variable) is bound.
    // `tag` above is exempt from this -- Nextflow always evaluates it per task.
    publishDir(path: { "${params.outdir}/sketches/${seq_type}" }, mode: 'copy', pattern: "*.sig.zip")
    publishDir "${params.outdir}/ledger", mode: 'copy', pattern: "*.ledger.csv"

    input:
    tuple val(accession), val(seq_type), path(fasta)

    output:
    tuple val(accession), val(seq_type), path("${accession}.${seq_type}.k${params.sourmash_ksize}.sig.zip"), emit: sig
    path("${accession}.${seq_type}.sketch.ledger.csv"), emit: ledger

    script:
    """
    set -euo pipefail

    sourmash sketch dna -f "${fasta}" \\
        -p k=${params.sourmash_ksize},scaled=${params.sourmash_scale},abund \\
        --name "${accession}" \\
        -o "${accession}.${seq_type}.k${params.sourmash_ksize}.sig.zip"

    write_ledger_row.sh "${accession}.${seq_type}.sketch.ledger.csv" "${accession}" "${seq_type}" sketch DONE
    """
}
