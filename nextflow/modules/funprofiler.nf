// Functional profiling against KEGG Orthology sketches. The shared decompression
// process supplies the FASTA used by this process and SOURMASH_SKETCH.

process FUNPROFILER {
    tag "${accession}.${seq_type}"
    // path must be a closure -- see the NOTE in modules/sketch.nf.
    publishDir(path: { "${params.outdir}/ko_profiles/${seq_type}" }, mode: 'copy', pattern: "*_ko_profiles.csv")
    publishDir(path: { "${params.outdir}/ko_profiles/${seq_type}" }, mode: 'copy', pattern: "*_prefetch_out.csv")
    publishDir "${params.outdir}/ledger", mode: 'copy', pattern: "*.ledger.csv"

    input:
    tuple val(accession), val(seq_type), path(fasta)

    output:
    tuple val(accession), val(seq_type), path("${accession}.${seq_type}_ko_profiles.csv"), emit: ko_csv
    path("${accession}.${seq_type}_prefetch_out.csv"), emit: prefetch
    path("${accession}.${seq_type}.funprofiler.ledger.csv"), emit: ledger

    script:
    """
    set -euo pipefail

    funcprofiler "${fasta}" "${params.ko_sig}" \\
        ${params.funprofiler_ksize} ${params.funprofiler_scale} \\
        "${accession}.${seq_type}_ko_profiles.csv" \\
        -p "${accession}.${seq_type}_prefetch_out.csv"

    # For a near-empty contig/unitig file, the underlying `sourmash prefetch` call can
    # bail out with "no query hashes!? exiting" before ever writing its output csv, even
    # though funcprofiler treats that as a valid zero-KO-match result and still writes
    # ko_profiles.csv. Make sure the declared prefetch output exists either way.
    : >> "${accession}.${seq_type}_prefetch_out.csv"

    write_ledger_row.sh "${accession}.${seq_type}.funprofiler.ledger.csv" "${accession}" "${seq_type}" funprofiler DONE
    """
}
