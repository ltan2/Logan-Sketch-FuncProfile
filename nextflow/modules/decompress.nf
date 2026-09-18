// Decompress each validated Logan assembly once, then share the resulting FASTA
// with all downstream analyses.

process DECOMPRESS_FASTA {
    tag "${accession}.${seq_type}"

    input:
    tuple val(accession), val(seq_type), path(zst_fasta)

    output:
    tuple val(accession), val(seq_type), path("${accession}.${seq_type}.fa"), emit: fasta

    script:
    """
    set -euo pipefail

    # -f is required because Nextflow stages input files as symlinks.
    zstd -d -f "${zst_fasta}" -o "${accession}.${seq_type}.fa"

    # FETCH_LOGAN's compressed output is only consumed here. Delete its real file now
    # (resolving through the symlink Nextflow staged it as) instead of leaving it in
    # -work-dir until the whole run finishes -- see modules/cleanup.nf for the equivalent
    # treatment of this process's own FASTA output, and benchmark/run_benchmark.log for
    # the "No space left on device" failures this avoids.
    rm -f "\$(readlink -f "${zst_fasta}")"
    """
}
