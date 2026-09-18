// Deletes the shared decompressed FASTA once both of its consumers (SOURMASH_SKETCH and
// FUNPROFILER) have finished with it, instead of leaving it in -work-dir until the whole
// run completes (nextflow.config's `cleanup = true` only fires at end-of-run). See
// ../subworkflows/analyze.nf for how this is wired up, and benchmark/run_benchmark.log for
// the "No space left on device" failures at n=5000 this is meant to avoid.

process CLEANUP_FASTA {
    tag "${accession}.${seq_type}"

    input:
    tuple val(accession), val(seq_type), path(fasta)

    script:
    """
    set -euo pipefail

    # 'fasta' is staged here as a symlink into DECOMPRESS_FASTA's original work dir;
    # resolve it and delete the real file so its disk space is actually freed.
    rm -f "\$(readlink -f "${fasta}")"
    """
}
