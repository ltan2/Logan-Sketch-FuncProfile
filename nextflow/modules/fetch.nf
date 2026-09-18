// Stage 3 (fetch + empty gate) from ../../design/latest_design-doc.md section 3.3.
// Commands are taken directly from logan-subworkflow.sh and reference-logan-analysis.sh.

process FETCH_LOGAN {
    tag "$accession"

    // Transient S3/network failures shouldn't fail a whole run; other stages keep the
    // default 'terminate' behaviour so a real bug surfaces instead of being retried away.
    // Once retries are exhausted, params.fetch_exhausted_strategy decides whether that
    // (now clearly non-transient) failure aborts the run or is skipped like any other
    // dropped accession -- see the 'benchmark' profile in nextflow.config.
    errorStrategy { task.attempt <= task.maxRetries ? 'retry' : params.fetch_exhausted_strategy }
    maxRetries 3

    input:
    val accession

    output:
    tuple val(accession), val('unitigs'), path("${accession}.unitigs.fa.zst"), env('U_STATUS'), emit: unitigs
    tuple val(accession), val('contigs'), path("${accession}.contigs.fa.zst"), env('C_STATUS'), emit: contigs

    script:
    """
    set -euo pipefail

    # Logan doesn't always publish both files for a given accession (e.g. an accession
    # with no assembled contigs has no c/ key at all, not just a small one). That's a
    # permanent 404, not the transient S3 error this process's errorStrategy retries --
    # so treat "key does not exist" as an empty download (caught by the size gate below)
    # instead of a retryable failure. Any other `aws s3 cp` error still fails the task.
    fetch_or_empty() {
        local uri="\$1" out="\$2"
        if aws s3 cp "\$uri" "\$out" --no-sign-request 2>fetch.err; then
            return 0
        elif grep -q "does not exist" fetch.err; then
            : > "\$out"
            return 0
        else
            cat fetch.err >&2
            return 1
        fi
    }

    fetch_or_empty "s3://logan-pub/u/${accession}/${accession}.unitigs.fa.zst" "${accession}.unitigs.fa.zst"
    fetch_or_empty "s3://logan-pub/c/${accession}/${accession}.contigs.fa.zst" "${accession}.contigs.fa.zst"

    # Logan can publish near-empty unitig/contig files for accessions with no assembled
    # sequence, and zstd hangs decompressing such files (reference-logan-analysis.sh).
    # Gate on the *compressed* size before ever calling zstd -d, and record the outcome
    # explicitly (OK / SKIPPED_EMPTY) rather than silently dropping the accession.
    u_size=\$(stat -c%s "${accession}.unitigs.fa.zst")
    c_size=\$(stat -c%s "${accession}.contigs.fa.zst")

    # NOTE: written as if/fi rather than `[ ... ] && VAR=...` because under `set -e` a
    # false test on the left of `&&` exits the script instead of just skipping the RHS.
    U_STATUS=OK
    if [ "\$u_size" -lt ${params.min_compressed_bytes} ]; then
        U_STATUS=SKIPPED_EMPTY
    fi

    C_STATUS=OK
    if [ "\$c_size" -lt ${params.min_compressed_bytes} ]; then
        C_STATUS=SKIPPED_EMPTY
    fi
    """
}
