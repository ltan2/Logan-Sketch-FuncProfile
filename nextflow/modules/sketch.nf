// Build sourmash sketches from the FASTA produced by the shared decompression process:
// a DNA sketch, and a protein sketch of the same sequence translated in six frames. 

process SOURMASH_SKETCH {
    tag "${accession}.${seq_type}"
    // NOTE: path must be a closure here, not a plain interpolated string -- publishDir's
    // path is evaluated once at process-definition time unless wrapped in `{ ... }` to
    // defer it until a task actually runs and `seq_type` (an input variable) is bound.
    publishDir(path: { "${params.outdir}/sketches/${seq_type}" }, mode: 'copy', pattern: "*.sig.zip")
    publishDir "${params.outdir}/ledger", mode: 'copy', pattern: "*.ledger.csv"

    input:
    tuple val(accession), val(seq_type), path(fasta)

    output:
    tuple val(accession), val(seq_type), path("${accession}.${seq_type}.k${params.sourmash_ksize}.sig.zip"), emit: sig
    // Nothing downstream consumes this one; it is published for later protein-space searches.
    tuple val(accession), val(seq_type), path("${accession}.${seq_type}.protein.k${params.sourmash_protein_ksize}.sig.zip"), emit: protein_sig
    path("${accession}.${seq_type}.sketch.ledger.csv"), emit: ledger

    script:
    // `sourmash sketch translate` is the DNA-input counterpart of `sketch protein`: it reads
    // nucleotide sequence and emits an amino-acid sketch. Its k is the amino-acid k-mer size,
    // which is why 11 (not 31) is the right number here -- it has to match the KO collection
    // this sketch is meant to be compared against, e.g.
    //   sourmash prefetch <acc>.<seq>.protein.k11.sig.zip <ko_sig> \
    //       --protein -k 11 --scaled 1000 --threshold-bp 1000 -o prefetch.csv
    """
    set -euo pipefail

    sourmash sketch dna -f "${fasta}" \\
        -p k=${params.sourmash_ksize},scaled=${params.sourmash_scale},abund \\
        --name "${accession}" \\
        -o "${accession}.${seq_type}.k${params.sourmash_ksize}.sig.zip"

    sourmash sketch translate -f "${fasta}" \\
        -p k=${params.sourmash_protein_ksize},scaled=${params.sourmash_protein_scale},abund \\
        --name "${accession}" \\
        -o "${accession}.${seq_type}.protein.k${params.sourmash_protein_ksize}.sig.zip"

    write_ledger_row.sh "${accession}.${seq_type}.sketch.ledger.csv" "${accession}" "${seq_type}" sketch DONE
    """
}
