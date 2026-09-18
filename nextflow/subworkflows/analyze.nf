// Stage 4 (analysis fan-out) from ../design/latest_design-doc.md section 3.4.
// Each validated compressed input is decompressed once. The resulting FASTA feeds both
// funprofiler (for the sequence types in params.funprofiler_seq_types) and sourmash.

include { DECOMPRESS_FASTA } from '../modules/decompress.nf'
include { SOURMASH_SKETCH }  from '../modules/sketch.nf'
include { FUNPROFILER }      from '../modules/funprofiler.nf'
include { CLEANUP_FASTA }    from '../modules/cleanup.nf'

workflow ANALYZE {
    take:
    zst_fasta_ch   // tuple(accession, seq_type, zst_fasta_path)

    main:
    funprofiler_seq_types = params.funprofiler_seq_types.tokenize(',')*.trim()
    decompressed = DECOMPRESS_FASTA(zst_fasta_ch)
    sketch_out   = SOURMASH_SKETCH(decompressed.fasta)
    ko_out       = FUNPROFILER(decompressed.fasta.filter { accession, seq_type, fasta -> seq_type in funprofiler_seq_types })

    // The decompressed FASTA (often multi-GB per accession) is only needed by the two
    // processes above. Delete it as soon as both have consumed it, keyed on
    // accession+seq_type, rather than waiting for end-of-run cleanup. If either consumer
    // fails (errorStrategy 'ignore' in the benchmark profile), its join key never
    // completes and that one FASTA is left for end-of-run cleanup instead -- a much
    // smaller leftover than "every FASTA, for the whole run".
    fasta_keyed = decompressed.fasta.map { accession, seq_type, fasta ->
        tuple("${accession}.${seq_type}", accession, seq_type, fasta)
    }
    sketch_done = sketch_out.sig.map { accession, seq_type, sig -> tuple("${accession}.${seq_type}", true) }
    // A seq type FUNPROFILER skips has no ko_csv to wait for, so mark it done as soon as it is
    // decompressed -- otherwise its join key never completes and the FASTA survives until
    // end-of-run cleanup. The sketch is still required before the FASTA can be deleted.
    ko_done = ko_out.ko_csv.map { accession, seq_type, csv -> tuple("${accession}.${seq_type}", true) }
        .mix(decompressed.fasta
            .filter { accession, seq_type, fasta -> !(seq_type in funprofiler_seq_types) }
            .map { accession, seq_type, fasta -> tuple("${accession}.${seq_type}", true) })

    CLEANUP_FASTA(
        fasta_keyed.join(sketch_done).join(ko_done)
            .map { key, accession, seq_type, fasta, _s, _k -> tuple(accession, seq_type, fasta) }
    )

    emit:
    sig    = sketch_out.sig
    ko     = ko_out.ko_csv
    ledger = sketch_out.ledger.mix(ko_out.ledger)
}
