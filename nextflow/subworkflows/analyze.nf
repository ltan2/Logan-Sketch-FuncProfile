// Stage 4 (analysis fan-out) from ../design/latest_design-doc.md section 3.4.
// Each validated compressed input is decompressed once. The resulting FASTA feeds both
// funprofiler (for the sequence types in params.funprofiler_seq_types) and sourmash.

include { DECOMPRESS_FASTA } from '../modules/decompress.nf'
include { CLEANUP_FASTA }    from '../modules/cleanup.nf'
// Aliased into a second copy of each analysis process so the big ones can be given their own
// maxForks. See the size branch below for why -- a process's maxForks applies to the whole
// process, so "cap only the big tasks" has to be expressed as a separate process.
include { SOURMASH_SKETCH; SOURMASH_SKETCH as SOURMASH_SKETCH_BIG } from '../modules/sketch.nf'
include { FUNPROFILER;     FUNPROFILER as FUNPROFILER_BIG }         from '../modules/funprofiler.nf'

workflow ANALYZE {
    take:
    zst_fasta_ch   // tuple(accession, seq_type, zst_fasta_path)

    main:
    funprofiler_seq_types = params.funprofiler_seq_types.tokenize(',')*.trim()
    decompressed = DECOMPRESS_FASTA(zst_fasta_ch)

    // Route the few huge inputs down a parallel, capped lane.
    //
    // Task size in this manifest spans ~734x (median 52 MB compressed, largest 38.5 GB), and
    // the top 5% of inputs carry ~49% of the bytes. Slot occupancy is count x duration, not
    // count, so a rare multi-hour task holds its slot far longer than the 64% of inputs that
    // clear in under a minute -- and the executor pool converts to big tasks by accumulation.
    // Measured 10 h into a 24,718-accession shard: 307 of 383 slots were held by tasks over an
    // hour old, 255 of them over three, while completions had fallen to ~33 accessions/hour.
    //
    // Capping the big lane is NOT the same as de-prioritising it. The makespan tail is bounded
    // by the longest single task (~9.2 h here), which is satisfied as long as the largest
    // inputs *start* early -- they still do, since the cap only queues big task N+1, never the
    // first N. What the cap buys is that the remaining slots stay free for short work instead
    // of being swallowed. Size the cap by core-time share, not intuition: too small and the
    // big lane becomes the tail it was meant to prevent (at ~2,800 big task-hours per shard,
    // 60 concurrent => ~47 h, against a ~37 h work-bound; ~190 => ~15 h, comfortably inside).
    by_size = decompressed.fasta.branch { accession, seq_type, fasta ->
        big:   fasta.size() >= params.big_fasta_bytes
        small: true
    }

    sketch_small = SOURMASH_SKETCH(by_size.small)
    sketch_big   = SOURMASH_SKETCH_BIG(by_size.big)
    ko_small = FUNPROFILER(by_size.small.filter { accession, seq_type, fasta -> seq_type in funprofiler_seq_types })
    ko_big   = FUNPROFILER_BIG(by_size.big.filter { accession, seq_type, fasta -> seq_type in funprofiler_seq_types })

    // Recombined immediately: nothing downstream (the CLEANUP join, the emitted channels)
    // should have to know which lane a task took.
    sketch_out = [sig:    sketch_small.sig.mix(sketch_big.sig),
                  ledger: sketch_small.ledger.mix(sketch_big.ledger)]
    ko_out     = [ko_csv: ko_small.ko_csv.mix(ko_big.ko_csv),
                  ledger: ko_small.ledger.mix(ko_big.ledger)]

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
