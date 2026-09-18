#!/usr/bin/env python3
"""Drop accessions that a previous run already finished, so a resumed run only does what's left.

Why this exists instead of `nextflow -resume`: DECOMPRESS_FASTA deletes FETCH_LOGAN's .fa.zst
and CLEANUP_FASTA deletes the decompressed FASTA (see nextflow/modules/cleanup.nf) as soon as
they are consumed, to keep a large run inside its disk budget. Nextflow invalidates a cached
task whose declared output files no longer exist, so on `-resume` FETCH_LOGAN is re-run,
which re-times-stamps its outputs, which cascades: every downstream task re-runs too. A
measured re-run of the 5-accession test profile with `-resume` reported `cached=0`, i.e. it
re-downloaded and re-analyzed everything.

Published results, on the other hand, are durable. So resume is done at the accession level:
an accession whose results are all on disk is removed from the manifest before launch.

An accession counts as complete when, for every sequence type, either
  * the fetch ledger recorded a skip                     <outdir>/ledger/<acc>.<seq>.fetch.ledger.csv
    (Logan publishes no such file for it, or it is below params.min_compressed_bytes), or
  * its sketch                                           <outdir>/sketches/<seq>/<acc>.<seq>.k<K>.sig.zip
    and, for the sequence types FUNPROFILER runs on, its KO profile
                                                         <outdir>/ko_profiles/<seq>/<acc>.<seq>_ko_profiles.csv
    are both published.

The protein sketch (`<acc>.<seq>.protein.k11.sig.zip`) is deliberately not part of this test:
it is published by the same task as the DNA sketch, so one implies the other. The exception is
a results tree that was built before protein sketching was turned on -- those accessions count
as complete and are not re-run, which is the intended behaviour, since re-making their protein
sketch means fetching and decompressing the whole assembly again. Run them as their own
manifest if you want the sketch added.

Only these exact paths are stat'ed -- the results directories of a full run hold millions of
files, so nothing here lists a directory or reads the ledger CSVs.

Usage:
  filter_completed.py --accessions manifest.txt --outdir results [--out remaining.txt]

Exit status is 0 whether or not anything remains; the caller decides what an empty remainder
means (nextflow/run_full.sh treats it as "this shard is already done").
"""
import argparse
import os
import sys


def parse_accessions(path):
    """Same convention as main.nf's channel: trimmed lines, no blanks, no '#' comments."""
    with open(path) as handle:
        for line in handle:
            accession = line.strip()
            if accession and not accession.startswith("#"):
                yield accession


def seq_type_complete(outdir, accession, seq_type, ksize, run_funprofiler):
    skipped = os.path.join(outdir, "ledger", f"{accession}.{seq_type}.fetch.ledger.csv")
    if os.path.exists(skipped):
        return True
    sketch = os.path.join(outdir, "sketches", seq_type, f"{accession}.{seq_type}.k{ksize}.sig.zip")
    if not os.path.exists(sketch):
        return False
    if not run_funprofiler:
        return True
    ko = os.path.join(outdir, "ko_profiles", seq_type, f"{accession}.{seq_type}_ko_profiles.csv")
    return os.path.exists(ko)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--accessions", required=True, help="Manifest to filter")
    parser.add_argument("--outdir", required=True, help="Published results directory (params.outdir)")
    parser.add_argument("--out", default="-", help="Where to write the remaining accessions (default: stdout)")
    parser.add_argument("--seq-types", default="unitigs,contigs",
                        help="Sequence types the pipeline produces (default: unitigs,contigs)")
    parser.add_argument("--funprofiler-seq-types", default="unitigs,contigs",
                        help="params.funprofiler_seq_types of the run being resumed (default: unitigs,contigs)")
    parser.add_argument("--ksize", type=int, default=31, help="params.sourmash_ksize of the run (default: 31)")
    parser.add_argument("--count-only", action="store_true",
                        help="Only print 'complete<TAB>remaining' to stdout; write no manifest")
    args = parser.parse_args()

    seq_types = [s.strip() for s in args.seq_types.split(",") if s.strip()]
    funprofiler_types = {s.strip() for s in args.funprofiler_seq_types.split(",") if s.strip()}

    remaining, complete = [], 0
    for accession in parse_accessions(args.accessions):
        if all(seq_type_complete(args.outdir, accession, seq_type, args.ksize, seq_type in funprofiler_types)
               for seq_type in seq_types):
            complete += 1
        else:
            remaining.append(accession)

    if args.count_only:
        print(f"{complete}\t{len(remaining)}")
        return

    handle = sys.stdout if args.out == "-" else open(args.out, "w")
    try:
        for accession in remaining:
            print(accession, file=handle)
    finally:
        if handle is not sys.stdout:
            handle.close()

    print(f"{complete} of {complete + len(remaining)} accessions already complete; "
          f"{len(remaining)} left to run", file=sys.stderr)


if __name__ == "__main__":
    main()
