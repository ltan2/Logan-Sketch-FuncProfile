#!/usr/bin/env python3
"""Randomly sample N accessions from the full manifest, for pipeline benchmarking.

The sample is deterministic for a given (source, n, seed) so a benchmark run can be
reproduced exactly. The output order is itself randomized (not sorted by input order),
so run_benchmark.sh can build nested batches by sampling once at the largest size and
taking a prefix of that list for every smaller size -- see run_benchmark.sh for why.

--max-compressed-gb optionally excludes accessions whose unitigs or contigs .fa.zst
exceeds that size on S3 (per benchmark/query_accession_sizes.py's survey), before
sampling. Without it, a single outlier -- e.g. SRR3507924, a 10.4GB-compressed/45GB-
decompressed accession that stalled a real n=10/50/100 sanity-check sweep for hours --
can land in the sample by pure luck and dominate the whole batch's wall-clock time,
since the analysis stages are single-threaded and their cost scales with input size.
"""
import argparse
import csv
import os
import random


def load_accessions(source_path):
    with open(source_path) as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]


def load_oversized(size_csv, max_bytes):
    """Returns the set of accessions whose unitigs or contigs compressed size (per
    size_csv, from query_accession_sizes.py) exceeds max_bytes. An accession missing
    from size_csv (e.g. it predates the survey, or --source doesn't match what was
    surveyed) is not flagged -- there's no evidence it's an outlier."""
    oversized = set()
    with open(size_csv, newline="") as f:
        for row in csv.DictReader(f):
            raw = row["compressed_bytes"]
            if raw and int(raw) > max_bytes:
                oversized.add(row["accession"])
    return oversized


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", required=True, help="Full accession manifest, one accession per line")
    parser.add_argument("--n", type=int, required=True, help="Number of accessions to sample")
    parser.add_argument("--seed", type=int, required=True, help="RNG seed -- fixes which accessions are picked")
    parser.add_argument("--out", required=True, help="Output path for the sampled accession list")
    parser.add_argument(
        "--max-compressed-gb", type=float, default=None,
        help="Exclude accessions whose unitigs or contigs .fa.zst exceeds this size on S3 "
             "(GB, per --size-csv). Default: no filtering.",
    )
    parser.add_argument(
        "--size-csv",
        default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "accession_size_analysis", "accession_sizes.csv"),
        help="Output of query_accession_sizes.py; only read when --max-compressed-gb is set.",
    )
    args = parser.parse_args()

    accessions = load_accessions(args.source)

    if args.max_compressed_gb is not None:
        max_bytes = int(args.max_compressed_gb * 1e9)
        oversized = load_oversized(args.size_csv, max_bytes)
        before = len(accessions)
        accessions = [a for a in accessions if a not in oversized]
        print(f"Excluded {before - len(accessions)} accessions with a compressed unitigs/contigs "
              f"file over {args.max_compressed_gb:.2f} GB (per {args.size_csv})")

    if args.n > len(accessions):
        raise SystemExit(
            f"Requested {args.n} accessions but {args.source} only has {len(accessions)} "
            "(after size filtering, if any)"
        )

    sample = random.Random(args.seed).sample(accessions, args.n)

    with open(args.out, "w") as f:
        f.write("\n".join(sample) + "\n")

    print(f"Wrote {len(sample)} accessions (seed={args.seed}) to {args.out}")


if __name__ == "__main__":
    main()
