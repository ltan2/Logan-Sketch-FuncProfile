#!/usr/bin/env python3
"""Query the compressed (.fa.zst) object size of every accession's unitigs/contigs file
on S3, without downloading them, via HeadObject. Used to characterize the size
distribution of the accession pool (benchmark/plot_accession_sizes.py plots it) -- e.g.
a 45GB-decompressed outlier accession was found to stall a benchmark run for hours; its
compressed S3 size (10.4GB) would have flagged it in seconds, with no download needed.

Usage:
  python3 benchmark/query_accession_sizes.py [--source PATH] [--out PATH] [--workers N]
"""
import argparse
import csv
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import boto3
from botocore import UNSIGNED
from botocore.config import Config
from botocore.exceptions import ClientError

BUCKET = "logan-pub"
_PREFIX = {"unitigs": "u", "contigs": "c"}


def load_accessions(path):
    with open(path) as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]


def head_size(s3, accession, seq_type, retries=3):
    """Returns the compressed size in bytes, or None if genuinely not published (404).
    Retries transient errors (throttling, connection resets) a few times before giving
    up, so a sustained multi-hour run doesn't misreport those as "not published"."""
    key = f"{_PREFIX[seq_type]}/{accession}/{accession}.{seq_type}.fa.zst"
    for attempt in range(retries):
        try:
            resp = s3.head_object(Bucket=BUCKET, Key=key)
            return resp["ContentLength"]
        except ClientError as e:
            if e.response.get("Error", {}).get("Code") in ("404", "NoSuchKey"):
                return None
            time.sleep(0.5 * (attempt + 1))
        except Exception:
            time.sleep(0.5 * (attempt + 1))
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--source",
        default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "database", "wgs_metagenome_accessions.txt"),
        help="Accession manifest, one per line (default: database/wgs_metagenome_accessions.txt)",
    )
    parser.add_argument(
        "--out",
        default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "accession_size_analysis", "accession_sizes.csv"),
        # On this machine /scratch is meant for bulk pipeline I/O (run_benchmark.sh's
        # multi-GB work dirs) -- this dataset is small (~85MB) and durability matters more
        # than /scratch's capacity, so it defaults under the repo (gitignored) on /home instead.
        help="Output CSV path (default: accession_size_analysis/accession_sizes.csv)",
    )
    parser.add_argument("--workers", type=int, default=128, help="Concurrent HeadObject requests (default: 128)")
    args = parser.parse_args()

    accessions = load_accessions(args.source)
    all_jobs = [(a, seq_type) for a in accessions for seq_type in ("unitigs", "contigs")]

    # Resume support: an (accession, seq_type) pair already on a prior partial run's
    # output is skipped rather than re-queried. The csv module can leave a truncated
    # final row if a previous run was killed mid-write; skip any row that doesn't parse
    # as exactly 3 fields instead of treating that as fatal.
    already_done = set()
    file_exists = os.path.exists(args.out)
    if file_exists:
        with open(args.out, newline="") as f:
            reader = csv.reader(f)
            next(reader, None)  # header
            for row in reader:
                if len(row) == 3:
                    already_done.add((row[0], row[1]))
    jobs = [j for j in all_jobs if j not in already_done]
    if already_done:
        print(f"Resuming: {len(already_done)} of {len(all_jobs)} objects already recorded in {args.out}", flush=True)
    print(f"Querying {len(jobs)} remaining objects ({len(accessions)} accessions x 2 seq types) with {args.workers} workers...", flush=True)

    s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED, max_pool_connections=args.workers * 2))
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    t0 = time.time()
    done = 0
    with open(args.out, "a" if file_exists else "w", newline="") as f, ThreadPoolExecutor(max_workers=args.workers) as ex:
        writer = csv.writer(f)
        if not file_exists:
            writer.writerow(["accession", "seq_type", "compressed_bytes"])
        futures = {ex.submit(head_size, s3, a, seq_type): (a, seq_type) for a, seq_type in jobs}
        for future in as_completed(futures):
            accession, seq_type = futures[future]
            size = future.result()
            writer.writerow([accession, seq_type, size if size is not None else ""])
            done += 1
            if done % 50_000 == 0 or done == len(jobs):
                elapsed = time.time() - t0
                rate = done / elapsed
                eta_min = (len(jobs) - done) / rate / 60 if rate else float("nan")
                print(f"{done}/{len(jobs)} done, {rate:.0f} req/s, ETA {eta_min:.1f} min", flush=True)

    print(f"Wrote {args.out} in {(time.time() - t0) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
