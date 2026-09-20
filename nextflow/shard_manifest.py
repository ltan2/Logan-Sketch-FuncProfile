#!/usr/bin/env python3
"""Split an accession manifest into the shards nextflow/run_full.sh runs one at a time,
balancing their total work and putting the heaviest accessions first inside each shard.

Why not just `split -l`: a shard's wall clock is its longest accession, not its average one.
Measured on a 1,000-accession benchmark, 90% of accessions finished in 28 minutes and the
remaining 10% took another 68 -- 71% of the run -- because a handful of multi-GB accessions
happened to be dispatched late. FUNPROFILER on the worst one ran 88.8 minutes pinned to a
single core, with the machine otherwise idle.

Two things fix that, both done here when --size-csv is given:

  * balance   accessions are dealt to shards largest-first, each going to the shard with the
              least work so far (greedy LPT bin packing), so no shard collects a wildly
              heavier set than its neighbours;
  * order     within a shard, the heavy accessions are spread through the list at a fixed
              stride (every --big-stride-th position is the next-heaviest, the gaps are filled
              from the light end), so the long-running tasks are submitted early and overlap
              with the small ones instead of stranding the machine at the end.

Why a stride rather than plain descending order: the decompressed FASTA of every accession
being analysed sits on disk until both its sketch and its KO profile are done (modules/
cleanup.nf), and a multi-GB accession occupies that disk for hours. With ~350 analysis tasks
running at once, strict descending order puts a shard's 350 largest accessions on disk
simultaneously -- measured at 2.0 TB compressed, ~7.2 TB decompressed at the observed 3.6x
ratio, on every shard. At the default stride the same first wave holds ~1.4 TB, while the ten
heaviest accessions still start within it (positions 0, 12, 24, ... of 350).

Size comes from benchmark/query_accession_sizes.py's CSV (accession,seq_type,compressed_bytes);
an accession's weight is the sum over its sequence types. Compressed size is a good proxy:
every downstream stage's cost scales with how much sequence it has to read. Accessions missing
from the CSV get the median weight, so they are spread evenly rather than all landing together.

Without --size-csv this degrades to plain sequential chunks, i.e. exactly what `split -l` did.

Usage:
  shard_manifest.py --accessions manifest.txt --out-dir shards/ --shard-size 10000 \
      [--size-csv accession_sizes.csv] [--prefix shard_]

Deterministic: the same inputs always produce the same shards, which is what lets a resumed
run trust the shard numbering recorded in its state directory.
"""
import argparse
import csv
import heapq
import os
import statistics
import sys


def parse_accessions(path):
    """Same convention as main.nf's channel: trimmed lines, no blanks, no '#' comments."""
    with open(path) as handle:
        seen, accessions = set(), []
        for line in handle:
            accession = line.strip()
            if accession and not accession.startswith("#") and accession not in seen:
                seen.add(accession)
                accessions.append(accession)
        return accessions


def load_weights(size_csv, accessions):
    """accession -> total compressed bytes across its sequence types, for the accessions we
    care about. Missing ones are filled with the median so they sort into the middle."""
    wanted = set(accessions)
    weights = {}
    with open(size_csv, newline="") as handle:
        for row in csv.DictReader(handle):
            accession = row.get("accession")
            if accession in wanted:
                try:
                    weights[accession] = weights.get(accession, 0) + int(row["compressed_bytes"] or 0)
                except ValueError:
                    continue
    if not weights:
        raise SystemExit(f"No accessions from the manifest found in {size_csv}")
    median = int(statistics.median(weights.values()))
    missing = 0
    for accession in accessions:
        if accession not in weights:
            weights[accession] = median
            missing += 1
    print(f"sizes: {len(accessions) - missing:,} accessions from {size_csv}, "
          f"{missing:,} missing (given the median, {median / 1e9:.2f} GB)", file=sys.stderr)
    return weights


def stride_order(ordered, stride):
    """Take a heaviest-first list and lay it out so every stride-th position is the next
    heaviest and the gaps are filled from the light end. Keeps the long poles early without
    making their FASTAs co-resident on disk (see the module docstring)."""
    if stride <= 1:
        return list(ordered)
    out, i, j = [], 0, len(ordered) - 1
    while i <= j:
        out.append(ordered[i])
        i += 1
        for _ in range(stride - 1):
            if i > j:
                break
            out.append(ordered[j])
            j -= 1
    return out


def build_shards(accessions, shard_size, weights, stride=12):
    """Greedy LPT: walk the accessions heaviest-first, and put each one in the least-loaded
    shard that still has room. Then lay each shard out with stride_order()."""
    n_shards = max(1, -(-len(accessions) // shard_size))
    if weights is None:
        return [accessions[i:i + shard_size] for i in range(0, len(accessions), shard_size)]

    ordered = sorted(accessions, key=lambda a: (-weights[a], a))
    shards = [[] for _ in range(n_shards)]
    # (total weight, count, index) -- only shards with room stay in the heap.
    heap = [(0, 0, i) for i in range(n_shards)]
    heapq.heapify(heap)
    for accession in ordered:
        total, count, index = heapq.heappop(heap)
        shards[index].append(accession)
        if count + 1 < shard_size:
            heapq.heappush(heap, (total + weights[accession], count + 1, index))
    for index, shard in enumerate(shards):
        shard.sort(key=lambda a: (-weights[a], a))
        shards[index] = stride_order(shard, stride)
    return shards


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--accessions", required=True, help="Manifest to split")
    parser.add_argument("--out-dir", required=True, help="Directory to write the shard files into")
    parser.add_argument("--shard-size", type=int, default=10000, help="Accessions per shard (default: 10000)")
    parser.add_argument("--size-csv", default=None,
                        help="benchmark/query_accession_sizes.py output; enables size balancing and ordering")
    parser.add_argument("--big-stride", type=int, default=12,
                        help="Spacing between heavy accessions inside a shard (default: 12). "
                             "1 means strict heaviest-first, which maximises disk residency.")
    parser.add_argument("--prefix", default="shard_", help="Shard filename prefix (default: shard_)")
    args = parser.parse_args()

    accessions = parse_accessions(args.accessions)
    if not accessions:
        raise SystemExit(f"No accessions in {args.accessions}")

    weights = load_weights(args.size_csv, accessions) if args.size_csv else None
    shards = build_shards(accessions, args.shard_size, weights, args.big_stride)

    os.makedirs(args.out_dir, exist_ok=True)
    for i, shard in enumerate(shards):
        with open(os.path.join(args.out_dir, f"{args.prefix}{i:05d}.txt"), "w") as handle:
            handle.write("".join(f"{a}\n" for a in shard))

    print(f"wrote {len(shards)} shards of up to {args.shard_size} accessions "
          f"({len(accessions):,} total) into {args.out_dir}", file=sys.stderr)
    if weights:
        totals = [sum(weights[a] for a in shard) / 1e12 for shard in shards]
        heaviest = [max(weights[a] for a in shard) / 1e9 for shard in shards]
        print(f"shard weight (compressed TB): min {min(totals):.2f}, max {max(totals):.2f} "
              f"-- spread {100 * (max(totals) - min(totals)) / max(totals):.1f}%", file=sys.stderr)
        print(f"largest accession in a shard (compressed GB): min {min(heaviest):.1f}, "
              f"max {max(heaviest):.1f}; each is that shard's first line", file=sys.stderr)
        wave = min(350, args.shard_size)
        first = [sum(weights[a] for a in shard[:wave]) / 1e12 for shard in shards]
        print(f"first {wave} accessions of a shard (what is on disk at once early on): "
              f"{max(first):.2f} TB compressed, ~{max(first) * 3.6:.1f} TB decompressed",
              file=sys.stderr)


if __name__ == "__main__":
    main()
