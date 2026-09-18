#!/usr/bin/env python3
"""Compare the size distribution of a benchmark sample against the full accession population.

benchmark/plot_accession_sizes.py answers "how big is each accession?" across all ~1.2M of
them. This answers the follow-up: "is the batch I actually benchmarked representative of
that?" -- benchmark/sample_accessions.py draws uniformly at random, so a 1,000-accession
batch should track the population, but a random draw of 1,000 from a distribution whose top
percentiles span two orders of magnitude can easily miss (or over-sample) the tail that
dominates wall-clock. This plots both so that's visible instead of assumed.

Reads:
  --size-csv   accession,seq_type,compressed_bytes -- the population (query_accession_sizes.py);
               a blank size means Logan doesn't publish that file for that accession
  --sample     one or more accession lists (one accession per line), e.g. the
               accessions_<N>.txt files run_benchmark.sh writes into its BENCH_DIR

Writes a histogram + ECDF figure and a stats CSV, and prints the same stats.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

POPULATION_COLOR = "#9aa0a6"
SAMPLE_COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf"]


def load_sizes(path):
    df = pd.read_csv(path)
    df["compressed_bytes"] = pd.to_numeric(df["compressed_bytes"], errors="coerce")
    return df


def load_sample(path):
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def describe(sizes, n_listed, label, seq_type):
    """One row of summary stats for a set of accessions of one seq_type."""
    gb = sizes.dropna() / 1e9
    row = {
        "set": label,
        "seq_type": seq_type,
        "accessions": n_listed,
        "published": len(gb),
        "published_pct": 100 * len(gb) / n_listed if n_listed else np.nan,
    }
    if len(gb):
        row.update({
            "mean_gb": gb.mean(), "median_gb": gb.median(),
            "p90_gb": gb.quantile(.90), "p95_gb": gb.quantile(.95), "p99_gb": gb.quantile(.99),
            "max_gb": gb.max(),
            "over_2gb": int((gb > 2).sum()), "over_2gb_pct": 100 * (gb > 2).mean(),
            "total_tb": gb.sum() / 1000,
        })
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--size-csv", default="accession_size_analysis/accession_sizes.csv")
    parser.add_argument("--sample", nargs="+", required=True,
                        help="accession list file(s); labelled by filename stem")
    parser.add_argument("--out", default="accession_size_analysis/sample_distribution.png")
    args = parser.parse_args()

    pop = load_sizes(args.size_csv)
    seq_types = sorted(pop["seq_type"].unique())
    samples = {os.path.splitext(os.path.basename(p))[0]: load_sample(p) for p in args.sample}

    rows = []
    fig, axes = plt.subplots(3, len(seq_types), figsize=(7 * len(seq_types), 13), squeeze=False)

    for col, seq_type in enumerate(seq_types):
        sub = pop[pop["seq_type"] == seq_type]
        pop_sizes = sub["compressed_bytes"]
        by_accession = sub.set_index("accession")["compressed_bytes"]
        rows.append(describe(pop_sizes, len(sub), "population", seq_type))

        pop_gb = pop_sizes.dropna() / 1e9
        # Log-spaced bins: sizes span ~8 orders of magnitude, so linear bins would put
        # everything in the first bucket. Bars are weighted to the fraction of the set in each
        # bin rather than a density -- a density divides by bin width, and log bins down at
        # 1e-6 GB are so narrow that they'd tower over the 0.01-10 GB range that drives runtime.
        bins = np.logspace(np.log10(max(pop_gb.min(), 1e-7)), np.log10(pop_gb.max()), 50)

        ax_hist, ax_ecdf, ax_tail = axes[0][col], axes[1][col], axes[2][col]
        ax_hist.hist(pop_gb, bins=bins, weights=np.full(len(pop_gb), 1 / len(pop_gb)),
                     color=POPULATION_COLOR, label=f"population (n={len(pop_gb):,})")
        x = np.sort(pop_gb)
        frac_below = np.arange(1, len(x) + 1) / len(x)
        ax_ecdf.step(x, frac_below, where="post", color=POPULATION_COLOR, lw=2, label="population")
        ax_tail.step(x, 1 - frac_below, where="post", color=POPULATION_COLOR, lw=2, label="population")

        for i, (label, accessions) in enumerate(samples.items()):
            color = SAMPLE_COLORS[i % len(SAMPLE_COLORS)]
            sizes = by_accession.reindex(accessions)
            rows.append(describe(sizes, len(accessions), label, seq_type))
            gb = sizes.dropna() / 1e9
            if gb.empty:
                continue
            ax_hist.hist(gb, bins=bins, weights=np.full(len(gb), 1 / len(gb)),
                         histtype="step", lw=1.8, color=color, label=f"{label} (n={len(gb):,})")
            xs = np.sort(gb)
            below = np.arange(1, len(xs) + 1) / len(xs)
            ax_ecdf.step(xs, below, where="post", lw=1.6, color=color, label=label)
            ax_tail.step(xs, 1 - below, where="post", lw=1.6, color=color, label=label)
            # Two-sample KS: how far the sample's ECDF strays from the population's. A large
            # p-value means the draw is consistent with being uniform from the population.
            ks = stats.ks_2samp(gb, pop_gb)
            rows[-1]["ks_stat"] = ks.statistic
            rows[-1]["ks_pvalue"] = ks.pvalue

        for ax in (ax_hist, ax_ecdf, ax_tail):
            ax.set_xscale("log")
            ax.set_xlabel("compressed size (GB, log scale)")
            ax.grid(alpha=.3)
            ax.legend(fontsize=8)
        ax_hist.set_ylabel("fraction of accessions in bin")
        ax_hist.set_title(f"{seq_type}: size distribution, sample vs population")
        ax_ecdf.set_ylabel("fraction of accessions at or below")
        ax_ecdf.set_title(f"{seq_type}: cumulative distribution")
        # The tail on log-log: a handful of multi-GB accessions decide a batch's wall-clock
        # (one 10 GB accession stalled an early run for 3+ hours), and they are invisible
        # anywhere else on this figure.
        ax_tail.set_yscale("log")
        ax_tail.set_ylabel("fraction of accessions above")
        ax_tail.set_title(f"{seq_type}: tail (fraction larger than x)")
        ax_tail.axvline(2, color="black", ls=":", lw=1.2)
        ax_tail.text(2, ax_tail.get_ylim()[1], " 2 GB", va="top", fontsize=8)

    fig.suptitle("Benchmark samples vs the full accession population", fontsize=13)
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    fig.savefig(args.out, dpi=150)

    table = pd.DataFrame(rows)
    csv_path = os.path.splitext(args.out)[0] + ".csv"
    table.to_csv(csv_path, index=False)

    with pd.option_context("display.width", 200, "display.max_columns", 30,
                           "display.float_format", lambda v: f"{v:.3f}"):
        print(table.to_string(index=False))
    print(f"\nWrote {args.out}\nWrote {csv_path}")


if __name__ == "__main__":
    main()
