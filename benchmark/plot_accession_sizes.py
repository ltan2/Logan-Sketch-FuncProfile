#!/usr/bin/env python3
"""Summarize the accession-size survey produced by benchmark/query_accession_sizes.py.

Reads:
  <csv> -- accession,seq_type,compressed_bytes (one row per accession x seq_type;
           compressed_bytes blank means Logan doesn't publish that file for that
           accession -- see fetch.nf)

Prints summary statistics (count, missing/published rate, mean/median/p90/p95/p99/max)
per seq_type, and writes a log-scale histogram to <out>.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def print_summary(df):
    for seq_type, sub in df.groupby("seq_type"):
        total = len(sub)
        sizes = sub["compressed_bytes"].dropna()
        missing = total - len(sizes)
        print(f"\n=== {seq_type} ===")
        print(f"accessions checked:  {total:,}")
        print(f"published on S3:     {len(sizes):,} ({100 * len(sizes) / total:.1f}%)")
        print(f"not published (404): {missing:,} ({100 * missing / total:.1f}%)")
        if sizes.empty:
            continue
        gb = sizes / 1e9
        print(f"compressed size (GB): mean={gb.mean():.3f}  median={gb.median():.3f}  "
              f"p90={gb.quantile(.90):.3f}  p95={gb.quantile(.95):.3f}  p99={gb.quantile(.99):.3f}  "
              f"max={gb.max():.3f}  min={gb.min() * 1000:.3f} MB")


def plot_histogram(df, out_path):
    seq_types = sorted(df["seq_type"].unique())
    fig, axes = plt.subplots(1, len(seq_types), figsize=(6 * len(seq_types), 5), sharey=True)
    if len(seq_types) == 1:
        axes = [axes]

    for ax, seq_type in zip(axes, seq_types):
        sizes = df.loc[df["seq_type"] == seq_type, "compressed_bytes"].dropna()
        if sizes.empty:
            continue
        bins = np.logspace(np.log10(max(sizes.min(), 1)), np.log10(sizes.max()), 60)
        ax.hist(sizes, bins=bins, color="tab:blue", alpha=0.75)
        ax.set_xscale("log")
        for q, style in [(0.50, "-"), (0.90, "--"), (0.95, ":"), (0.99, "-.")]:
            val = sizes.quantile(q)
            ax.axvline(val, color="tab:red", linestyle=style, linewidth=1.2, label=f"p{int(q * 100)} = {val / 1e9:.2f} GB")
        ax.set_title(f"{seq_type} (n={len(sizes):,} published)")
        ax.set_xlabel("Compressed size on S3 (bytes, log scale)")
        ax.legend(fontsize=8)
        ax.grid(True, which="both", linestyle="--", alpha=0.3)

    axes[0].set_ylabel("Number of accessions")
    fig.suptitle("Logan accession .fa.zst compressed-size distribution")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    _survey_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "accession_size_analysis")
    parser.add_argument(
        "--csv",
        default=os.path.join(_survey_dir, "accession_sizes.csv"),
        help="Input CSV from query_accession_sizes.py",
    )
    parser.add_argument(
        "--out",
        default=os.path.join(_survey_dir, "accession_size_distribution.png"),
        help="Output histogram PNG path",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.csv, dtype={"accession": str, "seq_type": str})
    df["compressed_bytes"] = pd.to_numeric(df["compressed_bytes"], errors="coerce")

    print_summary(df)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    plot_histogram(df, args.out)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
