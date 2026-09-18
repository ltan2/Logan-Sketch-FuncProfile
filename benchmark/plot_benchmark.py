#!/usr/bin/env python3
"""Plot benchmark results produced by run_benchmark.sh.

Reads:
  <bench-dir>/wall_clock_summary.csv           -- total wall-clock time per accession count
  <bench-dir>/results_<N>/pipeline_trace.tsv   -- per-task timing/resource-usage per accession
                                                   count (written by Nextflow itself; see the
                                                   trace{} block in nextflow/nextflow.config --
                                                   this script assumes trace.raw = true there)
  <bench-dir>/results_<N>/ledger/*.csv         -- per-task-outcome ledger rows (one file per
                                                   task; see nextflow/bin/write_ledger_row.sh
                                                   and design/latest_design-doc.md section 5)

Each results_<N> directory is one independent benchmark batch of N accessions.

Writes into <bench-dir>/plots/:
  01_total_wall_clock.png            total pipeline wall-clock time vs. accession count
  02_process_total_time.png          wall-clock time attributed to each process vs. accession count
  03_process_mean_task_time.png      mean per-task time per process vs. accession count
  04_cpu_utilization/*.png           per batch: CPUs allocated vs. actually used over time,
                                      plus active task count, and below it RAM reserved vs.
                                      actually used and swap -- are we filling server capacity,
                                      and if not, is memory what's holding it back?
  05_task_timeline/*.png             per batch: Gantt timeline for a sample of accessions,
                                      colored by stage, split into unitigs/contigs panels --
                                      what overlaps, and where are the gaps?
  06_runtime_distribution/*.png      per batch: task realtime distribution by stage x seq type
                                      (median, p90, p95, not just the mean) -- which subtasks
                                      are slow or unpredictable?
  07_throughput/*.png                per batch: accessions completed per time bin, and
                                      cumulative completions -- what sustained throughput are
                                      we achieving?
  08_throughput_vs_concurrency.png   all batches pooled: throughput vs. concurrency, colored by
                                      memory use and sized by disk I/O -- when does adding
                                      parallel work stop helping?
  09_task_resources/*.png, *.csv     per batch: per-stage table of task time, cores reserved vs.
                                      actually used, and RAM reserved vs. peak -- which stage is
                                      the resource bottleneck, and which over-reserve?
"""
import argparse
import csv
import getpass
import glob
import os
import re

# Matches run_benchmark.sh's own default: BENCH_DIR env var if set, else /scratch (not the
# repo's own disk -- a sweep can run into hundreds of GB), scoped by username.
_DEFAULT_BENCH_DIR = os.environ.get(
    "BENCH_DIR", f"/scratch/{getpass.getuser()}/logan_sketch_funcprofile_benchmark_runs"
)

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_wall_clock(bench_dir):
    path = os.path.join(bench_dir, "wall_clock_summary.csv")
    df = pd.read_csv(path)
    return df.sort_values("n_accessions")


def _numeric(series):
    return pd.to_numeric(series.astype(str).str.strip().replace({"-": None, "": None, "nan": None}), errors="coerce")


def _parse_utc_datetime(series):
    """Nextflow's trace submit/start/complete columns are raw epoch milliseconds under
    trace.raw = true (not formatted date strings) -- an unambiguous absolute instant, so
    this needs no timezone guessing to line up with write_ledger_row.sh's UTC timestamps."""
    return pd.to_datetime(_numeric(series), unit="ms", errors="coerce")


def _split_tag(process, tag):
    """'SOURMASH_SKETCH', 'DRR000836.unitigs' -> ('DRR000836', 'unitigs').
    FETCH_LOGAN's tag is just the accession (one task fetches both seq types), so its
    seq_type comes back NA -- callers that split by seq_type should treat NA as "applies
    to both"."""
    tag = "" if pd.isna(tag) else str(tag)
    if not tag:
        return pd.NA, pd.NA
    parts = tag.split(".")
    if parts[-1] in ("unitigs", "contigs"):
        return ".".join(parts[:-1]), parts[-1]
    if process == "WRITE_LEDGER_ROW" and len(parts) >= 2:
        return parts[0], parts[1]
    return tag, pd.NA


# Nextflow's "GB" (e.g. memory = '2400 GB') is 1024-based, as is `free -g`.
_GIB = 1024 ** 3

_REPORT_RECORD_RE = re.compile(r'"process":"([^"]*)".*?"tag":"([^"]*)".*?"memory":"(\d*)"')


def _reserved_memory_from_report(result_dir, df):
    """Traces written before `memory` was added to trace.fields don't record each task's
    reserved memory, but Nextflow's HTML report embeds every task record (not parseable as
    JSON -- scripts carry invalid escapes). Recover it record by record, matched on
    (process, tag)."""
    path = os.path.join(result_dir, "pipeline_report.html")
    if not os.path.exists(path):
        return pd.NA
    with open(path) as f:
        records = f.read().split('{"task_id":')[1:]
    lookup = {}
    for rec in records:
        m = _REPORT_RECORD_RE.search(rec)
        if m:
            lookup[(m.group(1), m.group(2))] = m.group(3)
    return [lookup.get((p, t), pd.NA) for p, t in zip(df["process"], df["tag"])]


def load_system_memory(bench_dir):
    """results_<N>/system_memory.tsv, sampled from /proc/meminfo by run_benchmark.sh while
    each batch ran. Whole-server numbers (includes other users and the OS)."""
    rows = []
    for path in sorted(glob.glob(os.path.join(bench_dir, "results_*", "system_memory.tsv"))):
        m = re.search(r"results_(\d+)", path)
        if not m:
            continue
        df = pd.read_csv(path, sep="\t")
        df["n_accessions"] = int(m.group(1))
        rows.append(df)
    if not rows:
        return pd.DataFrame(columns=["epoch_ms", "mem_total_bytes", "mem_used_bytes", "swap_used_bytes",
                                     "n_accessions", "ts"])
    mem = pd.concat(rows, ignore_index=True)
    mem["ts"] = pd.to_datetime(mem["epoch_ms"], unit="ms")
    return mem.sort_values("ts")


def load_trace(bench_dir):
    rows = []
    for trace_path in sorted(glob.glob(os.path.join(bench_dir, "results_*", "pipeline_trace.tsv"))):
        m = re.search(r"results_(\d+)", trace_path)
        if not m:
            continue
        df = pd.read_csv(trace_path, sep="\t", dtype=str)
        if "memory" not in df.columns:
            df["memory"] = _reserved_memory_from_report(os.path.dirname(trace_path), df)
        df["n_accessions"] = int(m.group(1))
        rows.append(df)
    if not rows:
        raise SystemExit(f"No pipeline_trace.tsv files found under {bench_dir}/results_*/")
    trace = pd.concat(rows, ignore_index=True)

    # Nextflow prefixes a process's trace name with its enclosing subworkflow (e.g.
    # "ANALYZE:SOURMASH_SKETCH") -- keep only the process itself.
    trace["process"] = trace["process"].str.split(":").str[-1]

    trace["cpus"] = _numeric(trace["cpus"])
    trace["pct_cpu"] = _numeric(trace["%cpu"])
    trace["realtime_s"] = _numeric(trace["realtime"]) / 1000.0
    trace["memory_bytes"] = _numeric(trace["memory"])
    trace["peak_rss_bytes"] = _numeric(trace["peak_rss"])
    trace["read_bytes"] = _numeric(trace["read_bytes"])
    trace["write_bytes"] = _numeric(trace["write_bytes"])
    trace["submit_ts"] = _parse_utc_datetime(trace["submit"])
    trace["start_ts"] = _parse_utc_datetime(trace["start"])
    trace["complete_ts"] = _parse_utc_datetime(trace["complete"])

    split = trace.apply(lambda r: _split_tag(r["process"], r["tag"]), axis=1, result_type="expand")
    trace[["accession", "seq_type"]] = split

    return trace


def load_ledger(bench_dir):
    """Union every per-task ledger CSV under results_*/ledger/ -- the same union-by-name
    approach design/latest_design-doc.md's own duckdb query uses. Plain csv parsing
    (rather than one pandas.read_csv call per file) because a large sweep produces tens
    of thousands of one-row files."""
    columns = ["accession", "seq_type", "stage", "status", "shard_id", "updated_at", "n_accessions"]
    records = []
    for result_dir in sorted(glob.glob(os.path.join(bench_dir, "results_*"))):
        m = re.search(r"results_(\d+)$", result_dir)
        if not m:
            continue
        n_accessions = int(m.group(1))
        for path in glob.glob(os.path.join(result_dir, "ledger", "*.csv")):
            with open(path, newline="") as f:
                reader = csv.reader(f)
                header = next(reader, None)
                row = next(reader, None)
            if not header or not row:
                continue
            record = dict(zip(header, row))
            record["n_accessions"] = n_accessions
            records.append(record)
    if not records:
        return pd.DataFrame(columns=columns)
    ledger = pd.DataFrame.from_records(records)
    ledger["updated_at"] = pd.to_datetime(ledger["updated_at"], errors="coerce", utc=True).dt.tz_localize(None)
    return ledger


def stage_color_map(processes):
    cmap = plt.get_cmap("tab10")
    return {p: cmap(i % 10) for i, p in enumerate(sorted(processes))}


def _step_series(starts, ends, values):
    """Event-based step function: each task contributes `values` for its whole [start,
    complete) interval. Returns (event_times, cumulative_value) sorted ascending -- plot
    directly with ax.step(..., where='post'), or sample at arbitrary times with
    _sample_step."""
    starts = pd.to_datetime(starts).to_numpy()
    ends = pd.to_datetime(ends).to_numpy()
    values = np.asarray(values, dtype=float)
    ev_time = np.concatenate([starts, ends])
    ev_delta = np.concatenate([values, -values])
    order = np.argsort(ev_time, kind="stable")
    ev_time = ev_time[order]
    cum = np.cumsum(ev_delta[order])
    return ev_time, cum


def _sample_step(ev_time, cum, query_times):
    query_times = np.asarray(pd.to_datetime(pd.Series(query_times)).to_numpy())
    idx = np.searchsorted(ev_time, query_times, side="right") - 1
    return np.where(idx >= 0, cum[np.clip(idx, 0, len(cum) - 1)], 0.0)


def plot_total_wall_clock(wall_df, out_path):
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(wall_df["n_accessions"], wall_df["wall_clock_seconds"], marker="o", linewidth=2)
    ax.set_xscale("log")
    ax.set_xlabel("Number of accessions processed")
    ax.set_ylabel("Total wall-clock time (seconds)")
    ax.set_title("Logan metagenomics pipeline: end-to-end wall-clock time vs. accession count")
    ax.grid(True, which="both", linestyle="--", alpha=0.4)
    for x, y in zip(wall_df["n_accessions"], wall_df["wall_clock_seconds"]):
        ax.annotate(f"{y:,.0f}s", (x, y), textcoords="offset points", xytext=(0, 8), ha="center")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_process_total_time(completed_df, wall_df, out_path):
    """Grouped bars per batch size (one bar per process), where each process's bar is
    that process's summed task time rescaled so that, summed across processes, they'd
    add up to plot_total_wall_clock's wall-clock value at that n_accessions. This
    answers "of the 1hr this batch took, how much was fetch vs. sketch vs. funprofiler." Note
    this is an attribution, not a measured per-process wall-clock: because tasks of
    different processes run concurrently, wall-clock time isn't actually "spent on"
    one process at a time, so each bar is that process's summed task time rescaled
    from the (much larger) sum of all concurrently-running task time down to the real
    wall-clock total."""
    agg = completed_df.groupby(["n_accessions", "process"])["realtime_s"].sum().reset_index()
    wall_clock = wall_df.set_index("n_accessions")["wall_clock_seconds"]

    # Only batch sizes with a recorded wall-clock total can be attributed -- a size
    # present in the trace files but missing here means that run errored out before
    # run_benchmark.sh appended a row to wall_clock_summary.csv (e.g. a partial/failed
    # sweep), so there's no real total to split across processes.
    all_sizes = sorted(agg["n_accessions"].unique())
    sizes = [n for n in all_sizes if n in wall_clock.index]
    skipped = [n for n in all_sizes if n not in wall_clock.index]
    if skipped:
        print(f"plot_process_total_time: skipping n_accessions without a wall-clock total: {skipped}")
    processes = sorted(agg["process"].unique())

    fig, ax = plt.subplots(figsize=(8, 5))
    x = list(range(len(sizes)))
    n_proc = len(processes)
    width = 0.8 / max(n_proc, 1)
    for i, proc in enumerate(processes):
        values = []
        for n in sizes:
            proc_time = agg[(agg["n_accessions"] == n) & (agg["process"] == proc)]["realtime_s"].sum()
            total_time = agg[agg["n_accessions"] == n]["realtime_s"].sum()
            scale = wall_clock[n] / total_time if total_time else 0.0
            values.append(proc_time * scale)
        offsets = [xi + (i - (n_proc - 1) / 2) * width for xi in x]
        ax.bar(offsets, values, width=width, label=proc)

    ax.set_xticks(x)
    ax.set_xticklabels([f"{n:,}" for n in sizes])
    ax.set_xlabel("Number of accessions processed")
    ax.set_ylabel("Wall-clock time attributed to that process (seconds)")
    ax.set_title("Logan metagenomics pipeline: wall-clock time by process vs. accession count")
    ax.legend(title="Process", bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.grid(True, axis="y", linestyle="--", alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_process_mean_task_time(completed_df, out_path):
    agg = completed_df.groupby(["n_accessions", "process"])["realtime_s"].mean().reset_index()
    processes = sorted(agg["process"].unique())

    fig, ax = plt.subplots(figsize=(8, 5))
    for proc in processes:
        sub = agg[agg["process"] == proc].sort_values("n_accessions")
        ax.plot(sub["n_accessions"], sub["realtime_s"], marker="o", label=proc)

    ax.set_xscale("log")
    ax.set_xlabel("Number of accessions processed")
    ax.set_ylabel("Mean single-task time (seconds)")
    ax.set_title(
        "Logan metagenomics pipeline: mean per-task time by process\n"
        "(flat lines confirm each accession's task cost is independent of run size)"
    )
    ax.legend(title="Process")
    ax.grid(True, which="both", linestyle="--", alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_cpu_utilization(trace_df, sysmem_df, out_dir, total_cpus=None, memory_limit_gb=None):
    """One figure per batch. Top: allocated vs. actually-used CPU cores over wall-clock time,
    plus the number of concurrently active tasks (right axis) -- "are we filling the
    available server capacity?" Bottom: RAM reserved by running tasks (what Nextflow counts
    against its memory limit) vs. RAM and swap actually used on the server -- explains
    stretches where CPU usage stays low because memory, not cores, is what's full. CPU and
    reserved-RAM lines are exact event-based step functions over each task's [start,
    complete) interval; actual RAM/swap are periodic samples from system_memory.tsv."""
    os.makedirs(out_dir, exist_ok=True)

    # Memory logs without a per-benchmark column can't tell this benchmark's RAM from other
    # users'. Use RAM in use within an hour before the sweep's first task or after its last one
    # (i.e. while no batch was running) as the "other users + OS" level.
    baseline_gb = None
    if not sysmem_df.empty:
        hour = pd.Timedelta(hours=1)
        first, last = trace_df["start_ts"].min(), trace_df["complete_ts"].max()
        quiet = sysmem_df[((sysmem_df["ts"] >= first - hour) & (sysmem_df["ts"] < first))
                          | ((sysmem_df["ts"] > last) & (sysmem_df["ts"] <= last + hour))]
        if not quiet.empty:
            baseline_gb = quiet["mem_used_bytes"].min() / _GIB

    for n, sub in trace_df.groupby("n_accessions"):
        sub = sub.dropna(subset=["start_ts", "complete_ts", "cpus"])
        if sub.empty:
            print(f"plot_cpu_utilization: no timestamped tasks for n={n}, skipping")
            continue
        pct_cpu = sub["pct_cpu"].fillna(0.0)

        ev_time, active = _step_series(sub["start_ts"], sub["complete_ts"], pd.Series(1.0, index=sub.index))
        _, cpu_alloc = _step_series(sub["start_ts"], sub["complete_ts"], sub["cpus"])
        _, cpu_actual = _step_series(sub["start_ts"], sub["complete_ts"], pct_cpu / 100.0)

        fig, (ax1, ax3) = plt.subplots(2, 1, figsize=(11, 9), sharex=True)
        ax1.step(ev_time, cpu_alloc, where="post", label="CPUs allocated", color="tab:blue")
        ax1.step(ev_time, cpu_actual, where="post", label="CPUs actually used (measured %cpu)", color="tab:orange")
        if total_cpus:
            ax1.axhline(total_cpus, color="tab:red", linestyle=":", linewidth=1,
                        label=f"Server capacity ({total_cpus} cpus)")
        ax1.set_ylabel("CPU cores")
        ax1.grid(True, linestyle="--", alpha=0.4)

        ax2 = ax1.twinx()
        ax2.step(ev_time, active, where="post", color="tab:green", alpha=0.6, label="Active tasks")
        ax2.set_ylabel("Active tasks", color="tab:green")
        ax2.tick_params(axis="y", labelcolor="tab:green")

        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left", fontsize=8)

        ax1.set_title(f"CPU allocation vs. actual usage over time (n={n:,} accessions)")

        if sub["memory_bytes"].notna().any():
            _, mem_reserved = _step_series(sub["start_ts"], sub["complete_ts"], sub["memory_bytes"].fillna(0.0))
            ax3.step(ev_time, mem_reserved / _GIB, where="post", color="tab:blue",
                     label="RAM reserved by running tasks (declared memory)")
        mem = sysmem_df[sysmem_df["n_accessions"] == n]
        pad = pd.Timedelta(minutes=10)
        mem = mem[(mem["ts"] >= sub["start_ts"].min() - pad) & (mem["ts"] <= sub["complete_ts"].max() + pad)]
        if not mem.empty:
            ts = mem["ts"]
            used = mem["mem_used_bytes"] / _GIB
            total_gb = mem["mem_total_bytes"].max() / _GIB
            # Plain lines, each measured from zero -- stacked areas made the thin "other users"
            # layer look as large as the benchmark underneath it.
            if "bench_pss_bytes" in mem and mem["bench_pss_bytes"].notna().any():
                bench = (mem["bench_pss_bytes"] / _GIB).clip(upper=used)
                ax3.plot(ts, bench, color="tab:orange", linewidth=1.5, label="Used by this benchmark")
                ax3.plot(ts, used - bench, color="tab:brown", linewidth=1.5, label="Used by other users + OS")
            else:
                # Older logs can't split the benchmark from other users -- total is all there is.
                ax3.plot(ts, used, color="tab:orange", linewidth=1.5,
                         label="Used (whole server: this benchmark + other users + OS)")
                if baseline_gb is not None:
                    ax3.axhline(baseline_gb, color="tab:brown", linestyle="--", linewidth=1.2,
                                label=f"Used while benchmark not running (other users + OS): {baseline_gb:,.0f} GB")
            ax3.fill_between(ts, used, total_gb, color="tab:green", alpha=0.12, label="Available (unused)")
            ax3.plot(ts, mem["swap_used_bytes"] / _GIB, color="tab:purple", label="Swap used (whole server)")
            ax3.axhline(total_gb, color="tab:red", linestyle=":", linewidth=1,
                        label=f"Server RAM ({total_gb:,.0f} GB)")
        else:
            print(f"plot_cpu_utilization: no system_memory.tsv for n={n}; RAM panel shows reserved memory only")
        if memory_limit_gb:
            ax3.axhline(memory_limit_gb, color="gray", linestyle="--", linewidth=1.2,
                        label=f"Nextflow memory limit ({memory_limit_gb:,.0f} GB)")
        ax3.set_ylabel("Memory (GB)")
        ax3.set_xlabel("Wall-clock time (UTC)")
        ax3.set_ylim(bottom=0)
        ax3.grid(True, linestyle="--", alpha=0.4)
        ax3.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=8)
        ax3.set_title("Memory: reserved by tasks vs. used vs. available")
        ax3.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))

        fig.autofmt_xdate()
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f"cpu_utilization_{n}.png"), dpi=150, bbox_inches="tight")
        plt.close(fig)


def plot_task_timeline(trace_df, out_dir, colors, sample_size=15, accessions=None):
    """One figure per batch: a Gantt-style timeline of individual tasks, one row per
    accession (a small deterministic sample -- plotting every accession in a
    multi-thousand-accession batch is illegible), split into unitigs/contigs panels and
    colored by pipeline stage. Answers "what overlaps, and where are there gaps?"
    FETCH_LOGAN fetches both seq types in one task, so its bar is drawn in both panels."""
    os.makedirs(out_dir, exist_ok=True)
    for n, sub in trace_df.groupby("n_accessions"):
        sub = sub.dropna(subset=["start_ts", "complete_ts", "accession"])
        if sub.empty:
            print(f"plot_task_timeline: no timestamped tasks for n={n}, skipping")
            continue

        if accessions:
            wanted = set(accessions)
            chosen = [a for a in sub["accession"].unique() if a in wanted]
        else:
            first_start = sub.groupby("accession")["start_ts"].min().sort_values()
            chosen = list(first_start.index[:sample_size])
        if not chosen:
            print(f"plot_task_timeline: none of the requested accessions found for n={n}, skipping")
            continue

        rows = sub[sub["accession"].isin(chosen)]
        order = {a: i for i, a in enumerate(chosen)}

        fig, axes = plt.subplots(1, 2, figsize=(14, max(3, 0.4 * len(chosen))), sharey=True)
        for ax, seq_type in zip(axes, ["unitigs", "contigs"]):
            panel_rows = rows[rows["seq_type"].isna() | (rows["seq_type"] == seq_type)]
            for _, row in panel_rows.iterrows():
                y = order[row["accession"]]
                start_num = mdates.date2num(row["start_ts"])
                width = max(mdates.date2num(row["complete_ts"]) - start_num, 1e-6)
                ax.broken_barh([(start_num, width)], (y - 0.4, 0.8), facecolors=colors.get(row["process"], "gray"))
            ax.xaxis_date()
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
            ax.set_title(seq_type)
            ax.set_xlabel("Wall-clock time (UTC)")
            ax.grid(True, axis="x", linestyle="--", alpha=0.4)

        axes[0].set_yticks(list(order.values()))
        axes[0].set_yticklabels(list(order.keys()))
        axes[0].set_ylabel("Accession")

        handles = [plt.Line2D([0], [0], color=c, lw=6, label=p) for p, c in sorted(colors.items())]
        fig.legend(handles=handles, loc="upper center", ncol=len(colors), fontsize=8, bbox_to_anchor=(0.5, 1.08))
        fig.suptitle(f"Task timeline, {len(chosen)} of {sub['accession'].nunique():,} accessions (n={n:,} batch)",
                     y=1.15)
        fig.autofmt_xdate()
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f"task_timeline_{n}.png"), dpi=150, bbox_inches="tight")
        plt.close(fig)


def plot_runtime_distribution(trace_df, out_dir):
    """One figure per batch: realtime distribution per stage, split into unitigs/contigs
    panels (same layout as plot_task_timeline) -- box is the IQR, the line through it is
    the median, whiskers extend to p5/p95, and a diamond marks p90. Answers "which
    subtasks are slow or unpredictable?": a wide box or a p90/p95 far above the median
    flags a stage whose cost an average alone would hide. FETCH_LOGAN fetches both seq
    types in one task, so its distribution is drawn in both panels."""
    os.makedirs(out_dir, exist_ok=True)
    for n, sub in trace_df.groupby("n_accessions"):
        sub = sub[sub["status"] == "COMPLETED"].dropna(subset=["realtime_s", "process"]).copy()
        if sub.empty:
            print(f"plot_runtime_distribution: no completed tasks for n={n}, skipping")
            continue

        processes = sorted(sub["process"].unique())

        fig, axes = plt.subplots(1, 2, figsize=(max(9, 1.4 * len(processes)), 5), sharey=True)
        for panel_i, (ax, seq_type) in enumerate(zip(axes, ["unitigs", "contigs"])):
            panel_rows = sub[sub["seq_type"].isna() | (sub["seq_type"] == seq_type)]
            positions, data, p90s = [], [], []
            for j, proc in enumerate(processes):
                values = panel_rows[panel_rows["process"] == proc]["realtime_s"]
                if values.empty:
                    continue
                positions.append(j)
                data.append(values.to_numpy())
                p90s.append(np.percentile(values, 90))

            if data:
                bp = ax.boxplot(data, positions=positions, widths=0.6, whis=(5, 95),
                                 patch_artist=True, showfliers=False)
                for box in bp["boxes"]:
                    box.set_facecolor("tab:blue")
                    box.set_alpha(0.6)
                ax.scatter(positions, p90s, marker="D", color="tab:orange", edgecolor="black", s=25,
                           zorder=5, label="p90" if panel_i == 0 else None)

            ax.set_xticks(range(len(processes)))
            ax.set_xticklabels(processes, rotation=30, ha="right")
            ax.set_title(seq_type)
            ax.grid(True, axis="y", which="both", linestyle="--", alpha=0.3)

        axes[0].set_yscale("log")
        axes[0].set_ylabel("Task realtime, seconds (log scale)")
        axes[0].legend(fontsize=8, loc="upper left")
        fig.suptitle(f"Runtime distribution by stage (n={n:,} accessions)\n"
                     "box = IQR, line = median, whiskers = p5/p95, diamond = p90")
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f"runtime_distribution_{n}.png"), dpi=150)
        plt.close(fig)


def plot_throughput(ledger_df, trace_df, out_dir, bin_minutes=None):
    """One figure per batch: accessions completed per time bin (bars) and cumulative
    completions (line). Answers "what sustained throughput are we achieving?" An
    accession counts as done at the latest ledger row recorded for it (sketch /
    funprofiler per seq type, or a skipped-empty fetch row) -- the ledger is the
    pipeline's own record of per-accession outcomes, so this doesn't need to guess from
    trace timestamps which task is "the last one" for an accession."""
    os.makedirs(out_dir, exist_ok=True)
    for n, sub in ledger_df.groupby("n_accessions"):
        completions = sub.groupby("accession")["updated_at"].max().dropna().sort_values()
        if completions.empty:
            print(f"plot_throughput: no ledger completions for n={n}, skipping")
            continue

        batch_start = trace_df.loc[trace_df["n_accessions"] == n, "submit_ts"].min()
        if pd.isna(batch_start):
            batch_start = completions.min()

        elapsed_hours = (completions - batch_start).dt.total_seconds() / 3600.0
        span_hours = max(elapsed_hours.max(), 1e-6)
        bin_h = (bin_minutes / 60.0) if bin_minutes else max(span_hours / 24.0, 1.0 / 60.0)
        edges = np.arange(0.0, span_hours + bin_h, bin_h)
        counts, _ = np.histogram(elapsed_hours, bins=edges)
        cumulative = np.cumsum(counts)
        centers = (edges[:-1] + edges[1:]) / 2

        bin_min = bin_h * 60
        bin_label = f"{bin_min:.1f}-min" if bin_min >= 1 else f"{bin_min * 60:.0f}-sec"

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
        ax1.bar(centers, counts, width=bin_h * 0.9, color="tab:blue")
        ax1.set_ylabel(f"Accessions completed / {bin_label} bin")
        ax1.set_title(f"Throughput (n={n:,} accessions)")
        ax1.grid(True, axis="y", linestyle="--", alpha=0.4)

        ax2.plot(centers, cumulative, marker=".", color="tab:purple")
        ax2.set_ylabel("Cumulative accessions completed")
        ax2.set_xlabel("Hours since batch start")
        ax2.grid(True, linestyle="--", alpha=0.4)

        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f"throughput_{n}.png"), dpi=150)
        plt.close(fig)


def plot_task_resources(trace_df, out_dir):
    """One table per batch (PNG + CSV): per pipeline stage, task count, time, and reserved vs.
    actual cores and RAM, as medians with p95. Actual cores = %cpu averaged over each task's run
    (bursts don't show); peak RAM = Nextflow's peak_rss, which sums a task's whole process tree
    and so overstates real use, compared against server-level memory."""
    os.makedirs(out_dir, exist_ok=True)
    for n, sub in trace_df[trace_df["status"] == "COMPLETED"].groupby("n_accessions"):
        sub = sub.assign(cores_used=sub["pct_cpu"] / 100,
                         mem_reserved_gb=sub["memory_bytes"] / _GIB,
                         peak_ram_gb=sub["peak_rss_bytes"] / _GIB,
                         task_hours=sub["realtime_s"] / 3600)
        p95 = lambda s: s.quantile(0.95)  # noqa: E731
        g = sub.groupby("process").agg(
            tasks=("process", "size"),
            total_task_hours=("task_hours", "sum"),
            median_task_min=("realtime_s", lambda s: s.median() / 60),
            p95_task_min=("realtime_s", lambda s: s.quantile(0.95) / 60),
            cores_reserved=("cpus", "median"),
            cores_used_median=("cores_used", "median"),
            cores_used_p95=("cores_used", p95),
            ram_reserved_gb=("mem_reserved_gb", "median"),
            peak_ram_gb_median=("peak_ram_gb", "median"),
            peak_ram_gb_p95=("peak_ram_gb", p95),
        ).sort_values("total_task_hours", ascending=False)
        g.round(3).to_csv(os.path.join(out_dir, f"task_resources_{n}.csv"))

        headers = ["Stage", "Tasks", "Task-hours\n(total)", "Task time\nmedian / p95",
                   "Cores\nreserved", "Cores used\nmedian / p95", "RAM\nreserved", "Peak RAM\nmedian / p95"]
        rows = [[proc, f"{int(r.tasks):,}", f"{r.total_task_hours:,.1f}",
                 f"{r.median_task_min:.1f} / {r.p95_task_min:.1f} min",
                 f"{r.cores_reserved:g}", f"{r.cores_used_median:.2f} / {r.cores_used_p95:.2f}",
                 "—" if pd.isna(r.ram_reserved_gb) else f"{r.ram_reserved_gb:g} GB",
                 f"{r.peak_ram_gb_median:.2f} / {r.peak_ram_gb_p95:.2f} GB"]
                for proc, r in g.iterrows()]

        fig, ax = plt.subplots(figsize=(14, 0.42 * (len(rows) + 1) + 0.8))
        ax.axis("off")
        table = ax.table(cellText=rows, colLabels=headers, cellLoc="center", bbox=[0, 0, 1, 1],
                         colWidths=[0.17, 0.08, 0.1, 0.13, 0.08, 0.14, 0.1, 0.2])
        table.auto_set_font_size(False)
        table.set_fontsize(9.5)
        for (row, col), cell in table.get_celld().items():
            cell.set_edgecolor("#d8dcd8")
            if row == 0:
                cell.set_facecolor("#eceee9")
                cell.set_text_props(weight="bold")
            if col == 0 and row > 0:
                cell.set_text_props(ha="left")
                cell.PAD = 0.04
        ax.set_title(f"Per-stage resource use, n={n:,} accessions (completed tasks)", fontsize=12, pad=12)
        fig.text(0.125, -0.02,
                 "Cores used = %cpu averaged over each task's run (short bursts don't show). Peak RAM = Nextflow "
                 "peak_rss, which counts shared memory more than once and overstates real use.",
                 fontsize=8, color="#59615c")
        fig.savefig(os.path.join(out_dir, f"task_resources_{n}.png"), dpi=150, bbox_inches="tight")
        plt.close(fig)


def plot_throughput_vs_concurrency(trace_df, ledger_df, out_path, bin_minutes=15):
    """One combined figure across every batch: throughput (completions/hour, from the
    ledger) vs. concurrency (active tasks, from the trace) in fixed time bins, colored by
    aggregate memory use and sized by aggregate disk I/O rate. Memory/disk are
    approximated by holding each task's one measured peak_rss / total read+write bytes
    constant over its [start, complete) interval -- the same convention
    plot_cpu_utilization uses for %cpu, since Nextflow's trace records one summary value
    per task rather than a continuous series. Answers "when does adding parallel work
    stop helping?" """
    bin_delta = pd.Timedelta(minutes=bin_minutes)
    points = []
    for n, sub in trace_df.groupby("n_accessions"):
        sub = sub.dropna(subset=["start_ts", "complete_ts"])
        if sub.empty:
            continue
        t0, t1 = sub["start_ts"].min(), sub["complete_ts"].max()
        if pd.isna(t0) or pd.isna(t1) or t0 >= t1:
            continue
        edges = pd.date_range(t0, t1 + bin_delta, freq=bin_delta)
        if len(edges) < 2:
            continue
        mids = edges[:-1] + bin_delta / 2

        active_ev, active_cum = _step_series(sub["start_ts"], sub["complete_ts"], pd.Series(1.0, index=sub.index))
        mem_ev, mem_cum = _step_series(sub["start_ts"], sub["complete_ts"], sub["peak_rss_bytes"].fillna(0.0))
        disk_rate = ((sub["read_bytes"].fillna(0.0) + sub["write_bytes"].fillna(0.0))
                     / sub["realtime_s"].replace(0, np.nan)).fillna(0.0)
        disk_ev, disk_cum = _step_series(sub["start_ts"], sub["complete_ts"], disk_rate)

        concurrency = _sample_step(active_ev, active_cum, mids)
        mem_gb = _sample_step(mem_ev, mem_cum, mids) / 1e9
        disk_mbps = _sample_step(disk_ev, disk_cum, mids) / 1e6

        n_completions = ledger_df.loc[ledger_df["n_accessions"] == n].groupby("accession")["updated_at"].max().dropna()
        comp_vals = n_completions.to_numpy()
        edge_vals = edges.to_numpy()
        bin_idx = np.searchsorted(edge_vals, comp_vals, side="right") - 1
        valid = (bin_idx >= 0) & (bin_idx < len(mids))
        comp_counts = np.bincount(bin_idx[valid], minlength=len(mids))
        throughput = comp_counts / (bin_minutes / 60.0)

        for i in range(len(mids)):
            points.append((n, concurrency[i], throughput[i], mem_gb[i], disk_mbps[i]))

    if not points:
        print("plot_throughput_vs_concurrency: no data, skipping")
        return
    df = pd.DataFrame(points, columns=["n_accessions", "concurrency", "throughput_per_hr", "mem_gb", "disk_mbps"])
    df = df[df["concurrency"] > 0]
    if df.empty:
        print("plot_throughput_vs_concurrency: no active-task bins, skipping")
        return

    fig, ax = plt.subplots(figsize=(9, 6))
    max_disk = max(df["disk_mbps"].max(), 1e-9)
    sizes = 20 + 200 * (df["disk_mbps"] / max_disk)
    sc = ax.scatter(df["concurrency"], df["throughput_per_hr"], c=df["mem_gb"], s=sizes,
                     cmap="viridis", alpha=0.7, edgecolor="black", linewidth=0.3)
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("Aggregate memory in use (GB, peak_rss held constant over each task)")

    # Median throughput per concurrency bucket, to show where the curve flattens.
    n_buckets = min(15, df["concurrency"].nunique())
    if n_buckets >= 2:
        df = df.assign(_bucket=pd.qcut(df["concurrency"], n_buckets, duplicates="drop"))
        trend = (df.groupby("_bucket", observed=True)
                   .agg(concurrency=("concurrency", "median"), throughput_per_hr=("throughput_per_hr", "median"))
                   .sort_values("concurrency"))
        ax.plot(trend["concurrency"], trend["throughput_per_hr"], color="tab:red", linewidth=2,
                marker="o", label="Median throughput per concurrency bucket")
        ax.legend(fontsize=8)

    ax.set_xlabel("Concurrency (active tasks)")
    ax.set_ylabel("Throughput (completed accessions / hour)")
    ax.set_title(f"Throughput vs. concurrency, {bin_minutes}-min bins, all batches pooled\n"
                 "point size = aggregate disk I/O rate (read+write bytes/s, held constant over each task)")
    ax.grid(True, linestyle="--", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--bench-dir",
        default=_DEFAULT_BENCH_DIR,
        help=f"Directory produced by run_benchmark.sh (default: {_DEFAULT_BENCH_DIR})",
    )
    parser.add_argument(
        "--total-cpus", type=int, default=None,
        help="Draw a reference line at this many cores on the CPU-utilization plot. This "
             "machine is shared with other users, so there's no "
             "single correct default -- pass whatever core count you actually had for the run.",
    )
    parser.add_argument(
        "--memory-limit-gb", type=float, default=None,
        help="Draw a reference line at Nextflow's memory limit (executor.memory in "
             "nextflow/nextflow.config) on the RAM panel of the CPU-utilization plot.",
    )
    parser.add_argument(
        "--timeline-sample-size", type=int, default=15,
        help="Number of accessions (earliest-starting) to show on the task-timeline plot (default: 15)",
    )
    parser.add_argument(
        "--timeline-accessions", default=None,
        help="Comma-separated accession list to show on the task-timeline plot instead of the "
             "default earliest-starting sample",
    )
    parser.add_argument(
        "--throughput-bin-minutes", type=float, default=None,
        help="Bin width for the per-batch completions/cumulative plot (default: batch span / 24)",
    )
    parser.add_argument(
        "--concurrency-bin-minutes", type=float, default=15,
        help="Bin width for the throughput-vs-concurrency plot (default: 15)",
    )
    args = parser.parse_args()

    wall_df = load_wall_clock(args.bench_dir)
    trace_df = load_trace(args.bench_dir)
    ledger_df = load_ledger(args.bench_dir)
    completed_df = trace_df[trace_df["status"] == "COMPLETED"]

    plots_dir = os.path.join(args.bench_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    plot_total_wall_clock(wall_df, os.path.join(plots_dir, "01_total_wall_clock.png"))
    plot_process_total_time(completed_df, wall_df, os.path.join(plots_dir, "02_process_total_time.png"))
    plot_process_mean_task_time(completed_df, os.path.join(plots_dir, "03_process_mean_task_time.png"))

    colors = stage_color_map(trace_df["process"].dropna().unique())
    timeline_accessions = args.timeline_accessions.split(",") if args.timeline_accessions else None

    plot_cpu_utilization(trace_df, load_system_memory(args.bench_dir),
                         os.path.join(plots_dir, "04_cpu_utilization"),
                         total_cpus=args.total_cpus, memory_limit_gb=args.memory_limit_gb)
    plot_task_timeline(trace_df, os.path.join(plots_dir, "05_task_timeline"), colors,
                        sample_size=args.timeline_sample_size, accessions=timeline_accessions)
    plot_runtime_distribution(trace_df, os.path.join(plots_dir, "06_runtime_distribution"))
    plot_throughput(ledger_df, trace_df, os.path.join(plots_dir, "07_throughput"),
                     bin_minutes=args.throughput_bin_minutes)
    plot_throughput_vs_concurrency(trace_df, ledger_df,
                                    os.path.join(plots_dir, "08_throughput_vs_concurrency.png"),
                                    bin_minutes=args.concurrency_bin_minutes)
    plot_task_resources(trace_df, os.path.join(plots_dir, "09_task_resources"))

    print(f"Wrote plots to {plots_dir}")


if __name__ == "__main__":
    main()
