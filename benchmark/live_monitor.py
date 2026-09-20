#!/usr/bin/env python3
"""Watch a run that is still going and redraw its plots every interval.

Same measurements benchmark/plot_benchmark.py makes after a benchmark sweep -- CPU
allocated vs. actually used, active tasks, memory, per-stage task time, per-stage resource
table -- plus the one a production run needs while it's happening: how much of the manifest
is done, how fast, and when it will finish.

Reads (all written while the run is in progress, nothing has to finish first):
  <run-dir>/traces/pipeline_trace.*.tsv  one per shard attempt, written by Nextflow as tasks
                                         complete (trace{} in nextflow/nextflow.config)
  <run-dir>/system_memory.tsv            server RAM/swap samples (benchmark/log_system_memory.sh)
  <run-dir>/manifest.txt                 the frozen manifest, i.e. the denominator
  <run-dir>/state/                       per-shard completion state from nextflow/run_full.sh

Writes into <run-dir>/plots_live/ every --interval seconds:
  index.html                 auto-refreshing dashboard tying the below together -- open this
  00_progress.png            accessions done vs. the whole manifest, completion rate and ETA,
                              completions per time bin, and per-stage task counts
  04_cpu_utilization.png     CPUs allocated vs. actually used, active tasks, RAM and swap
  05_task_timeline.png       Gantt timeline of the accessions currently moving through
  06_runtime_distribution.png  per-stage task runtimes
  09_task_resources.png/.csv  per-stage reserved vs. actual cores and RAM
  progress.json              the same numbers as data, for scripting or a status line

The resource plots cover the last --window-hours only: they're about what the machine is
doing now, and a multi-week trace doesn't fit in one readable figure. Progress, rate and ETA
always cover the whole run.

Usage:
  benchmark/live_monitor.py --run-dir /scratch/$USER/logan_full_run
  benchmark/live_monitor.py --run-dir ... --once          # one snapshot, then exit
  benchmark/live_monitor.py --trace 'results/pipeline_trace.tsv' --total 1000 --out-dir plots_live

nextflow/run_full.sh starts one of these automatically unless it's given --no-monitor.
"""
import argparse
import glob
import json
import os
import re
import sys
import time
from datetime import datetime, timedelta

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import plot_benchmark as pb

# The last task an accession runs for one sequence type: CLEANUP_FASTA after both analyses are
# done with the FASTA (subworkflows/analyze.nf), or the ledger row written for a sequence type
# that was skipped as empty (main.nf). Seeing one of these means that half of the accession is
# finished and published.
TERMINAL_PROCESSES = {"CLEANUP_FASTA", "WRITE_LEDGER_ROW"}


def read_trace_file(path):
    """One shard's trace. Tolerates the partial final line of a file Nextflow is writing."""
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, on_bad_lines="skip")
    except (pd.errors.EmptyDataError, FileNotFoundError):
        return None
    if df.empty or "process" not in df.columns:
        return None
    return df


class TraceReader:
    """Keeps the run's traces summarized instead of in memory.

    A full run produces millions of trace rows across hundreds of shard files, which is too
    much to re-parse (or hold) every interval. Each file is parsed once and reduced to what
    the dashboard needs; a file is only re-parsed when its size or mtime changes, which in
    practice means only the shard currently running. Full rows are kept only for files
    touched inside the plotting window, since that's all the resource plots draw."""

    def __init__(self, patterns, window_hours):
        self.patterns = patterns
        self.window_hours = window_hours
        self._cache = {}   # path -> (mtime, size, summary dict)

    def _cutoff(self):
        return pd.Timestamp.now("UTC").tz_convert(None) - pd.Timedelta(hours=self.window_hours)

    def _summarize(self, df):
        df = pb.prepare_trace(df.copy())
        terminal = df[df["process"].isin(TERMINAL_PROCESSES) & (df["status"] == "COMPLETED")]
        # One row per (accession, seq_type) that has finished, with when it finished.
        halves = (terminal.dropna(subset=["accession", "seq_type"])
                          .groupby(["accession", "seq_type"])["complete_ts"].max()
                          .reset_index())
        counts = (df.groupby(["process", "status"]).size()
                    .rename("tasks").reset_index())
        failed = sorted(set(df.loc[df["status"] == "FAILED", "accession"].dropna()))
        # Keep only the rows the resource plots can still use. Without this, first reading a
        # run that already has hundreds of shard traces would hold every one of their millions
        # of rows in memory at once, just to drop them again.
        cutoff = self._cutoff()
        rows = df[df["complete_ts"].isna() | (df["complete_ts"] >= cutoff)]
        return {"halves": halves, "counts": counts, "failed": failed, "rows": rows,
                "last_event": df["complete_ts"].max()}

    def load(self):
        paths = sorted({p for pattern in self.patterns for p in glob.glob(pattern)})
        for path in paths:
            try:
                stat = os.stat(path)
            except FileNotFoundError:
                continue
            key = (stat.st_mtime, stat.st_size)
            if self._cache.get(path, (None,))[0:2] != key:
                df = read_trace_file(path)
                if df is None:
                    continue
                self._cache[path] = (*key, self._summarize(df))

        summaries = [entry[2] for entry in self._cache.values()]
        if not summaries:
            return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), []

        halves = pd.concat([s["halves"] for s in summaries], ignore_index=True)
        counts = (pd.concat([s["counts"] for s in summaries], ignore_index=True)
                    .groupby(["process", "status"])["tasks"].sum().reset_index())
        failed = sorted({a for s in summaries for a in s["failed"]})

        # Rows for the resource plots: only what is still inside the window. A shard spans many
        # hours, so its file can be recent while most of its rows are not.
        cutoff = self._cutoff()
        recent = [s["rows"] for s in summaries if not s["rows"].empty]
        rows = pd.concat(recent, ignore_index=True) if recent else pd.DataFrame(columns=summaries[0]["rows"].columns)
        if not rows.empty:
            rows = rows[rows["complete_ts"].isna() | (rows["complete_ts"] >= cutoff)]
        return halves, counts, rows, failed

    def forget_stale(self):
        """Release the cached rows of files that have dropped out of the plotting window."""
        cutoff = self._cutoff()
        for path, (mtime, size, summary) in self._cache.items():
            last = summary.get("last_event")
            if summary.get("rows") is not None and pd.notna(last) and last < cutoff:
                summary["rows"] = summary["rows"].iloc[0:0]


def load_system_memory(path, n_accessions):
    """benchmark/log_system_memory.sh's TSV. It is appended to across pauses and resumes, so
    a header line can appear mid-file; those rows are dropped rather than parsed."""
    if not path or not os.path.exists(path):
        return pd.DataFrame(columns=["epoch_ms", "mem_total_bytes", "mem_used_bytes",
                                     "swap_used_bytes", "bench_pss_bytes", "n_accessions", "ts"])
    mem = pd.read_csv(path, sep="\t", dtype=str, on_bad_lines="skip")
    for column in ("epoch_ms", "mem_total_bytes", "mem_used_bytes", "swap_used_bytes", "bench_pss_bytes"):
        if column in mem.columns:
            mem[column] = pd.to_numeric(mem[column], errors="coerce")
    mem = mem.dropna(subset=["epoch_ms"])
    mem["ts"] = pd.to_datetime(mem["epoch_ms"], unit="ms")
    mem["n_accessions"] = n_accessions
    return mem.sort_values("ts")


def count_accessions(path):
    if not path or not os.path.exists(path):
        return 0
    with open(path) as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def shard_state(run_dir):
    if not run_dir:
        return {}
    state_dir = os.path.join(run_dir, "state")
    shards = glob.glob(os.path.join(run_dir, "shards", "shard_*.txt"))
    shards = [s for s in shards if not s.endswith(".todo.txt")]
    done = glob.glob(os.path.join(state_dir, "shard_*.done"))
    return {"shards_total": len(shards), "shards_done": len(done)}


def summarize_progress(halves, counts, failed, total, seq_types_expected, eta_window_hours):
    """Accession-level progress: an accession is done when every sequence type it has is
    finished (see TERMINAL_PROCESSES). Rate and ETA come from the completion times of the
    last eta_window_hours, so they describe current throughput rather than the run's
    average across every pause."""
    progress = {"total": total, "finished": 0, "failed": len(failed), "rate_per_hour": None,
                "eta_hours": None, "eta_time": None, "completions": pd.Series(dtype="datetime64[ns]")}
    if halves.empty:
        return progress

    per_accession = halves.groupby("accession").agg(seq_types=("seq_type", "nunique"),
                                                    finished_at=("complete_ts", "max"))
    finished = per_accession[per_accession["seq_types"] >= seq_types_expected]
    completions = finished["finished_at"].dropna().sort_values()
    progress["finished"] = int(len(finished))
    progress["completions"] = completions
    progress["in_progress"] = int(len(per_accession) - len(finished))

    if len(completions):
        window_start = completions.max() - pd.Timedelta(hours=eta_window_hours)
        recent = completions[completions >= window_start]
        span_hours = max((completions.max() - recent.min()).total_seconds() / 3600.0, 1e-9)
        # One completion in the window says nothing about a rate; fall back to the whole run.
        if len(recent) > 1 and span_hours > 1e-6:
            rate = len(recent) / span_hours
        else:
            whole = max((completions.max() - completions.min()).total_seconds() / 3600.0, 1e-9)
            rate = len(completions) / whole
        progress["rate_per_hour"] = rate
        remaining = max(total - progress["finished"], 0)
        if rate > 0 and total:
            progress["eta_hours"] = remaining / rate
            progress["eta_time"] = (datetime.utcnow() + timedelta(hours=remaining / rate)).strftime("%Y-%m-%d %H:%M UTC")
    return progress


def plot_progress(progress, counts, out_path, window_hours, rate_window_hours):
    """The plot a production run is watched on: how much of the manifest is done, how fast
    it is being done, and what the stages are doing right now."""
    completions = progress["completions"]
    total = progress["total"]

    fig = plt.figure(figsize=(12, 9))
    grid = fig.add_gridspec(3, 1, height_ratios=[1.3, 1.0, 1.0], hspace=0.45)
    ax1, ax2, ax3 = (fig.add_subplot(grid[i]) for i in range(3))

    # Cumulative completions against the size of the whole manifest.
    # A run this is watched on spans days; the first minutes of one span seconds. Pick the
    # time format from what is actually on screen so the ticks aren't all the same label.
    span_hours = 0.0
    if len(completions) > 1:
        span_hours = (completions.max() - completions.min()).total_seconds() / 3600.0
    time_format = mdates.DateFormatter("%H:%M:%S" if span_hours < 2 else "%m-%d %H:%M")

    if len(completions):
        cumulative = np.arange(1, len(completions) + 1)
        ax1.plot(completions.to_numpy(), cumulative, color="tab:purple", linewidth=2)
        ax1.xaxis.set_major_formatter(time_format)
    if total:
        ax1.axhline(total, color="tab:red", linestyle=":", linewidth=1.2,
                    label=f"manifest: {total:,} accessions")
        ax1.set_ylim(0, total * 1.05)
        ax1.legend(loc="upper left", fontsize=9)
    ax1.set_ylabel("Accessions complete")
    ax1.set_xlabel("Wall-clock time (UTC)")
    ax1.grid(True, linestyle="--", alpha=0.35)

    pct = (progress["finished"] / total * 100) if total else 0.0
    rate = progress["rate_per_hour"]
    rate_text = (f"{rate:,.0f} accessions/hour (last {rate_window_hours:g} h)"
                 if rate else "rate unknown yet")
    eta_text = "ETA unknown"
    if progress["eta_hours"] is not None:
        eta_text = f"ETA {progress['eta_hours'] / 24:,.1f} days ({progress['eta_time']})"
    ax1.set_title(f"{progress['finished']:,} of {total:,} accessions complete ({pct:.2f}%)   |   "
                  f"{rate_text}   |   {eta_text}\n"
                  f"{progress.get('in_progress', 0):,} accessions started but not finished, "
                  f"{progress['failed']:,} with a failed task", fontsize=11)

    # Completions per time bin over the recent window -- the shape of current throughput.
    if len(completions):
        cutoff = completions.max() - pd.Timedelta(hours=window_hours)
        recent = completions[completions >= cutoff]
        if len(recent) > 1:
            span_minutes = max((recent.max() - recent.min()).total_seconds() / 60.0, 1.0)
            bin_minutes = max(span_minutes / 60.0, 1.0)
            bin_width = pd.Timedelta(minutes=bin_minutes)
            edges = pd.date_range(recent.min(), recent.max() + bin_width, freq=bin_width)
            heights, _ = np.histogram(recent.to_numpy().astype("datetime64[ns]").astype("int64"),
                                      bins=edges.to_numpy().astype("datetime64[ns]").astype("int64"))
            centers = (edges[:-1] + bin_width / 2).to_pydatetime()
            ax2.bar(centers, heights, width=0.9 * bin_minutes / 1440.0, color="tab:blue")
            ax2.xaxis.set_major_formatter(time_format)
            ax2.set_xlabel("Wall-clock time (UTC)")
            ax2.set_ylabel(f"Completed / {bin_minutes:.0f} min")
    ax2.set_title(f"Accessions completed per time bin (last {window_hours:g} h)", fontsize=10)
    ax2.grid(True, axis="y", linestyle="--", alpha=0.35)

    # Per-stage task counts: what the pipeline is working on, and what has gone wrong.
    if not counts.empty:
        pivot = counts.pivot_table(index="process", columns="status", values="tasks",
                                   aggfunc="sum", fill_value=0)
        order = [s for s in ("COMPLETED", "RUNNING", "SUBMITTED", "CACHED", "FAILED", "ABORTED")
                 if s in pivot.columns]
        order += [c for c in pivot.columns if c not in order]
        pivot = pivot[order].sort_values(order[0], ascending=True)
        colors = {"COMPLETED": "tab:green", "RUNNING": "tab:blue", "SUBMITTED": "tab:cyan",
                  "CACHED": "tab:gray", "FAILED": "tab:red", "ABORTED": "tab:orange"}
        left = np.zeros(len(pivot))
        for status in pivot.columns:
            values = pivot[status].to_numpy(dtype=float)
            ax3.barh(pivot.index, values, left=left, label=status, color=colors.get(status, "tab:brown"))
            left += values
        ax3.legend(fontsize=8, ncol=len(pivot.columns), loc="lower right")
        ax3.set_xlabel("Tasks (whole run)")
    ax3.set_title("Tasks per stage, by status", fontsize=10)
    ax3.grid(True, axis="x", linestyle="--", alpha=0.35)

    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def render_with(function, produced_prefix, out_dir, target, *args, **kwargs):
    """plot_benchmark's per-batch plotters name their files <prefix>_<n_accessions>.<ext>.
    Run one and rename its output to a stable name the dashboard can link to."""
    before = set(os.listdir(out_dir)) if os.path.isdir(out_dir) else set()
    try:
        function(*args, **kwargs)
    except Exception as exc:                      # one unplottable figure must not stop the loop
        print(f"live_monitor: {function.__name__} failed: {exc}", file=sys.stderr)
        return
    for name in sorted(set(os.listdir(out_dir)) - before):
        if name.startswith(produced_prefix):
            extension = os.path.splitext(name)[1]
            os.replace(os.path.join(out_dir, name), os.path.join(out_dir, target + extension))


INDEX_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta http-equiv="refresh" content="{refresh}">
<title>Logan run progress</title>
<style>
  body {{ font-family: -apple-system, Segoe UI, Roboto, sans-serif; margin: 2rem auto; max-width: 1200px;
         color: #1d2220; background: #fbfbfa; }}
  h1 {{ font-size: 1.4rem; margin-bottom: .2rem; }}
  .meta {{ color: #59615c; font-size: .85rem; margin-bottom: 1.2rem; }}
  .stats {{ display: flex; flex-wrap: wrap; gap: 1.5rem; margin-bottom: 1.5rem; }}
  .stat {{ background: #fff; border: 1px solid #e3e6e3; border-radius: 8px; padding: .8rem 1.1rem; min-width: 9rem; }}
  .stat .value {{ font-size: 1.5rem; font-weight: 600; }}
  .stat .label {{ color: #59615c; font-size: .8rem; text-transform: uppercase; letter-spacing: .04em; }}
  .bar {{ height: 10px; background: #e3e6e3; border-radius: 5px; overflow: hidden; margin-bottom: 1.5rem; }}
  .bar > div {{ height: 100%; background: #6b5bd2; width: {pct:.3f}%; }}
  img {{ width: 100%; border: 1px solid #e3e6e3; border-radius: 8px; margin-bottom: 1.5rem; background: #fff; }}
  h2 {{ font-size: 1rem; margin: .3rem 0 .5rem; }}
  @media (prefers-color-scheme: dark) {{
    body {{ background: #14181a; color: #e7eae8; }}
    .stat {{ background: #1c2124; border-color: #2c3235; }}
    .bar {{ background: #2c3235; }}
    img {{ border-color: #2c3235; background: #fff; }}
    .meta, .stat .label {{ color: #9aa3a0; }}
  }}
</style>
</head>
<body>
<h1>Logan sketch / funcprofile &mdash; run in progress</h1>
<div class="meta">{run_dir} &middot; updated {updated} &middot; refreshes every {refresh}s</div>
<div class="stats">
  <div class="stat"><div class="value">{finished}</div><div class="label">accessions done</div></div>
  <div class="stat"><div class="value">{pct:.2f}%</div><div class="label">of {total}</div></div>
  <div class="stat"><div class="value">{rate}</div><div class="label">accessions / hour</div></div>
  <div class="stat"><div class="value">{eta_days}</div><div class="label">days remaining</div></div>
  <div class="stat"><div class="value">{shards}</div><div class="label">shards complete</div></div>
  <div class="stat"><div class="value">{failed}</div><div class="label">failed accessions</div></div>
</div>
<div class="bar"><div></div></div>
<h2>Progress, rate and per-stage task counts</h2>
<img src="00_progress.png?t={cache_bust}" alt="progress">
<h2>CPU and memory (last {window:g} h)</h2>
<img src="04_cpu_utilization.png?t={cache_bust}" alt="cpu utilization">
<h2>Task timeline (last {window:g} h)</h2>
<img src="05_task_timeline.png?t={cache_bust}" alt="task timeline">
<h2>Runtime distribution by stage (last {window:g} h)</h2>
<img src="06_runtime_distribution.png?t={cache_bust}" alt="runtime distribution">
<h2>Per-stage resource use (last {window:g} h)</h2>
<img src="09_task_resources.png?t={cache_bust}" alt="task resources">
</body>
</html>
"""


def write_index(out_dir, run_dir, progress, state, window_hours, refresh):
    total = progress["total"]
    pct = (progress["finished"] / total * 100) if total else 0.0
    rate = progress["rate_per_hour"]
    html = INDEX_TEMPLATE.format(
        refresh=refresh,
        run_dir=run_dir or out_dir,
        updated=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        finished=f"{progress['finished']:,}",
        total=f"{total:,}",
        pct=pct,
        rate=f"{rate:,.0f}" if rate else "&mdash;",
        eta_days=f"{progress['eta_hours'] / 24:,.1f}" if progress["eta_hours"] is not None else "&mdash;",
        shards=f"{state.get('shards_done', 0)} / {state.get('shards_total', 0)}" if state else "&mdash;",
        failed=f"{progress['failed']:,}",
        window=window_hours,
        cache_bust=int(time.time()),
    )
    with open(os.path.join(out_dir, "index.html"), "w") as handle:
        handle.write(html)


# Nextflow's own log is the only place a run in progress records tasks that have started but
# not finished: the trace gets a row only when a task completes, so every still-running task is
# invisible to it (see plot_benchmark.plot_cpu_utilization). Submission and completion are both
# timestamped and keyed by the work-directory hash, so the difference between the two sets is
# exactly what is in flight. Recovering it is what lets the CPU plot show true allocation right
# up to now, instead of a curve that sags toward zero precisely as long tasks pile up.
_SUBMITTED_RE = re.compile(r"^(\w{3}-\d{2} \d{2}:\d{2}:\d{2}\.\d{3})\s.*"
                           r"\[([0-9a-f]{2}/[0-9a-f]{6})\] Submitted process > (\S+)")
_COMPLETED_RE = re.compile(r"Task completed >.*workDir: (\S+)")

# How stale the log may be before unfinished submissions are read as killed rather than
# running: an interrupted shard leaves orphan submissions behind permanently.
_INFLIGHT_MAX_LOG_AGE_S = 300


def _log_ts_to_utc(stamps, ref_epoch):
    """Nextflow logs local wall-clock time with no year ('Sep-18 23:47:05.666'); the trace is
    UTC epoch milliseconds. Reconcile the two so they can share an axis. The year comes from
    the log's own mtime, stepping back one for any line that would otherwise land in the
    future (a run spanning New Year)."""
    ref = pd.Timestamp.fromtimestamp(ref_epoch)
    local = pd.to_datetime(f"{ref.year}-" + pd.Series(list(stamps), dtype="string"),
                           format="%Y-%b-%d %H:%M:%S.%f", errors="coerce")
    ahead = local > ref + pd.Timedelta(days=1)
    if ahead.any():
        local = local.where(~ahead, local - pd.DateOffset(years=1))
    return local - (pd.Timestamp.now() - pd.Timestamp.now("UTC").tz_convert(None)).round("s")


def read_inflight(run_dir, rows, now):
    """Tasks submitted but not yet completed, from the newest shard's Nextflow log, shaped like
    trace rows (start_ts from the submission, complete_ts = now) so they concatenate with it.
    Returns an empty frame when the log is missing, unparseable, or too stale to mean anything
    is still running -- callers then fall back to marking the region incomplete instead."""
    logs = sorted(glob.glob(os.path.join(run_dir, "logs", "nextflow.*.log")),
                  key=os.path.getmtime) if run_dir else []
    if not logs:
        return pd.DataFrame()
    path = logs[-1]
    mtime = os.path.getmtime(path)
    if time.time() - mtime > _INFLIGHT_MAX_LOG_AGE_S:
        return pd.DataFrame()

    submitted, completed = {}, set()
    try:
        with open(path, errors="replace") as handle:
            for line in handle:
                if "Submitted process > " in line:
                    match = _SUBMITTED_RE.match(line)
                    if match:
                        submitted[match.group(2)] = (match.group(1), match.group(3).split(":")[-1])
                elif "Task completed >" in line:
                    match = _COMPLETED_RE.search(line)
                    if match:
                        parts = match.group(1).rstrip("]").rstrip("/").split("/")
                        if len(parts) >= 2:
                            completed.add(f"{parts[-2]}/{parts[-1][:6]}")
    except OSError:
        return pd.DataFrame()

    live = [(stamp, proc) for key, (stamp, proc) in submitted.items() if key not in completed]
    if not live:
        return pd.DataFrame()

    df = pd.DataFrame(live, columns=["stamp", "process"])
    df["start_ts"] = _log_ts_to_utc(df["stamp"], mtime)
    df = df.dropna(subset=["start_ts"]).drop(columns=["stamp"])
    if df.empty:
        return df
    df["complete_ts"] = now
    # cpus/memory are per-process constants, so take them from what the trace has already
    # recorded for that process rather than re-reading the Nextflow config.
    declared = (rows.groupby("process")[["cpus", "memory_bytes"]].max()
                if not rows.empty else pd.DataFrame(columns=["cpus", "memory_bytes"]))
    df["cpus"] = df["process"].map(declared.get("cpus", pd.Series(dtype=float))).fillna(1.0)
    df["memory_bytes"] = df["process"].map(declared.get("memory_bytes", pd.Series(dtype=float))).fillna(0.0)
    # %cpu is only measured at completion. Leaving it NaN keeps these tasks out of the
    # "actually used" curve rather than inventing a figure for them -- which is why that curve
    # is labelled as covering completed tasks only.
    df["pct_cpu"] = float("nan")
    return df


def snapshot(args, reader):
    total = args.total or count_accessions(args.manifest)
    halves, counts, rows, failed = reader.load()

    seq_types_expected = 2
    if not halves.empty:
        seq_types_expected = max(int(halves["seq_type"].nunique()), 1)

    progress = summarize_progress(halves, counts, failed, total, seq_types_expected, args.eta_window_hours)
    state = shard_state(args.run_dir)

    os.makedirs(args.out_dir, exist_ok=True)
    plot_progress(progress, counts, os.path.join(args.out_dir, "00_progress.png"),
                  args.window_hours, args.eta_window_hours)

    if not rows.empty:
        rows = rows.copy()
        rows["n_accessions"] = total
        sysmem = load_system_memory(args.system_memory, total)
        now = pd.Timestamp.now("UTC").tz_convert(None)
        inflight = read_inflight(args.run_dir, rows, now)
        if inflight.empty:
            # Couldn't recover the running tasks, so the trace's own blind spot stands: every
            # moment within one longest-task of the newest completion may be missing tasks that
            # have not finished yet. Mark it rather than drawing it as measurement.
            longest = (rows["complete_ts"] - rows["start_ts"]).max()
            cpu_rows, incomplete_after = rows, (rows["complete_ts"].max() - longest
                                                if pd.notna(longest) else None)
        else:
            inflight["n_accessions"] = total
            cpu_rows, incomplete_after = pd.concat([rows, inflight], ignore_index=True), None
        render_with(pb.plot_cpu_utilization, "cpu_utilization", args.out_dir, "04_cpu_utilization",
                    cpu_rows, sysmem, args.out_dir, total_cpus=args.total_cpus,
                    memory_limit_gb=args.memory_limit_gb, incomplete_after=incomplete_after)
        colors = pb.stage_color_map(rows["process"].dropna().unique())
        # plot_benchmark's timeline samples the earliest-starting accessions, which for a run in
        # progress means the oldest thing on screen. Ask for the most recently started ones
        # instead: the timeline is here to show what the pipeline is working on now.
        latest = (rows.dropna(subset=["start_ts", "accession"])
                      .groupby("accession")["start_ts"].max()
                      .sort_values(ascending=False).head(args.timeline_sample_size).index.tolist())
        render_with(pb.plot_task_timeline, "task_timeline", args.out_dir, "05_task_timeline",
                    rows, args.out_dir, colors, sample_size=args.timeline_sample_size,
                    accessions=latest or None)
        render_with(pb.plot_runtime_distribution, "runtime_distribution", args.out_dir,
                    "06_runtime_distribution", rows, args.out_dir)
        render_with(pb.plot_task_resources, "task_resources", args.out_dir, "09_task_resources",
                    rows, args.out_dir)

    write_index(args.out_dir, args.run_dir, progress, state, args.window_hours, args.interval)

    payload = {
        "updated": datetime.now().isoformat(timespec="seconds"),
        "accessions_total": total,
        "accessions_finished": progress["finished"],
        "accessions_in_progress": progress.get("in_progress", 0),
        "accessions_failed": progress["failed"],
        "percent_complete": round((progress["finished"] / total * 100) if total else 0.0, 4),
        "rate_per_hour": progress["rate_per_hour"],
        "eta_hours": progress["eta_hours"],
        "eta_time_utc": progress["eta_time"],
        "shards_done": state.get("shards_done"),
        "shards_total": state.get("shards_total"),
        "failed_accessions_sample": failed[:50],
    }
    with open(os.path.join(args.out_dir, "progress.json"), "w") as handle:
        json.dump(payload, handle, indent=2)

    reader.forget_stale()
    rate = f"{progress['rate_per_hour']:,.0f}/h" if progress["rate_per_hour"] else "rate n/a"
    eta = f"ETA {progress['eta_hours'] / 24:,.1f} d" if progress["eta_hours"] is not None else "ETA n/a"
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {progress['finished']:,}/{total:,} accessions "
          f"({payload['percent_complete']:.2f}%), {rate}, {eta}", flush=True)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-dir", default=None,
                        help="A nextflow/run_full.sh run directory; supplies every path below")
    parser.add_argument("--trace", action="append", default=None,
                        help="Trace file or glob to read (repeatable). Default: <run-dir>/traces/*.tsv")
    parser.add_argument("--manifest", default=None,
                        help="Manifest whose accession count is the denominator "
                             "(default: <run-dir>/manifest.txt)")
    parser.add_argument("--total", type=int, default=None,
                        help="Total accessions, if there is no manifest to count")
    parser.add_argument("--system-memory", default=None,
                        help="system_memory.tsv from benchmark/log_system_memory.sh "
                             "(default: <run-dir>/system_memory.tsv)")
    parser.add_argument("--out-dir", default=None, help="Where to write plots (default: <run-dir>/plots_live)")
    parser.add_argument("--interval", type=float, default=60,
                        help="Seconds between refreshes, and the page's own refresh rate (default: 60)")
    parser.add_argument("--once", action="store_true", help="Render one snapshot and exit")
    parser.add_argument("--window-hours", type=float, default=24,
                        help="How much recent history the resource plots cover (default: 24)")
    parser.add_argument("--eta-window-hours", type=float, default=6,
                        help="Recent completions used for the rate and ETA (default: 6)")
    parser.add_argument("--timeline-sample-size", type=int, default=15,
                        help="Accessions shown on the task timeline (default: 15)")
    parser.add_argument("--total-cpus", type=int, default=None,
                        help="Reference line for server capacity on the CPU plot")
    parser.add_argument("--memory-limit-gb", type=float, default=None,
                        help="Reference line for executor.memory on the RAM panel")
    args = parser.parse_args()

    if args.run_dir:
        args.trace = args.trace or [os.path.join(args.run_dir, "traces", "*.tsv")]
        args.manifest = args.manifest or os.path.join(args.run_dir, "manifest.txt")
        args.system_memory = args.system_memory or os.path.join(args.run_dir, "system_memory.tsv")
        args.out_dir = args.out_dir or os.path.join(args.run_dir, "plots_live")
    if not args.trace:
        parser.error("pass --run-dir, or --trace with --total")
    args.out_dir = args.out_dir or "plots_live"

    os.makedirs(args.out_dir, exist_ok=True)
    reader = TraceReader(args.trace, args.window_hours)

    while True:
        try:
            snapshot(args, reader)
        except Exception as exc:                  # keep watching even if one cycle can't be drawn
            print(f"live_monitor: snapshot failed: {exc}", file=sys.stderr, flush=True)
        if args.once:
            return
        time.sleep(args.interval)


if __name__ == "__main__":
    main()
