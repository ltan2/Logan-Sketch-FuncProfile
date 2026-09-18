#!/usr/bin/env bash
# Samples server RAM/swap while a run is in progress, into a TSV that plot_benchmark.py
# and live_monitor.py read for the RAM panel of the CPU-utilization plot.
#
# The Nextflow trace only records per-task peak_rss, which can't show swapping or the
# memory other users of this shared server are holding -- that's what this adds.
#
# Columns:
#   epoch_ms          sample time
#   mem_total_bytes   whole server
#   mem_used_bytes    whole server (MemTotal - MemAvailable): this run + other users + OS
#   swap_used_bytes   whole server
#   bench_pss_bytes   this run's own share: summed PSS of every process in PGID's process
#                     group (Nextflow and all its tasks). PSS splits shared pages across
#                     processes, so unlike summed RSS it doesn't double-count.
#
# Usage:
#   log_system_memory.sh OUTFILE [PGID] [INTERVAL_SECONDS]
#
#   OUTFILE   TSV to write. Appended to (header written only when the file is new), so a
#             run that is paused and resumed keeps one continuous memory history.
#   PGID      process group whose PSS to attribute to this run (default: this script's own).
#   INTERVAL  seconds between samples (default: 15).
#
# Runs until killed. Callers are expected to kill it when the run ends, e.g.
#   ./log_system_memory.sh out.tsv "$(ps -o pgid= $$ | tr -d ' ')" & logger=$!
#   trap 'kill "$logger" 2>/dev/null' EXIT

# Tasks exit between pgrep and reading their smaps_rollup; that expected failure must not
# kill the logger, so this script deliberately does not use errexit/pipefail.
set +e +o pipefail

outfile=${1:?"Usage: $0 OUTFILE [PGID] [INTERVAL_SECONDS]"}
pgid=${2:-$(ps -o pgid= $$ | tr -d ' ')}
interval=${3:-15}

mkdir -p "$(dirname "$outfile")"
if [ ! -s "$outfile" ]; then
    printf 'epoch_ms\tmem_total_bytes\tmem_used_bytes\tswap_used_bytes\tbench_pss_bytes\n' > "$outfile"
fi

while true; do
    t=$(date +%s%3N)
    bench_kb=$(for p in $(pgrep -g "$pgid"); do cat "/proc/$p/smaps_rollup" 2>/dev/null; done \
        | awk '/^Pss:/{s+=$2} END{print s+0}')
    awk -v t="$t" -v b="$bench_kb" '/^MemTotal:/{mt=$2} /^MemAvailable:/{ma=$2}
        /^SwapTotal:/{st=$2} /^SwapFree:/{sf=$2}
        END{printf "%s\t%.0f\t%.0f\t%.0f\t%.0f\n", t, mt*1024, (mt-ma)*1024, (st-sf)*1024, b*1024}' /proc/meminfo >> "$outfile"
    sleep "$interval"
done
