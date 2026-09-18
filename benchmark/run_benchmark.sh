#!/usr/bin/env bash
# Benchmarks the workflow in nextflow/ at increasing accession-list sizes.
#
# Batches are nested: we draw ONE random sample at the largest requested size, and every
# smaller size is a prefix of that same sample (benchmark/sample_accessions.py, seeded by
# the largest size so it's reproducible). This guarantees size-100 accessions are a subset
# of size-500's, which are a subset of size-1000's, etc. Sampling each size independently
# would let a slow/oversized accession land in a small batch by pure luck (and not in a
# larger one), making a small batch look artificially slower than a bigger one -- nesting
# means every batch is exposed to exactly the same accessions the smaller batches were,
# plus more, so wall-clock time should scale monotonically with N.
#
# Before any batch, the pipeline unit test (nextflow/test/run_unit_test.sh) must pass; the benchmark
# stops if it doesn't.
#
# For each size N, this:
#   1. takes the first N accessions of the largest-size sample (see above)
#   2. runs the full pipeline against that sample, into its own --outdir
#   3. times the whole `nextflow run` invocation (this is the "entire workflow" number)
#   4. lets Nextflow's own trace (enabled in nextflow/nextflow.config) record per-task
#      timing for every process into <outdir>/pipeline_trace.tsv -- that TSV is the
#      "time for each process" record; benchmark/plot_benchmark.py parses it.
#
# Usage:
#   ./benchmark/run_benchmark.sh [SIZE...]
# Examples:
#   ./benchmark/run_benchmark.sh                      # default: 1000 10000 50000 100000
#   ./benchmark/run_benchmark.sh 1000 10000
# Optional environment:
#   MAX_COMPRESSED_GB=2   only sample accessions whose compressed unitigs/contigs are <= 2 GB
#   NF_EXTRA_ARGS='...'   extra `nextflow run` arguments for every batch and the unit test
#   SKIP_UNIT_TEST=1      skip the unit test gate
#

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(dirname "$script_dir")"


# Run inside a detached tmux session so the sweep survives that. Opt out with NO_TMUX=1
if [ -z "${NO_TMUX:-}" ] && [ -z "${TMUX:-}" ]; then
    if command -v tmux >/dev/null 2>&1; then
        session="run_benchmark_$(date +%Y%m%d_%H%M%S)"
        extra_env=""
        if [ -n "${MAX_COMPRESSED_GB:-}" ]; then
            extra_env="$extra_env MAX_COMPRESSED_GB=$(printf '%q' "$MAX_COMPRESSED_GB")"
        fi
        if [ -n "${SKIP_UNIT_TEST:-}" ]; then
            extra_env="$extra_env SKIP_UNIT_TEST=$(printf '%q' "$SKIP_UNIT_TEST")"
        fi
        if [ -n "${SIZE_CSV:-}" ]; then
            extra_env="$extra_env SIZE_CSV=$(printf '%q' "$SIZE_CSV")"
        fi
        # Extra `nextflow run` arguments, e.g. NF_EXTRA_ARGS='--funprofiler_seq_types unitigs'
        if [ -n "${NF_EXTRA_ARGS:-}" ]; then
            extra_env="$extra_env NF_EXTRA_ARGS=$(printf '%q' "$NF_EXTRA_ARGS")"
        fi
        tmux new-session -d -s "$session" \
            "NO_TMUX=1 PATH=$(printf '%q' "$PATH") ACCESSIONS_SOURCE=$(printf '%q' "${ACCESSIONS_SOURCE:-}") BENCH_DIR=$(printf '%q' "${BENCH_DIR:-}")$extra_env $(printf '%q ' "$0" "$@")"
        echo "Started in detached tmux session '$session' -- this will keep running if your terminal/SSH session drops."
        echo "Attach to watch progress:  tmux attach -t $session"
        echo "Detach without stopping it: Ctrl-b d"
        exit 0
    else
        echo "WARNING: tmux not found -- running in the foreground. A dropped terminal/SSH" >&2
        echo "connection will kill this run partway through a sweep. Install tmux, or run" >&2
        echo "this script under 'nohup ... &' / a batch scheduler instead." >&2
    fi
fi

if [ "$#" -gt 0 ]; then
    sizes=("$@")
else
    sizes=(100 500 1000 5000)
fi

# Default output lives on /scratch, not the repo's own (much smaller) disk -- a sweep at
# n=100000 can run into hundreds of GB of intermediate decompressed FASTA/sketch/work-dir
# data. Scoped by username so concurrent users of this shared box don't collide.
bench_dir="${BENCH_DIR:-/scratch/$(whoami)/logan_sketch_funcprofile_benchmark_runs}"
mkdir -p "$bench_dir"

# Guard against two invocations racing on the same BENCH_DIR: each nextflow run below is
# unthrottled enough (executor pool sized to the whole machine, no cleanup of intermediate
# task dirs) that two of them running concurrently against the same disk can exhaust it in
# minutes -- and both would stomp the same wall_clock_summary.csv while they're at it.
lockfile="$bench_dir/.run_benchmark.lock"
exec 9>"$lockfile"
if ! flock -n 9; then
    echo "ERROR: another run_benchmark.sh is already using BENCH_DIR=$bench_dir" >&2
    echo "(lock held on $lockfile). Wait for it to finish, or set BENCH_DIR to a different directory." >&2
    exit 1
fi

accessions_source="${ACCESSIONS_SOURCE:-$repo_root/wgs_metagenome_accessions.txt}"

# Gate: the full pipeline must reproduce known-good results for one accession before any benchmark
# batch runs (nextflow/test/run_unit_test.sh -- same NF_EXTRA_ARGS as the batches below). Skip only
# deliberately, with SKIP_UNIT_TEST=1.
if [ -n "${SKIP_UNIT_TEST:-}" ]; then
    echo "=== [$(date '+%Y-%m-%d %H:%M:%S')] SKIP_UNIT_TEST set -- skipping the pipeline unit test ==="
else
    echo "=== [$(date '+%Y-%m-%d %H:%M:%S')] running pipeline unit test ==="
    if ! "$repo_root/nextflow/test/run_unit_test.sh"; then
        echo "ERROR: pipeline unit test failed -- not running the benchmark." >&2
        exit 1
    fi
fi

summary_csv="$bench_dir/wall_clock_summary.csv"
echo "n_accessions,wall_clock_seconds,outdir,trace_file" > "$summary_csv"

# Batches sample from every accession by default, so benchmark timings reflect the real size mix
# (including the ~3% whose unitigs/contigs .fa.zst is over 2 GB on S3). Optionally set
# MAX_COMPRESSED_GB (e.g. 2) to exclude accessions above that size, per
# benchmark/query_accession_sizes.py's survey -- useful for a quick check that shouldn't be held
# up by a single outlier like SRR3507924 (10.4 GB compressed / 45 GB decompressed).
size_filter_args=()
if [ -n "${MAX_COMPRESSED_GB:-}" ]; then
    size_filter_args+=(--max-compressed-gb "$MAX_COMPRESSED_GB")
    if [ -n "${SIZE_CSV:-}" ]; then
        size_filter_args+=(--size-csv "$SIZE_CSV")
    fi
fi

# Draw the nested pool once, at the largest requested size. Every smaller size below is
# just a prefix of this same file (see the block comment above for why).
max_n="$(printf '%s\n' "${sizes[@]}" | sort -n | tail -1)"
pool_file="$bench_dir/accessions_pool_${max_n}.txt"
echo "=== [$(date '+%Y-%m-%d %H:%M:%S')] sampling nested pool of ${max_n} accessions ==="
python3 "$script_dir/sample_accessions.py" \
    --source "$accessions_source" \
    --n "$max_n" --seed "$max_n" --out "$pool_file" \
    "${size_filter_args[@]}"

for n in "${sizes[@]}"; do
    accession_file="$bench_dir/accessions_${n}.txt"
    outdir="$bench_dir/results_${n}"
    work_dir="$bench_dir/work_${n}"
    trace_file="${outdir}/pipeline_trace.tsv"   # written by nextflow/nextflow.config

    head -n "$n" "$pool_file" > "$accession_file"

    echo "=== [$(date '+%Y-%m-%d %H:%M:%S')] running pipeline for n=${n} ==="
    start_ts=$(date +%s.%N)

    # RAM/swap every 15s while this batch runs -- the Nextflow trace only has per-task peaks,
    # which can't show swapping or other users' memory (see the RAM panel in plot_benchmark.py).
    # mem_used/swap_used are whole-server; bench_pss is this benchmark's own share: the summed
    # PSS of every process in this script's process group (Nextflow and all its tasks). PSS
    # splits shared pages across processes, so unlike summed RSS it doesn't double-count.
    mkdir -p "$outdir"
    pgid=$(ps -o pgid= $$ | tr -d ' ')
    (
        # Tasks exit between pgrep and reading their smaps_rollup; that expected failure must
        # not kill the logger via the errexit/pipefail this subshell inherits from the script.
        set +e +o pipefail
        printf 'epoch_ms\tmem_total_bytes\tmem_used_bytes\tswap_used_bytes\tbench_pss_bytes\n'
        while true; do
            t=$(date +%s%3N)
            bench_kb=$(for p in $(pgrep -g "$pgid"); do cat "/proc/$p/smaps_rollup" 2>/dev/null; done \
                | awk '/^Pss:/{s+=$2} END{print s+0}')
            awk -v t="$t" -v b="$bench_kb" '/^MemTotal:/{mt=$2} /^MemAvailable:/{ma=$2}
                /^SwapTotal:/{st=$2} /^SwapFree:/{sf=$2}
                END{printf "%s\t%.0f\t%.0f\t%.0f\t%.0f\n", t, mt*1024, (mt-ma)*1024, (st-sf)*1024, b*1024}' /proc/meminfo
            sleep 15
        done
    ) > "$outdir/system_memory.tsv" &
    mem_logger_pid=$!
    trap 'kill "$mem_logger_pid" 2>/dev/null' EXIT

    # runs nextflow workflow that is calling main.nf, get unitig/contig, run analysis - sketch and funprofiler
    nextflow run "$repo_root/nextflow" \
        -profile benchmark \
        --accessions "$accession_file" \
        --outdir "$outdir" \
        -work-dir "$work_dir" \
        ${NF_EXTRA_ARGS:-}

    kill "$mem_logger_pid" 2>/dev/null
    end_ts=$(date +%s.%N)
    wall=$(awk -v a="$start_ts" -v b="$end_ts" 'BEGIN { printf "%.3f", b - a }')

    echo "${n},${wall},${outdir},${trace_file}" >> "$summary_csv"
    echo "=== n=${n} took ${wall}s (trace: ${trace_file}) ==="
done

echo
echo "Done. Wall-clock summary: $summary_csv"
echo "Next: python3 ${script_dir}/plot_benchmark.py --bench-dir ${bench_dir}"
