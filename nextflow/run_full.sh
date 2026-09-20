#!/usr/bin/env bash
# Runs the whole accession manifest through the pipeline, shard by shard, so the run can be
# paused and resumed without losing finished work. See ../how_to_run.md for the full
# procedure; this header documents the contract the script implements.
#
# Usage:
#   nextflow/run_full.sh [options]
#
# Options:
#   --accessions FILE   manifest to run          (default: <repo>/wgs_metagenome_accessions.txt)
#   --run-dir DIR       run state + work + results
#                                                (default: /scratch/$USER/logan_full_run)
#   --outdir DIR        published results        (default: <run-dir>/results)
#   --shard-size N      accessions per shard     (default: 10000)
#   --size-csv FILE     accession size survey (benchmark/query_accession_sizes.py output) used
#                       to balance shards and put the biggest accessions first within each one.
#                       Defaults to <repo>/accession_size_analysis/accession_sizes.csv if that
#                       exists; without it, shards are plain sequential chunks.
#   --shards N          run at most N shards this invocation, then stop (default: all)
#   --keep-work         keep each completed shard's work directory (default: delete it)
#   --no-monitor        don't start benchmark/live_monitor.py alongside the run
#   --dry-run           print what would run, change nothing
#   -h, --help          this message
#
# Environment:
#   SKIP_UNIT_TEST=1    skip the nextflow/test/run_unit_test.sh gate (not recommended)
#   NF_EXTRA_ARGS='...' extra `nextflow run` arguments, also passed to the unit test gate
#   NO_TMUX=1           run in the foreground instead of re-launching in a detached tmux session
#
# PAUSING
#   Graceful, at a shard boundary (preferred):
#       touch <run-dir>/PAUSE
#     The current shard runs to completion, then the driver stops. Nothing is lost or repeated.
#   Immediately, mid-shard:
#       Ctrl-C in the driver's tmux session, or  kill -INT <driver pid>
#     Nextflow stops and its in-flight tasks are discarded; every task that had already
#     finished has already published its results and is skipped on resume.
#   For a few minutes only (e.g. to hand the machine's cores to someone else):
#       kill -STOP -- -<driver pgid>   then   kill -CONT -- -<driver pgid>
#     This freezes the processes in place. Don't leave it stopped for long: tasks keep their
#     scratch files and open S3 connections, which can time out.
#
# RESUMING
#   Re-run this script with the same --run-dir (and remove <run-dir>/PAUSE if you created it).
#   Shards already marked complete in <run-dir>/state/ are skipped; the shard that was
#   interrupted is re-run with the accessions it already published filtered out of it
#   (nextflow/filter_completed.py).
#
#   Note it does NOT rely on `nextflow -resume`. DECOMPRESS_FASTA and CLEANUP_FASTA delete
#   their inputs as they go to stay inside the disk budget, and Nextflow invalidates any
#   cached task whose output files are gone -- so FETCH_LOGAN is invalidated and the whole
#   chain below it re-runs. Measured on the 5-accession test profile: a second run with
#   -resume reported `cached=0`. Published results are what makes this run resumable, not the
#   work cache.
#
# LAYOUT of --run-dir
#   manifest.txt                frozen copy of the manifest this run is sharded from
#   manifest.sha256             guards against the manifest changing under a resumed run
#   shards/shard_NNNNN.txt      the shards themselves (stable: never re-sharded), size-balanced
#                               and ordered heaviest-first when a size survey is available
#   shards/shard_NNNNN.todo.txt what the latest attempt at that shard actually ran
#   state/shard_NNNNN.done      shard finished (one line: when, and how many accessions)
#   state/shard_NNNNN.failed.txt  accessions with a FAILED task in that shard, from its trace
#   state/progress.tsv          one row per shard attempt: shard, accessions, start, end, status
#   logs/shard_NNNNN.log        console output of that shard's `nextflow run`
#   logs/nextflow.shard_NNNNN.log  Nextflow's own debug log for that shard
#   traces/pipeline_trace.shard_NNNNN[.aN].tsv  per-shard trace (input to the live plots)
#   work/shard_NNNNN/           Nextflow work dir, deleted once the shard completes
#   results/                    published results for every shard (one tree for the whole run)
#   system_memory.tsv           server RAM/swap samples for the whole run
#   plots_live/                 live plots, refreshed by benchmark/live_monitor.py
#   PAUSE                       create this file to stop after the current shard
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(dirname "$script_dir")"

accessions="$repo_root/wgs_metagenome_accessions.txt"
run_dir="/scratch/$(whoami)/logan_full_run"
outdir=""
shard_size=10000
size_csv=""
max_shards=0
keep_work=0
monitor=1
dry_run=0

# Keep a pristine copy: the parsing loop below consumes $@ with `shift`, and the tmux
# relaunch further down has to hand the session exactly what the user typed. Passing "$@"
# there silently relaunched with no arguments at all -- every flag reverted to its default.
original_args=("$@")

while [ "$#" -gt 0 ]; do
    case "$1" in
        --accessions) accessions="$2"; shift 2 ;;
        --run-dir)    run_dir="$2";    shift 2 ;;
        --outdir)     outdir="$2";     shift 2 ;;
        --shard-size) shard_size="$2"; shift 2 ;;
        --size-csv)   size_csv="$2";   shift 2 ;;
        --shards)     max_shards="$2"; shift 2 ;;
        --keep-work)  keep_work=1;     shift ;;
        --no-monitor) monitor=0;       shift ;;
        --dry-run)    dry_run=1;       shift ;;
        -h|--help)    sed -n '2,72p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "Unknown option: $1 (try --help)" >&2; exit 2 ;;
    esac
done

# A full run takes days to weeks, so it must not die with the terminal that started it.
# Re-launch into a detached tmux session unless we're already in one (or NO_TMUX is set).
if [ -z "${NO_TMUX:-}" ] && [ -z "${TMUX:-}" ] && [ "$dry_run" -eq 0 ]; then
    if command -v tmux >/dev/null 2>&1; then
        session="logan_full_run"
        if tmux has-session -t "$session" 2>/dev/null; then
            echo "ERROR: tmux session '$session' already exists -- attach with: tmux attach -t $session" >&2
            exit 1
        fi
        tmux new-session -d -s "$session" \
            "NO_TMUX=1 PATH=$(printf '%q' "$PATH") SKIP_UNIT_TEST=$(printf '%q' "${SKIP_UNIT_TEST:-}") NF_EXTRA_ARGS=$(printf '%q' "${NF_EXTRA_ARGS:-}") $(printf '%q ' "${BASH_SOURCE[0]}" ${original_args[@]+"${original_args[@]}"})"
        echo "Started in detached tmux session '$session' -- it keeps running if your SSH session drops."
        echo "  Watch it:            tmux attach -t $session      (detach again with Ctrl-b d)"
        echo "  Pause after shard:   touch $run_dir/PAUSE"
        echo "  Pause immediately:   attach and press Ctrl-C"
        exit 0
    fi
    echo "WARNING: tmux not found -- running in the foreground. A dropped SSH connection will" >&2
    echo "interrupt the run (which is resumable, but you'd have to notice). Install tmux, or" >&2
    echo "run this under nohup." >&2
fi

[ -f "$accessions" ] || { echo "ERROR: accession manifest not found: $accessions" >&2; exit 1; }
outdir="${outdir:-$run_dir/results}"

mkdir -p "$run_dir"/{shards,state,logs,traces,work} "$outdir"

# One driver per run directory. Two would race on the same shard state, the same results tree
# and the same disk.
exec 9>"$run_dir/.run_full.lock"
if ! flock -n 9; then
    echo "ERROR: another run_full.sh is already running against $run_dir" >&2
    exit 1
fi

log() { echo "=== [$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ---------------------------------------------------------------------------------------
# Freeze the manifest. Sharding has to stay stable across resumes: if the manifest changed,
# shard N would no longer mean the same accessions the state directory says are done.
# ---------------------------------------------------------------------------------------
frozen="$run_dir/manifest.txt"
if [ ! -f "$frozen" ]; then
    # Nothing here yet: a fresh run. Under --dry-run, report what this would become and stop.
    if [ "$dry_run" -eq 1 ]; then
        n_accessions="$(grep -v '^[[:space:]]*#' "$accessions" | awk 'NF' | wc -l)"
        log "dry run: would freeze $n_accessions accessions from $accessions into $frozen"
        log "dry run: would split them into $(( (n_accessions + shard_size - 1) / shard_size )) shards of up to $shard_size accessions"
        log "dry run: would publish results into $outdir"
        exit 0
    fi
    grep -v '^[[:space:]]*#' "$accessions" | awk 'NF' > "$frozen"
    sha256sum "$frozen" | awk '{print $1}' > "$run_dir/manifest.sha256"
    log "froze $(wc -l < "$frozen") accessions from $accessions into $frozen"
elif [ -f "$run_dir/manifest.sha256" ]; then
    if [ "$(sha256sum "$frozen" | awk '{print $1}')" != "$(cat "$run_dir/manifest.sha256")" ]; then
        echo "ERROR: $frozen changed since this run was sharded (checksum mismatch)." >&2
        echo "Shard numbering would no longer match $run_dir/state/. Start a new --run-dir." >&2
        exit 1
    fi
fi

# ---------------------------------------------------------------------------------------
# Shard once, then never again -- the shard files are this run's unit of progress.
# ---------------------------------------------------------------------------------------
if ! compgen -G "$run_dir/shards/shard_*.txt" > /dev/null; then
    # A shard's wall clock is its longest accession, so the shards are built by
    # nextflow/shard_manifest.py rather than `split`: balanced by compressed size, and ordered
    # heaviest-first inside each shard so the multi-hour accessions start at t=0 instead of
    # stranding the machine at the end. Falls back to sequential chunks with no size survey.
    if [ -z "$size_csv" ] && [ -f "$repo_root/accession_size_analysis/accession_sizes.csv" ]; then
        size_csv="$repo_root/accession_size_analysis/accession_sizes.csv"
    fi
    if [ -n "$size_csv" ]; then
        [ -f "$size_csv" ] || { echo "ERROR: --size-csv not found: $size_csv" >&2; exit 1; }
        log "sharding by size using $size_csv"
        python3 "$repo_root/nextflow/shard_manifest.py" --accessions "$frozen" \
            --out-dir "$run_dir/shards" --shard-size "$shard_size" --size-csv "$size_csv"
    else
        log "no accession size survey found -- shards will be sequential chunks, which lets a"
        log "  large accession land late in a shard and idle the machine (see how_to_run.md)."
        python3 "$repo_root/nextflow/shard_manifest.py" --accessions "$frozen" \
            --out-dir "$run_dir/shards" --shard-size "$shard_size"
    fi
    log "$(ls "$run_dir"/shards/shard_*.txt | grep -vc '\.todo\.txt$') shards of up to $shard_size accessions"
fi
shards=()
while IFS= read -r line; do shards+=("$line"); done < <(ls "$run_dir"/shards/shard_*.txt 2>/dev/null | grep -v '\.todo\.txt$' | sort)
[ "${#shards[@]}" -gt 0 ] || { echo "ERROR: no shards found under $run_dir/shards/" >&2; exit 1; }

done_count=0
for shard in "${shards[@]}"; do
    tag="$(basename "$shard" .txt)"
    [ -f "$run_dir/state/$tag.done" ] && done_count=$((done_count + 1))
done
log "$done_count of ${#shards[@]} shards already complete"

if [ "$dry_run" -eq 1 ]; then
    log "dry run: would run $(( ${#shards[@]} - done_count )) shards into $outdir"
    exit 0
fi

# ---------------------------------------------------------------------------------------
# Gate: the pipeline must still reproduce known-good results before it processes a manifest
# this size. Same gate benchmark/run_benchmark.sh uses, with the same NF_EXTRA_ARGS.
# ---------------------------------------------------------------------------------------
if [ -n "${SKIP_UNIT_TEST:-}" ]; then
    log "SKIP_UNIT_TEST set -- skipping the pipeline unit test"
else
    log "running pipeline unit test"
    if ! "$repo_root/nextflow/test/run_unit_test.sh"; then
        echo "ERROR: pipeline unit test failed -- not starting the full run." >&2
        exit 1
    fi
fi

# ---------------------------------------------------------------------------------------
# Background helpers: server memory sampling (for the RAM panel of the live CPU plot) and
# the live plots themselves. Both are stopped when this script exits, for any reason.
# ---------------------------------------------------------------------------------------
helpers=()
stop_helpers() { for pid in "${helpers[@]:-}"; do [ -n "$pid" ] && kill "$pid" 2>/dev/null; done; }
trap stop_helpers EXIT

"$repo_root/benchmark/log_system_memory.sh" "$run_dir/system_memory.tsv" \
    "$(ps -o pgid= $$ | tr -d ' ')" 15 &
helpers+=("$!")

if [ "$monitor" -eq 1 ]; then
    python3 "$repo_root/benchmark/live_monitor.py" --run-dir "$run_dir" \
        >> "$run_dir/logs/live_monitor.log" 2>&1 &
    helpers+=("$!")
    log "live plots: $run_dir/plots_live/index.html (log: $run_dir/logs/live_monitor.log)"
fi

# ---------------------------------------------------------------------------------------
# Pause handling. A signal stops the shard that's running; the PAUSE file stops the driver
# at the next shard boundary. Either way the next invocation picks up where this left off.
# ---------------------------------------------------------------------------------------
interrupted=0
nf_pid=""
on_signal() {
    [ "$interrupted" -eq 1 ] && return          # already shutting down; don't stack signals
    interrupted=1
    echo
    log "interrupt received -- stopping Nextflow. In-flight tasks are discarded; finished ones are kept."
    [ -n "$nf_pid" ] && kill -INT "$nf_pid" 2>/dev/null
}
trap on_signal INT TERM

# Waits out a Nextflow JVM's shutdown, escalating if it doesn't stop. `pattern` is matched
# against full command lines, so it must be something only this attempt's JVM carries.
await_nextflow() {
    local pattern="$1"
    local deadline=$(( $(date +%s) + 600 ))
    local escalated=0
    while pgrep -f "$pattern" > /dev/null 2>&1; do
        if [ "$(date +%s)" -ge "$deadline" ]; then
            if [ "$escalated" -eq 0 ]; then
                log "Nextflow still running 10 minutes after being asked to stop -- sending SIGTERM"
                pkill -TERM -f "$pattern" || true
                escalated=1
                deadline=$(( $(date +%s) + 120 ))
            else
                log "Nextflow still running after SIGTERM -- sending SIGKILL"
                pkill -KILL -f "$pattern" || true
                deadline=$(( $(date +%s) + 60 ))
            fi
        fi
        sleep 2
    done
}

progress_tsv="$run_dir/state/progress.tsv"
[ -s "$progress_tsv" ] || printf 'shard\taccessions_run\tstarted\tended\tseconds\tstatus\n' > "$progress_tsv"

# Which published files the resume filter should expect per accession depends on two run
# parameters. Take them from the same place the run itself will: NF_EXTRA_ARGS if it overrides
# them, otherwise the config -- the same precedence nextflow/test/run_unit_test.sh uses. Getting
# this wrong would make finished accessions look unfinished and re-run them on every resume.
nf_config="$( { cd "$run_dir" && nextflow config "$repo_root/nextflow" -profile full -flat 2>/dev/null; } || true)"
config_value() {
    printf '%s\n' "$nf_config" | awk -F" = " -v k="params.$1" '$1==k {gsub(/\x27/, "", $2); print $2}'
}
override() {
    printf '%s\n' "${NF_EXTRA_ARGS:-}" | { grep -oE -- "--$1[ =][^ ]+" || true; } | sed -E "s/--$1[ =]//" | tail -1
}
param() {
    local value; value="$(override "$1")"
    [ -n "$value" ] || value="$(config_value "$1")"
    printf '%s' "$value"
}
funprofiler_seq_types="$(param funprofiler_seq_types)"; funprofiler_seq_types="${funprofiler_seq_types:-unitigs,contigs}"
sourmash_ksize="$(param sourmash_ksize)";             sourmash_ksize="${sourmash_ksize:-31}"
log "resume filter uses funprofiler_seq_types=$funprofiler_seq_types, sourmash_ksize=$sourmash_ksize"

ran=0
for shard in "${shards[@]}"; do
    tag="$(basename "$shard" .txt)"

    [ -f "$run_dir/state/$tag.done" ] && continue

    if [ -f "$run_dir/PAUSE" ]; then
        log "PAUSE file present -- stopping before $tag. Resume with the same command once you remove it:"
        echo "    rm $run_dir/PAUSE && $0 --run-dir $run_dir"
        exit 0
    fi
    [ "$interrupted" -eq 1 ] && break
    if [ "$max_shards" -gt 0 ] && [ "$ran" -ge "$max_shards" ]; then
        log "--shards $max_shards reached -- stopping. Re-run the same command to continue."
        exit 0
    fi

    # Accession-level resume: whatever this shard already published is not run again.
    todo="$run_dir/shards/$tag.todo.txt"
    python3 "$repo_root/nextflow/filter_completed.py" \
        --accessions "$shard" --outdir "$outdir" --out "$todo" \
        --funprofiler-seq-types "$funprofiler_seq_types" --ksize "$sourmash_ksize"
    n_todo="$(wc -l < "$todo")"
    if [ "$n_todo" -eq 0 ]; then
        echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) 0 accessions left to run" > "$run_dir/state/$tag.done"
        log "$tag was already fully published -- marked complete"
        continue
    fi

    # Each attempt gets its own trace so the live plots keep the whole history of the run.
    # Counted by globbing rather than `ls | wc -l`: under `set -o pipefail` an `ls` that
    # matches nothing would take the whole driver down with it.
    attempt=1
    for previous_trace in "$run_dir/traces/pipeline_trace.$tag"*.tsv; do
        [ -e "$previous_trace" ] && attempt=$((attempt + 1))
    done
    run_tag="$tag"
    [ "$attempt" -gt 1 ] && run_tag="$tag.a$attempt"
    work_dir="$run_dir/work/$tag"
    # A previous attempt's work dir is not resumed from (see the header), so it is only taking
    # up disk. Its failure evidence -- the shard log and trace -- is kept.
    [ -d "$work_dir" ] && { log "removing $work_dir left by a previous attempt at $tag"; rm -rf "$work_dir"; }

    log "$tag: running $n_todo accessions (attempt $attempt) -- log: $run_dir/logs/$tag.log"
    start_ts=$(date +%s)
    started="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

    set +e
    (
        cd "$run_dir"
        # shellcheck disable=SC2086  # NF_EXTRA_ARGS is intentionally word-split
        nextflow -log "$run_dir/logs/nextflow.$run_tag.log" \
            run "$repo_root/nextflow" \
            -profile full \
            --accessions "$todo" \
            --outdir "$outdir" \
            --trace_dir "$run_dir/traces" \
            --run_tag "$run_tag" \
            -work-dir "$work_dir" \
            -ansi-log false \
            ${NF_EXTRA_ARGS:-}
    ) >> "$run_dir/logs/$tag.log" 2>&1 &
    nf_pid=$!
    wait "$nf_pid"
    rc=$?
    # Neither of the two things that just happened means "Nextflow has stopped": bash's `wait`
    # returns as soon as our signal handler runs, and the `nextflow` command is a launcher
    # script that forwards the signal to the JVM and exits before that JVM has finished killing
    # its tasks. Exiting here would release the driver's lock -- and invite a resume -- while
    # tasks are still writing into this run directory. So wait for the JVM itself, identified
    # by the -log path unique to this attempt.
    await_nextflow "$run_dir/logs/nextflow.$run_tag.log"
    set -e
    nf_pid=""

    end_ts=$(date +%s)
    ended="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

    # Accessions the shard couldn't finish (FETCH_LOGAN gave up under the 'full' profile, or a
    # task failed). Recovered from the trace so they're recorded rather than silently dropped.
    trace_file="$run_dir/traces/pipeline_trace.$run_tag.tsv"
    if [ -f "$trace_file" ]; then
        awk -F'\t' 'NR > 1 && $3 == "FAILED" {split($2, a, "."); print a[1]}' "$trace_file" \
            | sort -u > "$run_dir/state/$tag.failed.txt"
        [ -s "$run_dir/state/$tag.failed.txt" ] || rm -f "$run_dir/state/$tag.failed.txt"
    fi

    if [ "$interrupted" -eq 1 ]; then
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$run_tag" "$n_todo" "$started" "$ended" "$((end_ts - start_ts))" "interrupted" >> "$progress_tsv"
        log "$tag interrupted. Resume with:  $0 --run-dir $run_dir"
        exit 130
    fi

    if [ "$rc" -ne 0 ]; then
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$run_tag" "$n_todo" "$started" "$ended" "$((end_ts - start_ts))" "failed" >> "$progress_tsv"
        log "$tag FAILED (nextflow exit $rc). The driver stops here so the failure is visible."
        echo "  Shard log:     $run_dir/logs/$tag.log"
        echo "  Nextflow log:  $run_dir/logs/nextflow.$run_tag.log"
        echo "  Work dir kept: $work_dir"
        echo "  After fixing the cause, re-run:  $0 --run-dir $run_dir"
        exit "$rc"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$run_tag" "$n_todo" "$started" "$ended" "$((end_ts - start_ts))" "done" >> "$progress_tsv"
    echo "$ended $n_todo accessions run in $((end_ts - start_ts))s (attempt $attempt)" > "$run_dir/state/$tag.done"
    log "$tag complete in $((end_ts - start_ts))s"

    if [ "$keep_work" -eq 0 ]; then
        rm -rf "$work_dir"
    fi
    ran=$((ran + 1))
done

if [ "$interrupted" -eq 1 ]; then
    log "stopped by interrupt. Resume with:  $0 --run-dir $run_dir"
    exit 130
fi

remaining=0
for shard in "${shards[@]}"; do
    tag="$(basename "$shard" .txt)"
    [ -f "$run_dir/state/$tag.done" ] || remaining=$((remaining + 1))
done
if [ "$remaining" -eq 0 ]; then
    log "ALL ${#shards[@]} SHARDS COMPLETE -- results in $outdir"
else
    log "$remaining shards still outstanding -- re-run the same command to continue"
fi
