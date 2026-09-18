#!/usr/bin/env bash
# Launches run_benchmark.sh inside a detached tmux session so it keeps running
# after you disconnect. Any arguments are forwarded to run_benchmark.sh as sizes.
#
# Usage:
#   ./benchmark/launch_tmux_benchmark.sh [SIZE...]
# Examples:
#   ./benchmark/launch_tmux_benchmark.sh                # default sizes
#   ./benchmark/launch_tmux_benchmark.sh 1000 10000
#
# Reattach later with:  tmux attach -t benchmark
# Watch the log without attaching:  tail -f benchmark/run_benchmark.log

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(dirname "$script_dir")"

session="benchmark"
sizes="$*"

if tmux has-session -t "$session" 2>/dev/null; then
    echo "tmux session '$session' already exists -- attach with: tmux attach -t $session" >&2
    exit 1
fi

tmux new -s "$session" -d

# No ACCESSIONS_SOURCE here on purpose: run_benchmark.sh samples from the repo's own
# manifest, <repo>/wgs_metagenome_accessions.txt, the one `make get-accessions` writes.
# Set ACCESSIONS_SOURCE in the environment only to sample from somewhere else deliberately.
tmux send-keys -t "$session" "
cd '$repo_root'
conda activate logan
./benchmark/run_benchmark.sh $sizes 2>&1 | tee benchmark/run_benchmark.log
" C-m

echo "Started in tmux session '$session'."
echo "Reattach:      tmux attach -t $session"
echo "Watch the log: tail -f $script_dir/run_benchmark.log"
