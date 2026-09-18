#!/usr/bin/env bash
# Pipeline unit test: runs the full Nextflow pipeline on one accession with known-good results DRR001355_test_res and
# compares every published output against them (nextflow/test/compare_to_expected.py). Exits 0 only
# if everything matches, so it can gate a benchmark or production run:
#
#   nextflow/test/run_unit_test.sh && <real run>
#
# Environment (all optional):
#   UNIT_TEST_ACCESSION  accession to run                  (default DRR001355)
#   UNIT_TEST_EXPECTED   folder of known-good results       (default <repo>/DRR001355_test_res)
#   UNIT_TEST_DIR        scratch folder for the test run    (default /scratch/$USER/logan_sketch_funcprofile_unit_test)
#   NF_EXTRA_ARGS        extra `nextflow run` arguments, same as the real run gets, so the test
#                        exercises the same pipeline settings (e.g. --funprofiler_seq_types unitigs)
#
# Runs without -profile benchmark on purpose: a failing task should stop the test loudly instead of
# being ignored.
set -euo pipefail

test_dir_script="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$test_dir_script/../.." && pwd)"
accession="${UNIT_TEST_ACCESSION:-DRR001355}"
expected="${UNIT_TEST_EXPECTED:-$repo_root/DRR001355_test_res}"
test_dir="${UNIT_TEST_DIR:-/scratch/$(whoami)/logan_sketch_funcprofile_unit_test}"

[ -d "$expected" ] || { echo "UNIT TEST ERROR: expected results folder not found: $expected" >&2; exit 1; }

rm -rf "$test_dir"
mkdir -p "$test_dir"
echo "$accession" > "$test_dir/accessions.txt"

# The sequence types FUNPROFILER runs on decide which expected files apply: the pipeline
# config's defaults, unless NF_EXTRA_ARGS overrides them for this run. `nextflow config` runs inside
# the test folder because it writes .nextflow.log into the current directory.
config_value() {
    (cd "$test_dir" && nextflow config "$repo_root/nextflow" -flat 2>/dev/null) \
        | awk -F" = " -v k="params.$1" '$1==k {gsub(/\x27/, "", $2); print $2}'
}
override() {
    printf '%s\n' "${NF_EXTRA_ARGS:-}" | { grep -oE -- "--$1[ =][^ ]+" || true; } | sed -E "s/--$1[ =]//" | tail -1
}
funprofiler_seq_types="$(override funprofiler_seq_types)"; funprofiler_seq_types="${funprofiler_seq_types:-$(config_value funprofiler_seq_types)}"

# The translated protein sketch is always produced, so it is always checked; only its
# parameters vary -- same NF_EXTRA_ARGS-beats-config precedence as the sequence types above.
protein_ksize="$(override sourmash_protein_ksize)"; protein_ksize="${protein_ksize:-$(config_value sourmash_protein_ksize)}"
protein_scale="$(override sourmash_protein_scale)"; protein_scale="${protein_scale:-$(config_value sourmash_protein_scale)}"

echo "unit test settings: funprofiler_seq_types=$funprofiler_seq_types protein sketch k=${protein_ksize:-11}, scaled=${protein_scale:-1000}"
echo "=== [$(date '+%Y-%m-%d %H:%M:%S')] unit test: running pipeline on $accession ==="

# Nextflow writes .nextflow/ and its log into the working directory; keep those out of the repo.
(
    cd "$test_dir"
    # shellcheck disable=SC2086  # NF_EXTRA_ARGS is intentionally word-split into separate arguments
    nextflow run "$repo_root/nextflow" \
        --accessions "$test_dir/accessions.txt" \
        --outdir "$test_dir/results" \
        -work-dir "$test_dir/work" \
        ${NF_EXTRA_ARGS:-}
) > "$test_dir/nextflow_run.log" 2>&1 || {
    echo "UNIT TEST FAILED: the pipeline run itself failed -- see $test_dir/nextflow_run.log" >&2
    tail -20 "$test_dir/nextflow_run.log" >&2
    exit 1
}

echo "=== [$(date '+%Y-%m-%d %H:%M:%S')] unit test: comparing outputs ==="
python3 "$test_dir_script/compare_to_expected.py" \
    --expected "$expected" --results "$test_dir/results" --accession "$accession" \
    --funprofiler-seq-types "$funprofiler_seq_types" \
    --protein-ksize "${protein_ksize:-11}" --protein-scale "${protein_scale:-1000}"

rm -rf "$test_dir/work"   # passed: the intermediate task folders aren't needed; results stay for inspection
