# How to run the full Logan sketch / funcprofile pipeline

Everything needed to take the whole WGS-metagenome accession manifest through the pipeline:
the gate that must pass first, how to launch, how to watch it live, how to pause and resume
it, and what to do when something fails.

Run every command from the repository root with the `logan` environment active:

```bash
cd ~/Logan-Sketch-FuncProfile
conda activate logan
```

Contents:

1. [One-time setup](#1-one-time-setup)
2. [Build the accession manifest](#2-build-the-accession-manifest)
3. [Gate: the unit test must pass](#3-gate-the-unit-test-must-pass)
4. [Plan the run](#4-plan-the-run)
5. [Launch the full run](#5-launch-the-full-run)
6. [Watch it live](#6-watch-it-live)
7. [Pause](#7-pause)
8. [Resume](#8-resume)
9. [When something fails](#9-when-something-fails)
10. [Finishing up](#10-finishing-up)
11. [Reference](#11-reference)

---

## 1. One-time setup

```bash
make install          # creates/updates the `logan` conda env and downloads the KO signatures
conda activate logan
```

`make install` writes `KOs_sketched_scaled_1000.sig.zip` (98 MB) to the repository root, which
is where `params.ko_sig` looks for it. If you already have that file elsewhere, symlink it
instead of downloading it again:

```bash
ln -sfn /path/to/KOs_sketched_scaled_1000.sig.zip KOs_sketched_scaled_1000.sig.zip
```

Check the tools are the pinned ones:

```bash
nextflow -v          # nextflow version 26.04.6
sourmash --version   # sourmash 4.9.4
```

## 2. Build the accession manifest

This synchronizes the public SRA metadata, builds the DuckDB database, and writes
`wgs_metagenome_accessions.txt` — one accession per line, the input to the whole run:

```bash
make get-accessions
```

The metadata snapshot and database are large; put them somewhere with room if the repo disk
is small:

```bash
make get-accessions \
  METADATA_DIR=/scratch/$USER/aws_sra_metadata \
  METADATA_DB=/scratch/$USER/sra_metadata.db
```

Then confirm what you are about to run:

```bash
wc -l wgs_metagenome_accessions.txt
head -3 wgs_metagenome_accessions.txt
```

## 3. Gate: the unit test must pass

`nextflow/test/run_unit_test.sh` runs the complete pipeline on one accession (`DRR001355`) and
compares every published output — sketches, KO profiles, prefetch tables — against the
known-good results in `DRR001355_test_res/`. It exits nonzero on any mismatch.

```bash
nextflow/test/run_unit_test.sh
```

Expected (last verified 2026-09-18, passing):

```text
unit test settings: funprofiler_seq_types=unitigs,contigs
=== [...] unit test: running pipeline on DRR001355 ===
=== [...] unit test: comparing outputs ===
Comparing DRR001355: produced /scratch/$USER/logan_sketch_funcprofile_unit_test/results vs expected .../DRR001355_test_res
  PASS  unitigs sketch: 1600 hashes identical
  NOTE  unitigs sketch: abundance tracking differs (expected False, produced True); hashes compared only
  PASS  unitigs KO profile: 98 rows identical
  PASS  unitigs prefetch: 98 rows identical
  PASS  contigs sketch: 17 hashes identical
  NOTE  contigs sketch: abundance tracking differs (expected False, produced True); hashes compared only
  PASS  contigs KO profile: 0 rows identical
  PASS  contigs prefetch: 0 rows identical
UNIT TEST PASSED
```

The two `NOTE` lines are expected: the reference sketches were made without abundance
tracking, so only their hashes can be compared. A `FAIL` line, or a nonzero exit, means do not
start the full run.

Run it by hand with the same options the real run will use, if you plan to pass any:

```bash
NF_EXTRA_ARGS='--funprofiler_seq_types unitigs' nextflow/test/run_unit_test.sh
```

You do not have to remember to run it before launching: **`nextflow/run_full.sh` runs this same
test first and refuses to start if it fails** (override only deliberately, with
`SKIP_UNIT_TEST=1`).

## 4. Plan the run

Print what the driver would do. It creates the run directory, but runs, shards and freezes nothing:

```bash
nextflow/run_full.sh \
  --accessions wgs_metagenome_accessions.txt \
  --run-dir /scratch/$USER/logan_full_run \
  --shard-size 5000 \
  --dry-run
```

What to decide first:

| Decision | Guidance |
|---|---|
| **Where** | Put `--run-dir` on `/scratch`, not the repo disk. Results, work dirs and state all live under it. |
| **Shard size** | `--shard-size 5000` is the default. Smaller shards = finer pause granularity and less work lost to an interrupt, but more Nextflow startups and more end-of-shard drain time where the machine isn't full. |
| **Concurrency** | Set in `nextflow/nextflow.config`: `executor.cpus` (700), `executor.memory` (2700 GB), and `FETCH_LOGAN.maxForks` (50 concurrent S3 downloads). Agree these with the other users of this shared server before starting. |

One more caveat worth knowing before you start: a full run publishes millions of small files
into `results/sketches/*/` and `results/ko_profiles/*/`. That is fine to write and to open by
name, but don't expect `ls` in those directories to be quick. Everything below looks files up
by exact path for that reason.

## 5. Launch the full run

```bash
nextflow/run_full.sh \
  --accessions wgs_metagenome_accessions.txt \
  --run-dir /scratch/$USER/logan_full_run \
  --shard-size 5000
```

That single command:

1. re-launches itself in a **detached tmux session** (`logan_full_run`) so the run survives a
   dropped SSH connection — pass `NO_TMUX=1` to stay in the foreground;
2. freezes the manifest into `<run-dir>/manifest.txt` and splits it into shards of 5,000
   accessions, recording a checksum so a resumed run can't silently re-shard a changed
   manifest;
3. runs the **unit test gate** from section 3 and stops if it fails;
4. starts the **server memory logger** and the **live plots** (section 6);
5. runs the shards one at a time with `-profile full`, all publishing into one
   `<run-dir>/results` tree, marking each shard complete in `<run-dir>/state/` and deleting
   that shard's work directory as it goes.

Watch it, or detach again with `Ctrl-b d`:

```bash
tmux attach -t logan_full_run
```

Useful variants:

```bash
# Run just two shards and stop -- a good way to see real throughput before committing.
nextflow/run_full.sh --run-dir /scratch/$USER/logan_full_run --shards 2

# Keep every shard's work directory (for debugging; costs a lot of disk).
nextflow/run_full.sh --run-dir /scratch/$USER/logan_full_run --keep-work

# Extra Nextflow arguments, applied to the gate and to every shard.
NF_EXTRA_ARGS='--funprofiler_seq_types unitigs' \
  nextflow/run_full.sh --run-dir /scratch/$USER/logan_full_run

# Without the live monitor.
nextflow/run_full.sh --run-dir /scratch/$USER/logan_full_run --no-monitor

nextflow/run_full.sh --help
```

## 6. Watch it live

`benchmark/live_monitor.py` runs alongside the pipeline (started automatically by
`run_full.sh`) and redraws, every 60 s from the traces Nextflow writes as tasks complete:

| File in `<run-dir>/plots_live/` | What it shows |
|---|---|
| `index.html` | **Open this.** Auto-refreshing dashboard: accessions done, % of manifest, accessions/hour, days remaining, shards complete, failures — with every plot below |
| `00_progress.png` | Cumulative accessions complete against the whole manifest, completions per time bin, and per-stage task counts by status |
| `04_cpu_utilization.png` | CPUs allocated vs. actually used, active task count, RAM reserved vs. used, swap |
| `05_task_timeline.png` | Gantt timeline of accessions currently moving through the stages |
| `06_runtime_distribution.png` | Per-stage task runtime distribution (median, IQR, p90, p95) |
| `09_task_resources.png/.csv` | Per-stage table: task time, cores reserved vs. used, RAM reserved vs. peak |
| `progress.json` | The same numbers as data |

The resource plots cover the last 24 h (`--window-hours`); progress, rate and ETA always cover
the whole run.

One thing to know when reading them: Nextflow writes a trace row when a task *finishes*, so
the plots are complete up to the last finished task and a task in flight is invisible until it
ends. With a median task around a minute and a p95 around ten, the dashboard trails the real
state by minutes, not hours — but a sudden flat spot at the right edge usually means "tasks
still running", not "nothing happening".

**In the terminal:**

```bash
# the current numbers, refreshed every minute
watch -n 60 jq . /scratch/$USER/logan_full_run/plots_live/progress.json

# the monitor's own log, one line per refresh
tail -f /scratch/$USER/logan_full_run/logs/live_monitor.log

# shard-by-shard timing
column -t /scratch/$USER/logan_full_run/state/progress.tsv
```

**Watching a plain `nextflow run`** (no driver) works too — point the monitor at its trace:

```bash
python3 benchmark/live_monitor.py \
  --trace results/pipeline_trace.tsv \
  --total "$(grep -vc '^#' wgs_metagenome_accessions.txt)" \
  --out-dir plots_live
```

Add `--once` for a single snapshot instead of a loop (this is also how you regenerate the
plots after the run is over), `--total-cpus 700 --memory-limit-gb 2700` to draw the capacity
reference lines, and `--interval 30` to refresh faster.

## 7. Pause

Three ways to stop, from gentlest to most immediate. All of them are safe: every task that has
finished has already published its results, and section 8 picks up from there.

**a. Pause at the next shard boundary — preferred**

```bash
touch /scratch/$USER/logan_full_run/PAUSE
```

The shard that is running finishes normally, then the driver stops and prints how to resume.
Nothing is repeated, nothing is thrown away. This is the pause to use for anything planned
(freeing the machine overnight, a maintenance window). With the default shard size it takes
effect within a few hours; check with:

```bash
tail -5 /scratch/$USER/logan_full_run/state/progress.tsv
```

**b. Stop now, mid-shard**

```bash
tmux attach -t logan_full_run     # then press Ctrl-C
```

Or without attaching, by sending the driver the same interrupt:

```bash
pgrep -af 'run_full\.sh'          # the driver is the `bash .../run_full.sh --accessions ...` line
kill -INT <pid>
```

(Match the pid carefully: `pgrep -f` also matches any shell whose command line merely mentions
the script.)

The driver forwards the interrupt to Nextflow, which stops and discards its in-flight tasks. It
then waits for Nextflow to actually exit — a few tens of seconds while it kills those tasks —
before printing a resume command and exiting, so let the driver exit before you resume. The
accessions that had already finished in that shard are published and will be skipped on resume;
the ones in flight are re-run from the start. With a 5,000-accession shard you lose at most the
partial work of that shard, not the run.

**c. Freeze for a few minutes**

```bash
driver=$(pgrep -f 'run_full\.sh' | head -1)
pgid=$(ps -o pgid= -p "$driver" | tr -d ' ')

kill -STOP -- -"$pgid"     # freeze the driver, Nextflow and every running task
kill -CONT -- -"$pgid"     # let them continue
```

Use this only to hand the machine's cores to someone else briefly. Frozen tasks keep their
scratch files and their open S3 connections, and connections can time out if you leave them
stopped for long. For anything longer than a coffee break, use (a) or (b).

After (a) or (b), the live dashboard stops refreshing along with the driver — the plots it
already wrote stay where they are, and the resumed run starts a new one. To redraw them while
the run is paused:

```bash
python3 benchmark/live_monitor.py --run-dir /scratch/$USER/logan_full_run --once
```

## 8. Resume

Re-run the same command. That is the whole procedure:

```bash
rm -f /scratch/$USER/logan_full_run/PAUSE          # only if you created it

nextflow/run_full.sh \
  --accessions wgs_metagenome_accessions.txt \
  --run-dir /scratch/$USER/logan_full_run \
  --shard-size 5000
```

What happens, in order:

1. shards already marked complete in `<run-dir>/state/` are skipped outright;
2. the interrupted shard is re-run, but first `nextflow/filter_completed.py` removes from it
   every accession whose results are already published, so only the unfinished remainder runs;
3. each re-attempt writes its own trace (`pipeline_trace.shard_00042.a2.tsv`), so the live
   plots keep the history of the whole run across pauses.

**Why this does not use `nextflow -resume`.** The pipeline deletes its own intermediates as it
goes to stay inside the disk budget: `DECOMPRESS_FASTA` deletes the `.fa.zst` it consumed, and
`CLEANUP_FASTA` deletes the decompressed FASTA once the sketch and KO profile are done with
it. Nextflow invalidates any cached task whose declared outputs are missing, so on `-resume`
`FETCH_LOGAN` is invalidated, re-runs, re-timestamps its outputs, and every task below it
re-runs too. Measured on the 5-accession test profile, a second run with `-resume` reported:

```text
[SUCCESS] completed=45 failed=0 cached=0
```

— i.e. it re-downloaded and re-analyzed all five accessions. Published results, not the work
cache, are what make this pipeline resumable, which is what `filter_completed.py` uses.

You can ask what is left at any time (it only stats the exact expected paths, so it is fine to
run against a results tree with millions of files, though a full 1.2M-accession manifest still
takes a few minutes):

```bash
python3 nextflow/filter_completed.py \
  --accessions /scratch/$USER/logan_full_run/manifest.txt \
  --outdir /scratch/$USER/logan_full_run/results \
  --count-only
# prints: <complete>	<remaining>
```

## 9. When something fails

**A shard failed and the driver stopped.** That is deliberate — a stage failing on real data is
a signal, not noise. Look at, in this order:

```bash
RUN=/scratch/$USER/logan_full_run
tail -40 $RUN/logs/shard_00042.log             # the shard's console output
grep -iE 'error|caused by' $RUN/logs/nextflow.shard_00042.log | head
ls $RUN/work/shard_00042                       # kept on failure: the failing task's directory
```

The failing task's directory holds `.command.sh`, `.command.err` and `.command.log` — run
`.command.sh` by hand there to reproduce it. Once the cause is fixed, re-run the same
`run_full.sh` command; the shard restarts with its already-published accessions filtered out.

**Individual accessions that could not be fetched.** Under `-profile full`, `FETCH_LOGAN` gives
up on an accession after its retries instead of aborting the run (an SRA record can be removed
or embargoed). Those are recorded per shard rather than dropped silently:

```bash
cat /scratch/$USER/logan_full_run/state/*.failed.txt | sort -u | wc -l
```

Re-run them later as their own manifest, publishing into the same results tree:

```bash
cat /scratch/$USER/logan_full_run/state/*.failed.txt | sort -u > /tmp/retry_accessions.txt

nextflow/run_full.sh \
  --accessions /tmp/retry_accessions.txt \
  --run-dir /scratch/$USER/logan_full_run_retry \
  --outdir /scratch/$USER/logan_full_run/results \
  --shard-size 500
```

**The driver refuses to start.**

| Message | Meaning |
|---|---|
| `another run_full.sh is already running against <run-dir>` | A driver already holds that run directory's lock. Don't start a second one; attach to the tmux session instead. |
| `tmux session 'logan_full_run' already exists` | The previous run's session is still there. `tmux attach -t logan_full_run`, and `tmux kill-session -t logan_full_run` once you're sure it has exited. |
| `manifest.txt changed since this run was sharded` | The manifest differs from the one this run directory was sharded from, so shard numbering no longer means the same accessions. Use a new `--run-dir` (results can still go to the same `--outdir`). |
| `pipeline unit test failed` | Fix the pipeline first. Section 3. |

## 10. Finishing up

When every shard is complete the driver prints `ALL <n> SHARDS COMPLETE`. Then:

```bash
RUN=/scratch/$USER/logan_full_run

# 1. Confirm nothing is outstanding (expect "<total>	0").
python3 nextflow/filter_completed.py --accessions $RUN/manifest.txt --outdir $RUN/results --count-only

# 2. Union the per-task ledger CSVs into one queryable table.
nextflow/unify_ledger.sh $RUN/results $RUN/results/ledger.duckdb
duckdb $RUN/results/ledger.duckdb -c "SELECT stage, status, count(*) FROM ledger GROUP BY 1,2 ORDER BY 1,2;"

# 3. Final plots over the whole run instead of the last 24 h.
python3 benchmark/live_monitor.py --run-dir $RUN --once --window-hours 100000 \
  --total-cpus 700 --memory-limit-gb 2700
```

Results are in `$RUN/results/{sketches,ko_profiles,ledger}/`. Each sequence type's sketch
directory holds two sketches per accession: `<acc>.<seq>.k31.sig.zip` (DNA) and
`<acc>.<seq>.protein.k11.sig.zip` (translated protein, built with the KO collection's k and
scaled so it can be searched against it directly):

```bash
sourmash prefetch $RUN/results/sketches/unitigs/DRR001355.unitigs.protein.k11.sig.zip \
  KOs_sketched_scaled_1000.sig.zip \
  --protein -k 11 --scaled 1000 --threshold-bp 1000 -o prefetch_out.csv
```

That specific query is already published per accession as
`results/ko_profiles/<seq>/<acc>.<seq>_prefetch_out.csv` by FUNPROFILER; the sketch is what
lets you run any *other* protein-space query later, once the FASTA is gone.

## 11. Reference

### Run-directory layout

```text
<run-dir>/
  manifest.txt                  frozen manifest this run is sharded from
  manifest.sha256               guards against it changing under a resume
  shards/shard_NNNNN.txt        the shards (stable; never re-sharded)
  shards/shard_NNNNN.todo.txt   what the latest attempt at that shard actually ran
  state/shard_NNNNN.done        shard complete
  state/shard_NNNNN.failed.txt  accessions with a FAILED task in that shard
  state/progress.tsv            one row per shard attempt: shard, accessions, start, end, status
  logs/shard_NNNNN.log          that shard's console output
  logs/nextflow.*.log           Nextflow's own debug logs
  logs/live_monitor.log         one line per dashboard refresh
  traces/pipeline_trace.*.tsv   per-shard-attempt trace -- what the live plots read
  work/shard_NNNNN/             Nextflow work dir; deleted when the shard completes
  results/                      published results for the whole run
  system_memory.tsv             server RAM/swap samples
  plots_live/                   live dashboard and plots
  PAUSE                         create to stop after the current shard
```

### Command summary

| Task | Command |
|---|---|
| Set up | `make install && conda activate logan` |
| Build the manifest | `make get-accessions` |
| Unit test | `nextflow/test/run_unit_test.sh` |
| Plan | `nextflow/run_full.sh --accessions ... --run-dir ... --dry-run` |
| Launch | `nextflow/run_full.sh --accessions ... --run-dir ... --shard-size 5000` |
| Attach | `tmux attach -t logan_full_run` |
| Live dashboard | open `<run-dir>/plots_live/index.html` (or serve it with `python3 -m http.server`) |
| Pause after this shard | `touch <run-dir>/PAUSE` |
| Pause now | attach, `Ctrl-C` |
| Resume | re-run the launch command (remove `PAUSE` first) |
| What's left | `python3 nextflow/filter_completed.py --accessions <run-dir>/manifest.txt --outdir <run-dir>/results --count-only` |
| Ledger table | `nextflow/unify_ledger.sh <run-dir>/results` |
| Benchmark sweep | `benchmark/run_benchmark.sh 1000 10000` then `python3 benchmark/plot_benchmark.py --bench-dir ...` |

### Files this procedure uses

| Path | Role |
|---|---|
| `nextflow/run_full.sh` | Shard-by-shard driver: gate, sharding, pause/resume, state, live monitor |
| `nextflow/filter_completed.py` | Removes already-published accessions from a manifest (accession-level resume) |
| `benchmark/live_monitor.py` | Live plots and dashboard from the traces of a run in progress |
| `benchmark/log_system_memory.sh` | Server RAM/swap sampler, shared by the full run and the benchmark |
| `benchmark/plot_benchmark.py` | After-the-fact plots for a benchmark sweep across several batch sizes |
| `nextflow/test/run_unit_test.sh` | The gate: full pipeline on one accession vs. known-good results |
| `nextflow/unify_ledger.sh` | Unions the per-task ledger CSVs into one DuckDB table |
| `nextflow/nextflow.config` | Parameters, resources, and the `standard` / `test` / `full` / `benchmark` profiles |
