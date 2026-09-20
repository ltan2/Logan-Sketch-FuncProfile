# Logan Sketch / FuncProfile Pipeline

This repository selects WGS metagenomic accessions from SRA metadata, retrieves the corresponding Logan unitig and contig assemblies, and analyzes them with sourmash and FuncProfiler.

See the [design document](design/latest_design-doc.md) for the architecture, file-by-file Nextflow structure, failure behavior, output contract, and future extensions.

## Workflow

```text
SRA metadata -> DuckDB -> accession manifest -> Logan unitigs/contigs
                                                  `-> decompress once
                                                       |-> sourmash sketch
                                                       `-> FuncProfiler
```

Metadata preparation and accession selection are orchestrated by Make. Sequence retrieval and analysis are orchestrated by Nextflow.

## Requirements

- GNU Make
- Miniconda, Miniforge, or Mamba
- access to the public SRA and Logan S3 buckets
- sufficient local space for SRA metadata, the DuckDB database, and Nextflow work files

## Setup

Create the `logan` environment and download the KO signature collection:

```bash
make install
conda activate logan
```

`make install` is equivalent to running the `env` and `data` targets.

## Prepare the accession manifest

Run the complete metadata-to-manifest workflow with one command:

```bash
make get-accessions
```

This command:

1. synchronizes the public SRA metadata into `aws_sra_metadata/`;
2. creates or replaces the `sra_metadata` table in `sra_metadata.db` using `database/build_metadata_db.sql`; and
3. runs `database/wgs_metagenome.sql` to write `wgs_metagenome_accessions.txt`.

To build the database without generating a new accession manifest, run:

```bash
make create-db
```

The default paths can be overridden with Make variables. For example:

```bash
make get-accessions \
  METADATA_DIR=/path/to/sra_metadata \
  METADATA_DB=/path/to/sra_metadata.db
```

## Run the pipeline

Run all accessions in the default manifest:

```bash
nextflow run nextflow -resume
```

Run a smaller manifest or production shard:

```bash
nextflow run nextflow \
  --accessions path/to/accession_shard.txt \
  -resume
```

Run the repository smoke-test subset:

```bash
make run-test
```

Nextflow uses `nextflow/nextflow.config` for input paths, analysis parameters, process resources, retry behavior, output paths, and execution profiles.

Run the entire manifest with the shard-by-shard driver, which gates on the unit test, can be
paused and resumed, and plots progress live:

```bash
nextflow/run_full.sh \
  --accessions wgs_metagenome_accessions.txt \
  --run-dir /scratch/$USER/logan_full_run
```

[how_to_run.md](how_to_run.md) is the full procedure: setup, launch, live monitoring, pausing,
resuming, and recovering failures.

## Test and benchmark

`nextflow/test/run_unit_test.sh` runs the full pipeline on one accession and compares every published output against the known-good results in `DRR001355_test_res/`. It exits nonzero on any mismatch, so it can gate a real run:

```bash
nextflow/test/run_unit_test.sh && nextflow run nextflow -resume
```

`benchmark/run_benchmark.sh` times the workflow at increasing accession-list sizes (the unit test above gates every sweep), and `benchmark/plot_benchmark.py` turns the resulting Nextflow traces into per-stage timing, CPU, memory, and throughput plots.

## Repository structure

| Path | Role |
|---|---|
| `Makefile` | User-facing commands for setup, metadata preparation, accession selection, and testing |
| `database/aws_sra_metadata.sh` | Manual helper for synchronizing the public SRA metadata bucket; Make performs the same sync directly |
| `database/build_metadata_db.sql` | Builds the `sra_metadata` table from synchronized SRA Parquet files |
| `database/wgs_metagenome.sql` | Selects WGS metagenomic records and writes the accession manifest |
| `nextflow/main.nf` | Top-level Nextflow workflow and channel routing |
| `nextflow/nextflow.config` | Parameters, resources, profiles, retries, trace, and report settings |
| `nextflow/modules/fetch.nf` | Retrieves Logan unitigs and contigs and rejects missing or empty inputs |
| `nextflow/modules/decompress.nf` | Decompresses each valid Logan assembly once for both analysis branches |
| `nextflow/modules/sketch.nf` | Produces sourmash FracMinHash signatures |
| `nextflow/modules/funprofiler.nf` | Produces KEGG Orthology profiles |
| `nextflow/modules/cleanup.nf` | Deletes the shared decompressed FASTA once both consumers are done |
| `nextflow/modules/ledger.nf` | Records skip outcomes that bypass an analysis module |
| `nextflow/subworkflows/analyze.nf` | Fans valid sequences out to sourmash and FuncProfiler |
| `nextflow/bin/write_ledger_row.sh` | Writes one CSV progress record per accession, sequence type, and stage |
| `nextflow/unify_ledger.sh` | Unions the per-task ledger CSVs into one queryable DuckDB table |
| `nextflow/run_full.sh` | Shard-by-shard driver for a full run: unit-test gate, pause/resume, live plots |
| `nextflow/filter_completed.py` | Removes already-published accessions from a manifest (accession-level resume) |
| `nextflow/shard_manifest.py` | Splits the manifest into size-balanced shards, heaviest accession first |
| `nextflow/test/accessions_smoke.txt` | Small manifest used by the Nextflow `test` profile |
| `nextflow/test/run_unit_test.sh` | End-to-end pipeline test against known-good results |
| `benchmark/` | Throughput benchmarking harness and plotting scripts |
| `benchmark/live_monitor.py` | Live progress/resource plots for a run in progress |
| `benchmark/log_system_memory.sh` | Server RAM/swap sampler shared by the full run and the benchmark |
| `how_to_run.md` | End-to-end procedure for the full run |
| `DRR001355_test_res/` | Known-good outputs the unit test compares against |

## Outputs

Results are published under `results/` by default:

```text
results/
  sketches/{unitigs,contigs}/        # <acc>.<seq>.k31.sig.zip (DNA) and
                                     # <acc>.<seq>.protein.k11.sig.zip (translated protein)
  ko_profiles/{unitigs,contigs}/
  ledger/
  pipeline_trace.tsv
  pipeline_report.html
```

Results can be queried by running `nextflow/unify_ledger.sh [RESULTS_DIR] [OUTPUT_DB]`.

The ledger distinguishes completed work from valid skip conditions such as missing Logan objects. Use the same launch directory and `-resume` to reuse completed Nextflow tasks.
