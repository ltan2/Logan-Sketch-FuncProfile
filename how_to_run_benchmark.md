cd ~/Logan-Sketch-FuncProfile
conda activate logan
                                                                                                                    
# 1. Build the accession manifest into this repo.
#    Writes aws_sra_metadata/ (~13 GB), sra_metadata.db (~36 GB) and
#    wgs_metagenome_accessions.txt here. All three are gitignored.
make get-accessions

# 2. Sanity-check what you got
wc -l wgs_metagenome_accessions.txt
head -3 wgs_metagenome_accessions.txt
    
# 3. Run the 1,000-accession benchmark (samples randomly from the file above,
#    runs the unit-test gate first, then times the whole pipeline)
BENCH_DIR=/scratch/$USER/logan_bench_1k ./benchmark/run_benchmark.sh 1000
    
# 4. Attach to watch it (the script prints the session name it created)
tmux attach -t $(tmux ls -F '#S' | grep run_benchmark | tail -1)     # Ctrl-b d to detach
                                                                                                                    
# 5. When it finishes: the timing, then the plots                                                                                                                              
cat /scratch/$USER/logan_bench_1k/wall_clock_summary.csv
python3 benchmark/plot_benchmark.py --bench-dir /scratch/$USER/logan_bench_1k \
--total-cpus 700 --memory-limit-gb 2700

Optional, between 3 and 5 — live plots while it runs:

B=/scratch/$USER/logan_bench_1k
tmux new -d -s bench_monitor "~/miniconda3/envs/logan/bin/python3 \
~/Logan-Sketch-FuncProfile/benchmark/live_monitor.py \
--trace $B/results_1000/pipeline_trace.tsv --total 1000 \
--system-memory $B/results_1000/system_memory.tsv \
--out-dir $B/plots_live --total-cpus 700 --memory-limit-gb 2700"
# then open $B/plots_live/index.html, or: watch -n 60 jq . $B/plots_live/progress.json

