"""Annotate and plot the most shared k-mer hashes found by 02_hash_counts.py.

- protein k11 hashes -> KOs whose FuncProfiler KO sketch (k11, scaled 1000) contains the hash
- DNA k31 hashes     -> GTDB r232 genomes (YACHT reference sketches, k31, scaled 1000)
  that contain the hash. Unmatched hashes are not from any GTDB representative genome
  (candidates: host DNA, phiX/spike-ins, vectors, viruses, eukaryotes).
- prevalence of the top hashes per environment, from the per-sample downsampled cache.
"""
import csv
import gzip
import json
import zipfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import BLUE, COHORT, EDA, GTDB_SIGS, HASH_CACHE, KO_SIG, ORANGE, TAB, barh, con, save, write_csv

TOP_ANNOT = 5000
N_PLOT = 30
N_ENV = 12
WORKERS = 32
GTDB_TAX = "/scratch/akn5655/gtdb232/gtdb_genomes_reps_r232"


def ko_hits(top):
    with zipfile.ZipFile(KO_SIG) as z:
        man = [r for r in csv.DictReader(line for line in z.read("SOURMASH-MANIFEST.csv").decode().splitlines()
                                         if not line.startswith("#")) if r["ksize"] == "11"]
        hits = defaultdict(set)
        for r in man:
            data = z.read(r["internal_location"])
            sig = json.loads(gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data)
            for s in (sig if isinstance(sig, list) else [sig]):
                for mh in s["signatures"]:
                    for h in top.intersection(mh["mins"]):
                        hits[h].add(r["name"].replace("ko:", ""))
    return hits


def _scan_gtdb(args):
    paths, top = args
    out = []
    for p in paths:
        with open(p, "rb") as f:
            data = f.read()
        if data[:2] == b"\x1f\x8b":
            data = gzip.decompress(data)
        for s in json.loads(data):
            for mh in s["signatures"]:
                for h in top.intersection(mh["mins"]):
                    out.append((h, s.get("name", "")))
    return out


def gtdb_hits(top):
    paths = sorted(str(p) for p in (GTDB_SIGS / "signatures").iterdir())
    chunks = [(paths[i:i + 500], top) for i in range(0, len(paths), 500)]
    hits = defaultdict(set)
    with ProcessPoolExecutor(WORKERS) as ex:
        for i, res in enumerate(ex.map(_scan_gtdb, chunks), 1):
            for h, name in res:
                hits[h].add(name.split()[0].replace("_genomic", ""))
            if i % 50 == 0:
                print(f"  gtdb scan {i}/{len(chunks)} chunks", flush=True)
    return hits


def gtdb_genus():
    rows = []
    for f in ("bac120_taxonomy_r232.tsv", "ar53_taxonomy_r232.tsv"):
        with open(f"{GTDB_TAX}/{f}") as fh:
            for line in fh:
                acc, lin = line.rstrip("\n").split("\t")
                r = dict(x.split("__", 1) for x in lin.split(";"))
                rows.append((acc[7:].split(".")[0], r.get("p", ""), r.get("g", ""), r.get("s", "")))
    return {a: (p, g, s) for a, p, g, s in rows}


def summarize_gtdb(names, tax):
    if not names:
        return "", "", 0
    lin = [tax.get(n[4:].split(".")[0]) for n in names]
    lin = [x for x in lin if x]
    genera = pd.Series([x[1] for x in lin]).value_counts()
    phyla = pd.Series([x[0] for x in lin]).value_counts()
    return "; ".join(f"{g} ({n})" for g, n in genera.head(3).items()), \
        "; ".join(f"{p} ({n})" for p, n in phyla.head(3).items()), len(names)


def env_prevalence(view, hashes):
    coh = con().execute(f"SELECT accession, environment FROM '{COHORT}'").df()
    top_env = coh.environment.value_counts().head(N_ENV).index
    coh = coh[coh.environment.isin(top_env)]
    want = np.array(sorted(hashes), dtype=np.uint64)
    hits = np.zeros((len(coh), len(want)), dtype=bool)
    for i, a in enumerate(coh.accession):
        h = np.load(HASH_CACHE / view / f"{a}.npy")
        hits[i] = np.isin(want, h, assume_unique=True)
    df = pd.DataFrame(hits, columns=want.astype(str))
    df["environment"] = coh.environment.values
    return df.groupby("environment").mean().loc[top_env].T


def main():
    cohort_n = con().execute(f"SELECT count(*) FROM '{COHORT}'").fetchone()[0]
    names = pd.read_csv(EDA / "reference/kegg_ko_list.tsv", sep="\t", header=None, names=["ko", "definition"])
    sym = names.set_index("ko").definition.str.split(";").str[0].str.split(",").str[0]

    # --- prevalence distributions
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.3))
    for ax, view, label in zip(axes, ["dna", "prot"], ["DNA k=31", "protein k=11"]):
        h = pd.read_csv(TAB / f"hash_prevalence_hist_{view}.csv")
        total = h.n_hashes.sum()
        ax.loglog(h.n_samples, h.n_hashes, ".", color=BLUE, ms=3)
        single = h.n_hashes[h.n_samples == 1].sum() / total
        ge1 = h.n_hashes[h.n_samples >= 0.01 * cohort_n].sum()
        ax.set_xlabel("number of samples containing the hash")
        ax.set_ylabel("number of hashes")
        ax.set_title(f"{label}: {total:,} distinct hashes (scaled=20,000 subset)\n"
                     f"{single:.0%} seen in only 1 sample; {ge1:,} in ≥1% of samples", loc="left")
    fig.suptitle("How widely shared are k-mers? (hash prevalence across the cohort)", x=0.01, ha="left")
    fig.tight_layout()
    save(fig, "hash_01_prevalence_distribution")

    # --- annotate top hashes
    top_p = pd.read_csv(TAB / "top_hashes_prot.csv", dtype={"hash": str}).head(TOP_ANNOT)
    hits = ko_hits({int(h) for h in top_p.hash})
    top_p["kos"] = [";".join(sorted(hits.get(int(h), []))) for h in top_p.hash]
    top_p["ko_symbols"] = [";".join(str(sym.get(k, k)) for k in sorted(hits.get(int(h), []))) for h in top_p.hash]
    write_csv(top_p, "top_hashes_prot_annotated")
    print(f"protein: {(top_p.kos != '').mean():.0%} of top {TOP_ANNOT} hashes are in a KO sketch")

    top_d = pd.read_csv(TAB / "top_hashes_dna.csv", dtype={"hash": str}).head(TOP_ANNOT)
    ghits = gtdb_hits({int(h) for h in top_d.hash})
    tax = gtdb_genus()
    ann = [summarize_gtdb(ghits.get(int(h), set()), tax) for h in top_d.hash]
    top_d["top_genera"], top_d["top_phyla"], top_d["n_gtdb_genomes"] = zip(*ann)
    write_csv(top_d, "top_hashes_dna_annotated")
    print(f"dna: {(top_d.n_gtdb_genomes > 0).mean():.0%} of top {TOP_ANNOT} hashes are in a GTDB genome")

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    t = top_d.head(N_PLOT)
    lab = [f"{h[-6:]}  {g.split(' (')[0] if g else 'not in GTDB reps'}"
           + (f" +{n - 1}" if n > 1 else "") for h, g, n in zip(t.hash, t.top_genera, t.n_gtdb_genomes)]
    barh(axes[0], lab, t.frac_samples * 100, fmt="{:.0f}%")
    axes[0].set_title(f"Top {N_PLOT} DNA k31 hashes (label: hash tail, top GTDB genus, +other genomes)",
                      loc="left")
    t = top_p.head(N_PLOT)
    lab = [f"{h[-6:]}  {s.split(';')[0][:30] if s else 'no KO'}" for h, s in zip(t.hash, t.ko_symbols)]
    barh(axes[1], lab, t.frac_samples * 100, color=ORANGE, fmt="{:.0f}%")
    axes[1].set_title(f"Top {N_PLOT} protein k11 hashes (label: hash tail, KO containing it)", loc="left")
    for ax in axes:
        ax.set_xlabel("% of cohort samples")
        ax.tick_params(axis="y", labelsize=7)
    fig.tight_layout()
    save(fig, "hash_02_top_hashes")

    # --- per-environment prevalence of the top hashes
    for view, t, color in [("dna", top_d, BLUE), ("prot", top_p, ORANGE)]:
        m = env_prevalence(view, set(int(h) for h in t.hash.head(N_PLOT)))
        m = m.loc[t.hash.head(N_PLOT)]
        write_csv(m.reset_index(names="hash"), f"top_hashes_{view}_by_environment")
    fig, axes = plt.subplots(1, 2, figsize=(16, 9))
    for ax, view, t, lab_col in [(axes[0], "dna", top_d, "top_genera"), (axes[1], "prot", top_p, "ko_symbols")]:
        m = pd.read_csv(TAB / f"top_hashes_{view}_by_environment.csv", dtype={"hash": str}).set_index("hash")
        im = ax.imshow(m.values * 100, cmap="viridis", vmin=0, vmax=100, aspect="auto")
        labs = t.set_index("hash")[lab_col].fillna("").reindex(m.index)
        ax.set_yticks(range(len(m)), [f"{h[-6:]} {str(l).split(';')[0].split(' (')[0][:22]}" for h, l in labs.items()],
                      fontsize=7)
        ax.set_xticks(range(m.shape[1]), m.columns, rotation=40, ha="right", fontsize=7)
        ax.set_title(f"Top {N_PLOT} {'DNA' if view == 'dna' else 'protein'} hashes: % of samples per environment",
                     loc="left")
    fig.colorbar(im, ax=axes, shrink=0.5, label="% of samples in environment")
    save(fig, "hash_03_top_hashes_by_environment")


if __name__ == "__main__":
    main()
