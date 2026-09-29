"""Count, for every k-mer hash, how many cohort samples contain it (DNA k31 and protein k11).

Counting all scaled=1000 hashes for ~53k samples is ~20 billion values, so each sketch is
first downsampled to scaled=DOWNSAMPLE (keep hash < 2**64 / DOWNSAMPLE). That is itself a
FracMinHash sketch: a uniform random 1/(DOWNSAMPLE/1000) subset of the hash space, so the
per-hash prevalence distribution is unbiased and the top hashes are real, just fewer.

Writes:
  /scratch/lbt5343/early_eda_cache/{dna,prot}/<acc>.npy   downsampled sorted hashes per sample
  tables/hash_prevalence_hist_{dna,prot}.csv   n_samples -> number of hashes seen in that many
  tables/top_hashes_{dna,prot}.csv             the TOP_N most prevalent hashes
"""
import gzip
import json
import sys
import time
import zipfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from common import COHORT, HASH_CACHE, con, write_csv

DOWNSAMPLE = 20000
MAX_HASH = np.uint64(2**64 // DOWNSAMPLE)
TOP_N = 5000
WORKERS = 32
VIEWS = {"dna": "unitigs_dna_path", "prot": "unitigs_prot_path"}


def load_one(args):
    src, dst = args
    dst = Path(dst)
    if dst.exists():
        return str(dst), None
    try:
        with zipfile.ZipFile(src) as z:
            name = next(n for n in z.namelist() if n.startswith("signatures/"))
            data = z.read(name)
        if name.endswith(".gz"):
            data = gzip.decompress(data)
        sig = json.loads(data)
        sigs = sig if isinstance(sig, list) else [sig]
        mins = sigs[0]["signatures"][0]["mins"]
        h = np.asarray(mins, dtype=np.uint64)
        h = np.sort(h[h < MAX_HASH])
        tmp = dst.with_name(dst.stem + ".tmp.npy")
        np.save(tmp, h)
        tmp.rename(dst)
        return str(dst), None
    except Exception as exc:
        return str(dst), f"{type(exc).__name__}: {exc}"


def main():
    views = sys.argv[1:] or list(VIEWS)
    cohort = con().execute(f"SELECT accession, unitigs_dna_path, unitigs_prot_path FROM '{COHORT}'").df()
    for view in views:
        t = time.time()
        out = HASH_CACHE / view
        out.mkdir(parents=True, exist_ok=True)
        jobs = [(p, str(out / f"{a}.npy")) for a, p in zip(cohort.accession, cohort[VIEWS[view]])]
        errors = []
        with ProcessPoolExecutor(WORKERS) as ex:
            for i, (dst, err) in enumerate(ex.map(load_one, jobs, chunksize=16), 1):
                if err:
                    errors.append((dst, err))
                if i % 5000 == 0:
                    print(f"[{view}] cached {i:,}/{len(jobs):,}  {time.time() - t:.0f}s", flush=True)
        print(f"[{view}] cache done, {len(errors)} errors", flush=True)
        for e in errors[:10]:
            print("   ", e)

        arrays = [np.load(d) for _, d in jobs if Path(d).exists()]
        n_samples = len(arrays)
        allh = np.concatenate(arrays)
        del arrays
        print(f"[{view}] {n_samples:,} samples, {allh.size:,} downsampled hashes; counting", flush=True)
        uniq, counts = np.unique(allh, return_counts=True)
        del allh
        hist = pd.Series(counts).value_counts().sort_index()
        write_csv(pd.DataFrame({"n_samples": hist.index, "n_hashes": hist.values}),
                  f"hash_prevalence_hist_{view}")
        order = np.argsort(-counts, kind="stable")[:TOP_N]
        top = pd.DataFrame({"hash": uniq[order].astype(str), "n_samples": counts[order],
                            "frac_samples": counts[order] / n_samples})
        write_csv(top, f"top_hashes_{view}")
        print(f"[{view}] {uniq.size:,} distinct hashes; {time.time() - t:.0f}s total", flush=True)


if __name__ == "__main__":
    main()
