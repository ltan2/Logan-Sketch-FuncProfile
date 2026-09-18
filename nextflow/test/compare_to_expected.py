#!/usr/bin/env python3
"""Compare one accession's pipeline outputs against a folder of known-good results.

Checks, per sequence type (unitigs, contigs):
  sketch        results/sketches/<seq>/<acc>.<seq>.k31.sig.zip
                -> same ksize, scaled and hashes (abundances too, if both files track them)
  KO profile    results/ko_profiles/<seq>/<acc>.<seq>_ko_profiles.csv
                -> same KOs with the same relative abundances
  prefetch      results/ko_profiles/<seq>/<acc>.<seq>_prefetch_out.csv
                -> same matches with the same values (file-path columns ignored, so expected
                   results made on another machine still compare)
FUNPROFILER outputs are only checked for the sequence types the run was configured to produce
(--funprofiler-seq-types); expected files for other types are skipped.

Expected-results folder layout (as produced by the original logan-subworkflow.sh run):
  <acc>.<seq>.k31.sig.zip, <acc>_ko_profiles_<unitig|contig>, <acc>_<unitig|contig>_prefetch_out

Exit status 0 = everything matches, 1 = at least one mismatch or missing output.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import sourmash

SINGULAR = {"unitigs": "unitig", "contigs": "contig"}
TOL = 1e-9


class Report:
    def __init__(self):
        self.failures = 0

    def ok(self, what, detail=""):
        print(f"  PASS  {what}{': ' + detail if detail else ''}")

    def fail(self, what, detail):
        self.failures += 1
        print(f"  FAIL  {what}: {detail}")

    def note(self, what, detail):
        print(f"  NOTE  {what}: {detail}")


def frames_equal(a, b, key, ignore=()):
    """Row-order-independent comparison on shared columns; numbers within TOL, text exact.
    Returns a description of the first difference, or None."""
    if len(a) != len(b):
        return f"{len(a)} rows expected, {len(b)} produced"
    if len(a) == 0:
        return None
    if set(a[key]) != set(b[key]):
        missing, extra = set(a[key]) - set(b[key]), set(b[key]) - set(a[key])
        return f"{key} differs (missing {sorted(missing)[:5]}, unexpected {sorted(extra)[:5]})"
    cols = [c for c in a.columns if c in b.columns and not c.endswith("_filename") and c not in ignore]
    a = a.sort_values(key).reset_index(drop=True)[cols]
    b = b.sort_values(key).reset_index(drop=True)[cols]
    for c in cols:
        if pd.api.types.is_numeric_dtype(a[c]) and pd.api.types.is_numeric_dtype(b[c]):
            if not np.allclose(a[c].astype(float), b[c].astype(float), rtol=0, atol=TOL, equal_nan=True):
                return f"column {c} values differ"
        elif not a[c].astype(str).equals(b[c].astype(str)):
            return f"column {c} values differ"
    return None


def read_csv_or_empty(path):
    return pd.read_csv(path) if os.path.getsize(path) else pd.DataFrame()


def check_sketch(rep, expected, produced, label):
    """Compare two sourmash sketches: same k, scaled and hashes (abundances too, when both
    files track them)."""
    if not os.path.exists(produced):
        rep.fail(label, f"missing {produced}")
        return
    e = list(sourmash.load_file_as_signatures(expected))[0].minhash
    p = list(sourmash.load_file_as_signatures(produced))[0].minhash
    if (e.ksize, e.scaled) != (p.ksize, p.scaled):
        rep.fail(label, f"k/scaled {e.ksize}/{e.scaled} expected, {p.ksize}/{p.scaled} produced")
        return
    if set(e.hashes) != set(p.hashes):
        rep.fail(label, f"hashes differ ({len(e.hashes)} expected, {len(p.hashes)} produced)")
        return
    if e.track_abundance and p.track_abundance and dict(e.hashes) != dict(p.hashes):
        rep.fail(label, "same hashes but abundances differ")
        return
    rep.ok(label, f"{len(p.hashes)} hashes identical")
    if e.track_abundance != p.track_abundance:
        rep.note(label, f"abundance tracking differs (expected {e.track_abundance}, produced {p.track_abundance}); "
                        "hashes compared only")


def check_table(rep, expected, produced, key, label):
    if not os.path.exists(produced):
        return rep.fail(label, f"missing {produced}")
    diff = frames_equal(read_csv_or_empty(expected), read_csv_or_empty(produced), key)
    rep.fail(label, diff) if diff else rep.ok(label, f"{len(read_csv_or_empty(produced))} rows identical")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--expected", required=True, help="Folder of known-good results")
    parser.add_argument("--results", required=True, help="Pipeline --outdir for the test run")
    parser.add_argument("--accession", required=True)
    parser.add_argument("--funprofiler-seq-types", default="unitigs,contigs")
    args = parser.parse_args()

    acc, E, R = args.accession, args.expected, args.results
    funprofiler_types = {s.strip() for s in args.funprofiler_seq_types.split(",")}
    rep = Report()
    print(f"Comparing {acc}: produced {R} vs expected {E}")

    for seq, one in SINGULAR.items():
        expected_sketch = f"{E}/{acc}.{seq}.k31.sig.zip"
        if os.path.exists(expected_sketch):
            check_sketch(rep, expected_sketch, f"{R}/sketches/{seq}/{acc}.{seq}.k31.sig.zip", f"{seq} sketch")

        for expected, produced, key, label in (
                (f"{E}/{acc}_ko_profiles_{one}", f"{R}/ko_profiles/{seq}/{acc}.{seq}_ko_profiles.csv", "ko_id", f"{seq} KO profile"),
                (f"{E}/{acc}_{one}_prefetch_out", f"{R}/ko_profiles/{seq}/{acc}.{seq}_prefetch_out.csv", "match_name", f"{seq} prefetch")):
            if not os.path.exists(expected):
                continue
            if seq not in funprofiler_types:
                rep.note(label, "skipped (FUNPROFILER not run on this sequence type)")
                continue
            check_table(rep, expected, produced, key, label)

    if rep.failures:
        print(f"UNIT TEST FAILED: {rep.failures} check(s) did not match")
        sys.exit(1)
    print("UNIT TEST PASSED")


if __name__ == "__main__":
    main()
