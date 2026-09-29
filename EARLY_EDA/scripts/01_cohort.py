"""Define the EDA cohort: unitig samples whose DNA sketch, protein sketch and KO profile all
succeeded (`unitigs_usable_with_protein` and `unitigs_status_ko = 'ok'` in
analysis/environment_discordance/samples.parquet).

Writes cohort.parquet and tables/cohort_funnel.csv.
"""
import pandas as pd

from common import COHORT, SAMPLES, con, write_csv


def main():
    c = con()
    c.execute(f"CREATE TABLE s AS SELECT * FROM '{SAMPLES}'")
    steps = [
        ("processed (>= 1 published output)", "true"),
        ("unitig DNA sketch ok", "unitigs_status_dna = 'ok'"),
        ("+ unitig protein sketch ok", "unitigs_usable_with_protein"),
        ("+ KO profile with >= 1 KO", "unitigs_usable_with_protein AND unitigs_status_ko = 'ok'"),
    ]
    rows = []
    for label, where in steps:
        n, studies = c.execute(f"SELECT count(*), count(DISTINCT study) FROM s WHERE {where}").fetchone()
        rows.append({"step": label, "samples": n, "studies": studies})
    cohort_where = steps[-1][1]
    n_qc = c.execute(f"SELECT count(*) FROM s WHERE {cohort_where} AND qc_pass").fetchone()[0]
    rows.append({"step": "(of which qc_pass)", "samples": n_qc, "studies": None})
    write_csv(pd.DataFrame(rows), "cohort_funnel")
    c.execute(f"COPY (SELECT * FROM s WHERE {cohort_where}) TO '{COHORT}' (FORMAT parquet)")
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
