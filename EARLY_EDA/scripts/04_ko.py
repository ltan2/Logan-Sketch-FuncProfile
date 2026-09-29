"""Functional (KO) EDA on the cohort's unitig FuncProfiler profiles.

KO `abundance` is FuncProfiler's share of summed per-KO sketch containment, not gene
abundance (see analysis/environment_discordance/data_dictionary.md). Prevalence (fraction of
samples in which a KO is detected) is the more robust summary and is used for most plots.

KEGG names and the BRITE ko00001 hierarchy are read from reference/ (downloaded from
rest.kegg.jp; see README).
"""
import csv
import json
import zipfile
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import BLUE, COHORT, EDA, KO_PROFILES, KO_SIG, ORANGE, barh, con, save, write_csv

REF = EDA / "reference"
N_ENV = 12
# BRITE level-A groups that are pathways; 09180 (BRITE hierarchies) and 09190 (not in
# pathway/BRITE) are summarized as their own buckets.
PATHWAY_A = {"09100", "09120", "09130", "09140", "09150", "09160"}


def kegg_tables():
    names = pd.read_csv(REF / "kegg_ko_list.tsv", sep="\t", header=None, names=["ko", "definition"])
    names["symbol"] = names.definition.str.split(";").str[0].str.split(",").str[0]
    brite = json.load(open(REF / "kegg_brite_ko00001.json"))
    rows = []
    for a in brite["children"]:
        for b in a["children"]:
            for c in b.get("children", []):
                for k in c.get("children", []):
                    rows.append({"ko": k["name"].split()[0], "level_a": a["name"], "level_b": b["name"],
                                 "pathway": c["name"]})
    return names, pd.DataFrame(rows)


def ko_universe():
    with zipfile.ZipFile(KO_SIG) as z:
        man = z.read("SOURMASH-MANIFEST.csv").decode().splitlines()
    # One signature per KO per ksize (5/7/11/15); count distinct KO names.
    rows = csv.DictReader(line for line in man if not line.startswith("#"))
    return len({r["name"] for r in rows})


def main():
    names, brite = kegg_tables()
    n_db = ko_universe()
    c = con()
    c.execute(f"CREATE TABLE coh AS SELECT accession, environment, study FROM '{COHORT}'")
    c.execute(f"""CREATE TABLE ko AS SELECT p.accession, replace(p.ko_id, 'ko:', '') AS ko, p.abundance
                  FROM '{KO_PROFILES}' p JOIN coh USING (accession) WHERE p.seq_type = 'unitigs'""")
    n = c.execute("SELECT count(DISTINCT accession) FROM ko").fetchone()[0]
    print(f"{n:,} samples with KO rows; {n_db:,} KOs in the FuncProfiler KO database")

    # --- per-KO prevalence and abundance
    per_ko = c.execute(f"""
        SELECT ko, count(*) AS n_samples, count(*) / {n} AS prevalence,
               avg(abundance) AS mean_abund_when_present, sum(abundance) / {n} AS mean_abund_overall,
               count(DISTINCT study) AS n_studies
        FROM ko JOIN coh USING (accession) GROUP BY ko""").df()
    per_ko = per_ko.merge(names, on="ko", how="left").sort_values("prevalence", ascending=False)
    write_csv(per_ko, "ko_prevalence")

    # --- per-sample richness
    rich = c.execute("""SELECT accession, environment, count(*) AS n_kos FROM ko JOIN coh USING (accession)
                        GROUP BY ALL""").df()
    write_csv(rich.n_kos.describe(percentiles=[.05, .25, .5, .75, .95]).round(1).reset_index()
              .rename(columns={"index": "stat", "n_kos": "kos_per_sample"}), "ko_richness_summary")

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), gridspec_kw={"width_ratios": [1, 1, 1.4]})
    axes[0].hist(rich.n_kos, bins=60, color=BLUE)
    axes[0].set_xlabel("KOs detected per sample")
    axes[0].set_ylabel("samples")
    axes[0].set_title(f"KO richness (median {int(rich.n_kos.median()):,}; "
                      f"database has {n_db:,} KOs)", loc="left")
    axes[1].hist(per_ko.prevalence * 100, bins=50, color=BLUE)
    axes[1].set_yscale("log")
    axes[1].set_xlabel("% of samples in which the KO is detected")
    axes[1].set_ylabel("KOs (log)")
    core = int((per_ko.prevalence >= 0.9).sum())
    rare = int((per_ko.prevalence < 0.01).sum())
    axes[1].set_title(f"KO prevalence: {len(per_ko):,} KOs seen; {core:,} in ≥90% of samples,\n"
                      f"{rare:,} in <1%", loc="left")
    top_env = rich.environment.value_counts().head(N_ENV).index
    data = [rich.n_kos[rich.environment == e] for e in top_env]
    axes[2].boxplot(data, vert=False, showfliers=False, widths=0.6,
                    medianprops={"color": ORANGE}, boxprops={"color": BLUE})
    axes[2].set_yticks(range(1, len(top_env) + 1), [f"{e} ({len(d):,})" for e, d in zip(top_env, data)])
    axes[2].invert_yaxis()
    axes[2].set_xlabel("KOs detected per sample")
    axes[2].set_title(f"KO richness by environment (top {N_ENV})", loc="left")
    fig.tight_layout()
    save(fig, "ko_01_richness_prevalence")

    # --- top KOs by prevalence and by mean abundance
    fig, axes = plt.subplots(1, 2, figsize=(15, 8))
    lab = lambda d: [f"{k} {str(s)[:28]}" for k, s in zip(d.ko, d.symbol.fillna(""))]  # noqa: E731
    t = per_ko.head(30)
    barh(axes[0], lab(t), t.prevalence * 100, fmt="{:.1f}%")
    axes[0].set_title("Top 30 KOs by prevalence (% samples detected)", loc="left")
    t = per_ko.sort_values("mean_abund_overall", ascending=False).head(30)
    barh(axes[1], lab(t), t.mean_abund_overall * 100, color=ORANGE, fmt="{:.2f}%")
    axes[1].set_title("Top 30 KOs by mean relative abundance (all samples)", loc="left")
    for ax in axes:
        ax.tick_params(axis="y", labelsize=7)
    fig.tight_layout()
    save(fig, "ko_02_top_kos")

    # --- BRITE level-B composition by environment (abundance split evenly across a KO's categories)
    b = brite.copy()
    b["cat"] = np.where(b.level_a.str[:5].isin(PATHWAY_A), b.level_b.str[6:],
                        np.where(b.level_a.str.startswith("09180"), "BRITE hierarchies only", "not in pathway/BRITE"))
    kb = b[["ko", "cat", "level_a"]].drop_duplicates(["ko", "cat"])
    # Prefer pathway categories: a KO that is in a pathway is not also counted as "BRITE only".
    in_path = set(kb.ko[kb.level_a.str[:5].isin(PATHWAY_A)])
    kb = kb[~(kb.ko.isin(in_path) & ~kb.level_a.str[:5].isin(PATHWAY_A))]
    kb["w"] = 1 / kb.groupby("ko").ko.transform("size")
    c.register("kb", kb[["ko", "cat", "w"]])
    comp = c.execute("""
        SELECT coh.environment, coalesce(kb.cat, 'not in KEGG BRITE') AS cat,
               sum(ko.abundance * coalesce(kb.w, 1)) / count(DISTINCT ko.accession) AS share
        FROM ko JOIN coh USING (accession) LEFT JOIN kb USING (ko) GROUP BY ALL""").df()
    allcomp = c.execute("""
        SELECT coalesce(kb.cat, 'not in KEGG BRITE') AS cat,
               sum(ko.abundance * coalesce(kb.w, 1)) / count(DISTINCT ko.accession) AS share
        FROM ko LEFT JOIN kb USING (ko) GROUP BY ALL ORDER BY share DESC""").df()
    write_csv(allcomp, "ko_brite_category_share_all")
    wide = comp.pivot(index="environment", columns="cat", values="share").loc[top_env]
    write_csv(wide.reset_index(), "ko_brite_category_share_by_environment")
    cats = allcomp.cat.head(16).tolist()
    w = wide[cats]
    fig, ax = plt.subplots(figsize=(13, 6.5))
    im = ax.imshow((w / w.mean()).values, cmap="RdBu_r", vmin=0.5, vmax=1.5, aspect="auto")
    ax.set_xticks(range(len(cats)), [f"{c_} ({allcomp.share[i] * 100:.1f}%)" for i, c_ in enumerate(cats)],
                  rotation=40, ha="right", fontsize=7)
    ax.set_yticks(range(len(w)), w.index, fontsize=8)
    for i in range(w.shape[0]):
        for j in range(w.shape[1]):
            ax.text(j, i, f"{w.values[i, j] * 100:.1f}", ha="center", va="center", fontsize=6)
    fig.colorbar(im, ax=ax, label="share relative to mean across environments", shrink=0.8)
    ax.set_title("KEGG BRITE category share of KO abundance (%) by environment\n"
                 "(columns: 16 largest categories, overall share in parentheses; colour = over/under-representation)",
                 loc="left")
    save(fig, "ko_03_brite_by_environment")

    # --- KOs that distinguish environments: prevalence per environment for the most variable KOs
    prev_env = c.execute(f"""
        WITH e AS (SELECT environment, count(*) AS n_env FROM coh
                   WHERE environment IN ({",".join("'" + e.replace("'", "''") + "'" for e in top_env)})
                   GROUP BY 1)
        SELECT e.environment, ko.ko, count(*) / any_value(e.n_env) AS prev
        FROM ko JOIN coh USING (accession) JOIN e USING (environment) GROUP BY ALL""").df()
    pw = prev_env.pivot(index="ko", columns="environment", values="prev").fillna(0)[list(top_env)]
    pw = pw[pw.max(axis=1) >= 0.3]
    sel = pw.std(axis=1).sort_values(ascending=False).head(40).index
    out = pw.loc[sel].reset_index().merge(names[["ko", "symbol", "definition"]], on="ko", how="left")
    write_csv(out, "ko_most_variable_by_environment")
    h = pw.loc[sel]
    fig, ax = plt.subplots(figsize=(11, 11))
    im = ax.imshow(h.values * 100, cmap="viridis", vmin=0, vmax=100, aspect="auto")
    sym = names.set_index("ko").symbol
    ax.set_yticks(range(len(h)), [f"{k} {str(sym.get(k, ''))[:30]}" for k in h.index], fontsize=7)
    ax.set_xticks(range(h.shape[1]), h.columns, rotation=40, ha="right", fontsize=8)
    fig.colorbar(im, ax=ax, label="% of samples in environment with KO detected", shrink=0.6)
    ax.set_title("40 KOs whose prevalence differs most across the top environments", loc="left")
    save(fig, "ko_04_variable_kos_by_environment")


if __name__ == "__main__":
    main()
