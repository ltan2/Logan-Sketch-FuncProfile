"""Taxonomy EDA from YACHT (GTDB r232 reps, ANI 0.95) on the cohort's unitigs.

YACHT reports genomes it judges present (`in_sample_est`) at each `min_coverage` setting; it is
a presence call on the assembly, not cell abundance, so everything here is prevalence
(fraction of samples in which a taxon is detected). Only ~60% of the cohort has a YACHT
profile (it comes from the lab's shared Logan YACHT database, not from this pipeline).

Headline setting: min_coverage = HEADLINE (the multiview pipeline's choice). Lineages come
from the GTDB r232 taxonomy files (same release as the YACHT reference); YACHT genomes
that are not r232 representatives fall back to their NCBI taxid lineage (taxonkit), marked
lineage_source = 'ncbi_taxid' in taxa_genome_prevalence.csv.
"""
from fractions import Fraction
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import BLUE, COHORT, GREY, ORANGE, TAXA_PROFILES, TAXDUMP, barh, con, save, write_csv

HEADLINE = 0.0625
N_ENV = 12
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
GTDB_TAX = Path("/scratch/akn5655/gtdb232/gtdb_genomes_reps_r232")


def lineages():
    """GTDB r232 lineage per genome accession keyed by accession number (GCA/GCF prefix and version dropped)."""
    parts = [pd.read_csv(GTDB_TAX / f, sep="\t", header=None, names=["organism_id", "lin"])
             for f in ("bac120_taxonomy_r232.tsv", "ar53_taxonomy_r232.tsv")]
    df = pd.concat(parts)
    df["organism_id"] = df.organism_id.str[3:]
    ranks = df.lin.str.split(";", expand=True).apply(lambda col: col.str[3:])
    ranks.columns = RANKS
    # Join on the accession number without GCA/GCF prefix and version, so a GenBank/RefSeq or
    # version mismatch between YACHT's genome and the GTDB record still matches.
    out = pd.concat([df.organism_id.str[4:].str.split(".").str[0].rename("acc_core"),
                     ranks.replace("", np.nan)], axis=1)
    return out.drop_duplicates("acc_core")


def ncbi_lineages(taxids):
    """Fallback for YACHT genomes that are not GTDB r232 representatives: NCBI lineage via taxonkit."""
    fmt = "\t".join("{" + r + "}" for r in RANKS).replace("{domain}", "{domain|superkingdom}")
    out = subprocess.run(["taxonkit", "reformat2", "--data-dir", str(TAXDUMP), "-I", "1", "-f", fmt],
                         input="\n".join(map(str, taxids)), capture_output=True, text=True, check=True).stdout
    df = pd.DataFrame([line.split("\t") for line in out.splitlines() if line], columns=["tax_id"] + RANKS)
    df["tax_id"] = df.tax_id.astype("int64")
    return df.replace("", np.nan)


def main():
    c = con()
    c.execute(f"CREATE TABLE coh AS SELECT accession, environment, study FROM '{COHORT}'")
    n_coh = c.execute("SELECT count(*) FROM coh").fetchone()[0]
    c.execute(f"""CREATE TABLE tx AS SELECT t.accession, t.min_coverage, t.organism_id, t.organism_name,
                  t.tax_id, t.containment FROM '{TAXA_PROFILES}' t JOIN coh USING (accession)
                  WHERE t.seq_type = 'unitigs'""")

    # --- coverage-setting comparison
    cov = c.execute(f"""
        WITH per AS (SELECT min_coverage, accession, count(*) AS n FROM tx GROUP BY ALL)
        SELECT min_coverage, count(*) AS samples_with_taxa, round(100.0 * count(*) / {n_coh}, 1) AS pct_cohort,
               median(n) AS median_genomes_per_sample, quantile_cont(n, 0.9) AS p90_genomes_per_sample,
               (SELECT count(DISTINCT organism_id) FROM tx t2 WHERE t2.min_coverage = per.min_coverage) AS distinct_genomes
        FROM per GROUP BY min_coverage ORDER BY min_coverage""").df()
    write_csv(cov, "taxa_coverage_settings")
    per = c.execute("SELECT min_coverage, accession, count(*) AS n FROM tx GROUP BY ALL").df()
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.3))
    x = np.arange(len(cov))
    labels = [str(Fraction(v).limit_denominator(1024)) for v in cov.min_coverage]
    axes[0].bar(x, cov.samples_with_taxa, color=BLUE)
    axes[0].axhline(n_coh, color=GREY, ls="--", lw=1)
    axes[0].text(0, n_coh, f" cohort = {n_coh:,}", va="bottom", fontsize=7, color=GREY)
    axes[0].set_title("Cohort samples with ≥1 YACHT genome", loc="left")
    axes[1].boxplot([per.n[per.min_coverage == v] for v in cov.min_coverage], showfliers=False,
                    medianprops={"color": ORANGE})
    axes[1].set_xticks(x + 1, labels)
    axes[1].set_yscale("log")
    axes[1].set_title("Genomes detected per sample", loc="left")
    axes[2].bar(x, cov.distinct_genomes, color=BLUE)
    axes[2].set_title("Distinct GTDB genomes detected (all samples)", loc="left")
    for ax in (axes[0], axes[2]):
        ax.set_xticks(x, labels)
    for ax in axes:
        ax.set_xlabel("YACHT min_coverage")
    fig.suptitle(f"How the YACHT coverage setting changes detections (headline below uses {HEADLINE:g})",
                 x=0.01, ha="left")
    fig.tight_layout()
    save(fig, "tax_01_coverage_settings")

    # --- headline setting: lineages and prevalence
    c.execute(f"CREATE TABLE h AS SELECT * FROM tx WHERE min_coverage = {HEADLINE}")
    n_h = c.execute("SELECT count(DISTINCT accession) FROM h").fetchone()[0]
    lin = lineages()
    c.register("lin", lin)
    c.execute("""CREATE TABLE hl AS SELECT h.*, lin.* EXCLUDE (acc_core) FROM h LEFT JOIN lin
                 ON lin.acc_core = split_part(substr(h.organism_id, 5), '.', 1)""")
    miss = c.execute("SELECT count(DISTINCT organism_id) FILTER (WHERE phylum IS NULL), count(DISTINCT organism_id) FROM hl").fetchone()
    print(f"{miss[0]:,}/{miss[1]:,} detected genomes without a GTDB r232 lineage; using NCBI taxid for those")
    ids = [r[0] for r in c.execute("SELECT DISTINCT tax_id FROM hl WHERE phylum IS NULL AND tax_id IS NOT NULL").fetchall()]
    c.register("ncbi", ncbi_lineages(ids))
    sets = ", ".join(f"\"{r}\" = ncbi.\"{r}\"" for r in RANKS)
    c.execute("ALTER TABLE hl ADD COLUMN lineage_source VARCHAR")
    c.execute("UPDATE hl SET lineage_source = 'gtdb_r232' WHERE phylum IS NOT NULL")
    c.execute(f"UPDATE hl SET {sets}, lineage_source = 'ncbi_taxid' FROM ncbi WHERE hl.phylum IS NULL AND hl.tax_id = ncbi.tax_id")
    print(c.execute("SELECT lineage_source, count(DISTINCT organism_id), count(*) FROM hl GROUP BY 1").fetchall())
    rows = []
    fig, axes = plt.subplots(1, 3, figsize=(18, 8))
    for ax, rank, k in zip(axes, ["phylum", "genus", "species"], [20, 30, 30]):
        t = c.execute(f"""
            SELECT coalesce({rank}, '(unassigned)') AS taxon, count(DISTINCT accession) AS n_samples,
                   count(DISTINCT accession) / {n_h} AS prevalence, count(DISTINCT organism_id) AS n_genomes
            FROM hl GROUP BY 1 ORDER BY n_samples DESC""").df()
        t.insert(0, "rank", rank)
        rows.append(t)
        unassigned = t.prevalence[t.taxon == "(unassigned)"].sum()
        top = t[t.taxon != "(unassigned)"].head(k)
        barh(ax, [s[:38] for s in top.taxon], top.prevalence * 100, fmt="{:.0f}%")
        ax.set_title(f"Top {k} {rank} by prevalence ({len(t) - 1:,} detected)\n"
                     f"{unassigned:.0%} of samples also have a genome with no {rank} name", loc="left")
        ax.set_xlabel(f"% of {n_h:,} profiled samples")
        ax.tick_params(axis="y", labelsize=7)
    fig.suptitle(f"Most represented taxa, GTDB r232 names with NCBI fallback (YACHT min_coverage={HEADLINE:g})", x=0.01, ha="left")
    fig.tight_layout()
    save(fig, "tax_02_top_taxa")
    write_csv(pd.concat(rows), "taxa_prevalence_by_rank")

    genomes = c.execute(f"""
        SELECT organism_id, any_value(organism_name) AS organism_name, any_value(tax_id) AS tax_id,
               any_value(phylum) AS phylum, any_value(genus) AS genus, any_value(species) AS species,
               any_value(lineage_source) AS lineage_source,
               count(DISTINCT accession) AS n_samples, count(DISTINCT accession) / {n_h} AS prevalence,
               median(containment) AS median_containment
        FROM hl GROUP BY organism_id ORDER BY n_samples DESC""").df()
    write_csv(genomes, "taxa_genome_prevalence")

    # --- phylum composition by environment (share of detected genomes, averaged over samples)
    comp = c.execute("""
        WITH s AS (SELECT accession, coalesce(phylum, '(unassigned)') AS phylum, count(*) AS n FROM hl GROUP BY ALL),
             f AS (SELECT accession, phylum, n / sum(n) OVER (PARTITION BY accession) AS frac FROM s)
        SELECT coh.environment, phylum, sum(frac) AS s,
               count(DISTINCT accession) AS dummy FROM f JOIN coh USING (accession) GROUP BY ALL""").df()
    env_n = c.execute("SELECT environment, count(DISTINCT accession) AS n FROM hl JOIN coh USING (accession) GROUP BY 1").df()
    comp = comp.merge(env_n, on="environment")
    comp["share"] = comp.s / comp.n
    top_env = env_n.sort_values("n", ascending=False).environment.head(N_ENV).tolist()
    wide = comp.pivot(index="environment", columns="phylum", values="share").fillna(0).loc[top_env]
    top_phy = wide.mean().sort_values(ascending=False).head(10).index
    wide = wide[top_phy].assign(other=1 - wide[top_phy].sum(axis=1))
    write_csv(wide.reset_index(), "taxa_phylum_share_by_environment")
    fig, ax = plt.subplots(figsize=(12, 6))
    cmap = plt.get_cmap("tab10")
    left = np.zeros(len(wide))
    y = np.arange(len(wide))[::-1]
    for i, p in enumerate(wide.columns):
        ax.barh(y, wide[p], left=left, color=GREY if p == "other" else cmap(i), label=p)
        left += wide[p].values
    ns = env_n.set_index("environment").n
    ax.set_yticks(y, [f"{e} ({ns[e]:,})" for e in wide.index], fontsize=8)
    ax.set_xlim(0, 1)
    ax.set_xlabel("mean share of detected genomes per sample")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False, fontsize=8)
    ax.set_title(f"Phylum make-up of YACHT detections by environment (min_coverage={HEADLINE:g})", loc="left")
    save(fig, "tax_03_phylum_by_environment")

    # --- genera that distinguish environments
    pg = c.execute(f"""
        WITH e AS (SELECT environment, count(DISTINCT accession) AS n_env FROM hl JOIN coh USING (accession) GROUP BY 1)
        SELECT coh.environment, genus, count(DISTINCT accession) / any_value(e.n_env) AS prev
        FROM hl JOIN coh USING (accession) JOIN e USING (environment)
        WHERE genus IS NOT NULL GROUP BY ALL""").df()
    pw = pg.pivot(index="genus", columns="environment", values="prev").fillna(0)[top_env]
    pw = pw[pw.max(axis=1) >= 0.2]
    sel = pw.std(axis=1).sort_values(ascending=False).head(40).index
    write_csv(pw.loc[sel].reset_index(), "taxa_genus_prevalence_by_environment")
    h = pw.loc[sel]
    fig, ax = plt.subplots(figsize=(11, 11))
    im = ax.imshow(h.values * 100, cmap="viridis", vmin=0, vmax=100, aspect="auto")
    ax.set_yticks(range(len(h)), h.index, fontsize=7)
    ax.set_xticks(range(h.shape[1]), h.columns, rotation=40, ha="right", fontsize=8)
    fig.colorbar(im, ax=ax, label="% of profiled samples in environment with genus detected", shrink=0.6)
    ax.set_title("40 genera whose prevalence differs most across the top environments", loc="left")
    save(fig, "tax_04_variable_genera_by_environment")

    # --- which cohort samples have taxonomy at all
    cov_env = c.execute(f"""
        SELECT environment, count(*) AS n_cohort,
               count(*) FILTER (WHERE accession IN (SELECT accession FROM h)) AS n_with_yacht
        FROM coh GROUP BY 1 ORDER BY n_cohort DESC""").df()
    cov_env["pct_with_yacht"] = (cov_env.n_with_yacht / cov_env.n_cohort * 100).round(1)
    write_csv(cov_env, "taxa_profile_coverage_by_environment")
    print(f"{n_h:,}/{n_coh:,} cohort samples have a YACHT profile at min_coverage={HEADLINE:g}")


if __name__ == "__main__":
    main()
