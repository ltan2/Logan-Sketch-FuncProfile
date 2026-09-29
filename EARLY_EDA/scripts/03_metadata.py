"""Metadata EDA on the cohort.

1. Which metadata is present: completeness of each harmonized field, and the raw BioSample
   attribute keys submitters actually filled in.
2. For each field, what the values are (top values after light cleanup).
3. Numeric fields, collection year, sequencing platform, and a lat/lon map.
"""
import json
import re
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import BLUE, COHORT, FIELDS, GREY, ORANGE, VALUE_COL, barh, con, save, write_csv

STATUS_ORDER = ["value", "explicit_negative", "unknown", "not_applicable", "absent"]
STATUS_COLOR = {"value": BLUE, "explicit_negative": "#7fb3e6", "unknown": ORANGE,
                "not_applicable": "#e8c07d", "absent": "#e3e6ea"}

# Raw attribute keys that are identifiers or bookkeeping, not sample descriptions.
ID_KEY = re.compile(
    r"^(id_|sample_name$|title$|sample_title$|description$|sample_comment$|alias$|"
    r"host_subject_id$|subject_id$|submitter_id$|biosample|bioproject|external_id|"
    r"ena_|insdc_|sra_|gold_stamp_id|primary_search$|bases$|bytes$|run_file_|"
    r"lat_lon$|.*_accession$|.*_id$)")

# Light, documented cleanup so trivial spelling variants are counted together.
SYNONYMS = {
    "stool": "feces", "faeces": "feces", "fecal material": "feces", "faecal": "feces",
    "fecal": "feces", "fecal sample": "feces", "stool sample": "feces", "faecal material": "feces",
    "human-gut": "human gut", "human gut microbiome": "human gut",
    "human_stool": "feces", "human stool": "feces",
    "uc": "ulcerative colitis", "cd": "crohn's disease", "crohn disease": "crohn's disease",
    "crohn": "crohn's disease", "crohns disease": "crohn's disease", "crohn's": "crohn's disease",
    "not aplicable": "not applicable",
    "ibd": "inflammatory bowel disease",
}
CATEGORICAL = [
    ("Environment", ["environment", "site_category", "sample_material", "isolation_source",
                     "env_biome", "env_feature", "env_material"]),
    ("Host and clinical", ["host", "disease", "antibiotic", "treatment", "diet", "sex"]),
    ("Geography and sequencing", ["country", "sra_geo_loc_name_country_continent_calc",
                                  "sra_platform", "sra_instrument", "sra_library_layout",
                                  "samp_mat_process"]),
]
PRETTY = {"sra_geo_loc_name_country_continent_calc": "continent (SRA calc)"}
NUMERIC = [("age", "age (years)", (0, 100)), ("bmi", "BMI", (10, 60)), ("ph", "pH", (0, 14)),
           ("temperature", "temperature (°C)", (-5, 100)), ("depth", "depth", None),
           ("sra_mbases", "sequenced Mbases (log10)", None)]


def clean(v):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return None
    v = re.sub(r"\s*\[?(envo|uberon|efo)[:_]\d+\]?", "", str(v).lower()).strip(" ;,_")
    return SYNONYMS.get(v, v) or None


def completeness(df):
    rows = []
    for f in FIELDS:
        st = df[f"{f}_status"].fillna("absent").value_counts(normalize=True) * 100
        rows.append({"field": f, **{s: round(st.get(s, 0.0), 2) for s in STATUS_ORDER},
                     "n_with_value": int((df[f"{f}_status"] == "value").sum())})
    t = pd.DataFrame(rows).sort_values("value", ascending=False)
    write_csv(t, "metadata_field_completeness")

    fig, ax = plt.subplots(figsize=(8, 7))
    left = np.zeros(len(t))
    y = np.arange(len(t))[::-1]
    for s in STATUS_ORDER:
        ax.barh(y, t[s], left=left, color=STATUS_COLOR[s], label=s.replace("_", " "))
        left += t[s].values
    ax.set_yticks(y, [f"{f}  ({n:,})" for f, n in zip(t.field, t.n_with_value)])
    ax.set_xlabel(f"% of {len(df):,} samples")
    ax.set_xlim(0, 100)
    ax.legend(ncol=5, loc="lower center", bbox_to_anchor=(0.4, 1.0), frameon=False, fontsize=8)
    ax.set_title("Harmonized metadata fields: how often each is filled in\n"
                 "(count = samples with a real value)", loc="left", pad=28)
    save(fig, "meta_01_field_completeness")


def raw_keys(df):
    """Raw attribute keys from the BioSample record (falls back to the SRA jattr)."""
    from importlib import util
    spec = util.spec_from_file_location(
        "fields", str(COHORT.parent.parent / "analysis/environment_discordance/scripts/fields.py"))
    fields = util.module_from_spec(spec)
    spec.loader.exec_module(fields)
    cnt = Counter()
    for bs, jattr in zip(df.biosample_attributes_json, df.sra_jattr):
        keys = set()
        for blob in (bs, jattr):
            if isinstance(blob, str):
                keys |= {fields.canon_key(k) for k in json.loads(blob)}
        cnt.update(k for k in keys if k and not ID_KEY.match(k))
    t = pd.DataFrame(cnt.most_common(), columns=["attribute_key", "n_samples"])
    t["pct_samples"] = (t.n_samples / len(df) * 100).round(2)
    write_csv(t, "metadata_raw_attribute_keys")
    top = t.head(40)
    fig, ax = plt.subplots(figsize=(7, 9))
    barh(ax, top.attribute_key, top.pct_samples, fmt="{:.0f}%")
    ax.set_xlabel("% of samples with this attribute key")
    ax.set_title(f"Top 40 raw BioSample/SRA attribute keys (IDs excluded)\n"
                 f"{len(t):,} distinct keys in total", loc="left")
    save(fig, "meta_02_raw_attribute_keys")


def top_values(df):
    rows = []
    for group, cols in CATEGORICAL:
        fig, axes = plt.subplots(int(np.ceil(len(cols) / 2)), 2, figsize=(12, 3.3 * np.ceil(len(cols) / 2)))
        for ax, col in zip(axes.flat, cols):
            vals = df[col].map(clean).dropna()
            vc = vals.value_counts()
            for v, n in vc.items():
                rows.append({"group": group, "field": PRETTY.get(col, col), "value": v, "n_samples": n,
                             "pct_of_samples_with_field": round(n / len(vals) * 100, 2)})
            top = vc.head(15)
            barh(ax, [s[:45] for s in top.index], top.values)
            ax.set_title(f"{PRETTY.get(col, col)}: {len(vals):,} samples "
                         f"({len(vals) / len(df):.0%}), {len(vc):,} distinct", loc="left", fontsize=9)
            ax.tick_params(axis="y", labelsize=7)
        for ax in list(axes.flat)[len(cols):]:
            ax.axis("off")
        fig.suptitle(f"{group}: top 15 values per field", x=0.01, ha="left", fontsize=11)
        fig.tight_layout()
        save(fig, f"meta_03_values_{group.split()[0].lower()}")
    write_csv(pd.DataFrame(rows), "metadata_field_values")


def disease_detail(df):
    """Disease is sparse, so show it by environment too."""
    d = df.assign(disease_c=df.disease.map(clean))
    t = (d[d.disease_status.isin(["value", "explicit_negative"])]
         .assign(disease_c=lambda x: x.disease_c.where(x.disease_status == "value", "healthy / none"))
         .groupby(["disease_c", "environment"]).size().reset_index(name="n_samples")
         .sort_values("n_samples", ascending=False))
    write_csv(t, "metadata_disease_by_environment")
    by_study = (d[d.disease_status == "value"].groupby("disease_c")
                .agg(n_samples=("accession", "size"), n_studies=("study", "nunique"))
                .sort_values("n_samples", ascending=False).reset_index())
    write_csv(by_study, "metadata_disease_values")
    top = by_study.head(25)
    fig, ax = plt.subplots(figsize=(7, 7))
    barh(ax, [f"{v[:40]}  ({s} studies)" for v, s in zip(top.disease_c, top.n_studies)], top.n_samples, color=ORANGE)
    n_neg = int((df.disease_status == "explicit_negative").sum())
    ax.set_title(f"Disease: {int((df.disease_status == 'value').sum()):,} samples name a disease, "
                 f"{n_neg:,} say healthy/none\n(top 25; most disease labels come from 1–2 studies)",
                 loc="left")
    save(fig, "meta_04_disease")


def numeric(df):
    fig, axes = plt.subplots(2, 4, figsize=(14, 6))
    axes = axes.flat
    rows = []
    for ax, (col, label, rng) in zip(axes, NUMERIC):
        v = pd.to_numeric(df[col], errors="coerce").dropna()
        if col == "sra_mbases":
            v = np.log10(v[v > 0])
        if rng:
            v = v[(v >= rng[0]) & (v <= rng[1])]
        elif col == "depth":
            v = v[(v >= 0) & (v <= v.quantile(0.99))]
        rows.append({"field": col, "n": len(v), **v.describe(percentiles=[.05, .25, .5, .75, .95]).round(2).to_dict()})
        ax.hist(v, bins=40, color=BLUE)
        ax.set_title(f"{label}  (n={len(v):,})", loc="left")
    year = pd.to_numeric(df.collection_date.str[:4], errors="coerce")
    year = year[(year >= 2000) & (year <= 2026)]
    axes[6].hist(year, bins=np.arange(2000, 2027) - 0.5, color=BLUE)
    axes[6].set_title(f"collection year  (n={len(year):,})", loc="left")
    rel = df.sra_releasedate.dt.year.dropna()
    axes[7].hist(rel, bins=np.arange(rel.min(), rel.max() + 2) - 0.5, color=GREY)
    axes[7].set_title(f"SRA release year  (n={len(rel):,})", loc="left")
    fig.suptitle("Numeric metadata (out-of-range values dropped: age 0–100, BMI 10–60, pH 0–14)",
                 x=0.01, ha="left")
    fig.tight_layout()
    save(fig, "meta_05_numeric")
    write_csv(pd.DataFrame(rows), "metadata_numeric_summary")


def richness_and_map(df):
    n_fields = (df[[f"{f}_status" for f in FIELDS]] == "value").sum(axis=1)
    broad = df.environment.fillna("unknown").str.split(r"[_:]", n=1).str[0]
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.5), gridspec_kw={"width_ratios": [1, 1.6]})
    axes[0].hist(n_fields, bins=np.arange(0, len(FIELDS) + 2) - 0.5, color=BLUE)
    axes[0].set_xlabel(f"# harmonized fields with a value (of {len(FIELDS)})")
    axes[0].set_ylabel("samples")
    axes[0].set_title(f"Metadata richness per sample (median {int(n_fields.median())})", loc="left")
    colors = {"human": BLUE, "animal": ORANGE, "other": "#3a9d5d", "conflict": GREY, "unknown": GREY}
    ll = df.assign(broad=broad).dropna(subset=["lat", "lon"])
    for b, g in ll.groupby("broad"):
        axes[1].scatter(g.lon, g.lat, s=4, alpha=0.4, color=colors.get(b, GREY),
                        label=f"{'environmental/other' if b == 'other' else b} ({len(g):,})")
    axes[1].set_xlim(-180, 180)
    axes[1].set_ylim(-90, 90)
    axes[1].set_xlabel("longitude")
    axes[1].set_ylabel("latitude")
    axes[1].axhline(0, color="#ddd", lw=0.5, zorder=0)
    axes[1].legend(markerscale=4, fontsize=7, frameon=False, loc="lower left")
    axes[1].set_title(f"Sampling locations with lat/lon ({len(ll):,} samples)", loc="left")
    fig.tight_layout()
    save(fig, "meta_06_richness_and_locations")

    env = (df.assign(n_fields=n_fields).groupby("environment")
           .agg(n_samples=("accession", "size"), n_studies=("study", "nunique"),
                median_fields=("n_fields", "median"),
                pct_country=("geo_loc_name_status", lambda s: round((s == "value").mean() * 100, 1)),
                pct_disease=("disease_status", lambda s: round(s.isin(["value", "explicit_negative"]).mean() * 100, 1)),
                pct_age=("age_status", lambda s: round((s == "value").mean() * 100, 1)),
                pct_env_biome=("env_biome_status", lambda s: round((s == "value").mean() * 100, 1)))
           .sort_values("n_samples", ascending=False).reset_index())
    write_csv(env, "metadata_by_environment")


def main():
    df = con().execute(f"SELECT * FROM '{COHORT}'").df()
    print(f"{len(df):,} cohort samples")
    completeness(df)
    raw_keys(df)
    top_values(df)
    disease_detail(df)
    numeric(df)
    richness_and_map(df)


if __name__ == "__main__":
    main()
