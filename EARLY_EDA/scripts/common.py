"""Shared paths and helpers for EARLY_EDA."""
from pathlib import Path

import duckdb
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

EDA = Path(__file__).resolve().parent.parent
REPO = EDA.parent
FIG = EDA / "figures"
TAB = EDA / "tables"
# Harmonized per-sample metadata + output status built by analysis/environment_discordance.
ED = REPO / "analysis" / "environment_discordance"
SAMPLES = ED / "samples.parquet"
KO_PROFILES = ED / "ko_profiles.parquet"
TAXA_PROFILES = ED / "taxa_profiles.parquet"
COHORT = EDA / "cohort.parquet"
# Downsampled per-sample hash sets (large, so on /scratch).
HASH_CACHE = Path("/scratch/lbt5343/early_eda_cache")
KO_SIG = Path("/home/grads/lbt5343/KOs_sketched_scaled_1000.sig.zip")
GTDB_SIGS = Path("/scratch/shared_data/YACHT_pretrained_GTDB_r232/gtdb_r232_ani_thresh_0.95_intermediate_files")
TAXDUMP = Path("/scratch/dmk333/kraken2_db")

# Metadata fields harmonized upstream (see analysis/environment_discordance/data_dictionary.md).
# ID-like fields (subject_id, timepoint, host_taxid, lat_lon) are left out on purpose.
FIELDS = [
    "host", "body_site", "sample_material", "isolation_source", "env_biome", "env_feature",
    "env_material", "disease", "antibiotic", "treatment", "diet", "age", "sex", "bmi",
    "geo_loc_name", "collection_date", "ph", "temperature", "salinity", "depth",
    "samp_mat_process", "samp_store_temp", "nucleic_acid_extraction",
]
# Column holding the normalized value, where it differs from the field name.
VALUE_COL = {"geo_loc_name": "country"}

plt.rcParams.update({
    "figure.dpi": 120, "savefig.dpi": 150, "savefig.bbox": "tight",
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 9,
})
BLUE, ORANGE, GREY = "#2f6db3", "#d9822b", "#9aa3ad"


def con():
    c = duckdb.connect()
    c.execute("SET enable_progress_bar=false")
    c.execute("SET threads=32")
    return c


def save(fig, name):
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / f"{name}.png")
    plt.close(fig)
    print(f"  wrote figures/{name}.png")


def write_csv(df, name):
    TAB.mkdir(exist_ok=True)
    df.to_csv(TAB / f"{name}.csv", index=False)
    print(f"  wrote tables/{name}.csv ({len(df):,} rows)")


def barh(ax, labels, values, color=BLUE, fmt="{:,}"):
    """Horizontal bar chart, largest at top, value printed at bar end."""
    y = range(len(labels))[::-1]
    ax.barh(list(y), values, color=color)
    ax.set_yticks(list(y), labels)
    vmax = max(values) if len(values) else 1
    for yi, v in zip(y, values):
        ax.text(v + vmax * 0.01, yi, fmt.format(v), va="center", fontsize=7)
    ax.set_xlim(0, vmax * 1.15)
