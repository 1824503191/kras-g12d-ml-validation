from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/myriadfs/home/ucapmge/Scratch/omm_pipeline_screening")
TABLE_DIR = ROOT / "analysis_batch0092" / "tables"
FIG_DIR = ROOT / "analysis_batch0092" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

matched = pd.read_parquet(TABLE_DIR / "batch0092_features_with_barriers.parquet")

# Define barrier groups
def assign_group(x):
    if x < 28:
        return "low"
    elif x >= 29:
        return "high"
    else:
        return "middle"

matched["barrier_group"] = matched["barrier_energy"].apply(assign_group)

print("Barrier group counts:")
print(matched["barrier_group"].value_counts())

matched[["ligand_id", "barrier_energy", "barrier_group"]].to_csv(
    TABLE_DIR / "batch0092_ligand_groups.tsv",
    sep="\t",
    index=False
)

compare = matched[matched["barrier_group"].isin(["low", "high"])].copy()

print("Low/high comparison rows:", len(compare))
print(compare["barrier_group"].value_counts())

feature_file = TABLE_DIR / "model_used_ligand_related_features.txt"

with open(feature_file) as f:
    features = [line.strip() for line in f if line.strip()]

features = [f for f in features if f in compare.columns]

print("Features to compare:", len(features))

rows = []

for feat in features:
    low_vals = pd.to_numeric(
        compare.loc[compare["barrier_group"] == "low", feat],
        errors="coerce"
    ).dropna()

    high_vals = pd.to_numeric(
        compare.loc[compare["barrier_group"] == "high", feat],
        errors="coerce"
    ).dropna()

    if len(low_vals) < 2 or len(high_vals) < 2:
        continue

    low_mean = low_vals.mean()
    high_mean = high_vals.mean()
    low_median = low_vals.median()
    high_median = high_vals.median()

    diff = high_mean - low_mean

    pooled_sd = np.sqrt(
        ((len(low_vals) - 1) * low_vals.var(ddof=1) +
         (len(high_vals) - 1) * high_vals.var(ddof=1))
        / (len(low_vals) + len(high_vals) - 2)
    )

    cohens_d = diff / pooled_sd if pooled_sd > 0 else np.nan

    if diff < 0:
        interpretation = "high group has smaller distance; high-barrier ligands are closer to this residue"
    else:
        interpretation = "high group has larger distance; high-barrier ligands are farther from this residue"

    rows.append({
        "feature": feat,
        "low_n": len(low_vals),
        "high_n": len(high_vals),
        "low_mean": low_mean,
        "high_mean": high_mean,
        "low_median": low_median,
        "high_median": high_median,
        "high_minus_low_mean": diff,
        "cohens_d": cohens_d,
        "abs_cohens_d": abs(cohens_d) if pd.notna(cohens_d) else np.nan,
        "interpretation": interpretation,
    })

summary = pd.DataFrame(rows).sort_values("abs_cohens_d", ascending=False)

summary.to_csv(
    TABLE_DIR / "batch0092_model_used_ligand_features_summary.tsv",
    sep="\t",
    index=False
)

print("Top features:")
print(summary.head(20))

# Plot top 20 Cohen's d
top = summary.head(20).copy()

plt.figure(figsize=(9, 7))
plt.barh(top["feature"][::-1], top["cohens_d"][::-1])
plt.axvline(0, color="black", linewidth=1)
plt.xlabel("Cohen's d, high minus low")
plt.ylabel("Model-used PRO-LIG feature")
plt.title("Batch0092: low vs high barrier feature differences")
plt.tight_layout()

plt.savefig(FIG_DIR / "batch0092_model_used_ligand_features_cohens_d.png", dpi=300)
plt.savefig(FIG_DIR / "batch0092_model_used_ligand_features_cohens_d.svg")

# Boxplots for top 8 features
for feat in top["feature"].head(8):
    data_low = pd.to_numeric(
        compare.loc[compare["barrier_group"] == "low", feat],
        errors="coerce"
    ).dropna()

    data_high = pd.to_numeric(
        compare.loc[compare["barrier_group"] == "high", feat],
        errors="coerce"
    ).dropna()

    plt.figure(figsize=(5, 4))
    plt.boxplot([data_low, data_high], labels=["low", "high"])
    plt.ylabel("Distance feature value")
    plt.title(feat)
    plt.tight_layout()

    safe_feat = feat.replace("/", "_").replace(" ", "_")
    plt.savefig(FIG_DIR / f"boxplot_batch0092_{safe_feat}.png", dpi=300)
    plt.savefig(FIG_DIR / f"boxplot_batch0092_{safe_feat}.svg")
    plt.close()

print("Saved summary and figures.")
