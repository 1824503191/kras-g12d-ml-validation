from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from scipy.stats import ttest_ind
    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False

ROOT = Path("/myriadfs/home/ucapmge/Scratch/omm_pipeline_screening")
TABLE_DIR = ROOT / "analysis_batch0092" / "tables"
FIG_DIR = ROOT / "analysis_batch0092" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

LOW_CUTOFF = 28.0
HIGH_CUTOFF = 29.0

MATCHED_FILE = TABLE_DIR / "batch0092_features_with_barriers.parquet"
SELECTED_FILE = TABLE_DIR / "batch0092_PROPRO_selected_feature_frequency.tsv"

OUT_FILE = TABLE_DIR / "batch0092_top_PROPRO_feature_effects.tsv"

print("=== Analyse top PRO-PRO feature effects ===")

df = pd.read_parquet(MATCHED_FILE)
selected = pd.read_csv(SELECTED_FILE, sep="\t")

# Use the k=20 selection result, because PRO-PRO_select20 performed best
selected_k20 = selected[selected["k"] == 20].copy()

selected_k20 = selected_k20.sort_values(
    ["selected_count_out_of_5", "mean_f_score_across_folds"],
    ascending=[False, False]
)

# keep top 30 for interpretation
top = selected_k20.head(30).copy()

print("Top PRO-PRO features to analyse:", len(top))

# Exclude middle band
df = df[
    (df["barrier_energy"] < LOW_CUTOFF)
    | (df["barrier_energy"] >= HIGH_CUTOFF)
].copy()

low = df[df["barrier_energy"] < LOW_CUTOFF]
high = df[df["barrier_energy"] >= HIGH_CUTOFF]

print("Low samples:", len(low))
print("High samples:", len(high))

rows = []

for _, r in top.iterrows():
    f = r["feature"]

    if f not in df.columns:
        print("Missing feature in matched table:", f)
        continue

    x_low = pd.to_numeric(low[f], errors="coerce").dropna().values
    x_high = pd.to_numeric(high[f], errors="coerce").dropna().values

    n_low = len(x_low)
    n_high = len(x_high)

    mean_low = np.mean(x_low)
    mean_high = np.mean(x_high)
    median_low = np.median(x_low)
    median_high = np.median(x_high)

    diff = mean_high - mean_low

    sd_low = np.std(x_low, ddof=1)
    sd_high = np.std(x_high, ddof=1)

    pooled_sd = np.sqrt(
        ((n_low - 1) * sd_low**2 + (n_high - 1) * sd_high**2)
        / (n_low + n_high - 2)
    )

    cohens_d = diff / pooled_sd if pooled_sd > 0 else np.nan

    # approximate 95% CI for Cohen's d
    se_d = np.sqrt(
        (n_low + n_high) / (n_low * n_high)
        + (cohens_d**2) / (2 * (n_low + n_high - 2))
    )

    ci_low = cohens_d - 1.96 * se_d
    ci_high = cohens_d + 1.96 * se_d

    if HAS_SCIPY:
        _, p_value = ttest_ind(x_high, x_low, equal_var=False)
    else:
        p_value = np.nan

    rows.append({
        "feature": f,
        "residue_1_name": r.get("residue_1_name", ""),
        "residue_2_name": r.get("residue_2_name", ""),
        "selected_count_out_of_5": r["selected_count_out_of_5"],
        "selection_frequency": r["selection_frequency"],
        "mean_f_score_across_folds": r["mean_f_score_across_folds"],
        "n_low": n_low,
        "n_high": n_high,
        "mean_low": mean_low,
        "mean_high": mean_high,
        "median_low": median_low,
        "median_high": median_high,
        "high_minus_low_mean": diff,
        "cohens_d": cohens_d,
        "cohens_d_ci_low": ci_low,
        "cohens_d_ci_high": ci_high,
        "ci_crosses_zero": ci_low <= 0 <= ci_high,
        "p_value_welch": p_value,
        "direction": "larger distance in high" if diff > 0 else "smaller distance in high"
    })

res = pd.DataFrame(rows)

# Benjamini-Hochberg FDR correction if p-values exist
if HAS_SCIPY and res["p_value_welch"].notna().any():
    res = res.sort_values("p_value_welch").reset_index(drop=True)
    m = len(res)
    ranks = np.arange(1, m + 1)
    raw = res["p_value_welch"].values * m / ranks
    # monotonic BH adjustment
    bh = np.minimum.accumulate(raw[::-1])[::-1]
    res["p_value_bh_fdr"] = np.clip(bh, 0, 1)
else:
    res["p_value_bh_fdr"] = np.nan

# Sort back by selection stability and F-score
res = res.sort_values(
    ["selected_count_out_of_5", "mean_f_score_across_folds"],
    ascending=[False, False]
)

res.to_csv(OUT_FILE, sep="\t", index=False)

print("\nTop PRO-PRO effects:")
show_cols = [
    "feature",
    "residue_1_name",
    "residue_2_name",
    "selected_count_out_of_5",
    "mean_low",
    "mean_high",
    "high_minus_low_mean",
    "cohens_d",
    "cohens_d_ci_low",
    "cohens_d_ci_high",
    "ci_crosses_zero",
    "direction"
]
print(res[show_cols].to_string(index=False))
print("\nSaved:", OUT_FILE)

# Plot Cohen's d with 95% CI for top 20
plot_df = res.head(20).iloc[::-1].copy()

labels = [
    f"{row['residue_1_name']}–{row['residue_2_name']}"
    for _, row in plot_df.iterrows()
]

ypos = np.arange(len(plot_df))
x = plot_df["cohens_d"].values
xerr_low = x - plot_df["cohens_d_ci_low"].values
xerr_high = plot_df["cohens_d_ci_high"].values - x

fig, ax = plt.subplots(figsize=(9, max(6, 0.38 * len(plot_df))))
ax.errorbar(x, ypos, xerr=[xerr_low, xerr_high], fmt="o", capsize=3)
ax.axvline(0, linestyle="--", linewidth=1)
ax.set_yticks(ypos)
ax.set_yticklabels(labels)
ax.set_xlabel("Cohen's d with 95% CI")
ax.set_title("Top selected PRO-PRO distance differences: high vs low RF-predicted barrier")
plt.tight_layout()

for ext in ["png", "svg", "pdf"]:
    fig.savefig(FIG_DIR / f"batch0092_top_PROPRO_cohens_d_with_CI.{ext}", dpi=300)

print("Saved figures:")
print(FIG_DIR / "batch0092_top_PROPRO_cohens_d_with_CI.png")
print(FIG_DIR / "batch0092_top_PROPRO_cohens_d_with_CI.svg")
print(FIG_DIR / "batch0092_top_PROPRO_cohens_d_with_CI.pdf")
print("Done.")
