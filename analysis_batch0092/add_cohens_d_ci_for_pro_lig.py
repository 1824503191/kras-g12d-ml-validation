from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path("/myriadfs/home/ucapmge/Scratch/omm_pipeline_screening")
TABLE_DIR = ROOT / "analysis_batch0092" / "tables"
FIG_DIR = ROOT / "analysis_batch0092" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

LOW_CUTOFF = 28.0
HIGH_CUTOFF = 29.0

df = pd.read_parquet(TABLE_DIR / "batch0092_features_with_barriers.parquet")

feature_file = TABLE_DIR / "model_used_ligand_related_features.txt"
features = [x.strip() for x in feature_file.read_text().splitlines() if x.strip()]
features = [f for f in features if f in df.columns]

df = df[(df["barrier_energy"] < LOW_CUTOFF) | (df["barrier_energy"] >= HIGH_CUTOFF)].copy()

low = df[df["barrier_energy"] < LOW_CUTOFF]
high = df[df["barrier_energy"] >= HIGH_CUTOFF]

rows = []

for f in features:
    x_low = pd.to_numeric(low[f], errors="coerce").dropna().values
    x_high = pd.to_numeric(high[f], errors="coerce").dropna().values

    n_low = len(x_low)
    n_high = len(x_high)

    mean_low = np.mean(x_low)
    mean_high = np.mean(x_high)
    diff = mean_high - mean_low

    sd_low = np.std(x_low, ddof=1)
    sd_high = np.std(x_high, ddof=1)

    pooled_sd = np.sqrt(
        ((n_low - 1) * sd_low**2 + (n_high - 1) * sd_high**2)
        / (n_low + n_high - 2)
    )

    d = diff / pooled_sd if pooled_sd > 0 else np.nan

    # approximate SE and 95% CI for Cohen's d
    se_d = np.sqrt((n_low + n_high) / (n_low * n_high) + (d**2) / (2 * (n_low + n_high - 2)))
    ci_low = d - 1.96 * se_d
    ci_high = d + 1.96 * se_d

    rows.append({
        "feature": f,
        "n_low": n_low,
        "n_high": n_high,
        "mean_low": mean_low,
        "mean_high": mean_high,
        "high_minus_low_mean": diff,
        "cohens_d": d,
        "cohens_d_ci_low": ci_low,
        "cohens_d_ci_high": ci_high,
        "ci_crosses_zero": ci_low <= 0 <= ci_high,
        "interpretation": "high farther" if diff > 0 else "high closer"
    })

res = pd.DataFrame(rows)
res = res.sort_values("cohens_d", key=lambda s: s.abs(), ascending=False)

out = TABLE_DIR / "batch0092_PRO_LIG_cohens_d_with_CI.tsv"
res.to_csv(out, sep="\t", index=False)

print(res.to_string(index=False))
print("Saved:", out)

# plot
plot_df = res.iloc[::-1].copy()
ypos = np.arange(len(plot_df))
x = plot_df["cohens_d"].values
xerr_low = x - plot_df["cohens_d_ci_low"].values
xerr_high = plot_df["cohens_d_ci_high"].values - x

fig, ax = plt.subplots(figsize=(8, max(5, 0.35 * len(plot_df))))
ax.errorbar(x, ypos, xerr=[xerr_low, xerr_high], fmt="o", capsize=3)
ax.axvline(0, linestyle="--", linewidth=1)
ax.set_yticks(ypos)
ax.set_yticklabels(plot_df["feature"])
ax.set_xlabel("Cohen's d with 95% CI")
ax.set_title("PRO-LIG differences: high vs low RF-predicted barrier ligands")
plt.tight_layout()

for ext in ["png", "svg", "pdf"]:
    fig.savefig(FIG_DIR / f"batch0092_PRO_LIG_cohens_d_with_CI.{ext}", dpi=300)

print("Saved figures to:", FIG_DIR)
print("Done.")
