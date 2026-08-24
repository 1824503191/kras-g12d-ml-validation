from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

WORK = Path("/data/student/mingsong/qmmm_validation_47_complete")
PRED = WORK / "predictions"
FIG = WORK / "figures"
FIG.mkdir(exist_ok=True)

df = pd.read_csv(PRED / "old_rf_vs_qmmm_complete_current.csv")

x = df["computed_barrier_kcal_per_mol"].to_numpy(float)
y = df["old_rf_predicted_barrier"].to_numpy(float)
offset = y - x

def boot_ci(x, y, func, n_boot=10000, seed=123):
    rng = np.random.default_rng(seed)
    vals = []
    n = len(x)

    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        try:
            val = func(x[idx], y[idx])
            if not np.isnan(val):
                vals.append(val)
        except Exception:
            pass

    return np.percentile(vals, [2.5, 97.5])

pearson = pearsonr(x, y).statistic
spearman = spearmanr(x, y).statistic

pearson_ci = boot_ci(x, y, lambda a, b: pearsonr(a, b).statistic)
spearman_ci = boot_ci(x, y, lambda a, b: spearmanr(a, b).statistic)

rng = np.random.default_rng(456)
offset_boot = []
for _ in range(10000):
    idx = rng.integers(0, len(offset), len(offset))
    offset_boot.append(np.mean(offset[idx]))

offset_ci = np.percentile(offset_boot, [2.5, 97.5])

summary = {
    "n": len(df),
    "pearson_r": pearson,
    "pearson_ci_low": pearson_ci[0],
    "pearson_ci_high": pearson_ci[1],
    "spearman_rho": spearman,
    "spearman_ci_low": spearman_ci[0],
    "spearman_ci_high": spearman_ci[1],
    "mean_offset_pred_minus_computed": np.mean(offset),
    "offset_ci_low": offset_ci[0],
    "offset_ci_high": offset_ci[1],
    "mae": np.mean(np.abs(offset)),
    "rmse": np.sqrt(np.mean(offset ** 2)),
    "reference_mean": np.mean(x),
    "reference_sd": np.std(x, ddof=1),
    "reference_min": np.min(x),
    "reference_max": np.max(x),
    "prediction_mean": np.mean(y),
    "prediction_sd": np.std(y, ddof=1),
    "prediction_min": np.min(y),
    "prediction_max": np.max(y),
    "prediction_sd_as_percent_of_reference_sd": 100 * np.std(y, ddof=1) / np.std(x, ddof=1),
    "prediction_range_as_percent_of_reference_range": 100 * (np.max(y) - np.min(y)) / (np.max(x) - np.min(x)),
}

summary_df = pd.DataFrame([summary])
summary_df.to_csv(PRED / "summary_stats_complete_current.csv", index=False)

with open(PRED / "summary_stats_complete_current.txt", "w") as f:
    for k, v in summary.items():
        f.write(f"{k}: {v}\n")

print(summary_df.T)

lo = min(x.min(), y.min()) - 0.5
hi = max(x.max(), y.max()) + 0.5

plt.figure(figsize=(6, 5))
plt.scatter(x, y)
plt.plot([lo, hi], [lo, hi], linestyle="--")
plt.xlabel("Computed ligand-in-QM barrier (kcal/mol)")
plt.ylabel("Old RF predicted barrier (kcal/mol)")
plt.title("Old RF prediction vs complete QM/MM barriers")
plt.tight_layout()
plt.savefig(FIG / "old_rf_vs_qmmm_complete_current_scatter_raw.png", dpi=300)
plt.savefig(FIG / "old_rf_vs_qmmm_complete_current_scatter_raw.pdf")
plt.close()

mean_offset = np.mean(offset)
y_corr = y - mean_offset

lo2 = min(x.min(), y_corr.min()) - 0.5
hi2 = max(x.max(), y_corr.max()) + 0.5

plt.figure(figsize=(6, 5))
plt.scatter(x, y_corr)
plt.plot([lo2, hi2], [lo2, hi2], linestyle="--")
plt.xlabel("Computed ligand-in-QM barrier (kcal/mol)")
plt.ylabel("Old RF prediction minus mean offset (kcal/mol)")
plt.title("Offset-corrected old RF prediction")
plt.tight_layout()
plt.savefig(FIG / "old_rf_vs_qmmm_complete_current_offset_corrected.png", dpi=300)
plt.savefig(FIG / "old_rf_vs_qmmm_complete_current_offset_corrected.pdf")
plt.close()

plt.figure(figsize=(6, 5))
plt.scatter(x, offset)
plt.axhline(0, linestyle="--")
plt.xlabel("Computed ligand-in-QM barrier (kcal/mol)")
plt.ylabel("Prediction − computed (kcal/mol)")
plt.title("Old RF residuals")
plt.tight_layout()
plt.savefig(FIG / "old_rf_vs_qmmm_complete_current_residuals.png", dpi=300)
plt.savefig(FIG / "old_rf_vs_qmmm_complete_current_residuals.pdf")
plt.close()

print("\nSaved summary and figures.")
