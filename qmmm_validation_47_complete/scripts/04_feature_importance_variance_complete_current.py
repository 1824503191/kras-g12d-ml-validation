from pathlib import Path
import joblib
import numpy as np
import pandas as pd

WORK = Path("/data/student/mingsong/qmmm_validation_47_complete")
FEATURE_DIR = WORK / "features"
PRED = WORK / "predictions"

MODEL_PATH = WORK / "model_features/rf_model_dist_below_10_new.pkl"
FEATURE_LIST_PATH = WORK / "model_features/dist_features_below_10_ang_list_new.npy"

required = np.load(FEATURE_LIST_PATH, allow_pickle=True).astype(str).tolist()
model = joblib.load(MODEL_PATH)

dfs = []
for p in sorted(FEATURE_DIR.glob("*_required_features.csv.gz")):
    df = pd.read_csv(p)
    molid = df["Molecule"].iloc[0]
    row = df[required].copy()
    row.insert(0, "ligand_id", molid)
    dfs.append(row)

all_df = pd.concat(dfs, ignore_index=True)
X = all_df[required]

sd = X.std(axis=0, ddof=0)
mean = X.mean(axis=0)
minv = X.min(axis=0)
maxv = X.max(axis=0)

importances = getattr(model, "feature_importances_", None)
if importances is None:
    raise RuntimeError("Model has no feature_importances_ attribute.")

out = pd.DataFrame({
    "feature": required,
    "importance": importances,
    "mean_across_structures": mean.values,
    "sd_across_structures": sd.values,
    "min_across_structures": minv.values,
    "max_across_structures": maxv.values,
})

out["is_fixed_sd_le_1e_6"] = out["sd_across_structures"] <= 1e-6

out["feature_group"] = np.where(
    out["feature"].str.contains("LIG"), "PRO-LIG",
    np.where(out["feature"].str.contains("GTP"), "PRO-GTP", "PRO-PRO")
)

total_imp = out["importance"].sum()

fixed_imp = out.loc[out["is_fixed_sd_le_1e_6"], "importance"].sum()
propro_imp = out.loc[out["feature_group"] == "PRO-PRO", "importance"].sum()
progtp_imp = out.loc[out["feature_group"] == "PRO-GTP", "importance"].sum()
prolig_imp = out.loc[out["feature_group"] == "PRO-LIG", "importance"].sum()

fixed_propro_imp = out.loc[
    (out["feature_group"] == "PRO-PRO") & (out["is_fixed_sd_le_1e_6"]),
    "importance"
].sum()

out = out.sort_values("importance", ascending=False)
out.to_csv(PRED / "feature_importance_variance_complete_current.csv", index=False)

summary_lines = []
summary_lines.append(f"Number of structures: {len(all_df)}")
summary_lines.append(f"Number of required features: {len(required)}")
summary_lines.append(f"Total importance: {total_imp:.6f}")
summary_lines.append("")
summary_lines.append(f"Importance on fixed features SD<=1e-6: {fixed_imp:.6f} ({100*fixed_imp/total_imp:.2f}%)")
summary_lines.append(f"Importance on PRO-PRO features: {propro_imp:.6f} ({100*propro_imp/total_imp:.2f}%)")
summary_lines.append(f"Importance on PRO-GTP features: {progtp_imp:.6f} ({100*progtp_imp/total_imp:.2f}%)")
summary_lines.append(f"Importance on PRO-LIG features: {prolig_imp:.6f} ({100*prolig_imp/total_imp:.2f}%)")
summary_lines.append(f"Importance on fixed PRO-PRO features: {fixed_propro_imp:.6f} ({100*fixed_propro_imp/total_imp:.2f}%)")
summary_lines.append("")
summary_lines.append("Top 30 features by old RF importance:")
summary_lines.append(out.head(30).to_string(index=False))

with open(PRED / "feature_importance_variance_summary_complete_current.txt", "w") as f:
    f.write("\n".join(summary_lines))

print("\n".join(summary_lines))
print("\nSaved:", PRED / "feature_importance_variance_complete_current.csv")
