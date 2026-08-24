from pathlib import Path
import numpy as np
import pandas as pd

NEW_DIR = Path("/data/student/mingsong/qmmm_validation_47_complete")
WORK = Path("/data/student/mingsong/high_low_transfer_reduced_features")

FEATURE_DIR = NEW_DIR / "features"
PRED_DIR = NEW_DIR / "predictions"
OUT = WORK / "tables"
OUT.mkdir(parents=True, exist_ok=True)

feature_list_path = NEW_DIR / "model_features/dist_features_below_10_ang_list_new.npy"
required = np.load(feature_list_path, allow_pickle=True).astype(str).tolist()

dfs = []
for p in sorted(FEATURE_DIR.glob("*_required_features.csv.gz")):
    df = pd.read_csv(p)
    molid = df["Molecule"].iloc[0]
    row = df[required].copy()
    row.insert(0, "ligand_id", molid)
    dfs.append(row)

new_X = pd.concat(dfs, ignore_index=True)
X = new_X[required]

sd = X.std(axis=0, ddof=0)
mean = X.mean(axis=0)
minv = X.min(axis=0)
maxv = X.max(axis=0)

feature_info = pd.DataFrame({
    "feature": required,
    "mean_across_new": mean.values,
    "sd_across_new": sd.values,
    "min_across_new": minv.values,
    "max_across_new": maxv.values,
})

feature_info["feature_group"] = np.where(
    feature_info["feature"].str.contains("LIG"), "PRO-LIG",
    np.where(feature_info["feature"].str.contains("GTP"), "PRO-GTP", "PRO-PRO")
)

feature_info["is_varying"] = feature_info["sd_across_new"] > 1e-6

reduced = feature_info.loc[feature_info["is_varying"], "feature"].tolist()

feature_info.to_csv(OUT / "new_structure_feature_variance_all_4086.csv", index=False)

with open(OUT / "reduced_varying_feature_list.txt", "w") as f:
    for feat in reduced:
        f.write(feat + "\n")

pred = pd.read_csv(PRED_DIR / "old_rf_vs_qmmm_complete_current.csv")
pred = pred[["ligand_id", "computed_barrier_kcal_per_mol"]].copy()

pred_sorted = pred.sort_values("computed_barrier_kcal_per_mol").reset_index(drop=True)

n = len(pred_sorted)
half = n // 2

labels = []
for i in range(n):
    if i < half:
        labels.append("low")
    elif i >= n - half:
        labels.append("high")
    else:
        labels.append("middle")

pred_sorted["new_label"] = labels
new_labels = pred_sorted[pred_sorted["new_label"].isin(["low", "high"])].copy()
new_labels["y"] = (new_labels["new_label"] == "high").astype(int)

new_labels.to_csv(OUT / "new_ligand_in_QM_high_low_labels.csv", index=False)

print("New structures:", len(new_X))
print("Original required features:", len(required))
print("Reduced varying features:", len(reduced))
print("")
print("Reduced feature groups:")
print(feature_info[feature_info["is_varying"]]["feature_group"].value_counts())
print("")
print("New high/low labels:")
print(new_labels["new_label"].value_counts())
print("")
print("Saved reduced feature list:")
print(OUT / "reduced_varying_feature_list.txt")
print("Saved new labels:")
print(OUT / "new_ligand_in_QM_high_low_labels.csv")
