from pathlib import Path
import numpy as np
import pandas as pd

from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

ROOT = Path("/myriadfs/home/ucapmge/Scratch/omm_pipeline_screening")
TABLE_DIR = ROOT / "analysis_batch0092" / "tables"
OUT_DIR = TABLE_DIR

matched_path = TABLE_DIR / "batch0092_features_with_barriers.parquet"
category_path = TABLE_DIR / "batch0092_model_features_with_categories.tsv"

matched = pd.read_parquet(matched_path)
categories = pd.read_csv(category_path, sep="\t")

# Define low / middle / high groups
def assign_group(x):
    if x < 28:
        return "low"
    elif x >= 29:
        return "high"
    else:
        return "middle"

matched["barrier_group"] = matched["barrier_energy"].apply(assign_group)

# Only compare low vs high
df = matched[matched["barrier_group"].isin(["low", "high"])].copy()

# y: low = 0, high = 1
y = (df["barrier_group"] == "high").astype(int)

print("Low/high sample size:")
print(df["barrier_group"].value_counts())

feature_sets = {}

# 1. Ligand-residue features
lig_res = categories.loc[
    categories["category"] == "ligand_residue_PRO_LIG", "feature"
].tolist()

# 2. Residue-residue features
res_res = categories.loc[
    categories["category"] == "residue_residue_PRO_PRO", "feature"
].tolist()

# 3. GTP-residue features
gtp_res = categories.loc[
    categories["category"] == "GTP_residue_PRO_GTP", "feature"
].tolist()

# 4. All model-used features available
all_model = categories["feature"].tolist()

# Keep only columns actually present in matched table
feature_sets["ligand_residue_PRO_LIG"] = [f for f in lig_res if f in df.columns]
feature_sets["residue_residue_PRO_PRO"] = [f for f in res_res if f in df.columns]
feature_sets["GTP_residue_PRO_GTP"] = [f for f in gtp_res if f in df.columns]
feature_sets["all_model_used_features"] = [f for f in all_model if f in df.columns]

print("\nFeature set sizes:")
for name, feats in feature_sets.items():
    print(name, len(feats))

models = {
    "logistic_regression_balanced": Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(
            max_iter=5000,
            class_weight="balanced",
            solver="liblinear"
        ))
    ]),
    "random_forest_balanced": Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("clf", RandomForestClassifier(
            n_estimators=300,
            random_state=42,
            class_weight="balanced",
            max_depth=None
        ))
    ])
}

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

scoring = {
    "accuracy": "accuracy",
    "balanced_accuracy": "balanced_accuracy",
    "f1": "f1",
    "roc_auc": "roc_auc"
}

result_rows = []

for fs_name, feats in feature_sets.items():
    if len(feats) == 0:
        print(f"Skipping {fs_name}: no features")
        continue

    X = df[feats].apply(pd.to_numeric, errors="coerce")

    for model_name, model in models.items():
        print(f"\nRunning {model_name} with {fs_name} ({len(feats)} features)")

        scores = cross_validate(
            model,
            X,
            y,
            cv=cv,
            scoring=scoring,
            return_train_score=False
        )

        row = {
            "feature_set": fs_name,
            "n_features": len(feats),
            "model": model_name,
            "n_samples": len(df),
            "n_low": int((df["barrier_group"] == "low").sum()),
            "n_high": int((df["barrier_group"] == "high").sum()),
        }

        for metric in scoring:
            values = scores[f"test_{metric}"]
            row[f"{metric}_mean"] = values.mean()
            row[f"{metric}_std"] = values.std()

        result_rows.append(row)

results = pd.DataFrame(result_rows)
results.to_csv(OUT_DIR / "batch0092_classifier_results_by_feature_set.tsv", sep="\t", index=False)

print("\nClassifier results:")
print(results.to_string(index=False))

print("\nSaved:")
print(OUT_DIR / "batch0092_classifier_results_by_feature_set.tsv")
