from pathlib import Path
import warnings

import numpy as np
import pandas as pd

from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    roc_auc_score,
    confusion_matrix,
)

warnings.filterwarnings("ignore")

ROOT = Path("/myriadfs/home/ucapmge/Scratch/omm_pipeline_screening")
TABLE_DIR = ROOT / "analysis_batch0092" / "tables"
OUT_FILE = TABLE_DIR / "batch0092_classifier_results_full_with_baseline.tsv"

LOW_CUTOFF = 28.0
HIGH_CUTOFF = 29.0

print("=== Full batch0092 classifier analysis ===", flush=True)

# --------------------------------------------------
# Load full matched table
# --------------------------------------------------
matched_file = TABLE_DIR / "batch0092_features_with_barriers.parquet"

if not matched_file.exists():
    raise FileNotFoundError(f"Missing matched table: {matched_file}")

df = pd.read_parquet(matched_file)

print("Original matched rows:", len(df), flush=True)
print("Original matched ligands:", df["ligand_id"].nunique(), flush=True)

# --------------------------------------------------
# Exclude middle band
# --------------------------------------------------
df = df[
    (df["barrier_energy"] < LOW_CUTOFF)
    | (df["barrier_energy"] >= HIGH_CUTOFF)
].copy()

df["class"] = np.where(df["barrier_energy"] >= HIGH_CUTOFF, 1, 0)

y = df["class"].values

n_samples = len(y)
n_low = int((y == 0).sum())
n_high = int((y == 1).sum())

print("After excluding middle band:", flush=True)
print("n_samples:", n_samples, flush=True)
print("n_low:", n_low, flush=True)
print("n_high:", n_high, flush=True)

if n_low < 5 or n_high < 5:
    raise RuntimeError("Too few low or high samples for 5-fold CV.")

# --------------------------------------------------
# Load RF model-used features
# --------------------------------------------------
feature_list_file = TABLE_DIR / "model_used_features_available_in_batch0092.txt"

if not feature_list_file.exists():
    raise FileNotFoundError(f"Missing feature list: {feature_list_file}")

features = [
    x.strip()
    for x in feature_list_file.read_text().splitlines()
    if x.strip()
]

features = [f for f in features if f in df.columns]

print("Available RF model-used features in table:", len(features), flush=True)

# --------------------------------------------------
# Categorise features
# --------------------------------------------------
def feature_category(f):
    parts = f.split("-")
    has_lig = "LIG" in parts
    has_gtp = "GTP" in parts
    n_pro = sum(p.startswith("PRO") for p in parts)

    if has_lig and n_pro >= 1:
        return "PRO-LIG"
    elif has_gtp and n_pro >= 1 and not has_lig:
        return "PRO-GTP"
    elif n_pro >= 2 and not has_lig and not has_gtp:
        return "PRO-PRO"
    else:
        return "other"

pro_lig = [f for f in features if feature_category(f) == "PRO-LIG"]
pro_gtp = [f for f in features if feature_category(f) == "PRO-GTP"]
pro_pro = [f for f in features if feature_category(f) == "PRO-PRO"]

print("PRO-LIG features:", len(pro_lig), flush=True)
print("PRO-GTP features:", len(pro_gtp), flush=True)
print("PRO-PRO features:", len(pro_pro), flush=True)
print("All model-used features:", len(features), flush=True)

# --------------------------------------------------
# Define feature sets
# Feature selection is inside CV pipeline to avoid leakage
# --------------------------------------------------
feature_sets = [
    ("PRO-LIG_all", pro_lig, None),
    ("PRO-GTP_all", pro_gtp, None),

    # Keep residue-residue features, but reduce dimensionality
    ("PRO-PRO_select20", pro_pro, 20),
    ("PRO-PRO_select50", pro_pro, 50),
    ("PRO-PRO_select100", pro_pro, 100),
    ("PRO-PRO_select200", pro_pro, 200),

    # All RF model-used features with selection
    ("ALL_select50", features, 50),
    ("ALL_select100", features, 100),
    ("ALL_select200", features, 200),
]

# --------------------------------------------------
# Cross-validation setup
# --------------------------------------------------
n_splits = 5
cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)

rows = []

# --------------------------------------------------
# Majority baseline
# --------------------------------------------------
majority_class = 1 if n_high >= n_low else 0
y_base = np.full_like(y, majority_class)

tn, fp, fn, tp = confusion_matrix(y, y_base, labels=[0, 1]).ravel()

rows.append({
    "feature_set": "majority_baseline",
    "n_features_total": 0,
    "n_features_used": 0,
    "model": "always_predict_majority",
    "n_samples": n_samples,
    "n_low": n_low,
    "n_high": n_high,
    "accuracy_mean": accuracy_score(y, y_base),
    "accuracy_std": 0.0,
    "balanced_accuracy_mean": balanced_accuracy_score(y, y_base),
    "balanced_accuracy_std": 0.0,
    "f1_mean": f1_score(y, y_base, zero_division=0),
    "f1_std": 0.0,
    "roc_auc_mean": 0.5,
    "roc_auc_std": 0.0,
    "mean_predicted_low": int((y_base == 0).sum()),
    "mean_predicted_high": int((y_base == 1).sum()),
    "mean_tn": tn,
    "mean_fp": fp,
    "mean_fn": fn,
    "mean_tp": tp,
})

# --------------------------------------------------
# Helper: create model pipeline
# --------------------------------------------------
def make_pipeline(model_name, k, n_total_features):
    steps = []
    steps.append(("imputer", SimpleImputer(strategy="median")))

    if k is not None:
        steps.append(("select", SelectKBest(score_func=f_classif, k=min(k, n_total_features))))

    if model_name == "logistic_regression_balanced":
        steps.append(("scaler", StandardScaler()))
        steps.append((
            "clf",
            LogisticRegression(
                class_weight="balanced",
                max_iter=5000,
                solver="liblinear",
                random_state=42,
            )
        ))

    elif model_name == "random_forest_balanced":
        steps.append((
            "clf",
            RandomForestClassifier(
                n_estimators=500,
                class_weight="balanced",
                min_samples_leaf=3,
                random_state=42,
                n_jobs=-1,
            )
        ))

    else:
        raise ValueError(model_name)

    return Pipeline(steps)

# --------------------------------------------------
# Run models
# --------------------------------------------------
for fs_name, fs_features, k in feature_sets:
    if len(fs_features) == 0:
        print("Skipping empty feature set:", fs_name, flush=True)
        continue

    print(f"\nRunning feature set: {fs_name} ({len(fs_features)} total features, k={k})", flush=True)

    X = df[fs_features].apply(pd.to_numeric, errors="coerce")
    n_used = len(fs_features) if k is None else min(k, len(fs_features))

    for model_name in ["logistic_regression_balanced", "random_forest_balanced"]:
        print(f"  Model: {model_name}", flush=True)

        accs = []
        baccs = []
        f1s = []
        aucs = []

        pred_lows = []
        pred_highs = []
        tns = []
        fps = []
        fns = []
        tps = []

        for fold, (train_idx, test_idx) in enumerate(cv.split(X, y), 1):
            X_train = X.iloc[train_idx]
            X_test = X.iloc[test_idx]
            y_train = y[train_idx]
            y_test = y[test_idx]

            pipe = make_pipeline(model_name, k, len(fs_features))
            pipe.fit(X_train, y_train)

            y_pred = pipe.predict(X_test)

            if hasattr(pipe.named_steps["clf"], "predict_proba"):
                y_score = pipe.predict_proba(X_test)[:, 1]
            else:
                y_score = pipe.decision_function(X_test)

            tn, fp, fn, tp = confusion_matrix(y_test, y_pred, labels=[0, 1]).ravel()

            accs.append(accuracy_score(y_test, y_pred))
            baccs.append(balanced_accuracy_score(y_test, y_pred))
            f1s.append(f1_score(y_test, y_pred, zero_division=0))
            aucs.append(roc_auc_score(y_test, y_score))

            pred_lows.append(int((y_pred == 0).sum()))
            pred_highs.append(int((y_pred == 1).sum()))
            tns.append(tn)
            fps.append(fp)
            fns.append(fn)
            tps.append(tp)

        rows.append({
            "feature_set": fs_name,
            "n_features_total": len(fs_features),
            "n_features_used": n_used,
            "model": model_name,
            "n_samples": n_samples,
            "n_low": n_low,
            "n_high": n_high,
            "accuracy_mean": np.mean(accs),
            "accuracy_std": np.std(accs),
            "balanced_accuracy_mean": np.mean(baccs),
            "balanced_accuracy_std": np.std(baccs),
            "f1_mean": np.mean(f1s),
            "f1_std": np.std(f1s),
            "roc_auc_mean": np.mean(aucs),
            "roc_auc_std": np.std(aucs),
            "mean_predicted_low": np.mean(pred_lows),
            "mean_predicted_high": np.mean(pred_highs),
            "mean_tn": np.mean(tns),
            "mean_fp": np.mean(fps),
            "mean_fn": np.mean(fns),
            "mean_tp": np.mean(tps),
        })

# --------------------------------------------------
# Save results
# --------------------------------------------------
results = pd.DataFrame(rows)

results.to_csv(OUT_FILE, sep="\t", index=False)

print("\n=== Classifier results ===", flush=True)
display_cols = [
    "feature_set",
    "n_features_used",
    "model",
    "n_samples",
    "n_low",
    "n_high",
    "accuracy_mean",
    "balanced_accuracy_mean",
    "roc_auc_mean",
    "mean_predicted_low",
    "mean_predicted_high",
]
print(results[display_cols].to_string(index=False), flush=True)

print("\nSaved:", OUT_FILE, flush=True)
print("Done.", flush=True)
