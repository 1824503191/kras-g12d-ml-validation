from pathlib import Path
import numpy as np
import pandas as pd
import joblib

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    roc_auc_score,
    f1_score,
    confusion_matrix
)

WORK = Path("/data/student/mingsong/high_low_transfer_reduced_features")
NEW_DIR = Path("/data/student/mingsong/qmmm_validation_47_complete")
STRUCTURE_ANALYSIS = Path("/data/student/mingsong/structure_analysis")

TABLES = WORK / "tables"
MODELS = WORK / "models"
MODELS.mkdir(exist_ok=True)

OLD_FEATURE_TABLE = WORK / "old_data/batch0092_features_with_barriers.parquet"
NEW_FEATURE_DIR = NEW_DIR / "features"

def read_id_list(path):
    ids = []
    if not path.exists():
        return ids
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s and not s.startswith("#"):
                ids.append(s.split()[0])
    return ids

def get_id_col(df):
    for c in ["ligand_id", "Molecule", "molid", "mol_id"]:
        if c in df.columns:
            return c
    return df.columns[0]

def find_barrier_col(df):
    preferred = [
        "barrier_kcal_per_mol",
        "barrier",
        "predicted_barrier",
        "activation_barrier",
        "RF_predicted_barrier",
        "old_rf_predicted_barrier"
    ]
    for c in preferred:
        if c in df.columns:
            return c
    candidates = [c for c in df.columns if "barrier" in c.lower()]
    if candidates:
        return candidates[0]
    return None

def bootstrap_ci(y_true, y_score, y_pred, metric, n_boot=5000, seed=123):
    rng = np.random.default_rng(seed)
    vals = []
    n = len(y_true)

    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        yt = y_true[idx]
        ys = y_score[idx]
        yp = y_pred[idx]

        if len(np.unique(yt)) < 2:
            continue

        try:
            if metric == "auc":
                vals.append(roc_auc_score(yt, ys))
            elif metric == "balanced_accuracy":
                vals.append(balanced_accuracy_score(yt, yp))
            elif metric == "accuracy":
                vals.append(accuracy_score(yt, yp))
            elif metric == "f1":
                vals.append(f1_score(yt, yp))
        except Exception:
            pass

    vals = np.array(vals)
    if len(vals) == 0:
        return [np.nan, np.nan]
    return np.percentile(vals, [2.5, 97.5])

def permutation_test(y_true, y_score, y_pred, observed_auc, observed_bacc, n_perm=10000, seed=456):
    rng = np.random.default_rng(seed)
    auc_vals = []
    bacc_vals = []

    for _ in range(n_perm):
        shuffled = rng.permutation(y_true)
        try:
            auc_vals.append(roc_auc_score(shuffled, y_score))
            bacc_vals.append(balanced_accuracy_score(shuffled, y_pred))
        except Exception:
            pass

    auc_vals = np.array(auc_vals)
    bacc_vals = np.array(bacc_vals)

    auc_p = (1 + np.sum(auc_vals >= observed_auc)) / (len(auc_vals) + 1)
    bacc_p = (1 + np.sum(bacc_vals >= observed_bacc)) / (len(bacc_vals) + 1)

    return auc_p, bacc_p

with open(TABLES / "reduced_varying_feature_list.txt") as f:
    reduced_features = [line.strip() for line in f if line.strip()]

print("Reduced features from new structures:", len(reduced_features))

if not OLD_FEATURE_TABLE.exists():
    raise FileNotFoundError(f"Old feature table not found: {OLD_FEATURE_TABLE}")

old_df = pd.read_parquet(OLD_FEATURE_TABLE)
old_id_col = get_id_col(old_df)
barrier_col = find_barrier_col(old_df)

print("Old feature table:", OLD_FEATURE_TABLE)
print("Old table shape:", old_df.shape)
print("Old ID column:", old_id_col)
print("Old barrier column:", barrier_col)

old_df = old_df.rename(columns={old_id_col: "ligand_id"})
old_df["ligand_id"] = old_df["ligand_id"].astype(str)

available_features = [f for f in reduced_features if f in old_df.columns]
missing_features = [f for f in reduced_features if f not in old_df.columns]

print("Reduced features available in old table:", len(available_features))
print("Reduced features missing in old table:", len(missing_features))

if len(available_features) < 5:
    raise RuntimeError("Too few reduced features available in old feature table.")

# 优先使用 structure_analysis 里的 old high/low labels
low_ids = read_id_list(STRUCTURE_ANALYSIS / "low_barrier_ligands.txt")
high_ids = read_id_list(STRUCTURE_ANALYSIS / "high_barrier_ligands.txt")

label_source = None

if len(low_ids) > 0 and len(high_ids) > 0:
    old_labels = pd.DataFrame({
        "ligand_id": low_ids + high_ids,
        "old_label": ["low"] * len(low_ids) + ["high"] * len(high_ids)
    })
    old_labels["ligand_id"] = old_labels["ligand_id"].astype(str)
    old_labels["y"] = (old_labels["old_label"] == "high").astype(int)

    train = old_labels.merge(old_df[["ligand_id"] + available_features], on="ligand_id", how="inner")

    if len(train) >= 20:
        label_source = "structure_analysis_high_low_txt"
    else:
        print("Warning: too few rows matched using high/low txt labels:", len(train))
        train = None

else:
    train = None

# 如果 high/low txt match 不上，则用 old table 里的 barrier median split
if train is None:
    if barrier_col is None:
        raise RuntimeError("No old high/low txt labels matched and no barrier column found.")
    
    tmp = old_df[["ligand_id", barrier_col] + available_features].copy()
    tmp = tmp.dropna(subset=[barrier_col])

    tmp_sorted = tmp.sort_values(barrier_col).reset_index(drop=True)
    n = len(tmp_sorted)
    half = n // 2

    old_labels = []
    for i in range(n):
        if i < half:
            old_labels.append("low")
        elif i >= n - half:
            old_labels.append("high")
        else:
            old_labels.append("middle")

    tmp_sorted["old_label"] = old_labels
    train = tmp_sorted[tmp_sorted["old_label"].isin(["low", "high"])].copy()
    train["y"] = (train["old_label"] == "high").astype(int)
    label_source = f"median_split_on_{barrier_col}"

print("Old label source:", label_source)
print("Old train rows:", len(train))
print(train["old_label"].value_counts())

X_train = train[available_features].copy()
y_train = train["y"].to_numpy()

new_labels = pd.read_csv(TABLES / "new_ligand_in_QM_high_low_labels.csv")

new_dfs = []
for p in sorted(NEW_FEATURE_DIR.glob("*_required_features.csv.gz")):
    df = pd.read_csv(p)
    molid = df["Molecule"].iloc[0]
    row = df[available_features].copy()
    row.insert(0, "ligand_id", molid)
    new_dfs.append(row)

new_feat = pd.concat(new_dfs, ignore_index=True)
test = new_labels.merge(new_feat, on="ligand_id", how="inner")

print("New test rows:", len(test))
print(test["new_label"].value_counts())

X_test = test[available_features].copy()
y_test = test["y"].to_numpy()

models = {
    "logistic_regression_balanced": Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(max_iter=5000, class_weight="balanced", solver="liblinear"))
    ]),
    "random_forest_balanced": Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("clf", RandomForestClassifier(
            n_estimators=500,
            random_state=42,
            class_weight="balanced",
            max_features="sqrt",
            min_samples_leaf=2,
            n_jobs=-1
        ))
    ])
}

summary_rows = []

for model_name, model in models.items():
    print("\n==============================")
    print("Model:", model_name)

    model.fit(X_train, y_train)

    y_score = model.predict_proba(X_test)[:, 1]
    y_pred = (y_score >= 0.5).astype(int)

    acc = accuracy_score(y_test, y_pred)
    bacc = balanced_accuracy_score(y_test, y_pred)
    auc = roc_auc_score(y_test, y_score)
    f1 = f1_score(y_test, y_pred)
    cm = confusion_matrix(y_test, y_pred)

    auc_ci = bootstrap_ci(y_test, y_score, y_pred, "auc")
    bacc_ci = bootstrap_ci(y_test, y_score, y_pred, "balanced_accuracy")
    acc_ci = bootstrap_ci(y_test, y_score, y_pred, "accuracy")
    f1_ci = bootstrap_ci(y_test, y_score, y_pred, "f1")

    auc_p, bacc_p = permutation_test(y_test, y_score, y_pred, auc, bacc)

    print("Accuracy:", acc)
    print("Balanced accuracy:", bacc)
    print("ROC AUC:", auc)
    print("F1:", f1)
    print("Confusion matrix:")
    print(cm)
    print("AUC 95% CI:", auc_ci)
    print("Balanced accuracy 95% CI:", bacc_ci)
    print("AUC permutation p:", auc_p)
    print("Balanced accuracy permutation p:", bacc_p)

    summary_rows.append({
        "model": model_name,
        "old_label_source": label_source,
        "n_train": len(train),
        "n_test": len(test),
        "n_features": len(available_features),
        "accuracy": acc,
        "accuracy_ci_low": acc_ci[0],
        "accuracy_ci_high": acc_ci[1],
        "balanced_accuracy": bacc,
        "balanced_accuracy_ci_low": bacc_ci[0],
        "balanced_accuracy_ci_high": bacc_ci[1],
        "roc_auc": auc,
        "roc_auc_ci_low": auc_ci[0],
        "roc_auc_ci_high": auc_ci[1],
        "f1": f1,
        "f1_ci_low": f1_ci[0],
        "f1_ci_high": f1_ci[1],
        "auc_permutation_p": auc_p,
        "balanced_accuracy_permutation_p": bacc_p,
        "tn": cm[0,0],
        "fp": cm[0,1],
        "fn": cm[1,0],
        "tp": cm[1,1],
    })

    pred_out = test[["ligand_id", "computed_barrier_kcal_per_mol", "new_label", "y"]].copy()
    pred_out["score_high_probability"] = y_score
    pred_out["predicted_label"] = np.where(y_pred == 1, "high", "low")
    pred_out.to_csv(TABLES / f"transfer_predictions_{model_name}.csv", index=False)

    joblib.dump(model, MODELS / f"{model_name}.joblib")

    if model_name == "random_forest_balanced":
        importances = model.named_steps["clf"].feature_importances_
    else:
        importances = np.abs(model.named_steps["clf"].coef_[0])

    imp_df = pd.DataFrame({
        "feature": available_features,
        "importance": importances,
    }).sort_values("importance", ascending=False)

    imp_df.to_csv(TABLES / f"transfer_feature_importance_{model_name}.csv", index=False)

summary = pd.DataFrame(summary_rows)
summary.to_csv(TABLES / "transfer_test_summary.csv", index=False)

with open(TABLES / "transfer_test_summary.txt", "w") as f:
    f.write(summary.to_string(index=False))
    f.write("\n")

print("\n==============================")
print("FINAL SUMMARY")
print(summary.to_string(index=False))
print("\nSaved:")
print(TABLES / "transfer_test_summary.csv")
print(TABLES / "transfer_test_summary.txt")
