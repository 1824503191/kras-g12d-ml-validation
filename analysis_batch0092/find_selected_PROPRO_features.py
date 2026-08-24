from pathlib import Path
import numpy as np
import pandas as pd

from sklearn.model_selection import StratifiedKFold
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import SelectKBest, f_classif

ROOT = Path("/myriadfs/home/ucapmge/Scratch/omm_pipeline_screening")
TABLE_DIR = ROOT / "analysis_batch0092" / "tables"

LOW_CUTOFF = 28.0
HIGH_CUTOFF = 29.0

df = pd.read_parquet(TABLE_DIR / "batch0092_features_with_barriers.parquet")

df = df[(df["barrier_energy"] < LOW_CUTOFF) | (df["barrier_energy"] >= HIGH_CUTOFF)].copy()
df["class"] = np.where(df["barrier_energy"] >= HIGH_CUTOFF, 1, 0)
y = df["class"].values

feature_file = TABLE_DIR / "model_used_features_available_in_batch0092.txt"
features = [x.strip() for x in feature_file.read_text().splitlines() if x.strip()]
features = [f for f in features if f in df.columns]

def is_pro_pro(f):
    parts = f.split("-")
    n_pro = sum(p.startswith("PRO") for p in parts)
    return n_pro >= 2 and "LIG" not in parts and "GTP" not in parts

pro_pro = [f for f in features if is_pro_pro(f)]

print("Low/high samples:", len(df))
print("Low:", int((y == 0).sum()))
print("High:", int((y == 1).sum()))
print("PRO-PRO features:", len(pro_pro))

X = df[pro_pro].apply(pd.to_numeric, errors="coerce")

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

records = []

for k in [20, 50, 100, 200]:
    print(f"Running SelectKBest for k={k}")

    selected_counts = {f: 0 for f in pro_pro}
    score_sums = {f: 0.0 for f in pro_pro}

    for fold, (train_idx, test_idx) in enumerate(cv.split(X, y), 1):
        X_train = X.iloc[train_idx]
        y_train = y[train_idx]

        imputer = SimpleImputer(strategy="median")
        X_train_imp = imputer.fit_transform(X_train)

        selector = SelectKBest(score_func=f_classif, k=min(k, len(pro_pro)))
        selector.fit(X_train_imp, y_train)

        scores = selector.scores_
        selected_mask = selector.get_support()

        for f, selected, score in zip(pro_pro, selected_mask, scores):
            if selected:
                selected_counts[f] += 1
            if np.isfinite(score):
                score_sums[f] += float(score)

    for f in pro_pro:
        records.append({
            "k": k,
            "feature": f,
            "selected_count_out_of_5": selected_counts[f],
            "selection_frequency": selected_counts[f] / 5,
            "mean_f_score_across_folds": score_sums[f] / 5,
        })

res = pd.DataFrame(records)

# Add residue numbers
def parse_pro_pair(f):
    parts = f.split("-")
    pros = [p for p in parts if p.startswith("PRO")]
    if len(pros) >= 2:
        r1 = pros[0].replace("PRO", "")
        r2 = pros[1].replace("PRO", "")
        return r1, r2
    return "", ""

pairs = res["feature"].apply(parse_pro_pair)
res["residue_1_number"] = [p[0] for p in pairs]
res["residue_2_number"] = [p[1] for p in pairs]

# Try to map residue names from PDB
pdb_file = ROOT / "manual_charmm_setup" / "step1_pdbreader.pdb"

standard_aa = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
}

residue_map = {}

if pdb_file.exists():
    with open(pdb_file, "r", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            resname = line[17:20].strip()
            chain = line[21].strip()
            resseq = line[22:26].strip()

            if resname not in standard_aa:
                continue

            if resseq not in residue_map:
                residue_map[resseq] = f"{chain}:{resname}{resseq}"

res["residue_1_name"] = res["residue_1_number"].map(residue_map).fillna("")
res["residue_2_name"] = res["residue_2_number"].map(residue_map).fillna("")

res = res.sort_values(
    ["k", "selected_count_out_of_5", "mean_f_score_across_folds"],
    ascending=[True, False, False]
)

out = TABLE_DIR / "batch0092_PROPRO_selected_feature_frequency.tsv"
res.to_csv(out, sep="\t", index=False)

print("Saved:", out)

for k in [20, 50, 100, 200]:
    print("\nTop selected PRO-PRO features for k =", k)
    sub = res[res["k"] == k].head(30)
    print(sub[[
        "feature",
        "selected_count_out_of_5",
        "selection_frequency",
        "mean_f_score_across_folds",
        "residue_1_name",
        "residue_2_name"
    ]].to_string(index=False))

print("Done.")
