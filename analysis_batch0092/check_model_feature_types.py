from pathlib import Path
import pandas as pd

ROOT = Path("/myriadfs/home/ucapmge/Scratch/omm_pipeline_screening")
TABLE_DIR = ROOT / "analysis_batch0092" / "tables"

feature_file = TABLE_DIR / "model_used_features_available_in_batch0092.txt"

features = [x.strip() for x in feature_file.read_text().splitlines() if x.strip()]

rows = []

for f in features:
    parts = f.split("-")

    has_lig = "LIG" in parts
    has_gtp = "GTP" in parts
    pro_parts = [p for p in parts if p.startswith("PRO")]
    n_pro = len(pro_parts)

    if has_lig and n_pro >= 1:
        category = "ligand_residue_PRO_LIG"
    elif n_pro >= 2 and not has_lig and not has_gtp:
        category = "residue_residue_PRO_PRO"
    elif has_gtp and n_pro >= 1 and not has_lig:
        category = "GTP_residue_PRO_GTP"
    elif has_lig and has_gtp:
        category = "ligand_GTP"
    elif n_pro >= 1:
        category = "protein_related_other"
    else:
        category = "other"

    rows.append({
        "feature": f,
        "category": category
    })

df = pd.DataFrame(rows)

summary = (
    df["category"]
    .value_counts()
    .reset_index()
)

summary.columns = ["category", "count"]

summary.to_csv(TABLE_DIR / "batch0092_model_feature_type_summary.tsv", sep="\t", index=False)
df.to_csv(TABLE_DIR / "batch0092_model_features_with_categories.tsv", sep="\t", index=False)

print("Feature type summary:")
print(summary.to_string(index=False))

print("\nSaved:")
print(TABLE_DIR / "batch0092_model_feature_type_summary.tsv")
print(TABLE_DIR / "batch0092_model_features_with_categories.tsv")
