import re
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/myriadfs/home/ucapmge/Scratch/omm_pipeline_screening")
BATCH = "batch0092"

BATCH_DIR = ROOT / "work_dir" / BATCH
FEATURE_DIR = BATCH_DIR / "features"
OUT_DIR = ROOT / "analysis_batch0092" / "tables"
OUT_DIR.mkdir(parents=True, exist_ok=True)

id_re = re.compile(r"PV-\d+")
num_re = re.compile(r"[-+]?\d*\.\d+|[-+]?\d+")

# -------------------------
# 1. Read barrier file
# -------------------------
barrier_file = BATCH_DIR / "ligand_barriers_JOB_1.txt"

if not barrier_file.exists():
    summary_files = sorted(BATCH_DIR.glob("ligand_summaries_*.txt"))
    if not summary_files:
        raise FileNotFoundError("No ligand_barriers_JOB_1.txt or ligand_summaries file found.")
    barrier_file = summary_files[-1]

print("Using barrier file:", barrier_file)

barrier_rows = []

with open(barrier_file, "r", errors="ignore") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        id_match = id_re.search(line)
        if not id_match:
            continue

        ligand_id = id_match.group(0)

        line_no_id = id_re.sub(" ", line)
        nums = [float(x) for x in num_re.findall(line_no_id)]
        candidates = [x for x in nums if 10 <= x <= 60]

        if not candidates:
            continue

        barrier = candidates[-1]

        barrier_rows.append({
            "ligand_id": ligand_id,
            "barrier_energy": barrier,
            "barrier_line": line,
        })

barriers = pd.DataFrame(barrier_rows).drop_duplicates("ligand_id", keep="last")

print("Barrier rows:", len(barriers))
print(barriers.head())

barriers.to_csv(OUT_DIR / "batch0092_barriers.tsv", sep="\t", index=False)

# -------------------------
# 2. Read feature parquet files
# -------------------------
parquet_files = sorted(FEATURE_DIR.glob("*.parquet"))

print("Feature parquet files:", len(parquet_files))

feature_rows = []

for p in parquet_files:
    m = id_re.search(p.name)
    ligand_id = m.group(0) if m else p.stem

    df = pd.read_parquet(p)

    if df.empty:
        continue

    row = df.iloc[0].to_dict()
    row["ligand_id"] = ligand_id
    feature_rows.append(row)

features = pd.DataFrame(feature_rows)

print("Feature table shape:", features.shape)
print(features[["ligand_id"]].head())

features.to_parquet(OUT_DIR / "batch0092_feature_table.parquet", index=False)

# -------------------------
# 3. Merge features and barriers
# -------------------------
matched = features.merge(barriers, on="ligand_id", how="inner")

print("Matched shape:", matched.shape)
print("Matched ligands:", matched["ligand_id"].nunique())
print(matched[["ligand_id", "barrier_energy"]].head())

matched.to_parquet(OUT_DIR / "batch0092_features_with_barriers.parquet", index=False)

matched[["ligand_id", "barrier_energy"]].to_csv(
    OUT_DIR / "batch0092_features_with_barriers_preview.tsv",
    sep="\t",
    index=False
)

# -------------------------
# 4. Load RF model-used features
# -------------------------
model_feature_path = ROOT / "ML_files" / "dist_features_below_10_ang_list_new.npy"

if model_feature_path.exists():
    model_features = [str(x) for x in np.load(model_feature_path, allow_pickle=True)]

    available_model_features = [f for f in model_features if f in matched.columns]
    missing_model_features = [f for f in model_features if f not in matched.columns]

    ligand_related = [
        f for f in available_model_features
        if "PRO" in f and "LIG" in f
    ]

    pd.Series(model_features).to_csv(
        OUT_DIR / "model_used_features.txt",
        index=False,
        header=False
    )

    pd.Series(available_model_features).to_csv(
        OUT_DIR / "model_used_features_available_in_batch0092.txt",
        index=False,
        header=False
    )

    pd.Series(missing_model_features).to_csv(
        OUT_DIR / "model_used_features_missing_from_batch0092.txt",
        index=False,
        header=False
    )

    pd.Series(ligand_related).to_csv(
        OUT_DIR / "model_used_ligand_related_features.txt",
        index=False,
        header=False
    )

    print("Model-used features:", len(model_features))
    print("Available model-used features:", len(available_model_features))
    print("Model-used PRO-LIG features:", len(ligand_related))
else:
    print("Model-used feature list not found:", model_feature_path)

print("Done.")
