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
num_re = re.compile(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?")

# -------------------------
# Read all barrier files
# -------------------------
barrier_files = sorted(BATCH_DIR.glob("ligand_barriers_JOB_*.txt"))
print("Barrier files found:", len(barrier_files))

if not barrier_files:
    raise FileNotFoundError("No ligand_barriers_JOB_*.txt files found")

barrier_rows = []

for barrier_file in barrier_files:
    job_match = re.search(r"JOB_(\d+)", barrier_file.name)
    job_id = job_match.group(1) if job_match else ""

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
                "job_id": job_id,
                "barrier_file": barrier_file.name,
                "barrier_line": line
            })

barriers = pd.DataFrame(barrier_rows).drop_duplicates("ligand_id", keep="last")
print("Unique barrier ligands:", barriers["ligand_id"].nunique())

barriers.to_csv(OUT_DIR / "batch0092_barriers_all_jobs.tsv", sep="\t", index=False)

# -------------------------
# Load model-used features
# -------------------------
model_feature_file = ROOT / "ML_files" / "dist_features_below_10_ang_list_new.npy"
model_features = np.load(model_feature_file, allow_pickle=True).astype(str).tolist()
print("Model-used features:", len(model_features))

# -------------------------
# Read only model-used columns from feature parquet files
# -------------------------
feature_files = sorted(FEATURE_DIR.glob("*.parquet"))
print("Feature parquet files found:", len(feature_files))

if not feature_files:
    raise FileNotFoundError("No feature parquet files found")

# use first readable parquet to identify available columns
sample_cols = None
for fp in feature_files:
    try:
        sample_cols = set(pd.read_parquet(fp).columns)
        break
    except Exception as e:
        print("Could not read sample:", fp.name, e)

if sample_cols is None:
    raise RuntimeError("Could not read any parquet file")

available_features = [f for f in model_features if f in sample_cols]
missing_features = [f for f in model_features if f not in sample_cols]

print("Available model-used features:", len(available_features))
print("Missing model-used features:", len(missing_features))

feature_rows = []

for i, fp in enumerate(feature_files, 1):
    if i % 50 == 0:
        print(f"Reading feature file {i}/{len(feature_files)}")

    m = id_re.search(fp.name)
    if not m:
        print("Could not identify ligand ID from filename:", fp.name)
        continue

    ligand_id = m.group(0)

    try:
        df = pd.read_parquet(fp, columns=available_features)
    except Exception as e:
        print("Failed to read:", fp.name, e)
        continue

    if df.empty:
        continue

    row = df.iloc[0].to_dict()
    row["ligand_id"] = ligand_id
    row["feature_file"] = fp.name
    feature_rows.append(row)

features = pd.DataFrame(feature_rows)
features = features.drop_duplicates("ligand_id", keep="last")

print("Feature rows:", len(features))
print("Unique feature ligands:", features["ligand_id"].nunique())

# -------------------------
# Merge
# -------------------------
matched = barriers.merge(features, on="ligand_id", how="inner")

print("Matched rows:", len(matched))
print("Matched unique ligands:", matched["ligand_id"].nunique())
print("Barrier min:", matched["barrier_energy"].min())
print("Barrier max:", matched["barrier_energy"].max())
print("Barrier mean:", matched["barrier_energy"].mean())

# overwrite the old pilot table with full batch table
features.to_parquet(OUT_DIR / "batch0092_feature_table_all_jobs_model_used_only.parquet", index=False)
matched.to_parquet(OUT_DIR / "batch0092_features_with_barriers.parquet", index=False)
matched.head(20).to_csv(OUT_DIR / "batch0092_features_with_barriers_preview.tsv", sep="\t", index=False)

(OUT_DIR / "model_used_features.txt").write_text("\n".join(model_features) + "\n")
(OUT_DIR / "model_used_features_available_in_batch0092.txt").write_text("\n".join(available_features) + "\n")
(OUT_DIR / "model_used_features_missing_from_batch0092.txt").write_text("\n".join(missing_features) + "\n")

lig_related = [f for f in available_features if "LIG" in f and "PRO" in f]
(OUT_DIR / "model_used_ligand_related_features.txt").write_text("\n".join(lig_related) + "\n")

print("Saved full matched table to:")
print(OUT_DIR / "batch0092_features_with_barriers.parquet")
print("Done.")