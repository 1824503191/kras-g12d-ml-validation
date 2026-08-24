import re
import sys
from pathlib import Path
import pandas as pd

if len(sys.argv) != 2:
    print("Usage: python map_batch0092_features_to_residues.py /path/to/pdb")
    sys.exit(1)

pdb_path = Path(sys.argv[1])
if not pdb_path.exists():
    raise FileNotFoundError(pdb_path)

ROOT = Path("/myriadfs/home/ucapmge/Scratch/omm_pipeline_screening")
SUMMARY = ROOT / "analysis_batch0092/tables/batch0092_model_used_ligand_features_summary.tsv"
OUT = ROOT / "analysis_batch0092/tables/batch0092_feature_residue_mapping.tsv"

AA3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL"
}

summary = pd.read_csv(SUMMARY, sep="\t")

feature_rows = []

for _, row in summary.iterrows():
    feat = row["feature"]
    m = re.search(r"PRO(\d+)-LIG", feat)

    if not m:
        continue

    pro_num = int(m.group(1))

    feature_rows.append({
        "feature": feat,
        "PRO_number": pro_num,
        "low_mean": row["low_mean"],
        "high_mean": row["high_mean"],
        "high_minus_low_mean": row["high_minus_low_mean"],
        "cohens_d": row["cohens_d"],
        "abs_cohens_d": row["abs_cohens_d"],
        "interpretation": row["interpretation"],
    })

feature_df = pd.DataFrame(feature_rows)

residue_rows = []

with open(pdb_path, "r", errors="ignore") as f:
    for line in f:
        if not line.startswith("ATOM"):
            continue

        resname = line[17:20].strip()
        chain = line[21].strip()
        resseq_raw = line[22:26].strip()

        if resname not in AA3:
            continue

        try:
            resseq = int(resseq_raw)
        except ValueError:
            continue

        residue_rows.append({
            "PRO_number": resseq,
            "chain": chain,
            "residue_number": resseq,
            "residue_name": resname,
        })

residue_df = pd.DataFrame(residue_rows).drop_duplicates()

mapped = feature_df.merge(residue_df, on="PRO_number", how="left")

mapped["mapping_status"] = mapped["residue_name"].apply(
    lambda x: "matched" if pd.notna(x) else "not_found"
)

mapped.to_csv(OUT, sep="\t", index=False)

print("PDB used:", pdb_path)
print("Saved:", OUT)
print(mapped.head(30).to_string(index=False))
