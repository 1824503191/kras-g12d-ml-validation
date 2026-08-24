import re
from pathlib import Path
import pandas as pd

ROOT = Path("/myriadfs/home/ucapmge/Scratch/omm_pipeline_screening")
TOP_FILE = ROOT / "analysis_batch0092/tables/batch0092_top10_PRO_LIG_features.txt"
OUT = ROOT / "analysis_batch0092/tables/pdb_mapping_candidates.tsv"

features = [x.strip() for x in TOP_FILE.read_text().splitlines() if x.strip()]

target_nums = []
for feat in features:
    m = re.search(r"PRO(\d+)-LIG", feat)
    if m:
        target_nums.append(int(m.group(1)))

target_nums = sorted(set(target_nums))

print("Target PRO numbers:")
print(target_nums)

rows = []

for pdb in ROOT.rglob("*.pdb"):
    found = set()
    matched_residues = {}

    try:
        with open(pdb, "r", errors="ignore") as f:
            for line in f:
                if not line.startswith(("ATOM", "HETATM")):
                    continue

                resname = line[17:20].strip()
                chain = line[21].strip()
                resseq_raw = line[22:26].strip()

                try:
                    resseq = int(resseq_raw)
                except:
                    continue

                if resseq in target_nums:
                    found.add(resseq)
                    matched_residues.setdefault(resseq, set()).add(
                        f"chain {chain}, {resname}{resseq}"
                    )

    except Exception:
        continue

    if found:
        rows.append({
            "pdb": str(pdb),
            "matched_count": len(found),
            "matched_PRO_numbers": ",".join(map(str, sorted(found))),
            "matched_residues": "; ".join(
                f"PRO{k}=" + ",".join(sorted(v))
                for k, v in sorted(matched_residues.items())
            )
        })

df = pd.DataFrame(rows)

if df.empty:
    print("No matching PDB found.")
else:
    df = df.sort_values("matched_count", ascending=False)
    df.to_csv(OUT, sep="\t", index=False)

    print("Top candidates:")
    print(df.head(20).to_string(index=False))
    print("Saved:", OUT)
