import os
import csv
import shutil

SOURCE_DIR = "/data/student/mingsong/work_dir"
TARGET_DIR = "/data/student/mingsong/structure_analysis"

SELECTED_FILE = os.path.join(TARGET_DIR, "selected_ligands.tsv")

HIGH_DIR = os.path.join(TARGET_DIR, "high_pdbqt")
LOW_DIR = os.path.join(TARGET_DIR, "low_pdbqt")
RECEPTOR_DIR = os.path.join(TARGET_DIR, "receptor")

os.makedirs(HIGH_DIR, exist_ok=True)
os.makedirs(LOW_DIR, exist_ok=True)
os.makedirs(RECEPTOR_DIR, exist_ok=True)

selected = []

with open(SELECTED_FILE, "r") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        selected.append((row["ligand_id"], row["group"]))

print("Selected ligands:", len(selected))

# Build pdbqt file index
pdbqt_files = []

for root, dirs, files in os.walk(SOURCE_DIR):
    # Skip very large intermediate folders if present
    dirs[:] = [d for d in dirs if d not in {"features"}]

    for file in files:
        if file.endswith(".pdbqt"):
            pdbqt_files.append(os.path.join(root, file))

print("Total pdbqt files found:", len(pdbqt_files))

missing = []

for ligand_id, group in selected:
    matches = [p for p in pdbqt_files if ligand_id in os.path.basename(p)]

    if not matches:
        missing.append(ligand_id)
        continue

    src = matches[0]

    if group == "high":
        dst_dir = HIGH_DIR
    else:
        dst_dir = LOW_DIR

    dst = os.path.join(dst_dir, os.path.basename(src))
    shutil.copy2(src, dst)

print("Copied high pdbqt:", len(os.listdir(HIGH_DIR)))
print("Copied low pdbqt:", len(os.listdir(LOW_DIR)))

if missing:
    with open("missing_pdbqt.txt", "w") as f:
        for ligand_id in missing:
            f.write(ligand_id + "\n")
    print("Missing pdbqt files:", len(missing))
    print("Saved missing_pdbqt.txt")
else:
    print("No missing pdbqt files.")

# Try to copy receptor files if names contain receptor
receptor_candidates = []

for root, dirs, files in os.walk(SOURCE_DIR):
    for file in files:
        lower = file.lower()
        if "receptor" in lower and (lower.endswith(".pdb") or lower.endswith(".pdbqt")):
            receptor_candidates.append(os.path.join(root, file))

for src in receptor_candidates[:5]:
    shutil.copy2(src, os.path.join(RECEPTOR_DIR, os.path.basename(src)))

print("Receptor candidates copied:", min(len(receptor_candidates), 5))
