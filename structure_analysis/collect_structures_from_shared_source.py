import os
import csv
import shutil
import glob

TARGET_DIR = "/data/student/mingsong/structure_analysis"
SELECTED_FILE = os.path.join(TARGET_DIR, "selected_ligands.tsv")

PDBQT_SOURCE = "/data/student/ligand_screening/top_docked_pdbqt"
RECEPTOR_SOURCE = "/data/student/ligand_screening/g12d_receptor.pdb"

HIGH_DIR = os.path.join(TARGET_DIR, "high_pdbqt")
LOW_DIR = os.path.join(TARGET_DIR, "low_pdbqt")
RECEPTOR_DIR = os.path.join(TARGET_DIR, "receptor")

os.makedirs(HIGH_DIR, exist_ok=True)
os.makedirs(LOW_DIR, exist_ok=True)
os.makedirs(RECEPTOR_DIR, exist_ok=True)

# 清空之前可能残留的 high/low 结构文件
for folder in [HIGH_DIR, LOW_DIR]:
    for f in glob.glob(os.path.join(folder, "*")):
        os.remove(f)

# 复制 receptor
if os.path.exists(RECEPTOR_SOURCE):
    shutil.copy2(RECEPTOR_SOURCE, os.path.join(RECEPTOR_DIR, "g12d_receptor.pdb"))
    print("Copied receptor.")
else:
    print("Warning: receptor file not found.")

selected = []

with open(SELECTED_FILE, "r") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        ligand_id = row["ligand_id"].strip()
        group = row["group"].strip().lower()
        selected.append((ligand_id, group))

print("Selected ligands:", len(selected))

missing = []

for ligand_id, group in selected:
    src = os.path.join(PDBQT_SOURCE, ligand_id + ".pdbqt")

    if not os.path.exists(src):
        missing.append(ligand_id)
        continue

    if group == "high":
        dst = os.path.join(HIGH_DIR, ligand_id + ".pdbqt")
    elif group == "low":
        dst = os.path.join(LOW_DIR, ligand_id + ".pdbqt")
    else:
        continue

    shutil.copy2(src, dst)

print("Copied high ligands:", len(os.listdir(HIGH_DIR)))
print("Copied low ligands:", len(os.listdir(LOW_DIR)))

if missing:
    with open(os.path.join(TARGET_DIR, "missing_pdbqt.txt"), "w") as f:
        for ligand_id in missing:
            f.write(ligand_id + "\n")
    print("Missing pdbqt files:", len(missing))
    print("Saved missing_pdbqt.txt")
else:
    print("No missing pdbqt files.")
