from pathlib import Path
import csv
import MDAnalysis as mda

WORK = Path("/data/student/mingsong/qmmm_validation_47_complete")
COMPLETE = Path("/home/edina/shared/qmmm_validation/complete_structures")

barriers_path = WORK / "barriers.csv"

rows = list(csv.DictReader(open(barriers_path)))
print("Columns:", list(rows[0].keys()))
print("Rows in barriers.csv:", len(rows))

id_col = "ligand_id" if "ligand_id" in rows[0] else list(rows[0].keys())[0]
print("Using ID column:", id_col)

ids = [r[id_col] for r in rows]

missing = []
for molid in ids:
    pdb = COMPLETE / f"{molid}.pdb"
    if not pdb.exists():
        missing.append(molid)

print("Missing complete structures:", len(missing))
for x in missing[:20]:
    print(x)

if missing:
    raise SystemExit("Some complete structures are missing. Stop here.")

# Inspect first structure
first = ids[0]
pdb = COMPLETE / f"{first}.pdb"
print("\nFirst structure:", pdb)

u = mda.Universe(str(pdb))
resnames = sorted(set(u.atoms.resnames))

print("Number of atoms:", u.atoms.n_atoms)
print("First 80 resnames:", resnames[:80])

for rn in ["GTP", "GDP", "LIG", "MG", "MG2", "MG1", "WAT", "HOH"]:
    sel = u.select_atoms(f"resname {rn}")
    print(f"{rn} atoms:", sel.n_atoms)

with open(WORK / "qmmm_validation_47_ids.txt", "w") as f:
    for molid in ids:
        f.write(molid + "\n")

print("\nSaved qmmm_validation_47_ids.txt")
