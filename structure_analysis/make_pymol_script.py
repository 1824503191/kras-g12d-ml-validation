import os

TARGET_DIR = "/data/student/mingsong/structure_analysis"

HIGH_DIR = os.path.join(TARGET_DIR, "high_pdbqt")
LOW_DIR = os.path.join(TARGET_DIR, "low_pdbqt")
RECEPTOR_DIR = os.path.join(TARGET_DIR, "receptor")

pml = []

pml.append("reinitialize")
pml.append("bg_color white")

# Load receptor if available
if os.path.isdir(RECEPTOR_DIR):
    for file in os.listdir(RECEPTOR_DIR):
        if file.endswith(".pdb") or file.endswith(".pdbqt"):
            path = os.path.join(RECEPTOR_DIR, file)
            pml.append(f'load "{path}", receptor')
            pml.append("hide everything, receptor")
            pml.append("show cartoon, receptor")
            pml.append("color grey70, receptor")
            break

# Load low barrier ligands
for i, file in enumerate(os.listdir(LOW_DIR)):
    if file.endswith(".pdbqt"):
        path = os.path.join(LOW_DIR, file)
        obj = f"low_{i+1}"
        pml.append(f'load "{path}", {obj}')
        pml.append(f"show sticks, {obj}")
        pml.append(f"color blue, {obj}")

# Load high barrier ligands
for i, file in enumerate(os.listdir(HIGH_DIR)):
    if file.endswith(".pdbqt"):
        path = os.path.join(HIGH_DIR, file)
        obj = f"high_{i+1}"
        pml.append(f'load "{path}", {obj}')
        pml.append(f"show sticks, {obj}")
        pml.append(f"color red, {obj}")

pml.append("set ray_opaque_background, off")
pml.append("set stick_radius, 0.18")
pml.append("zoom")
pml.append("orient")
pml.append("png low_high_ligands_overview.png, width=1600, height=1200, dpi=300, ray=1")
pml.append("save low_high_ligands.pse")

with open("pymol_low_high.pml", "w") as f:
    f.write("\n".join(pml))

print("Saved pymol_low_high.pml")
