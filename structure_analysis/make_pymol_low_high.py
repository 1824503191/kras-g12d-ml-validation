import os

TARGET_DIR = "/data/student/mingsong/structure_analysis"

HIGH_DIR = os.path.join(TARGET_DIR, "high_pdbqt")
LOW_DIR = os.path.join(TARGET_DIR, "low_pdbqt")
RECEPTOR_FILE = os.path.join(TARGET_DIR, "receptor", "g12d_receptor.pdb")

pml = []

pml.append("reinitialize")
pml.append("bg_color white")

# Receptor
if os.path.exists(RECEPTOR_FILE):
    pml.append(f'load "{RECEPTOR_FILE}", receptor')
    pml.append("hide everything, receptor")
    pml.append("show cartoon, receptor")
    pml.append("color grey70, receptor")

# Low barrier ligands: cyan
for i, file in enumerate(sorted(os.listdir(LOW_DIR)), start=1):
    if file.endswith(".pdbqt"):
        path = os.path.join(LOW_DIR, file)
        obj = f"low_{i}"
        pml.append(f'load "{path}", {obj}')
        pml.append(f"show sticks, {obj}")
        pml.append(f"color cyan, {obj}")

# High barrier ligands: orange
for i, file in enumerate(sorted(os.listdir(HIGH_DIR)), start=1):
    if file.endswith(".pdbqt"):
        path = os.path.join(HIGH_DIR, file)
        obj = f"high_{i}"
        pml.append(f'load "{path}", {obj}')
        pml.append(f"show sticks, {obj}")
        pml.append(f"color orange, {obj}")

pml.append("set stick_radius, 0.18")
pml.append("set ray_opaque_background, off")
pml.append("zoom")
pml.append("orient")
pml.append("png low_high_ligands_overview.png, width=1600, height=1200, dpi=300, ray=1")
pml.append("save low_high_ligands.pse")

with open(os.path.join(TARGET_DIR, "pymol_low_high.pml"), "w") as f:
    f.write("\n".join(pml))

print("Saved pymol_low_high.pml")
