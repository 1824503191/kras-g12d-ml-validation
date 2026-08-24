import re
import csv
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INPUT_FILE = "all_barriers.txt"
TOP_N = 20

id_re = re.compile(r"(PV-\d+)")
num_re = re.compile(r"[-+]?\d*\.\d+|[-+]?\d+")

rows = []

with open(INPUT_FILE, "r") as f:
    for i, line in enumerate(f, start=1):
        line = line.strip()
        if not line:
            continue

        id_match = id_re.search(line)
        if not id_match:
            continue

        ligand_id = id_match.group(1)

        # Remove ligand ID first, otherwise digits inside PV-xxx may be read as numbers.
        line_without_id = id_re.sub(" ", line)
        nums = [float(x) for x in num_re.findall(line_without_id)]

        candidates = [x for x in nums if 10 <= x <= 60]
        if not candidates:
            continue

        barrier = candidates[-1]
        rows.append((ligand_id, barrier, line))

if not rows:
    raise RuntimeError("No valid barrier values found. Check all_barriers.txt format.")

values = np.array([r[1] for r in rows])

print("Number of barrier values:", len(values))
print("Min:", values.min())
print("Max:", values.max())
print("Mean:", values.mean())
print("Std:", values.std())

with open("barrier_values_cleaned.tsv", "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["ligand_id", "barrier_energy", "original_line"])
    writer.writerows(rows)

# Remove duplicate ligand IDs if any
dedup = {}
for ligand_id, barrier, original_line in rows:
    if ligand_id not in dedup:
        dedup[ligand_id] = (barrier, original_line)

dedup_rows = [(ligand_id, v[0], v[1]) for ligand_id, v in dedup.items()]
dedup_rows_sorted = sorted(dedup_rows, key=lambda x: x[1])

low = dedup_rows_sorted[:TOP_N]
high = dedup_rows_sorted[-TOP_N:][::-1]

with open("low_barrier_ligands.txt", "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["ligand_id", "barrier_energy", "original_line"])
    writer.writerows(low)

with open("high_barrier_ligands.txt", "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["ligand_id", "barrier_energy", "original_line"])
    writer.writerows(high)

with open("selected_ligands.tsv", "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["ligand_id", "barrier_energy", "group", "original_line"])
    for ligand_id, barrier, original_line in low:
        writer.writerow([ligand_id, barrier, "low", original_line])
    for ligand_id, barrier, original_line in high:
        writer.writerow([ligand_id, barrier, "high", original_line])

# Plot histogram + smoothed density
plt.figure(figsize=(10, 6))
plt.hist(values, bins=60, density=True, alpha=0.65, label="Barrier energies")

x_grid = np.linspace(values.min(), values.max(), 500)
bandwidth = 1.06 * values.std() * (len(values) ** (-1 / 5))
if bandwidth <= 0:
    bandwidth = 0.25

density = np.zeros_like(x_grid)
chunk_size = 1000

for start in range(0, len(values), chunk_size):
    chunk = values[start:start + chunk_size]
    density += np.exp(-0.5 * ((x_grid[None, :] - chunk[:, None]) / bandwidth) ** 2).sum(axis=0)

density = density / (len(values) * bandwidth * np.sqrt(2 * np.pi))

plt.plot(x_grid, density, linewidth=2.5, label="Smoothed density")

plt.xlabel("Barrier energy (kcal/mol)")
plt.ylabel("Density")
plt.title("Barrier energies: bimodal distribution")
plt.legend()
plt.tight_layout()

plt.savefig("barrier_bimodal_distribution.png", dpi=300)
plt.savefig("barrier_bimodal_distribution.pdf")

print("Saved barrier plot and high/low ligand lists.")
