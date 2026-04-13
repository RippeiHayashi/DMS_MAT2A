#!/usr/bin/env python3

from pathlib import Path
import math
import pandas as pd

# -----------------------------
# User settings
# -----------------------------
workdir = Path("/scratch/lf10/rh1772/MAT2A/analysis/RNAseq/DMS_reactivity")
prefix = "MT-CO1"

dotbracket_file = workdir / f"{prefix}.folded.dotbracket.txt"
#reactivity_file = workdir / f"{prefix}.reactivity_25.txt"
reactivity_file = workdir / f"{prefix}.reactivity_05.txt"

# Clip reactivities for display in VARNA
# Values above this are set to this maximum for coloring
clip_max = 0.05

# -----------------------------
# Read dot-bracket file
# -----------------------------
with open(dotbracket_file) as fh:
    lines = [x.strip() for x in fh if x.strip()]

if len(lines) < 3:
    raise ValueError(f"Unexpected dot-bracket file format: {dotbracket_file}")

header = lines[0]
sequence = lines[1]
structure = lines[2]

if len(sequence) != len(structure):
    raise ValueError("Sequence and structure lengths do not match")

L = len(sequence)

# -----------------------------
# Read reactivity file
# -----------------------------
# Expected format: tx_pos \t reactivity
react = pd.read_csv(
    reactivity_file,
    sep="\t",
    header=None,
    names=["tx_pos", "reactivity"]
)

react["tx_pos"] = react["tx_pos"].astype(int)
react["reactivity"] = react["reactivity"].astype(float)

# Fill full-length vector with missing positions as NaN
full = pd.DataFrame({"tx_pos": range(1, L + 1)})
full = full.merge(react, on="tx_pos", how="left")

# Clip for visualization
#full["reactivity_clipped"] = full["reactivity"].clip(lower=0, upper=clip_max)
full["reactivity_clipped"] = full["reactivity"].clip(lower=0)

# -----------------------------
# Prepare VARNA values string
# -----------------------------
# VARNA accepts one value per nucleotide.
# For missing positions, use 0.0 so the nucleotide stays uncolored/low-colored.
varna_values = []
for val in full["reactivity_clipped"]:
    if pd.isna(val):
        varna_values.append("0.0")
    else:
        varna_values.append(f"{val:.6f}")

varna_values_str = ";".join(varna_values)

# -----------------------------
# Optional annotation labels
# -----------------------------
# Label only positions with strong reactivity
label_positions = full.loc[full["reactivity_clipped"] >= 0.02, "tx_pos"].tolist()
annotation_str = ";".join(str(x) for x in label_positions)

# -----------------------------
# Write outputs
# -----------------------------
seq_out = workdir / f"{prefix}.varna.sequence.txt"
struct_out = workdir / f"{prefix}.varna.structure.txt"
values_out = workdir / f"{prefix}.varna.values.txt"
#cmd_out = workdir / f"{prefix}.varna.command.sh"
table_out = workdir / f"{prefix}.varna.input_table.tsv"

with open(seq_out, "w") as fh:
    fh.write(sequence + "\n")

with open(struct_out, "w") as fh:
    fh.write(structure + "\n")

with open(values_out, "w") as fh:
    fh.write(varna_values_str + "\n")

full.to_csv(table_out, sep="\t", index=False)

print("Wrote:")
print(seq_out)
print(struct_out)
print(values_out)
print(table_out)
#print(cmd_out)

if annotation_str:
    print("\nHighly reactive positions (>=0.02):")
    print(annotation_str)
