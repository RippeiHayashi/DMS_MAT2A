#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# User settings
# -----------------------------
input_tsv = "/scratch/lf10/rh1772/MAT2A/analysis/RNAseq/DMS_reactivity/MT-CO1.mean_reactivity.tsv"
outdir = "/scratch/lf10/rh1772/MAT2A/analysis/RNAseq/DMS_reactivity"
prefix = "MT-CO1"

# Coverage filter:
# keep positions if BOTH untreated and treated mean depths are >= this threshold
min_depth = 20

# Clip high values for plotting so a few outliers do not dominate the color scale
heatmap_vmax = 0.05

# Rolling window for optional smoothed line plot
smooth_window = 15

# -----------------------------
# Read table
# -----------------------------
df = pd.read_csv(input_tsv, sep="\t")

required_cols = [
    "tx_pos", "ref",
    "depth_mean_0", "depth_mean_05", "depth_mean_25",
    "reactivity_05", "reactivity_25"
]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

df = df.copy()
df = df.sort_values("tx_pos")

# -----------------------------
# Coverage filtering
# -----------------------------
df["keep_05"] = (df["depth_mean_0"] >= min_depth) & (df["depth_mean_05"] >= min_depth)
df["keep_25"] = (df["depth_mean_0"] >= min_depth) & (df["depth_mean_25"] >= min_depth)

# Use NaN for low-depth positions so they appear blank in the heatmap
df["reactivity_05_plot"] = df["reactivity_05"].where(df["keep_05"], np.nan)
df["reactivity_25_plot"] = df["reactivity_25"].where(df["keep_25"], np.nan)

# Clip for plotting only
df["reactivity_05_clip"] = df["reactivity_05_plot"].clip(lower=0, upper=heatmap_vmax)
df["reactivity_25_clip"] = df["reactivity_25_plot"].clip(lower=0, upper=heatmap_vmax)

# -----------------------------
# Export compact matrix
# -----------------------------
matrix_cols = [
    "tx_pos", "ref",
    "depth_mean_0", "depth_mean_05", "depth_mean_25",
    "reactivity_05", "reactivity_25",
    "keep_05", "keep_25"
]
matrix_out = Path(outdir) / f"{prefix}.reactivity_matrix.tsv"
df[matrix_cols].to_csv(matrix_out, sep="\t", index=False)

# Export simple 2-column files useful later for structure tools
react05_out = Path(outdir) / f"{prefix}.reactivity_05.txt"
react25_out = Path(outdir) / f"{prefix}.reactivity_25.txt"

df.loc[df["keep_05"], ["tx_pos", "reactivity_05"]].to_csv(
    react05_out, sep="\t", index=False, header=False
)
df.loc[df["keep_25"], ["tx_pos", "reactivity_25"]].to_csv(
    react25_out, sep="\t", index=False, header=False
)

# -----------------------------
# Heat strip plot
# -----------------------------
heat = np.vstack([
    df["reactivity_05_clip"].to_numpy(),
    df["reactivity_25_clip"].to_numpy()
])

fig, ax = plt.subplots(figsize=(16, 2.8))
im = ax.imshow(
    heat,
    aspect="auto",
    interpolation="nearest"
)

ax.set_yticks([0, 1])
ax.set_yticklabels(["0.5% - 0%", "2.5% - 0%"])
ax.set_xlabel("MT-CO1 transcript position")
ax.set_title(f"{prefix} DMS reactivity heat strip")

# x ticks every 100 nt
max_pos = int(df["tx_pos"].max())
xticks = np.arange(0, max_pos, 100)
ax.set_xticks(xticks - 1)
ax.set_xticklabels([str(x) for x in xticks], rotation=0)

cbar = plt.colorbar(im, ax=ax, pad=0.02)
cbar.set_label(f"Reactivity (clipped at {heatmap_vmax})")

plt.tight_layout()
heat_png = Path(outdir) / f"{prefix}.reactivity_heatstrip.png"
heat_pdf = Path(outdir) / f"{prefix}.reactivity_heatstrip.pdf"
plt.savefig(heat_png, dpi=300)
plt.savefig(heat_pdf)
plt.close()

# -----------------------------
# Line plot
# -----------------------------
fig, ax = plt.subplots(figsize=(16, 4.5))

ax.plot(
    df["tx_pos"],
    df["reactivity_05_plot"],
    linewidth=0.8,
    alpha=0.5,
    label="0.5% - 0% raw"
)
ax.plot(
    df["tx_pos"],
    df["reactivity_25_plot"],
    linewidth=0.8,
    alpha=0.5,
    label="2.5% - 0% raw"
)

# smoothed lines
df["reactivity_05_smooth"] = df["reactivity_05_plot"].rolling(
    window=smooth_window, center=True, min_periods=1
).mean()
df["reactivity_25_smooth"] = df["reactivity_25_plot"].rolling(
    window=smooth_window, center=True, min_periods=1
).mean()

ax.plot(
    df["tx_pos"],
    df["reactivity_05_smooth"],
    linewidth=2.0,
    label=f"0.5% - 0% smoothed ({smooth_window} nt)"
)
ax.plot(
    df["tx_pos"],
    df["reactivity_25_smooth"],
    linewidth=2.0,
    label=f"2.5% - 0% smoothed ({smooth_window} nt)"
)

ax.set_xlabel("MT-CO1 transcript position")
ax.set_ylabel("Reactivity")
ax.set_title(f"{prefix} DMS reactivity profile")
ax.set_xlim(1, max_pos)
ax.set_ylim(bottom=0)
ax.grid(True, axis="y", alpha=0.3)
ax.legend(frameon=False)

plt.tight_layout()
line_png = Path(outdir) / f"{prefix}.reactivity_profile.png"
line_pdf = Path(outdir) / f"{prefix}.reactivity_profile.pdf"
plt.savefig(line_png, dpi=300)
plt.savefig(line_pdf)
plt.close()

print("Wrote:")
print(matrix_out)
print(react05_out)
print(react25_out)
print(heat_png)
print(heat_pdf)
print(line_png)
print(line_pdf)
