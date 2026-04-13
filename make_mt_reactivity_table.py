#!/usr/bin/env python3

import gzip
from pathlib import Path
import pandas as pd

analysis = "/scratch/lf10/rh1772/MAT2A/analysis/RNAseq"
outroot = Path(f"{analysis}/DMS_mutrates")
outdir = Path(f"{analysis}/DMS_reactivity")
outdir.mkdir(parents=True, exist_ok=True)

samples = {
    "0": [
        "HepG2_DMS_0_pA_rep1",
        "HepG2_DMS_0_pA_rep2",
    ],
    "05": [
        "HepG2_DMS_05_pA_rep1",
        "HepG2_DMS_05_pA_rep2",
    ],
    "25": [
        "HepG2_DMS_25_pA_rep1",
        "HepG2_DMS_25_pA_rep2",
    ],
}

# Optional gene region on MT.
# Set to None if you want all MT positions.
gene_name = "MT-CO1"
gene_region = ("MT", 5904, 7445, "+")   # Ensembl MT-CO1 locus on GRCh38
# Examples:
# gene_name = "MT-ND1";  gene_region = ("MT", 3307, 4262, "+")
# gene_name = "MT-ND4L"; gene_region = ("MT", 10470, 10763, "+")
# gene_name = None; gene_region = None

mut_cols = [
    "A>C", "A>G", "A>T",
    "C>A", "C>G", "C>T",
    "G>A", "G>C", "G>T",
    "T>A", "T>C", "T>G",
]

base_cols = ["chrom", "start", "end", "ref", "depth", "mismatches", "mut_rate"]
all_cols = base_cols + mut_cols

def read_one(sample_name: str) -> pd.DataFrame:
    path = outroot / f"{sample_name}.AorC.primary.q20.Q30.muttype.tsv.gz"
    # despite the filename, your current script includes A/C/G/T
    df = pd.read_csv(path, sep="\t", compression="gzip")
    df = df[df["chrom"] == "MT"].copy()

    if gene_region is not None:
        chrom, gstart1, gend1, strand = gene_region
        # muttype files use BED-like coords: start = pos1-1, end = pos1
        # keep positions whose 1-based coordinate lies within gene interval
        df_pos1 = df["end"]
        df = df[(df["chrom"] == chrom) & (df_pos1 >= gstart1) & (df_pos1 <= gend1)].copy()

        # add transcript coordinate if gene is on +
        if strand == "+":
            df["tx_pos"] = df["end"] - gstart1 + 1
        else:
            df["tx_pos"] = gend1 - df["end"] + 1
    else:
        df["tx_pos"] = pd.NA

    keep_cols = ["chrom", "start", "end", "ref", "tx_pos", "depth", "mismatches", "mut_rate"] + mut_cols
    return df[keep_cols]

def mean_by_condition(sample_list, condition_name):
    dfs = []
    for s in sample_list:
        df = read_one(s).copy()
        rename_map = {
            "depth": f"depth_{s}",
            "mismatches": f"mismatches_{s}",
            "mut_rate": f"mut_rate_{s}",
        }
        for m in mut_cols:
            rename_map[m] = f"{m}_{s}"
        df = df.rename(columns=rename_map)
        dfs.append(df)

    # merge replicates by genomic position
    merged = dfs[0]
    for df in dfs[1:]:
        merged = merged.merge(
            df,
            on=["chrom", "start", "end", "ref", "tx_pos"],
            how="outer"
        )

    # mean across replicates
    rep_names = sample_list
    merged[f"depth_mean_{condition_name}"] = merged[[f"depth_{r}" for r in rep_names]].mean(axis=1)
    merged[f"mismatches_mean_{condition_name}"] = merged[[f"mismatches_{r}" for r in rep_names]].mean(axis=1)
    merged[f"mut_rate_mean_{condition_name}"] = merged[[f"mut_rate_{r}" for r in rep_names]].mean(axis=1)

    for m in mut_cols:
        merged[f"{m}_mean_{condition_name}"] = merged[[f"{m}_{r}" for r in rep_names]].mean(axis=1)

    keep = ["chrom", "start", "end", "ref", "tx_pos",
            f"depth_mean_{condition_name}",
            f"mismatches_mean_{condition_name}",
            f"mut_rate_mean_{condition_name}"] + [f"{m}_mean_{condition_name}" for m in mut_cols]
    return merged[keep]

# Build mean tables
mean0 = mean_by_condition(samples["0"], "0")
mean05 = mean_by_condition(samples["05"], "05")
mean25 = mean_by_condition(samples["25"], "25")

# Merge conditions
df = mean0.merge(mean05, on=["chrom", "start", "end", "ref", "tx_pos"], how="outer")
df = df.merge(mean25, on=["chrom", "start", "end", "ref", "tx_pos"], how="outer")

# Reactivity columns from total mutation rate
df["reactivity_05_raw"] = df["mut_rate_mean_05"] - df["mut_rate_mean_0"]
df["reactivity_25_raw"] = df["mut_rate_mean_25"] - df["mut_rate_mean_0"]

# Clipped versions are often convenient for plotting
df["reactivity_05"] = df["reactivity_05_raw"].clip(lower=0)
df["reactivity_25"] = df["reactivity_25_raw"].clip(lower=0)

# Reactivity for each mutation type
for m in mut_cols:
    df[f"{m}_reactivity_05_raw"] = df[f"{m}_mean_05"] - df[f"{m}_mean_0"]
    df[f"{m}_reactivity_25_raw"] = df[f"{m}_mean_25"] - df[f"{m}_mean_0"]
    df[f"{m}_reactivity_05"] = df[f"{m}_reactivity_05_raw"].clip(lower=0)
    df[f"{m}_reactivity_25"] = df[f"{m}_reactivity_25_raw"].clip(lower=0)

# Sort by transcript position if available, else genomic position
sort_cols = ["tx_pos"] if gene_region is not None else ["chrom", "start", "end"]
df = df.sort_values(sort_cols)

# Write outputs
prefix = gene_name if gene_name is not None else "MT_all"
out_tsv = outdir / f"{prefix}.mean_reactivity.tsv"
out_csv = outdir / f"{prefix}.mean_reactivity.csv"

df.to_csv(out_tsv, sep="\t", index=False)
df.to_csv(out_csv, index=False)

print(f"Wrote: {out_tsv}")
print(f"Wrote: {out_csv}")
print(df.head())
