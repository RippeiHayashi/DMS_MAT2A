#!/usr/bin/env python3

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from pyfaidx import Fasta
except ImportError:
    raise SystemExit("ERROR: pyfaidx is required. Install with: pip install pyfaidx")


RCMAP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(RCMAP)[::-1]


def run_rnafold(seq: str):
    cmd = ["RNAfold", "--noPS"]
    p = subprocess.run(
        cmd,
        input=(seq + "\n").encode(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True
    )
    lines = p.stdout.decode().strip().splitlines()
    if len(lines) < 2:
        raise RuntimeError("Unexpected RNAfold output")
    folded_seq = lines[0].strip()
    parts = lines[1].strip().split()
    struct = parts[0]
    mfe = None
    if len(parts) > 1:
        mfe = parts[-1].strip("()")
    return folded_seq, struct, mfe


def dotbracket_pairs(struct: str):
    stack = []
    pairs = {}
    for i, ch in enumerate(struct, start=1):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if not stack:
                raise ValueError("Unbalanced dot-bracket structure")
            j = stack.pop()
            pairs[i] = j
            pairs[j] = i
    if stack:
        raise ValueError("Unbalanced dot-bracket structure")
    return pairs


def make_arc_plot(df, out_png, title, max_arc_pairs=800):
    """
    Arc diagram with nucleotides colored by reactivity.
    For long RNAs, plotting every base pair can get dense, so max_arc_pairs limits arcs.
    """
    L = len(df)
    x = np.arange(1, L + 1)

    fig, ax = plt.subplots(figsize=(18, 5))

    # nucleotide baseline colored by reactivity
    sc = ax.scatter(
        x,
        np.zeros(L),
        c=df["reactivity_25_plot"],
        s=18,
        cmap="viridis",
        vmin=0,
        vmax=np.nanquantile(df["reactivity_25_plot"], 0.95) if df["reactivity_25_plot"].notna().any() else 0.05
    )

    # draw arcs for base pairs
    pairs_drawn = 0
    seen = set()
    for _, row in df.iterrows():
        i = int(row["tx_pos"])
        j = row["pair_pos"]
        if pd.isna(j):
            continue
        j = int(j)
        if i >= j:
            continue
        if (i, j) in seen:
            continue
        seen.add((i, j))
        span = j - i
        xs = np.linspace(i, j, 100)
        # semicircle-like arc
        ys = np.sqrt(np.maximum(0, (span / 2.0) ** 2 - (xs - (i + j) / 2.0) ** 2))
        ax.plot(xs, ys, linewidth=0.5, alpha=0.35)
        pairs_drawn += 1
        if pairs_drawn >= max_arc_pairs:
            break

    ax.set_xlim(0, L + 1)
    ax.set_ylim(bottom=-0.2)
    ax.set_xlabel("MT-CO1 transcript position")
    ax.set_ylabel("Base-pair arc height")
    ax.set_title(title)
    ax.set_yticks([])
    cbar = plt.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("Reactivity (2.5% - 0%)")

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def main():
    analysis = Path("/scratch/lf10/rh1772/MAT2A/analysis/RNAseq/DMS_reactivity")
    genome_fa = "/g/data/lf10/as7425/genomes/human_genome/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
    input_tsv = analysis / "MT-CO1.mean_reactivity.tsv"

    # Ensembl / GRCh38 MT-CO1 coordinates
    chrom = "MT"
    gene_name = "MT-CO1"
    gene_start_1based = 5904
    gene_end_1based = 7445
    strand = "+"

    out_prefix = analysis / gene_name

    # read reactivity table
    df = pd.read_csv(input_tsv, sep="\t")
    if "tx_pos" not in df.columns:
        raise ValueError("Input table must contain tx_pos column")

    df = df.sort_values("tx_pos").copy()

    # fetch sequence from genome
    fasta = Fasta(genome_fa, as_raw=True, sequence_always_upper=True)
    seq = fasta[chrom][gene_start_1based - 1:gene_end_1based].upper()
    if strand == "-":
        seq = revcomp(seq)

    seq = seq.replace("T", "U")

    # run RNAfold
    folded_seq, struct, mfe = run_rnafold(seq)
    if len(folded_seq) != len(struct):
        raise RuntimeError("Sequence and structure lengths do not match")

    # map pairs
    pairs = dotbracket_pairs(struct)

    # make per-position table
    fold_df = pd.DataFrame({
        "tx_pos": np.arange(1, len(seq) + 1),
        "base": list(folded_seq),
        "dotbracket": list(struct),
        "pair_pos": [pairs.get(i, np.nan) for i in range(1, len(seq) + 1)],
    })

    # merge with reactivity
    merged = fold_df.merge(df, on="tx_pos", how="left")

    # plotting column
    merged["reactivity_25_plot"] = merged["reactivity_25"].clip(lower=0)

    # save outputs
    fasta_out = out_prefix.with_suffix(".folded.fa")
    db_out = out_prefix.with_suffix(".folded.dotbracket.txt")
    table_out = out_prefix.with_suffix(".structure_with_reactivity.tsv")

    with open(fasta_out, "w") as fh:
        fh.write(f">{gene_name}\n{folded_seq}\n")

    with open(db_out, "w") as fh:
        fh.write(f">{gene_name}\n")
        fh.write(f"{folded_seq}\n")
        fh.write(f"{struct}\n")
        if mfe is not None:
            fh.write(f"# MFE {mfe}\n")

    merged.to_csv(table_out, sep="\t", index=False)

    # --------------------------------------------------
    # Plot 1: linear reactivity with paired/unpaired bars
    # --------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(18, 7), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]}
    )

    ax1.plot(
        merged["tx_pos"],
        merged["reactivity_25_plot"],
        linewidth=1.2
    )
    ax1.set_ylabel("Reactivity (2.5% - 0%)")
    ttl = f"{gene_name} RNAfold structure with DMS reactivity"
    if mfe is not None:
        ttl += f" | MFE {mfe}"
    ax1.set_title(ttl)
    ax1.grid(True, axis="y", alpha=0.3)

    paired = merged["pair_pos"].notna().astype(int)
    ax2.fill_between(
        merged["tx_pos"],
        0,
        paired,
        step="mid",
        alpha=0.6
    )
    ax2.set_ylabel("Paired")
    ax2.set_xlabel("Transcript position")
    ax2.set_ylim(0, 1.1)
    ax2.set_yticks([0, 1])
    ax2.set_yticklabels(["no", "yes"])

    plt.tight_layout()
    plt.savefig(out_prefix.with_suffix(".reactivity_on_structure.png"), dpi=300)
    plt.savefig(out_prefix.with_suffix(".reactivity_on_structure.pdf"))
    plt.close()

    # --------------------------------------------------
    # Plot 2: arc diagram
    # --------------------------------------------------
    make_arc_plot(
        merged,
        out_prefix.with_suffix(".arc_reactivity.png"),
        f"{gene_name} arc diagram colored by DMS reactivity"
    )

    print("Wrote:")
    print(fasta_out)
    print(db_out)
    print(table_out)
    print(out_prefix.with_suffix(".reactivity_on_structure.png"))
    print(out_prefix.with_suffix(".reactivity_on_structure.pdf"))
    print(out_prefix.with_suffix(".arc_reactivity.png"))


if __name__ == "__main__":
    main()
