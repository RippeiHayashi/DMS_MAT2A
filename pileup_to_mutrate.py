#!/usr/bin/env python3

import sys

MUT_KEYS = [
    "A>C", "A>G", "A>T",
    "C>A", "C>G", "C>T",
    "G>A", "G>C", "G>T",
    "T>A", "T>C", "T>G",
]

def parse_indel(bases, i):
    """
    Skip over an insertion/deletion in the mpileup read-bases string.
    bases[i] is '+' or '-'
    """
    i += 1
    n = 0
    while i < len(bases) and bases[i].isdigit():
        n = n * 10 + int(bases[i])
        i += 1
    i += n
    return i

def count_from_pileup(ref, bases):
    """
    Parse the mpileup bases column.

    Counts:
      - depth (matches + mismatches; excludes deletions/ref-skips)
      - total mismatches
      - mutation-type counts (e.g. A>G, C>T)

    Ignores:
      - read start/end markers (^ and $)
      - indel annotations (+/-<len><seq>)
      - deletion placeholders (* and #)
      - reference skips (< and >), common with spliced RNA-seq
      - N/n and other non-ACGT symbols
    """
    ref = ref.upper()
    i = 0
    depth = 0
    mismatches = 0
    mut_counts = {k: 0 for k in MUT_KEYS}

    while i < len(bases):
        c = bases[i]

        # start of read segment; next char is mapping quality
        if c == '^':
            i += 2
            continue

        # end of read segment
        if c == '$':
            i += 1
            continue

        # insertion / deletion annotation
        if c == '+' or c == '-':
            i = parse_indel(bases, i)
            continue

        # deletion placeholder
        if c == '*' or c == '#':
            i += 1
            continue

        # reference skip (e.g. spliced N in CIGAR)
        if c == '<' or c == '>':
            i += 1
            continue

        # match to reference on forward/reverse strand
        if c == '.' or c == ',':
            depth += 1
            i += 1
            continue

        b = c.upper()

        # mismatch base
        if b in ("A", "C", "G", "T"):
            depth += 1
            if b != ref:
                mismatches += 1
                key = f"{ref}>{b}"
                if key in mut_counts:
                    mut_counts[key] += 1
            i += 1
            continue

        # skip N/n and any other symbols
        i += 1

    return depth, mismatches, mut_counts

def main():
    header = [
        "chrom", "start", "end", "ref",
        "depth", "mismatches", "mut_rate",
        "A>C", "A>G", "A>T",
        "C>A", "C>G", "C>T",
        "G>A", "G>C", "G>T",
        "T>A", "T>C", "T>G",
    ]
    sys.stdout.write("\t".join(header) + "\n")

    for line in sys.stdin:
        if not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            continue

        chrom = fields[0]
        pos1 = int(fields[1])   # 1-based
        ref = fields[2].upper()
        bases = fields[4]

        if ref not in ("A", "C", "G", "T"):
            continue

        depth, mismatches, mut_counts = count_from_pileup(ref, bases)
        if depth == 0:
            continue

        mut_rate = mismatches / depth
        pos0 = pos1 - 1

        row = [
            chrom, str(pos0), str(pos1), ref,
            str(depth), str(mismatches), f"{mut_rate:.6g}",
            str(mut_counts["A>C"]), str(mut_counts["A>G"]), str(mut_counts["A>T"]),
            str(mut_counts["C>A"]), str(mut_counts["C>G"]), str(mut_counts["C>T"]),
            str(mut_counts["G>A"]), str(mut_counts["G>C"]), str(mut_counts["G>T"]),
            str(mut_counts["T>A"]), str(mut_counts["T>C"]), str(mut_counts["T>G"]),
        ]
        sys.stdout.write("\t".join(row) + "\n")

if __name__ == "__main__":
    main()
