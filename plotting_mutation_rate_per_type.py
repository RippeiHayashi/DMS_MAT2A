#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'

# Input data from your awk output
data = {
    "HepG2_DMS_0_pA_rep1": {
        "A>C": 0.000494942, "A>G": 0.000786102, "A>T": 0.000316835,
        "C>A": 0.000307588, "C>G": 0.000187601, "C>T": 0.000747352,
        "G>A": 0.00058771,  "G>C": 0.0000552423, "G>T": 0.0000319064,
        "T>A": 0.000531682, "T>C": 0.00111672, "T>G": 0.000432454
    },
    "HepG2_DMS_0_pA_rep2": {
        "A>C": 0.000437136, "A>G": 0.000738417, "A>T": 0.000290154,
        "C>A": 0.000281143, "C>G": 0.000168516, "C>T": 0.000724836,
        "G>A": 0.000575909, "G>C": 0.0000489514, "G>T": 0.0000294893,
        "T>A": 0.000485664, "T>C": 0.000976828, "T>G": 0.000393684
    },
    "HepG2_DMS_05_pA_rep1": {
        "A>C": 0.000950493, "A>G": 0.000964448, "A>T": 0.00455207,
        "C>A": 0.00239646,  "C>G": 0.000787949, "C>T": 0.00632836,
        "G>A": 0.000452117, "G>C": 0.0000621862, "G>T": 0.0000787625,
        "T>A": 0.00078846,  "T>C": 0.000867957, "T>G": 0.000456524
    },
    "HepG2_DMS_05_pA_rep2": {
        "A>C": 0.00101897,  "A>G": 0.000997203, "A>T": 0.0045757,
        "C>A": 0.00232883,  "C>G": 0.000793643, "C>T": 0.00628798,
        "G>A": 0.000444078, "G>C": 0.0000665171, "G>T": 0.0000826707,
        "T>A": 0.000860092, "T>C": 0.000960401, "T>G": 0.00047989
    },
    "HepG2_DMS_25_pA_rep1": {
        "A>C": 0.00128254,  "A>G": 0.00126338, "A>T": 0.00655028,
        "C>A": 0.00300922,  "C>G": 0.00103882, "C>T": 0.00781817,
        "G>A": 0.00089583,  "G>C": 0.000157253, "G>T": 0.000291727,
        "T>A": 0.00282474,  "T>C": 0.00104108, "T>G": 0.000792393
    },
    "HepG2_DMS_25_pA_rep2": {
        "A>C": 0.00138009,  "A>G": 0.00141193, "A>T": 0.00672407,
        "C>A": 0.00320851,  "C>G": 0.00114832, "C>T": 0.00816062,
        "G>A": 0.000931041, "G>C": 0.000162684, "G>T": 0.00029166,
        "T>A": 0.00272783,  "T>C": 0.00125605, "T>G": 0.00088164
    }
}

mutation_order = [
    "A>C", "A>G", "A>T",
    "C>A", "C>G", "C>T",
    "G>A", "G>C", "G>T",
    "T>A", "T>C", "T>G"
]

sample_order = [
    "HepG2_DMS_0_pA_rep1",
    "HepG2_DMS_0_pA_rep2",
    "HepG2_DMS_05_pA_rep1",
    "HepG2_DMS_05_pA_rep2",
    "HepG2_DMS_25_pA_rep1",
    "HepG2_DMS_25_pA_rep2"
]

# Convert to tidy dataframe
rows = []
for sample in sample_order:
    for mut in mutation_order:
        rows.append({
            "sample": sample,
            "mutation_type": mut,
            "mutation_rate": data[sample][mut]
        })

df = pd.DataFrame(rows)

# Plot
fig, ax = plt.subplots(figsize=(14, 6))

for sample in sample_order:
    sub = df[df["sample"] == sample]
    ax.plot(
        sub["mutation_type"],
        sub["mutation_rate"],
        marker="o",
        linewidth=1.5,
        label=sample
    )

ax.set_xlabel("Mutation type")
ax.set_ylabel("Mutation rate")
ax.set_title("chrMT mutation spectrum across 6 libraries")
ax.tick_params(axis="x", rotation=45)
ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
ax.grid(True, axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig("chrMT_mutation_spectrum.png", dpi=300)
plt.savefig("chrMT_mutation_spectrum.pdf")
plt.show()
