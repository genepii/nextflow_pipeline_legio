#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt


# -------------------------
# Input / Output
# -------------------------

if len(sys.argv) != 4:
    sys.exit("Usage: blast_amplicons_plot_blast.py <input.tsv> <total_seq> <output.png>")

filtered_file = sys.argv[1]
total_seq = int(sys.argv[2])
output_file = sys.argv[3]


# -------------------------
# Load data
# -------------------------
df = pd.read_csv(
    filtered_file,
    sep="\t",
    header=None,
    names=["seq", "taxo", "count"]
)


# -------------------------
# reFormat
# -------------------------
# Parse taxonomy
taxo_split = df["taxo"].str.split("|", expand=True)

df["species"] = taxo_split[0]
#df["strain"] = taxo_split[1].fillna("strain=NA")
if 1 in taxo_split.columns:
    df["strain"] = taxo_split[1].fillna("strain=NA")
else:
    df["strain"] = "strain=NA"

# Aggregation (species + strain)
agg = (
    df.groupby(["species", "strain"], as_index=False)["count"]
    .sum()
)


# -------------------------
# Add "not assigned"
# -------------------------
assigned_sum = agg["count"].sum()
not_assigned = total_seq - assigned_sum

# No negative value
if not_assigned < 0:
    not_assigned = 0

# Add line if not_assigned > 0
if not_assigned > 0:
    agg = pd.concat(
        [
            agg,
            pd.DataFrame([{
                "species": "Not_assigned",
                "strain": "strain=NA",
                "count": not_assigned
            }])
        ],
        ignore_index=True
    )


# -------------------------
# Abundance and labels
# -------------------------
# Relative abundance (%)
agg["rel_abundance"] = (agg["count"] / total_seq) * 100


# Sorting + labels
agg = agg.sort_values("rel_abundance")

agg["label"] = agg["species"] + ", " + agg["strain"]


# -------------------------
# Colors
# -------------------------
## Default in deepskyblue3, "Not_assigned" in azure4 
## and "L-ambiguous" in deepskyblue4
colors = [
    "#838B8B" if sp == "Not_assigned"
    else "#00688B" if sp == "L-ambiguous"
    else "#009ACD"
    for sp in agg["species"]
]


# -------------------------
# Plots
# -------------------------
plt.style.use('ggplot')

# Size
fig, ax = plt.subplots(figsize=(10, 6))

# Content
bars = ax.barh(agg["label"], agg["rel_abundance"], color=colors)

# X axis limit
ax.set_xlim(0, 100)

# Axis label
ax.set_xlabel("Relative abundance (%)")
ax.set_ylabel("Species, Strain")

# Title + subtitle
ax.set_title(
    f"Taxonomic abundance (including not assigned)\nTotal sequences: {total_seq}",
    fontsize=12,
    pad=15
)

# Add labels inside bars
xmax = ax.get_xlim()[1]
for bar, count in zip(bars, agg["count"]):
    width = bar.get_width()
    y = bar.get_y() + bar.get_height() / 2

    # Position: inside if possible, otherwise outside
    if width > 5:
        x_text = min(width - 1, xmax * 0.98) #if bar = 100%
        ax.text(
            x_text,
            y,
            f"{int(count)}",
            va="center",
            ha="right",
            fontsize=8,
            color="white",
            fontweight="bold",
            clip_on=True
        )
    else:
        ax.text(
            width + 0.5,
            y,
            f"{int(count)}",
            va="center",
            ha="left",
            fontsize=8,
            color="black",
            fontweight="bold"
        )

plt.tight_layout()

# Save
plt.savefig(output_file, dpi=300, bbox_inches="tight")
plt.close()