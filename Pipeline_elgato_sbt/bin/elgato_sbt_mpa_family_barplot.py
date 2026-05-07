#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt

# -------------------------
# Input / Output
# -------------------------

if len(sys.argv) != 5:
    sys.exit("Usage: blast_amplicons_mpa_family_barplot.py <input_mpa.tsv> <total_reads> <output.tsv> <output.png>")

mpa_file = sys.argv[1]
total_reads = int(sys.argv[2])
tsv_out = sys.argv[3]
png_out = sys.argv[4]


# -------------------------
# Load data
# -------------------------
rows = []

with open(mpa_file, "r") as f:
    for line in f:
        cols = line.strip().split("\t")

        count = cols[0] if len(cols) > 0 else "NA"
        family = cols[5] if len(cols) > 5 else "NA"

        rows.append([count, family])

df = pd.DataFrame(rows, columns=["count", "family"])
df = df.fillna("NA")

# Ensure numeric counts
df["count"] = pd.to_numeric(df["count"], errors="coerce")
df = df.dropna(subset=["count"])


# -------------------------
# Sum count
# -------------------------
df = df.groupby("family", as_index=False)["count"].sum()

# Normalisation in percent
df["abundance"] = (df["count"] / total_reads) * 100

# Keep top 15 but ensure Legionellaceae is present
top = df.sort_values("abundance", ascending=False).head(15)

if "Legionellaceae" in df["family"].values and "Legionellaceae" not in top["family"].values:
    legion = df[df["family"] == "Legionellaceae"]
    top = pd.concat([top.iloc[:-1], legion])

df = top.sort_values("abundance", ascending=False)


# -------------------------
# Save TSV
# -------------------------
df.to_csv(tsv_out, sep="\t", index=False)


# -------------------------
# Colors
# -------------------------
## Default in dark grey/blue, "Not_assigned" in azure4 
## and "Legionellaceae" in deepskyblue3
colors = [
    "#838B8B" if fam == "Not_assigned"
    else "#009ACD" if fam == "Legionellaceae"
    else "#4A7CD870"
    for fam in df["family"]
]


# -------------------------
# Plots
# -------------------------
plt.style.use('ggplot')

# Size
fig, ax = plt.subplots(figsize=(10, 6))

# Content
bars = ax.barh(df["family"], df["abundance"], color=colors)
ax.invert_yaxis()

# X axis limit
ax.set_xlim(0, 100)

# Axis label
ax.set_xlabel("Relative abundance (%)")
ax.set_ylabel("Family")

# Title + subtitle
ax.set_title(
    f"Top families detected by Kraken2\nTotal reads: {total_reads} | Reads shown: {int(df['count'].sum())}",
    fontsize=12,
    pad=15
)

# Add labels inside bars
xmax = ax.get_xlim()[1]
for bar, count in zip(bars, df["count"]):
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
plt.savefig(png_out, dpi=300, bbox_inches="tight")
plt.close()