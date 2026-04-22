#!/usr/bin/env python3

import sys
from collections import defaultdict

# -------------------------
# Input / Output
# -------------------------

mpa_file = sys.argv[1]
out_file = sys.argv[2]


# -------------------------
# Load data
# -------------------------

def clean(x):
    # Remove taxonomic prefix (e.g., d__, p__, f__)
    return x.split("__")[-1]

# Store all taxonomic paths with their counts
paths = []

for line in open(mpa_file):
    tax, count = line.strip().split("\t")
    count = float(count)

    # Split taxonomy into hierarchical nodes
    nodes = [clean(n) for n in tax.split("|")]

    paths.append((nodes, count))


# -------------------------
# Reformat
# -------------------------

# Sort paths from deepest (most specific) to shallowest (most general)
paths.sort(key=lambda x: len(x[0]), reverse=True)

# Store exclusive counts per node
exclusive = defaultdict(float)

for nodes, count in paths:
    # Assign count to the most specific node
    exclusive[tuple(nodes)] += count

    # Subtract contribution from all parent nodes
    for i in range(len(nodes) - 1):
        parent = tuple(nodes[:i+1])
        exclusive[parent] -= count

# Remove negative values caused by over-subtraction
for k in list(exclusive.keys()):
    if exclusive[k] < 0:
        exclusive[k] = 0


# -------------------------
# Krona-compatible output
# -------------------------

# count <tab> taxonomic path
with open(out_file, "w") as f:
    for path, count in exclusive.items():
        if count > 0:
            f.write(str(int(count)) + "\t" + "\t".join(path) + "\n")