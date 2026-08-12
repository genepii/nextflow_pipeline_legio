#!/usr/bin/env python3

"""
Create a circular phylogenetic tree colored by metadata groups.

Inputs:
    - Newick tree file
    - Metadata TSV file containing strain identifiers and annotation groups
    - Metadata column used for coloring
    - Output prefix

Output:
    - SVG figure
"""

import os
import sys


# ----------------------------------------------------------------------
# Configure writable temporary directories
# ----------------------------------------------------------------------

tmp_dir = os.environ.get(
    "TMPDIR",
    "/tmp"
)

matplotlib_cache = os.path.join(
    tmp_dir,
    "matplotlib"
)

font_cache = os.path.join(
    tmp_dir,
    "fontconfig"
)

os.makedirs(
    matplotlib_cache,
    exist_ok=True
)

os.makedirs(
    font_cache,
    exist_ok=True
)

os.environ["QT_QPA_PLATFORM"] = "offscreen"
os.environ["MPLCONFIGDIR"] = matplotlib_cache
os.environ["XDG_CACHE_HOME"] = tmp_dir
os.environ["FONTCONFIG_CACHE"] = font_cache


# ----------------------------------------------------------------------
# Imports
# ----------------------------------------------------------------------

import pandas as pd
import matplotlib as mpl

from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace


# ----------------------------------------------------------------------
# Arguments
# ----------------------------------------------------------------------

if len(sys.argv) != 5:
    raise ValueError(
        "Usage: script.py <tree.nwk> <metadata.tsv> "
        "<annotation_column> <output_prefix>"
    )


tree_file = sys.argv[1]
metadata_file = sys.argv[2]
annotation_column = sys.argv[3]
output_prefix = sys.argv[4]


# ----------------------------------------------------------------------
# Load tree
# ----------------------------------------------------------------------

tree = Tree(
    tree_file,
    format=1
)


# ----------------------------------------------------------------------
# Load metadata
# ----------------------------------------------------------------------

metadata = pd.read_csv(
    metadata_file,
    sep="\t",
    dtype={"ID": str}
)

required_columns = {
    "ID",
    annotation_column
}

if not required_columns.issubset(metadata.columns):
    raise ValueError(
        "Metadata must contain columns: "
        + ", ".join(required_columns)
    )


# Keep identifiers as strings and preserve leading zeros
metadata["ID"] = (
    metadata["ID"]
    .astype(str)
    .str.strip()
)

metadata[annotation_column] = (
    metadata[annotation_column]
    .astype(str)
    .str.strip()
)


strain_to_group = dict(
    zip(
        metadata["ID"],
        metadata[annotation_column]
    )
)


# ----------------------------------------------------------------------
# Generate dynamic colors
# ----------------------------------------------------------------------

tree_strains = {
    leaf.name
    for leaf in tree.iter_leaves()
}


groups = sorted(
    {
        strain_to_group[strain]
        for strain in tree_strains
        if strain in strain_to_group
    }
)


if groups:

    colormap = mpl.colormaps.get_cmap(
        "tab20"
    ).resampled(
        len(groups)
    )

    group_colors = {
        group: "#%02x%02x%02x" % tuple(
            int(channel * 255)
            for channel in colormap(index)[:3]
        )
        for index, group in enumerate(groups)
    }

else:

    group_colors = {}


UNKNOWN_COLOR = "#000000"


# ----------------------------------------------------------------------
# Styling functions
# ----------------------------------------------------------------------

def apply_branch_style(node, color, width, size):
    """
    Apply branch and node colors.
    """

    style = NodeStyle()

    style["fgcolor"] = color
    style["hz_line_color"] = color
    style["vt_line_color"] = color

    style["hz_line_width"] = width
    style["vt_line_width"] = width

    style["size"] = size

    node.set_style(style)


def add_leaf_label(node, color, size):
    """
    Add colored leaf label.
    """

    face = AttrFace(
        "name",
        fgcolor=color,
        fsize=size
    )

    node.add_face(
        face,
        column=0,
        position="branch-right"
    )


# ----------------------------------------------------------------------
# Apply default style to all nodes and branches
# ----------------------------------------------------------------------

DEFAULT_COLOR = "#808080"
## Note: internal nodes color 

for node in tree.traverse():

    apply_branch_style(
        node,
        DEFAULT_COLOR,
        width=1,
        size=10
    )


# ----------------------------------------------------------------------
# Annotate terminal leaves only
# ----------------------------------------------------------------------

for leaf in tree.iter_leaves():

    strain = leaf.name

    if strain in strain_to_group:

        group = strain_to_group[strain]
        color = group_colors[group]

        # Color terminal branch and leaf node
        apply_branch_style(
            leaf,
            color,
            width=2,
            size=12
        )

        # Color leaf name
        add_leaf_label(
            leaf,
            color,
            size=12
        )

    else:

        # Unknown strains remain grey
        apply_branch_style(
            leaf,
            UNKNOWN_COLOR,
            width=1,
            size=10
        )

        add_leaf_label(
            leaf,
            UNKNOWN_COLOR,
            size=11
        )


# ----------------------------------------------------------------------
# Configure tree style
# ----------------------------------------------------------------------

tree_style = TreeStyle()

tree_style.mode = "c"
tree_style.show_leaf_name = False
tree_style.scale = 50


# Move tree and legend away from SVG borders.
# 8 px corresponds approximately to 0.2 cm at 96 DPI.
tree_style.margin_left = 8
tree_style.margin_right = 8
tree_style.margin_top = 8
tree_style.margin_bottom = 8


# ----------------------------------------------------------------------
# Legend
# ----------------------------------------------------------------------

legend_column = 0

for group, color in sorted(group_colors.items()):

    tree_style.legend.add_face(
        faces.CircleFace(
            radius=6,
            color=color
        ),
        column=legend_column
    )

    tree_style.legend.add_face(
        faces.TextFace(
            f" {group}",
            fsize=10
        ),
        column=legend_column + 1
    )

    legend_column += 2


tree_style.legend.add_face(
    faces.CircleFace(
        radius=6,
        color=UNKNOWN_COLOR
    ),
    column=legend_column
)

tree_style.legend.add_face(
    faces.TextFace(
        " Previous",
        fsize=10
    ),
    column=legend_column + 1
)

tree_style.legend_position = 1


# ----------------------------------------------------------------------
# Render SVG output
# ----------------------------------------------------------------------

svg_output = f"{output_prefix}.svg"


tree.render(
    svg_output,
    tree_style=tree_style,
    w=2000,
    units="px"
)
