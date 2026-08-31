#!/usr/bin/env python3
"""
assembly_mlst_summary_html.py

Generate HTML report from merged TSV.
Standalone version with inline Jinja2 template.
"""

import argparse
from pathlib import Path
from datetime import datetime

import pandas as pd
from jinja2 import Template


# Main summary columns displayed for each strain section
MAIN_COLUMNS = [
    "Sample_ID", "Total_bases", "Total_length", "Depth",
    "Number_of_contigs", "GC", "N50", "auN", "L90",
    "Legionella_pneumophila_percent",
    "Legionella_spp_percent", "Kraken2_results",
    "FastANI_strain", "FastANI_value",
    "ST", "flaA", "pilE", "asd", "mip",
    "mompS", "proA", "neuA", "lag",
    "lpeA", "lpeB"
]

# AMR summary columns
AMR_COLUMNS = [
    "Sample_ID",
    "AMR_Nb_Mutated_Genes",
    "AMR_Nb_with_Impact",
    "AMR_Nb_Non_Coding"
]

# Other summary columns
OTHER_COLUMNS = [
    "Sample_ID",
    "Other_Nb_Mutated_Genes",
    "Other_Nb_with_Impact",
    "Other_Nb_Non_Coding"
]

# Description of Main summary columns
COLUMN_DESC = {
    "Sample_ID": "Unique sample identifier",
    "Total_bases": "Total number of sequencing bases generated",
    "Total_length": "Assembly genome length",
    "Depth": "Sequencing depth (Total_bases / genome size)",
    "Number_of_contigs": "Number of contigs in the assembly (> threshold)",
    "GC": "GC content of the genome",
    "N50": "N50 contig statistic",
    "auN": "Area under N50 curve metric",
    "L90": "Number of contigs covering 90% of the genome",
    "Legionella_pneumophila_percent": "Proportion of reads assigned to L. pneumophila by Kraken2",
    "Legionella_spp_percent": "Proportion of reads assigned to Legionella spp. by Kraken2",
    "Kraken2_results": "Kraken2-based classification, Legionella pneumophila if present; Legionella spp. if only genus is detected; extended if only Legionellaceae is assigned; others are considered contamination",
    "FastANI_strain": "Closest reference strain determined by FastANI",
    "FastANI_value": "Average nucleotide identity with closest strain",
    "ST": "Sequence type (MLST)",
    "flaA": "Allele of flaA gene",
    "pilE": "Allele of pilE gene",
    "asd": "Allele of asd gene",
    "mip": "Allele of mip gene",
    "mompS": "Allele of mompS gene",
    "proA": "Allele of proA gene",
    "neuA": "Allele of neuA gene",
    "lag": "Allele of lag gene",
    "lpeA": "Allele of lpeA gene",
    "lpeB": "Allele of lpeB gene",
    "XX_Nb_Mutated_Genes": "Number of mutated XX genes detected",
    "XX_Nb_with_Impact": "Number of XX mutations with predicted functional impact",
    "XX_Nb_Non_Coding": "Number of non-coding XX-related variants"
}

# Values = missing value
MISSING_VALUES = [
    "",
    "nan",
    "NaN",
    "NA",
    "Na",
    "N/A",
    "<NA>",
    "null",
    "NULL",
    "Null",
]

# Columns containing float values
FLOAT_COLUMNS = [
    "Depth",
    "GC",
    "FastANI_value",
    "auN",
]

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--analyse_id", required=True)
    parser.add_argument("--sequencing_id", required=True)
    parser.add_argument("--software", required=True)
    return parser.parse_args()


def amr_warning(df):
    """Return True if at least one sample contains AMR mutations"""
    if "AMR_Nb_Mutated_Genes" not in df.columns:
        return False

    col = pd.to_numeric(df["AMR_Nb_Mutated_Genes"], errors="coerce")
    return (col.fillna(0) != 0).any()


def main():
    args = parse_args()

    # Load merged TSV
    df = pd.read_csv(args.input, sep="\t", dtype=str).sort_values("Sample_ID")

    # AMR table cleanup:
    # remove rows where every AMR field except Sample_ID is missing/nan
    amr_df = df.copy()

    amr_cols_without_sample = [
        c for c in AMR_COLUMNS
        if c != "Sample_ID" and c in amr_df.columns
    ]

    if amr_cols_without_sample:
        amr_tmp = amr_df[amr_cols_without_sample].apply(
            lambda s: s.astype(str).str.strip()
        )
        amr_df = amr_df[
            ~amr_tmp.replace(
                MISSING_VALUES,
                pd.NA
            )
            .isna()
            .all(axis=1)
        ]

    # create table only if OTHER_COLUMNS are present in input
    other_df = df.copy()
    other_available_columns = [
        c for c in OTHER_COLUMNS
        if c in other_df.columns
    ]

    if len(other_available_columns) == len(OTHER_COLUMNS):
        other_cols_without_sample = [
            c for c in OTHER_COLUMNS
            if c != "Sample_ID"
        ]
        other_tmp = other_df[other_cols_without_sample].apply(
            lambda s: s.astype(str).str.strip()
        )
        other_df = other_df[
            ~other_tmp.replace(
                MISSING_VALUES,
                pd.NA
            )
            .isna()
            .all(axis=1)
        ][OTHER_COLUMNS]
    else:
        other_df = pd.DataFrame()
        
    # Collect strain names for report sections
    strains = (
        sorted(df["FastANI_strain"].dropna().unique())
        if "FastANI_strain" in df.columns
        else []
    )

    # Read software/settings file
    try:
        software_text = Path(args.software).read_text(encoding="utf-8")
    except Exception:
        software_text = "Unable to read software file"

    template_str = r"""
<!DOCTYPE html>
<html>
<head>
<title>Report Seq. {{ sequencing_id }} - Analyse {{ analyse_id }}</title>

<style>
body {
    font-family: system-ui,-apple-system,"Segoe UI",Roboto,Arial,sans-serif;
    margin: 0;
    color: #263238;
}

h2 {
    text-decoration: underline;
}

.page {
    padding: 48px;
    padding-top: 4px; 
}

.titlebar {
    background: #263238;
    color: white;
    padding: 15px;
}

.report-title {
    font-size: 22px;
    font-weight: bold;
}

.report-subtitle {
    font-size: 13px;
    font-weight: bold;
    margin-top: 4px;
    margin-bottom: 4px;
    padding-left: calc(48px - 15px);
    color: white;
}

a {
    color: inherit;
    text-decoration: underline;
}

a:hover {
    color: color: #1C2529;
}

.info-block {
    margin-top: 4px;
    margin-bottom: 8px;
    font-size: 13px;
    line-height: 1.5;
}

.info-block a {
    font-weight: bold;
    color: #3949AB;
}

.info-block a:hover {
    color: #253179;
}

.section {
    margin-top: 40px;
}

.table-content {
    margin:10px 43px 0 43px;
}

.table-wrapper {
    width: calc(100% - 174px);
    margin-left: 87px;
    margin-right: 87px;
    margin-bottom: 10px;
}

.table-controls {
    width: 100%;
}

.table-scroll {
    overflow-x: auto;
    max-height: 600px;
    overflow-y: auto;
}

table {
    border-collapse: collapse;
    min-width: max-content;
}

th,td {
    border: 1px solid #CFD8DC;
    padding: 6px;
    font-size: 13px;
    white-space: nowrap;
}

th {
    position: sticky;
    top: 0;
    z-index 2;
    background: #455A64;
    color: white;
}

tr:hover td {
    background-color: #E8EAF6;
}

tr:hover td:first-child {
    background-color: #D5D8F0;
}

table th:first-child,
table td:first-child {
    position: sticky;
    left: 0;
    z-index: 1;
    background: white;
}

table th:first-child {
    top: 0;
    left: 0;
    z-index: 3;
    background: #455A64;
}

.table-description {
    width:calc(100% - 174px);
    margin:10px 87px 0 87px;
    font-size:13px;
    color: #607D8B;
    line-height:1.45;
}

.table-description ul {
    margin:4px 0 0 18px;
    padding:0;
}

.table-description li {
    margin:2px 0;
}

.amr-alert {
    background: #FFCDD2;
}

.amr-warning {
    background-color: #C62828;
    color: white;
    padding: 15px;
    font-weight: bold;
}

.amr-warning a {
    color: inherit;
    text-decoration: underline;
}

.amr-warning a:hover {
    color: #FFCDD2;
}

hr {
    border:none;
    border-top: 1px solid #e5e5e5;
    margin: 40px 0;
}

.hidden {
    display:none;
}

#settings {
    margin-left: 87px;
}

button {
    color: white;
    background-color: #546E7A;
    font-weight: 500;
    border-radius: 8px;
    font-size: 12px;
    line-height: 20px;
    width: 60px;
    height: 26px;
    padding: 3px 10px;
    cursor: pointer;
    text-align: center;
    margin-bottom: 6px;
    display: flex;
    align-items: center;
    justify-content: center;
    border: none;
}

button:hover {
    background-color: #37474F;
}

button.copy-button {
    position: relative;
    left: -44px;
}

button.show-button {
    margin-left: 43px;
}

.highlight {
    color: #c00000;
    font-weight: bold;
}
</style>

<script>
function toggleSettings(button) {
    const settings = document.getElementById("settings");

    settings.classList.toggle("hidden");

    button.innerText = settings.classList.contains("hidden")
        ? "Show"
        : "Hide";
}

async function copyTable(button) {
    const wrapper = button.closest(".table-wrapper");
    const table = wrapper?.querySelector("table");

    if (!table) {
        console.error("Table not found.");
        return;
    }

    const text = Array.from(table.rows)
        .map(row =>
            Array.from(row.cells)
                .map(cell => cell.innerText.trim())
                .join("\t")
        )
        .join("\n");

    try {
        await navigator.clipboard.writeText(text);

        const originalText = button.innerText;
        button.innerText = "Copied";

        setTimeout(() => {
            button.innerText = originalText;
        }, 1500);
    } catch (error) {
        console.error("Failed to copy table:", error);
    }
}
</script>

</head>

<body>

<div class="titlebar">

    <div class="report-title">
        Analysis Report of Sequencing {{ sequencing_id }} - Analyse {{ analyse_id }}
    </div>

    <div class="report-subtitle">
        <a href="https://github.com/genepii/nextflow_pipeline_legio/tree/main/Pipeline_assembly_mlst">
            Assembly + MLST Pipeline (Nextflow)
        </a>
    </div>

</div>

<div class="page">

<div class="info-block">
    <div>{{ date }}</div>
    <div>
        This pipeline assembles Illumina data and performs MLST-based strain characterization with mutation screening.
    </div>
    <br>
    <div>
        <a href="https://cohesive.izs.it/spread/">Link to SPREAD</a> for an interactive view of the MLST results; 
        drag and drop Newick trees (.nwk) and metadata onto it to view them. 
    </div>
</div>

<hr>

{% if amr_warning %}
<div class="amr-warning">
AMR variants detected : <a href="#AMR">see table</a>
</div>
<hr>
{% endif %}

<h2>Table of Contents</h2>

<div class="table-content">
<ol>
<li><a href="#Quality"><b>Quality indicators</b></a></li>
<li><a href="#Desc"><b>Column descriptions</b></a></li>

<li>
<a href="#Mlst"><b>Assembly &amp; MLST Summary</b></a>
    <ul>
    {% for s in strains %}
        <li>
            <a href="#{{ s }}">
                <i>{{ s.replace("_", " ") }}</i>
            </a>
        </li>
    {% endfor %}
    </ul>
</li>

<li><a href="#AMR"><b>AMR Summary</b></a></li>
<li><a href="#OTHER"><b>Other Summary</b></a></li>
<li><a href="#Settings"><b>Settings</b></a></li>
</ol>

</div>

<hr>

<div id="Quality">
<h2>Quality indicators</h2>

<div class="table-description">

Values highlighted in <span class="highlight">red and bold</span> indicate potential quality issues or results requiring attention.

<br>
<ul>
    <li><b>Depth</b> &lt; 30X</li>
    <li><b>Number_of_contigs</b> &gt; 200</li>
    <li><b>Total_length</b> &lt; 2,500,000 b or &gt; 4,500,000 b</li>
    <li><b>GC</b> &lt; 36% or &gt; 40%</li>
    <li>Any cell containing <b>Contamination</b></li>
    <li>Any cell containing <b>Potential new spp.</b></li>
    <li>In the AMR summary, rows with <b>AMR_Nb_with_Impact &gt; 0</b> are highlighted with a light red background.</li>
</ul>
<br>

These highlights are intended as visual quality-control indicators and should be interpreted together with the complete analysis results.

</div>
</div>

<hr>

<div id="Desc" >
<h2>Column descriptions</h2>

<div class="table-description">

<ul>
{% for col in column_desc %}
<li><b>{{ col }}</b>: {{ column_desc[col] }}</li>
{% endfor %}
</ul>

</div>

</div>

<hr>

<div id="Mlst">
<h2>Assembly &amp; MLST Summary</h2>

{% for s in strains %}
<div class="section" id="{{ s }}">

<h3><i>{{ s.replace("_", " ") }}</i></h3>

<div class="table-wrapper">

<div class="table-controls">
<button class="copy-button" onclick="copyTable(this)">Copy</button>
</div>

<div class="table-scroll">
<table>

<tr>
{% for col in main_columns %}
<th>{{ col }}</th>
{% endfor %}
</tr>

{% for _, row in df[df["FastANI_strain"] == s].iterrows() %}
<tr>
{% for col in main_columns %}
<td>

{% set value = row[col] if col in row else "" %}
{% set num = value|float if value not in missing_value else None %}

{% if col in ["Legionella_pneumophila_percent", "Legionella_spp_percent"] and num is not none %}
    {% set display_value = "%.2f"|format(num * 100) %}
{% elif col in float_columns and num is not none %}
    {% set display_value = "%.2f"|format(num) %}
{% else %}
    {% set display_value = value %}
{% endif %}

{% if value in missing_value %}
    <span style="color:#e5e5e5;">{{ display_value }}</span>

{% elif value == "Contamination" %}
    <span class="highlight">{{ display_value }}</span>

{% elif value == "Potential new spp." %}
    <span class="highlight">{{ display_value }}</span>

{% elif col == "Depth" and num is not none and num < 30 %}
    <span class="highlight">{{ display_value }}</span>

{% elif col == "Number_of_contigs" and num is not none and num > 200 %}
    <span class="highlight">{{ display_value }}</span>

{% elif col == "Total_length" and num is not none and (num < 2500000 or num > 4500000) %}
    <span class="highlight">{{ display_value }}</span>

{% elif col == "GC" and num is not none and (num < 36 or num > 40) %}
    <span class="highlight">{{ display_value }}</span>

{% else %}
    {{ display_value }}
{% endif %}

</td>
{% endfor %}
</tr>
{% endfor %}

</table>
</div>
</div>

</div>

{% endfor %}

</div>

<hr>

<div class="section" id="AMR">

<h2>AMR Summary</h2>

<div class="table-wrapper">

<div class="table-controls">
<button class="copy-button" onclick="copyTable(this)">Copy</button>
</div>

<div class="table-scroll">
<table>

<tr>
{% for col in amr_columns %}
<th>{{ col }}</th>
{% endfor %}
</tr>

{% for _, row in amr_df[amr_columns].iterrows() %}

{% set impact = row["AMR_Nb_with_Impact"]|default("0") %}

<tr
{% if impact|float != 0 %}
class="amr-alert"
{% endif %}
>

{% for col in amr_columns %}
<td>{{ row[col] }}</td>
{% endfor %}

</tr>
{% endfor %}

</table>
</div>
</div>

</div>

{% if not other_df.empty %}

<hr>

<div class="section" id="OTHER">

<h2>Other Summary</h2>

<div class="table-wrapper">

<div class="table-controls">
<button class="copy-button" onclick="copyTable(this)">Copy</button>
</div>

<div class="table-scroll">
<table>

<tr>
{% for col in other_columns %}
<th>{{ col }}</th>
{% endfor %}
</tr>

{% for _, row in other_df.iterrows() %}

{% set impact = row["Other_Nb_with_Impact"]|default("0") %}

<tr
{% if impact|float != 0 %}
class="amr-alert"
{% endif %}
>

{% for col in other_columns %}
<td>{{ row[col] }}</td>
{% endfor %}

</tr>

{% endfor %}

</table>
</div>
</div>

</div>

{% endif %}

<hr>

<div class="section" id="Settings">

<h2>Settings - Traceability of software and parameters</h2>

<button class="show-button" onclick="toggleSettings(this)">Show</button>

<pre id="settings" class="hidden">{{ software_text }}</pre>

</div>

</div>

</body>
</html>
"""

    template = Template(template_str)

    html = template.render(
        analyse_id=args.analyse_id,
        sequencing_id=args.sequencing_id,
        date=datetime.now().strftime("%A %d %B %Y, %H:%M:%S"),
        df=df,
        amr_df=amr_df,
        other_df=other_df,
        strains=strains,
        main_columns=MAIN_COLUMNS,
        column_desc=COLUMN_DESC,
        amr_columns=AMR_COLUMNS,
        other_columns=OTHER_COLUMNS,
        float_columns=FLOAT_COLUMNS,
        missing_value=MISSING_VALUES,
        software_text=software_text,
        amr_warning=amr_warning(df)
    )

    Path(args.output).write_text(html, encoding="utf-8")


if __name__ == "__main__":
    main()
