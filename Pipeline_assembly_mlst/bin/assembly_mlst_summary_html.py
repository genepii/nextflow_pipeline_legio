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
    "AMR_Nb_Mutated_Genes": "Number of mutated AMR genes detected",
    "AMR_Nb_with_Impact": "Number of AMR mutations with predicted functional impact",
    "AMR_Nb_Non_Coding": "Number of non-coding AMR-related variants"
}


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
                [
                    "nan",
                    "NaN",
                    "NA",
                    "Na",
                    "N/A",
                    "<NA>",
                    "null",
                    "NULL",
                    ""
                ],
                pd.NA
            )
            .isna()
            .all(axis=1)
        ]

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
body{
    font-family: system-ui,-apple-system,"Segoe UI",Roboto,Arial,sans-serif;
    margin: 0;
}

h2{
    text-decoration: underline;
}

.page{
    padding: 48px;
    padding-top: 4px; 
}

.titlebar{
    background: black;
    color: white;
    padding: 15px;
}

.report-title{
    font-size: 22px;
    font-weight: bold;
}

.report-subtitle{
    font-size: 13px;
    font-weight: bold;
    margin-top: 4px;
    margin-bottom: 4px;
    padding-left: calc(48px - 15px);
}

.report-subtitle a{
    color: white;
    text-decoration: underline;
}

.info-block{
    margin-top: 4px;
    margin-bottom: 8px;
    font-size: 13px;
    line-height: 1.5;
}

.info-block a{
    font-weight: bold;
    text-decoration: underline;
}

.section{
    margin-top: 40px;
}

.table-wrapper{
    width: calc(100% - 174px);
    margin-left: 87px;
    margin-right: 87px;
    overflow-x: auto;
    overflow-y: auto;
    max-height: 600px;
}

table{
    border-collapse: collapse;
    min-width: max-content;
}

th,td{
    border:1px solid #ccc;
    padding: 6px;
    font-size: 13px;
    white-space: nowrap;
}

.table-description{
    width:calc(100% - 174px);
    margin:10px 48px 0 48px;
    font-size:13px;
    color:#555;
    line-height:1.45;
}

.table-description ul{
    margin:4px 0 0 18px;
    padding:0;
}

.table-description li{
    margin:2px 0;
}

th{
    background: #4d4d4d;
    color: white;
}

.amr-alert{
    background: #ffdddd;
}

hr{
    border:none;
    border-top: 1px solid #e5e5e5;
    margin: 40px 0;
}

.hidden{
    display:none;
}
</style>

<script>
function toggleSettings(){
    document
        .getElementById("settings")
        .classList
        .toggle("hidden");
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
<div style="background:red;color:white;padding:15px;font-weight:bold;">
AMR variants detected : <a href="#AMR" style="color:white;">see table</a>
</div>
<hr>
{% endif %}

<h2>Table of Contents</h2>

<ol>
<li><a href="#Quality"><b>Quality indicators</b></a></li>
<li><a href="#Desc"><b>Column descriptions</b></a></li>

<li>
Assembly &amp; MLST Summary
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
<li><a href="#Settings"><b>Settings</b></a></li>
</ol>

<hr>

<h2>Quality indicators</h2>

<div id="Quality" class="table-description">

Values highlighted in <span style="color:#c00000; font-weight:bold;">red and bold</span> indicate potential quality issues or results requiring attention.

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

<h2>Assembly &amp; MLST Summary</h2>

{% for s in strains %}
<div class="section" id="{{ s }}">

<h3><i>{{ s.replace("_", " ") }}</i></h3>

<div class="table-wrapper">
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
{% set num = value|float if value not in ["", "nan", "NaN", None, "NA", "Na", "Null"] else None %}
{% set float_columns = ["Depth", "GC", "FastANI_value", "auN"] %}

{% if col in ["Legionella_pneumophila_percent", "Legionella_spp_percent"] and num is not none %}
    {% set display_value = "%.2f"|format(num * 100) %}
{% elif col in float_columns and num is not none %}
    {% set display_value = "%.2f"|format(num) %}
{% else %}
    {% set display_value = value %}
{% endif %}

{% if value in ["", "nan", "NaN", None, "NA", "Na", "Null"] %}
    <span style="color:#e5e5e5;">{{ display_value }}</span>

{% elif value == "Contamination" %}
    <span style="color:red;font-weight:bold;">{{ display_value }}</span>

{% elif value == "Potential new spp." %}
    <span style="color:red;font-weight:bold;">{{ display_value }}</span>

{% elif col == "Depth" and num is not none and num < 30 %}
    <span style="color:red;font-weight:bold;">{{ display_value }}</span>

{% elif col == "Number_of_contigs" and num is not none and num > 200 %}
    <span style="color:red;font-weight:bold;">{{ display_value }}</span>

{% elif col == "Total_length" and num is not none and (num < 2500000 or num > 4500000) %}
    <span style="color:red;font-weight:bold;">{{ display_value }}</span>

{% elif col == "GC" and num is not none and (num < 36 or num > 40) %}
    <span style="color:red;font-weight:bold;">{{ display_value }}</span>

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

{% endfor %}

<hr>

<div class="section" id="AMR">

<h2>AMR Summary</h2>

<div class="table-wrapper">
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

<hr>

<div class="section" id="Settings">

<h2>Settings - Traceability of software and parameters</h2>

<button onclick="toggleSettings()">Show</button>

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
        strains=strains,
        main_columns=MAIN_COLUMNS,
        column_desc=COLUMN_DESC,
        amr_columns=AMR_COLUMNS,
        software_text=software_text,
        amr_warning=amr_warning(df)
    )

    Path(args.output).write_text(html, encoding="utf-8")


if __name__ == "__main__":
    main()
