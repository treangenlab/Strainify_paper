import pandas as pd
import plotly.graph_objects as go
import plotly.subplots as sp
import numpy as np
import re
from functools import reduce
import os

# pantax writes the string 'no_strains_detected' in its JSD column when it finds
# nothing. Those samples get no bar (height 0) and are marked with a star instead.
PANTAX_UNDETECTED_MARKER = "*"
PANTAX_UNDETECTED_COL = "__pantax_undetected__"

# === Input Configuration ===
species_list = {
    "C. acnes": [
        "/home/Users/rl152/Strainify/30_strains/jsd/Strainify/cacnes_not_weighted.csv",
        "/home/Users/rl152/Strainify/30_strains/jsd/Strainify/cacnes_weighted.csv",
        "/home/Users/rl152/Strainify/30_strains/jsd/StrainScan/cacnes_default.csv",
        "/home/Users/rl152/Strainify/30_strains/jsd/StrainScan/cacnes_low_coverage.csv",
        "/dodo/rl152/pantax/30_strains/jsd/cacnes_pantax.csv",
        "/home/Users/rl152/Strainify_dev/StrainGE_benchmarks/30_strain_simulated/jsds/cacnes.csv"
    ],
    "C. difficile": [
        "/home/Users/rl152/Strainify/30_strains/jsd/Strainify/cdiff_not_weighted.csv",
        "/home/Users/rl152/Strainify/30_strains/jsd/Strainify/cdiff_weighted.csv",
        "/home/Users/rl152/Strainify/30_strains/jsd/StrainScan/cdiff_default.csv",
        "/home/Users/rl152/Strainify/30_strains/jsd/StrainScan/cdiff_low_coverage.csv",
        "/dodo/rl152/pantax/30_strains/jsd/cdiff_pantax.csv",
        "/home/Users/rl152/Strainify_dev/StrainGE_benchmarks/30_strain_simulated/jsds/cdiff.csv"
    ],
    "M. tuberculosis": [
        "/home/Users/rl152/Strainify/30_strains/jsd/Strainify/mtuberculosis_not_weighted.csv",
        "/home/Users/rl152/Strainify/30_strains/jsd/Strainify/mtuberculosis_weighted.csv",
        "/home/Users/rl152/Strainify/30_strains/jsd/StrainScan/mtuberculosis_default.csv",
        "/home/Users/rl152/Strainify/30_strains/jsd/StrainScan/mtuberculosis_low_coverage.csv",
        "/dodo/rl152/pantax/30_strains/jsd/mtuberculosis_pantax.csv",
        "/home/Users/rl152/Strainify_dev/StrainGE_benchmarks/30_strain_simulated/jsds/mtuberculosis.csv"
    ],
    "E. coli": [
        "/home/Users/rl152/Strainify/30_strains/jsd/Strainify/ecoli_not_weighted.csv",
        "/home/Users/rl152/Strainify/30_strains/jsd/Strainify/ecoli_weighted.csv",
        "/home/Users/rl152/Strainify/30_strains/jsd/StrainScan/ecoli_default.csv",
        "/home/Users/rl152/Strainify/30_strains/jsd/StrainScan/ecoli_low_coverage.csv",
        "/dodo/rl152/pantax/30_strains/jsd/ecoli_pantax.csv",
        "/home/Users/rl152/Strainify_dev/StrainGE_benchmarks/30_strain_simulated/jsds/ecoli.csv"
    ],
    "S. epidermidis": [
        "/home/Users/rl152/Strainify/30_strains/jsd/Strainify/sepidermidis_not_weighted.csv",
        "/home/Users/rl152/Strainify/30_strains/jsd/Strainify/sepidermidis_weighted.csv",
        "/home/Users/rl152/Strainify/30_strains/jsd/StrainScan/sepidermidis_default.csv",
        "/home/Users/rl152/Strainify/30_strains/jsd/StrainScan/sepidermidis_low_coverage.csv",
        "/dodo/rl152/pantax/30_strains/jsd/sepidermidis_pantax.csv",
        "/home/Users/rl152/Strainify_dev/StrainGE_benchmarks/30_strain_simulated/jsds/sepidermidis.csv"
    ]
}

label_map = {
    "not_weighted": "Strainify (unweighted)",
    "weighted": "Strainify (weighted)",
    "default": "StrainScan (default)",
    "low_coverage": "StrainScan (low coverage)",
    "pantax": "PanTax",
    "strainge": "StrainGE"
}

ratio_name_map = {
    "ratio_1": "dominant",
    "ratio_2": "uniform",
    "ratio_3": "dirichlet"
}

colors = {
    "Strainify (unweighted)": "#66c2a5",
    "Strainify (weighted)": "#fc8d62",
    "StrainScan (default)": "#8da0cb",
    "StrainScan (low coverage)": "#e78ac3",
    "PanTax": "#a6d854",
    "StrainGE": "#ffd92f"
}

def extract_coverage(sample):
    match = re.match(r"(\d+)x", sample)
    return int(match.group(1)) if match else 0

def extract_ratio(sample):
    match = re.search(r"(ratio_\d+)", sample)
    return match.group(1) if match else sample


def extract_label_from_filename(path):
    filename = os.path.basename(path)
    for key, val in label_map.items():
        if key in filename:
            return val

    if "strainge" in path.lower():
        return "StrainGE"
    # Fail loudly if nothing matches
    raise ValueError(f"Filename {filename} did not match any known label")


def make_species_subplot(species_name, file_paths, showlegend=False):

    dfs = []
    for f in file_paths:
        df = pd.read_csv(f)
        label = extract_label_from_filename(f)
        if label is None:
            label = os.path.basename(f).split('.')[0]

        # JSD may be non-numeric (e.g. pantax writes 'no_strains_detected');
        # coerce to a number so the merge/plot doesn't choke on strings.
        df["JSD"] = pd.to_numeric(df["JSD"], errors="coerce")

        keep_cols = ["sample", label]
        if label == label_map["pantax"]:
            # Mark samples pantax couldn't resolve, then zero out the bar height
            # so they render as a star (added below) rather than a bar.
            df[PANTAX_UNDETECTED_COL] = df["JSD"].isna()
            df["JSD"] = df["JSD"].fillna(0.0)
            keep_cols.append(PANTAX_UNDETECTED_COL)

        df = df.rename(columns={"JSD": label})
        # Keep only the columns we need; pantax files carry extra columns
        # (truth_ratio, num_estimated_nonzero, ...) that would clutter the merge.
        df = df[keep_cols]
        dfs.append(df)

    df = reduce(lambda l, r: pd.merge(l, r, on="sample", how="outer"), dfs).fillna(0)
    df["coverage"] = df["sample"].apply(extract_coverage)
    df["ratio_label"] = df["sample"].apply(extract_ratio)
    df = df.sort_values(by=["coverage", "ratio_label"]).reset_index(drop=True)

    jsd_cols = [col for col in df.columns if col in colors]
    x = np.arange(len(df))
    x_labels = df["ratio_label"].map(ratio_name_map).fillna(df["ratio_label"])

    # Boolean per row: True where pantax detected no strains (gets a star, no bar).
    if PANTAX_UNDETECTED_COL in df.columns:
        undetected_mask = df[PANTAX_UNDETECTED_COL].fillna(False).astype(bool).values
    else:
        undetected_mask = np.zeros(len(df), dtype=bool)

    print(df[jsd_cols].head())
    print("Columns in merged df:", df.columns.tolist())
    print("Selected jsd_cols:", jsd_cols)

    bars = []
    for i, col in enumerate(jsd_cols):
        # For pantax, place a star above the (zero-height) bar wherever no strains
        # were detected. Text rides on the bar trace, so it inherits the group offset.
        if col == label_map["pantax"]:
            text_vals = np.where(undetected_mask, PANTAX_UNDETECTED_MARKER, "")
        else:
            text_vals = None
        bars.append(go.Bar(
            x=x,
            y=df[col],
            name=col,
            marker_color=colors[col],
            offsetgroup=i,
            text=text_vals,
            textposition="outside",
            textfont=dict(size=26, color="red"),
            cliponaxis=False,
            hovertext=df["sample"],
            hoverinfo="text+y",
            showlegend=showlegend
        ))

    shapes = []
    annotations = []
    coverage_groups = df["coverage"].ne(df["coverage"].shift()).cumsum()
    coverage_meta = df.groupby(coverage_groups).agg({"coverage": "first", "ratio_label": "count"}).reset_index(drop=True)

    start = 0
    for _, row in coverage_meta.iterrows():
        count = row["ratio_label"]
        coverage = row["coverage"]
        end = start + count - 1
        center = start + (count - 1) / 2

        shapes.append(dict(
            type="rect", x0=start - 0.5, x1=end + 0.5, y0=1.0, y1=1.15,
            fillcolor="lightgrey", line=dict(width=0), layer="below"
        ))
        annotations.append(dict(
            x=center, y=1.075, text=f"{coverage}x", showarrow=False,
            font=dict(size=15), align="center"
        ))
        if start != 0:
            shapes.append(dict(
                type="line", x0=start - 0.5, x1=start - 0.5, y0=0, y1=1.15,
                line=dict(color="gray", dash="dash")
            ))
        start += count

    return bars, shapes, annotations, x, x_labels


# --- Make subplots: 5 rows, 1 column ---
fig = sp.make_subplots(
    rows=5, cols=1,
    subplot_titles=[
        "C. acnes", "C. difficile", "M. tuberculosis",
        "E. coli", "S. epidermidis"
    ],
    vertical_spacing=0.1
)

# --- Species positions (stacked) ---
species_positions = {
    "C. acnes": (1, 1),
    "C. difficile": (2, 1),
    "M. tuberculosis": (3, 1),
    "E. coli": (4, 1),
    "S. epidermidis": (5, 1),
}

for i, (species, files) in enumerate(species_list.items()):
    row, col = species_positions[species]
    showlegend = (i == 0)  # show legend only for first plot
    bars, shapes, annotations, x_vals, x_labels = make_species_subplot(species, files, showlegend=showlegend)

    # Add bars
    for bar in bars:
        fig.add_trace(bar, row=row, col=col)

    # Axes
    fig.update_xaxes(
        tickmode='array',
        tickvals=x_vals,
        ticktext=x_labels,
        tickangle=90,
        #title_text="Abundance Ratio",
        title_text="Abundance Ratio" if row == len(species_list) else None,
        title_font=dict(family="Arial", size=28, color="black"), 
        row=row, col=col,
        tickfont=dict(family="Arial", size=23, color="black"),  
    )
    fig.update_yaxes(
        title_text="JSD", range=[0, 1.15],
        title_font=dict(family="Arial", size=28, color="black"),
        row=row, col=col,
        tickfont=dict(family="Arial", size=23, color="black"),
    )

    # Shapes
    for s in shapes:
        fig.add_shape(s, row=row, col=col)

    # Correct annotation xref/yref for stacked layout
    subplot_index = row  # because each row is one subplot
    xref = f"x{subplot_index}" if subplot_index > 1 else "x"
    yref = f"y{subplot_index}" if subplot_index > 1 else "y"

    for a in annotations:
        fig.add_annotation(dict(**a, xref=xref, yref=yref))

# Adjust subplot title font size
for annotation in fig.layout.annotations:
    annotation.font.size = 32
    annotation.font.family = "Arial"
    annotation.font.color = "black"

num_species = len(species_list)  # 5 in your example

# Bold only the subplot titles (first 5 annotations)
for annotation in fig.layout.annotations[:num_species]:
    annotation.text = f"<b><i>{annotation.text}</i></b>"
    annotation.font.size = 32
    annotation.font.family = "Arial"
    annotation.font.color = "black"



fig.update_layout(
    font=dict(
        family="Arial",  # or "Helvetica", etc.
        size=20,  # default font size
        color="black"
    ),
    height=2000, width=1500,
    barmode='group',
    legend=dict(
        orientation="h", yanchor="top", y=-0.1,
        xanchor="center", x=0.5, font=dict(size=30)
    ),
    plot_bgcolor="white", paper_bgcolor="white",
    margin=dict(t=100, b=190)
)


#fig.write_image("/home/Users/rl152/Strainify/plot_scripts/jsd_all_species_grid_one_column_v5.png", width=1500, height=1950, scale=3)
fig.write_image("/home/Users/rl152/Strainify/plot_scripts/jsd_all_species_grid_one_column_v5.pdf", width=1500, height=1950, scale=3)