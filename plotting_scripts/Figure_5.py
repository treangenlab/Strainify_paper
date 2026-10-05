import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# === Load CSVs ===
df1 = pd.read_csv("/home/Users/rl152/Strainify/plot_scripts/mock_community_comparison_v2.csv")
df2 = pd.read_csv("/home/Users/rl152/Strainify/plot_scripts/mock_community_extra_comparison_v2.csv")

# === Rename Strainify columns ===
df1 = df1.rename(columns={
    "Strainify (plasmid included)": "Strainify (a)",
    "Strainify (plasmid excluded entirely)": "Strainify (b)",
    "Strainify (plasmid included in genome alignment and read mapping)": "Strainify (c)"
})

# === Extract strains and methods ===
strains = df1["strain"].values

# === Methods in first CSV (excluding Ground truth) ===
methods_main = ["Ground truth"] + [col for col in df1.columns if col not in ["strain", "Ground truth"]]

# === All methods including StrainGE/EST ===
all_methods = methods_main + ["StrainGE", "StrainEst", "PHLAME", "PanTax"]


# === Desired strain orders ===
strain_order_main = ["H10407", "UTI89", "Sakai", "E24377A"]
strain_order_ge   = ["H10407", "UTI89", "ECP17-1298", "E24377A"]
strain_order_est  = ["RM14721", "UTI89", "109", "E24377A", "Other"]
phlame_strain_order = ["H10407", "UTI89", "Sakai", "E24377A", "Other"]
pantax_strain_order = ["H10407", "UTI89", "Sakai", "E24377A"]

# === Prepare data dictionary: method -> list of (strain, value) ===
method_data = {}

# First CSV methods (main)
for method in methods_main:
    vals = []
    for strain in strain_order_main:  # use desired order
        if strain in df1["strain"].values:
            val = df1.loc[df1["strain"]==strain, method].values[0]
            vals.append((strain, val))
    method_data[method] = vals

# StrainGE
df_ge = df2[["StrainGE_strains","StrainGE"]].dropna().rename(columns={"StrainGE_strains":"strain"})
vals = []
for strain in strain_order_ge:
    row = df_ge[df_ge["strain"]==strain]
    if not row.empty:
        vals.append((strain, row["StrainGE"].values[0]))
method_data["StrainGE"] = vals

# StrainEST
df_est = df2[["StrainEST_strains","StrainEst"]].dropna().rename(columns={"StrainEST_strains":"strain"})
vals = []
for strain in strain_order_est:
    row = df_est[df_est["strain"]==strain]
    if not row.empty:
        vals.append((strain, row["StrainEst"].values[0]))
method_data["StrainEst"] = vals

# PHLAME
df_phlame = df2[["PHLAME_strains","PHLAME"]].dropna().rename(columns={"PHLAME_strains":"strain"})
vals = []
for strain in phlame_strain_order:
    row = df_phlame[df_phlame["strain"]==strain]
    if not row.empty:
        vals.append((strain, row["PHLAME"].values[0]))
method_data["PHLAME"] = vals

# PanTax
df_pantax = pd.read_csv(
    "/dodo/rl152/pantax/4_strain_mock/pantax_results/SRR13355226_interleaved_strains_abundance.txt",
    sep="\t"
)
# Normalize predicted_abundance to sum to 100% (same scale as the other tools)
df_pantax["abundance_pct"] = (
    df_pantax["predicted_abundance"] / df_pantax["predicted_abundance"].sum() * 100
)
vals = []
for strain in pantax_strain_order:
    row = df_pantax[df_pantax["genome_ID"] == strain]
    if not row.empty:
        vals.append((strain, row["abundance_pct"].values[0]))
method_data["PanTax"] = vals


all_strains = (
    strain_order_main
    + [s for s in strain_order_ge if s not in strain_order_main]
    + [s for s in strain_order_est if s not in strain_order_main + strain_order_ge]
)

color_palette = [
    "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
    "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3",
    "#1b9e77", "#7570b3", "#d95f02", "#a6761d"
]

colors = {strain: color_palette[i % len(color_palette)] for i, strain in enumerate(all_strains)}

# === Plot ===
bar_width = 0.55
spacing = 4.1  # increase spacing between methods
x_base = np.arange(len(all_methods)) * spacing

#x_base = np.arange(len(all_methods))  # center positions for each method

fig, ax = plt.subplots(figsize=(28,10))

# Track legend handles to avoid duplicates
legend_handles = {}

for i, method in enumerate(all_methods):
    vals = method_data.get(method, [])
    n = len(vals)
    if n == 0:
        continue
    # calculate bar positions for this method
    #offsets = np.linspace(-bar_width*(n-1)/2, bar_width*(n-1)/2, n)
    inner_spacing = bar_width * 1.35  # >1 adds space between bars
    offsets = np.linspace(-inner_spacing*(n-1)/2, inner_spacing*(n-1)/2, n)
    for j, (strain, value) in enumerate(vals):
        pos = x_base[i] + offsets[j]
        bar = ax.bar(pos, value, width=bar_width, color=colors[strain])
        ax.text(pos, value+1, f"{value:.1f}", ha="center", va="bottom", fontsize=13)
        if strain not in legend_handles:
            legend_handles[strain] = bar

# X-axis
ax.set_xticks(x_base)
ax.set_xticklabels(all_methods, fontsize=26, rotation=0)
ax.set_ylabel("Abundance Estimate (%)", fontsize=30)
ax.tick_params(axis='y', labelsize=28)
ax.set_ylim(0, 85)

# Legend
ax.legend(legend_handles.values(), legend_handles.keys(), bbox_to_anchor=(0.5,-0.1), loc="upper center", ncol=4, fontsize=28)

plt.tight_layout()
#plt.savefig("/home/Users/rl152/Strainify/plot_scripts/mock_community_comparison_v10_w_ground_truth.png", dpi=600, bbox_inches="tight")
plt.savefig("/home/Users/rl152/Strainify/plot_scripts/mock_community_comparison_v10_w_ground_truth.pdf", dpi=600, bbox_inches="tight")