import pandas as pd

df_abundance = pd.read_csv("/home/Users/rl152/Strainify/UMB/UMB18_results_weighted/abundance_estimates_combined.csv")

df_mapping = pd.read_csv("/home/Users/rl152/Strainify/UMB/acc_list_unique.txt", delim_whitespace=True, header=None, names=["nucleotide_accession", "assembly_accession"])

df_sample_name_mapping = pd.read_csv("/home/Users/rl152/Strainify/UMB/UMB18_stool_SRRs.csv")

df_strain_info_mapping = pd.read_csv("/home/Users/rl152/Strainify/UMB/Source_Data_Figure_3_UMB18_strainGST.csv")

df_strain_info_mapping["T"] = df_strain_info_mapping["T"] - 405


################ plotting #####################
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ---- inputs ----
abundance_csv = "/home/Users/rl152/Strainify/UMB/UMB18_results_weighted/abundance_estimates_combined_mapped.csv"

out_png = "/home/Users/rl152/Strainify/UMB/UMB18_results_weighted/abundance_over_time_plot.pdf"

# ---- load abundance ----
df_abundance = pd.read_csv(abundance_csv)

strain_col = "strain name"
day_cols = [c for c in df_abundance.columns if c != strain_col]

def to_num(x):
    try:
        return float(str(x).strip())
    except Exception:
        return None

# --- sort day columns numerically if possible ---
day_nums = [to_num(c) for c in day_cols]
if all(v is not None for v in day_nums):
    order = np.argsort(day_nums)
    day_cols_sorted = [day_cols[i] for i in order]
    x_vals = np.array([day_nums[i] for i in order], dtype=float)
    x_label_vals = [str(int(v)) if float(v).is_integer() else str(v) for v in x_vals]
else:
    day_cols_sorted = day_cols
    x_vals = np.arange(len(day_cols_sorted))
    x_label_vals = [str(c) for c in day_cols_sorted]

# ---- merge phylogroup onto abundance table ----
df_strain_info_mapping.columns = df_strain_info_mapping.columns.str.strip()

required = {"ClusterName", "Phylogroup"}
missing = required - set(df_strain_info_mapping.columns)
if missing:
    raise ValueError(f"df_strain_info_mapping missing required columns: {sorted(missing)}")

df_abundance[strain_col] = df_abundance[strain_col].astype(str).str.strip()
df_strain_info_mapping["ClusterName"] = df_strain_info_mapping["ClusterName"].astype(str).str.strip()
df_strain_info_mapping["Phylogroup"] = df_strain_info_mapping["Phylogroup"].astype(str).str.strip()

# Ensure mapping is one-to-one: one row per ClusterName
map_df = (
    df_strain_info_mapping[["ClusterName", "Phylogroup"]]
    .dropna(subset=["ClusterName"])
    .drop_duplicates(subset=["ClusterName"], keep="first")
)

df = df_abundance.merge(
    map_df,
    left_on=strain_col,
    right_on="ClusterName",
    how="left",
)
df["Phylogroup"] = df["Phylogroup"].fillna("Unknown")

# ---- plotting setup ----
phylo_colors = {
    "A": "green",
    "B1": "lightskyblue",
    "B2": "red",
    "D": "#FFD700",
}
linestyles = ["-", "--", ":"]  # only 3 styles, cycled within each phylogroup
phylo_style_idx = {}           # linestyle counter per phylogroup

epsilon = 1e-6  # for log scale: replace 0s / negatives

fig, ax = plt.subplots(figsize=(25, 14))

# Optional: make ordering stable (so same strain always gets same style across runs)

df = df.sort_values(["Phylogroup", strain_col], kind="stable").reset_index(drop=True)

for _, row in df.iterrows():
    strain = row[strain_col]
    pg = row["Phylogroup"]
    color = phylo_colors.get(pg, "gray")

    if pg not in phylo_style_idx:
        phylo_style_idx[pg] = 0

    idx = phylo_style_idx[pg]
    ls = linestyles[idx % len(linestyles)]
    phylo_style_idx[pg] += 1

    y = pd.to_numeric(row[day_cols_sorted], errors="coerce").to_numpy(dtype=float)


    ax.plot(
        x_vals,
        y,
        color=color,
        linestyle=ls,
        linewidth=5,
        alpha=0.95,
        label=strain
    )
print(df["Phylogroup"].value_counts())

# ---- add vertical event lines ----
uti_days = [0, 41, 231, 293]
nitrofuran_days = [6, 299]
betalactam_days = [234, 278]
unknown_abx_days = [48]

# Draw lines (only label first in each group to avoid legend spam)
for j, d in enumerate(uti_days):
    ax.axvline(d, color="black", linestyle="-", linewidth=4, alpha=0.9,
               label="UTI" if j == 0 else "_nolegend_")

for j, d in enumerate(nitrofuran_days):
    ax.axvline(d, color="brown", linestyle="--", linewidth=4, alpha=0.9,
               label="Nitrofuran" if j == 0 else "_nolegend_")

for j, d in enumerate(betalactam_days):
    ax.axvline(d, color="purple", linestyle="--", linewidth=4, alpha=0.9,
               label="Beta lactam" if j == 0 else "_nolegend_")

for j, d in enumerate(unknown_abx_days):
    ax.axvline(d, color="gray", linestyle="--", linewidth=4, alpha=0.9,
               label="Unknown antibiotic" if j == 0 else "_nolegend_")

event_legend_elems = [
    Line2D([0], [0], color="black", lw=5, linestyle="-", label="UTI"),
    Line2D([0], [0], color="brown",   lw=5, linestyle="--", label="Nitrofuran"),
    Line2D([0], [0], color="purple", lw=5, linestyle="--", label="Beta lactam"),
    Line2D([0], [0], color="gray",  lw=5, linestyle="--", label="Unknown antibiotic"),
]


ax.set_xlabel("Sample time (days)", fontsize=30, labelpad=13, fontweight="bold")
ax.set_ylabel("Relative abundance (%)", fontsize=30, fontweight="bold")
#ax.set_title("Strain abundances over time")
ax.set_xticks(x_vals)
ax.set_xticklabels(x_label_vals, rotation=45, ha="right", fontsize=25)
ax.tick_params(axis="y", labelsize=25)
ax.set_yscale("log")
# Major grid (stronger)
ax.grid(True, which="major", linestyle="-", linewidth=1.5, alpha=0.5)

# Minor grid (lighter, especially useful for log scale)
ax.grid(True, which="minor", linestyle=":", linewidth=1.0, alpha=0.3)


# ---- phylogroup legend (colors) ----
phylo_present = [k for k in ["A", "B1", "B2", "D"] if k in set(df["Phylogroup"])]
phylo_legend_elems = [
    Line2D([0], [0], color=phylo_colors[k], lw=4, linestyle="-", label=k)
    for k in phylo_present
]
leg1 = ax.legend(
    handles=phylo_legend_elems,
    title="Phylogroup",
    title_fontsize=28,
    loc="upper center",
    bbox_to_anchor=(0.15, -0.14),
    #ncol=max(1, min(5, len(phylo_legend_elems))),
    ncol=2,
    frameon=True,
    fontsize=26,
)


# ---- strain legend (each line/strain) ----
handles, labels = ax.get_legend_handles_labels()

event_labels = {"UTI", "Nitrofuran", "Beta lactam", "Unknown antibiotic", "_nolegend_"}
uniq = {}
for h, lab in zip(handles, labels):
    if lab in event_labels:
        continue
    if lab not in uniq:
        uniq[lab] = h

leg2 = ax.legend(
    handles=list(uniq.values()),
    labels=list(uniq.keys()),
    title="Strain",
    title_fontsize=28,
    loc="upper center",
    bbox_to_anchor=(0.85, -0.14),
    #ncol=max(1, min(4, len(uniq))),
    ncol = 2,
    frameon=True,
    fontsize=26,
)

leg3 = ax.legend(
    handles=event_legend_elems,
    title="Event",
    title_fontsize=28,
    loc="upper center",
    bbox_to_anchor=(0.48, -0.14),
    ncol=2,
    frameon=True,
    fontsize=26,
)

ax.add_artist(leg1)
ax.add_artist(leg3)
ax.add_artist(leg2)

#plt.tight_layout(rect=[0, 0.17, 1, 1])
fig.subplots_adjust(bottom=0.32, top = 0.98, left=0.06, right=0.98)
plt.savefig(out_png, dpi=600)
plt.close(fig)

print(f"Saved: {out_png}")

