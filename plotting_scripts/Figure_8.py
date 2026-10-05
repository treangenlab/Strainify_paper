import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from matplotlib.legend_handler import HandlerTuple

# === Load CSVs ===
df1 = pd.read_csv(
    "/home/Users/rl152/Strainify_dev/Strainify_benchmarks/4_strain_mock_benchmark/benchmarks/pipeline_runtime_summary.csv",
    index_col=0
)
df2 = pd.read_csv(
    "/home/Users/rl152/Strainify_dev/StrainScan_benchmark/logs/4_strain_mock/combined_default.csv",
    index_col=0
)
df3 = pd.read_csv(
    "/home/Users/rl152/Strainify_dev/StrainR2_benchmarks/combined_metrics.csv",
    index_col=0
)
df4 = pd.read_csv(
    "/home/Users/rl152/Strainify/phlame/4_strain_mock_community/benchmarks/runtime_ram_summary.csv",
    index_col=0
)

# === Normalize Strainify column names ===
df1 = df1.T
df1.columns = [col.replace("build_variant_matrix", "build") for col in df1.columns]

idx = list(df1.index)
idx[0] = "Wall clock time (s)"
idx[1] = "Peak RAM (MB)"
idx[7] = "CPU time (s)"
df1.index = idx

# === Standardize RAM metric names ===
df2 = df2.rename(index={"Max RSS (MB)": "Peak RAM (MB)"})
df3 = df3.rename(index={"Max RSS (MB)": "Peak RAM (MB)"})
df4 = df4.rename(index={"Max RSS (MB)": "Peak RAM (MB)"})

metrics_to_keep = ["Wall clock time (s)", "Peak RAM (MB)"]

# === Fix orientation if needed ===
if set(metrics_to_keep).issubset(df2.columns):
    df2 = df2.T
if set(metrics_to_keep).issubset(df3.columns):
    df3 = df3.T
if set(metrics_to_keep).issubset(df4.columns):
    df4 = df4.T

# === Keep build + sample columns ===
df1 = df1.loc[metrics_to_keep, ["build", "SRR13355226"]]
df2 = df2.loc[metrics_to_keep, ["build", "SRR13355226"]]
df3 = df3.loc[metrics_to_keep, ["build", "SRR13355226"]]
df4 = df4.loc[metrics_to_keep, ["build", "SRR13355226"]]

# === Label columns by tool ===
df1.columns = [f"Strainify_{c}" for c in df1.columns]
df2.columns = [f"StrainScan_{c}" for c in df2.columns]
df3.columns = [f"StrainR2_{c}" for c in df3.columns]
df4.columns = [f"PHLAME_{c}" for c in df4.columns]

# === Merge all ===
merged = pd.concat([df1, df2, df3, df4], axis=1)

df_build = merged.filter(regex="_build$")
df_sample = merged.filter(regex="_SRR13355226$")

df_build.columns = ["Strainify", "StrainScan", "StrainR2", "PHLAME"]
df_sample.columns = ["Strainify", "StrainScan", "StrainR2", "PHLAME"]

# === PanTax: parse GNU /usr/bin/time logs ===
# PanTax doesn't have a summary CSV like the other tools, so we read the two
# `/usr/bin/time -v` logs directly (build database + compute abundances).
PANTAX_LOG_DIR = "/dodo/rl152/pantax/4_strain_mock/pantax_time_logs"
PANTAX_BUILD_LOG = os.path.join(PANTAX_LOG_DIR, "build_db.time.log")
PANTAX_SAMPLE_LOG = os.path.join(PANTAX_LOG_DIR, "SRR13355226_interleaved.time.log")


def parse_gnu_time_log(path):
    """Return (wall_clock_seconds, peak_ram_MB) from a `/usr/bin/time -v` log."""
    wall_seconds = None
    peak_ram_mb = None
    with open(path) as f:
        for line in f:
            if "Elapsed (wall clock) time" in line:
                # e.g. "6:38.92" (m:ss) or "1:06:38.92" (h:mm:ss)
                parts = [float(p) for p in line.split()[-1].split(":")]
                if len(parts) == 3:
                    wall_seconds = parts[0] * 3600 + parts[1] * 60 + parts[2]
                elif len(parts) == 2:
                    wall_seconds = parts[0] * 60 + parts[1]
                else:
                    wall_seconds = parts[0]
            elif "Maximum resident set size" in line:
                # reported in kbytes -> MB (divide by 1024, same MiB convention
                # used for the other tools' "Max RSS (MB)")
                peak_ram_mb = float(line.split()[-1]) / 1024.0
    return wall_seconds, peak_ram_mb


pantax_build_wall, pantax_build_ram = parse_gnu_time_log(PANTAX_BUILD_LOG)
pantax_sample_wall, pantax_sample_ram = parse_gnu_time_log(PANTAX_SAMPLE_LOG)

df_build.loc["Wall clock time (s)", "PanTax"] = pantax_build_wall
df_build.loc["Peak RAM (MB)", "PanTax"] = pantax_build_ram
df_sample.loc["Wall clock time (s)", "PanTax"] = pantax_sample_wall
df_sample.loc["Peak RAM (MB)", "PanTax"] = pantax_sample_ram

# === StrainGE: parse timing-summary CSVs ===
# Each CSV has rows `wall_clock_time_sec` and `peak_ram_gb`, with a single
# data column (build keyed by "mock", compute by "SRR13355226").
STRAINGE_DIR = "/home/Users/rl152/Strainify_dev/StrainGE_benchmarks/timing_summary"
STRAINGE_BUILD_CSV = os.path.join(STRAINGE_DIR, "4_strain_mock_strainge_build_timing.csv")
STRAINGE_SAMPLE_CSV = os.path.join(STRAINGE_DIR, "4_strain_mock__mock_strainge_timing.csv")


def parse_strainge_timing(path):
    """Return (wall_clock_seconds, peak_ram_MB) from a StrainGE timing CSV."""
    s = pd.read_csv(path, index_col=0).iloc[:, 0]          # first (only) data column
    wall_seconds = float(s["wall_clock_time_sec"])
    peak_ram_mb = float(s["peak_ram_gb"]) * 1024.0         # GiB -> MiB (matches others)
    return wall_seconds, peak_ram_mb


strainge_build_wall, strainge_build_ram = parse_strainge_timing(STRAINGE_BUILD_CSV)
strainge_sample_wall, strainge_sample_ram = parse_strainge_timing(STRAINGE_SAMPLE_CSV)

df_build.loc["Wall clock time (s)", "StrainGE"] = strainge_build_wall
df_build.loc["Peak RAM (MB)", "StrainGE"] = strainge_build_ram
df_sample.loc["Wall clock time (s)", "StrainGE"] = strainge_sample_wall
df_sample.loc["Peak RAM (MB)", "StrainGE"] = strainge_sample_ram

tool_colors = {
    "Strainify": ("#1f77b4", "#aec7e8"),
    "StrainScan": ("#ff7f0e", "#ffbb78"),
    "StrainR2": ("#2ca02c", "#98df8a"),
    "PHLAME": ("#d62728", "#ff9896"),
    "PanTax": ("#9467bd", "#c5b0d5"),
    "StrainGE": ("#8c564b", "#c49c94"),
}

# Tool sets per panel
time_tools = ["Strainify", "StrainScan", "StrainR2", "PanTax", "StrainGE"]
ram_tools  = ["Strainify", "StrainScan", "StrainR2", "PHLAME", "PanTax", "StrainGE"]

# === Subplots ===
fig, axes = plt.subplots(1, 2, figsize=(22, 6))

panel_specs = [
    ("Wall clock time (s)", axes[0], time_tools),
    ("Peak RAM (MB)", axes[1], ram_tools),
]

width = 0.34

for metric, ax, tools in panel_specs:
    x = np.arange(len(tools))

    build_vals = df_build.loc[metric, tools].astype(float).values
    sample_vals = df_sample.loc[metric, tools].astype(float).values

    for i, tool in enumerate(tools):
        dark, light = tool_colors[tool]

        ax.bar(
            x[i] - width / 2,
            build_vals[i],
            width=width,
            color=dark,
            zorder=3
        )

        ax.bar(
            x[i] + width / 2,
            sample_vals[i],
            width=width,
            color=light,
            zorder=3
        )

    max_val = max(build_vals.max(), sample_vals.max())

    ax.set_ylim(0, max_val * 1.18)
    ax.set_xticks(x)
    ax.set_xticklabels(tools, fontsize=20)
    ax.tick_params(axis="y", labelsize=20)
    ax.grid(axis="y", linestyle="--", alpha=0.6, zorder=0)
    ax.set_title(metric.replace(" (s)", "").replace(" (MB)", ""), fontsize=23, pad=12)

axes[0].set_ylabel("Second", fontsize=22)
axes[1].set_ylabel("MB", fontsize=22)

# === Legend: use the RAM panel's tool set (includes all tools) ===
legend_tools = ram_tools
build_triplet = tuple(Patch(facecolor=tool_colors[t][0]) for t in legend_tools)
sample_triplet = tuple(Patch(facecolor=tool_colors[t][1]) for t in legend_tools)

handles = [build_triplet, sample_triplet]
labels = ["Build database / variant matrix", "Compute abundances"]

fig.legend(
    handles,
    labels,
    handler_map={tuple: HandlerTuple(ndivide=None)},
    loc="lower center",
    bbox_to_anchor=(0.5, -0.02),
    ncol=2,
    frameon=False,
    fontsize=22,
    handlelength=4.8
)

plt.tight_layout(rect=[0, 0.1, 1, 1])
plt.savefig(
    "/home/Users/rl152/Strainify_dev/benchmark_plots/four_strain_community_separate_6_tools.pdf",
    dpi=600,
    bbox_inches="tight"
)