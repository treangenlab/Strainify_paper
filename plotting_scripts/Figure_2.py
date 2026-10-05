import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator
import math
import warnings


warnings.filterwarnings('ignore', message='This figure includes Axes that are not compatible with tight_layout')


OUT_PDF = "/home/Users/rl152/Strainify/4_strain_ecoli_simulated/plots/4_strain_simulated_combined.pdf"
TRUTH_CSV = "/home/Users/rl152/Strainify/4_strain_ecoli_simulated/4-strain_simulated_ratios.csv"

ratio_labels = {
    'ratio_1': 'log',
    'ratio_2': 'uniform',
    'ratio_3': 'asym 1',
    'ratio_4': 'asym 2'
}

colors = {
    'Strainify': '#1f77b4',
    'StrainScan': '#ff7f0e',
    'PHLAME': '#2ca02c',
    'StrainScan (low cov)': '#d62728',
    'PHLAME (low cov)': '#9467bd',
    'PanTax': '#8c564b',
    'StrainGE': '#17becf'
}


# ---- Shared helpers (unchanged from the originals) ----
def renormalize_columns_to_100(df):
    """Scale each column so its values sum to 100.

    - NaN entries (e.g. PanTax 'undetected') are ignored in the sum and kept
      as NaN, so downstream asterisk handling is unaffected.
    - Columns that sum to 0 are left unchanged (stay all-zero, not NaN).
    """
    out = df.copy()
    col_sums = out.sum(axis=0, skipna=True)
    nonzero = col_sums[col_sums != 0].index
    out[nonzero] = out[nonzero].divide(col_sums[nonzero], axis=1) * 100
    return out


def group_columns_by_coverage(df):
    groups = {}
    for col in df.columns:
        match = re.match(r'(\d+x)_ratio_(\d+)', col)
        if match:
            cov, ratio = match.groups()
            groups.setdefault(cov, []).append((int(ratio), col))
    for cov in groups:
        groups[cov] = [col for _, col in sorted(groups[cov])]
    return groups


def load_pantax(path):
    raw = pd.read_csv(path, index_col=0)
    undetected = raw.astype(str).apply(lambda col: col.str.lower().eq("undetected"))
    df = raw.replace("undetected", np.nan).apply(pd.to_numeric, errors="coerce")
    return df, undetected


def load_truth():
    df = pd.read_csv(TRUTH_CSV, index_col=0)
    df.index = df.index.str.replace(".fna", "", regex=False)
    return df * 100



def load_high_cov():
    STRAINGE_CSV = "/home/Users/rl152/Strainify_dev/StrainGE_benchmarks/combined_rapct/4_strain_simulated__ecoli.csv"

    t1 = pd.read_csv("/home/Users/rl152/Strainify/4_strain_ecoli_simulated/results_500_0/abundance_estimates_combined_final.csv", index_col=0)
    t2 = pd.read_csv("/home/Users/rl152/Strainify/4_strain_ecoli_simulated/StrainScan/combined_abundance.csv", index_col=0)
    t3 = pd.read_csv("/home/Users/rl152/Strainify/phlame/4_strain_ecoli_simulated/results_customed_clades/combined_relative_abundance_with_novel.csv", index_col=0)
    t4 = pd.read_csv("/home/Users/rl152/Strainify/4_strain_ecoli_simulated/StrainScan/low_coverage/combined_abundance.csv", index_col=0)
    t5 = pd.read_csv("/home/Users/rl152/Strainify/phlame/4_strain_ecoli_simulated/results_customed_clades_2/combined_relative_abundance_with_novel.csv", index_col=0)
    t6, t6_und = load_pantax("/dodo/rl152/pantax/4_strain_simulated/pantax_results/pantax_4strain_relative_abundances.csv")
    t7 = pd.read_csv(STRAINGE_CSV, index_col=0)
    truth = load_truth()

    t1 = t1.rename(columns={'strain name': 'Strain_Name'})
    t2.columns = [c.split('.txt')[0] for c in t2.columns]
    t4.columns = [c.split('.txt')[0] for c in t4.columns]

    # Remove novel strain row from PHLAME tables
    t3 = t3[~t3.index.str.contains("novel", case=False, na=False)]
    t5 = t5[~t5.index.str.contains("novel", case=False, na=False)]

    # Convert to percent (Strainify already %)
    t2, t3, t4, t5, t6, t7 = [d * 100 for d in (t2, t3, t4, t5, t6, t7)]

    common = sorted(set(t1.index) & set(t2.index) & set(t3.index) & set(t4.index)
                    & set(t5.index) & set(t6.index) & set(t7.index) & set(truth.index))
    if not common:
        raise ValueError("No common strains found (high coverage).")

    t1, t2, t3, t4, t5, t6, t7, t6_und, truth = [
        d.loc[common] for d in (t1, t2, t3, t4, t5, t6, t7, t6_und, truth)
    ]

    # Renormalize everything except PHLAME / PHLAME (low cov)
    t1 = renormalize_columns_to_100(t1)
    t2 = renormalize_columns_to_100(t2)
    t4 = renormalize_columns_to_100(t4)
    t6 = renormalize_columns_to_100(t6)
    t7 = renormalize_columns_to_100(t7)

    return dict(
        tool_names=['Strainify', 'StrainScan', 'PHLAME', 'StrainScan (low cov)',
                    'PHLAME (low cov)', 'PanTax', 'StrainGE'],
        dfs=[t1, t2, t3, t4, t5, t6, t7],
        # Strainify / StrainScan / StrainGE: "*" when prediction == 0;
        # PHLAME always plotted; PanTax "*" only when "undetected"
        always_plot=[False, False, True, False, True, True, False],
        pantax_undetected=t6_und,
        truth=truth,
        strains=common,
        coverages_to_plot=['10x', '20x', '50x'],
        width=0.16,
        legend_kw=dict(ncol=7, fontsize=110),
    )


def load_low_cov():
    STRAINGE_CSV = "/home/Users/rl152/Strainify_dev/StrainGE_benchmarks/combined_rapct/4_strain_1-2-5x__ecoli.csv"

    t1 = pd.read_csv("/home/Users/rl152/Strainify/4_strain_ecoli_simulated/results_1-2-5x/strainify_results/abundance_estimates_combined.csv", index_col=0)
    t6, t6_und = load_pantax("/dodo/rl152/pantax/4_strain_simulated/pantax_results_low_cov/pantax_4strain_relative_abundances_low_cov.csv")
    t7 = pd.read_csv(STRAINGE_CSV, index_col=0)
    truth = load_truth()

    t1 = t1.rename(columns={'strain name': 'Strain_Name'})

    t6 = t6 * 100
    t7 = t7 * 100

    common = sorted(set(t1.index) & set(t6.index) & set(t7.index) & set(truth.index))
    if not common:
        raise ValueError("No common strains found (low coverage).")

    t1, t6, t7, t6_und, truth = [d.loc[common] for d in (t1, t6, t7, t6_und, truth)]

    t1 = renormalize_columns_to_100(t1)
    t6 = renormalize_columns_to_100(t6)
    t7 = renormalize_columns_to_100(t7)

    return dict(
        tool_names=['Strainify', 'PanTax', 'StrainGE'],
        dfs=[t1, t6, t7],
        always_plot=[False, True, False],
        pantax_undetected=t6_und,
        truth=truth,
        strains=common,
        coverages_to_plot=['1x', '2x', '5x'],
        width=0.28,
        legend_kw=dict(ncol=3, fontsize=130, handlelength=2.2, handleheight=1.4,
                       handletextpad=0.6, columnspacing=3.0),
    )


def prepare_panel(cfg):
    """Coverage grouping, grid size and bar positions for one panel."""
    tool_names = cfg['tool_names']
    dfs = cfg['dfs']
    strains = cfg['strains']
    width = cfg['width']

    groups = [group_columns_by_coverage(d) for d in dfs]
    common_coverages = set(cfg['coverages_to_plot'])
    for g in groups:
        common_coverages &= set(g)
    common_coverages = sorted(common_coverages, key=lambda c: int(c.rstrip('x')))
    if not common_coverages:
        raise ValueError(f"No common coverage groups for {cfg['coverages_to_plot']}.")

    num_ratios = len(groups[0][common_coverages[0]])

    strain_spacing = 1.8
    x = np.arange(len(strains)) * strain_spacing
    offsets = (np.arange(len(tool_names)) - (len(tool_names) - 1) / 2) * width

    return groups, common_coverages, num_ratios, x, offsets


def fill_axes(axes, cfg, groups, common_coverages, x, offsets):
    tool_names = cfg['tool_names']
    dfs = cfg['dfs']
    always_plot = cfg['always_plot']
    truth_df = cfg['truth']
    strains = cfg['strains']
    width = cfg['width']
    pantax_idx = tool_names.index('PanTax')

    for i, (col, display_col) in enumerate(ratio_labels.items()):
        for j, cov in enumerate(common_coverages):
            ax = axes[i][j]

            cols_per_tool = [g[cov] for g in groups]
            true_cols = [f'ratio_{k+1}' for k in range(len(cols_per_tool[0]))]
            true = truth_df[true_cols].copy()

            preds, diffs = [], []
            for d, cols in zip(dfs, cols_per_tool):
                p = d[cols].copy()
                p.columns = true.columns
                preds.append(p[col])
                diffs.append((p[col] - true[col]).replace([np.inf, -np.inf], np.nan).fillna(0))

            undetected6 = cfg['pantax_undetected'][cols_per_tool[pantax_idx]].copy()
            undetected6.columns = true.columns

            visible_heights = []
            for idx in range(len(strains)):
                for k, tool_name in enumerate(tool_names):
                    pred_val = preds[k].iloc[idx]
                    diff_val = diffs[k].iloc[idx]
                    xpos = x[idx] + offsets[k]

                    is_pantax_undetected = (
                        tool_name == 'PanTax' and bool(undetected6[col].iloc[idx])
                    )

                    if is_pantax_undetected or (pred_val == 0 and not always_plot[k]):
                        ax.text(xpos, 1, '*', ha='center', va='bottom',
                                fontsize=85, color='red')
                    else:
                        ax.bar(xpos, diff_val, width=width, alpha=0.8,
                               color=colors[tool_name])
                        visible_heights.append(diff_val)

            if not visible_heights:
                ymin, ymax = -1, 1
            else:
                ymin = min(min(visible_heights), 0)
                ymax = max(max(visible_heights), 0)

            pad = 2
            ymin -= pad
            ymax += pad
            ymin = 10 * math.floor(ymin / 10)
            ymax = 10 * math.ceil(ymax / 10)
            if ymin == ymax:
                ymin -= 10
                ymax += 10

            ax.set_ylim(ymin, ymax)
            ax.yaxis.set_major_locator(MultipleLocator(10))
            ax.axhline(0, color='black', linewidth=1.5)

            if i == 0:
                ax.set_title(f'{cov}', fontsize=130, fontweight='bold')

            if j == 0:
                ax.text(-0.22, 0.5, display_col, transform=ax.transAxes,
                        fontsize=130, rotation=90, va='center', ha='right',
                        fontweight='bold')
                ax.set_ylabel('Delta', fontsize=125, rotation=90, labelpad=-2)

            if i == len(ratio_labels) - 1:
                ax.set_xticks(x)
                ax.set_xticklabels(strains, rotation=90, ha='center', fontsize=100)
            else:
                ax.set_xticks([])

            ax.tick_params(axis='both', labelsize=115)


# =============================================================================
# Build the combined figure
# =============================================================================
high = load_high_cov()
low = load_low_cov()

# Prepare both panels (grid dimensions)
groups_hi, covs_hi, nrat_hi, x_hi, off_hi = prepare_panel(high)
groups_lo, covs_lo, nrat_lo, x_lo, off_lo = prepare_panel(low)


h_hi = 15 * nrat_hi + 25
h_lo = 15 * nrat_lo + 25
fig_width = 36 * max(len(covs_hi), len(covs_lo)) + 10
fig_height = h_hi + h_lo
fig = plt.figure(figsize=(fig_width, fig_height))

# Fraction of the canvas each half occupies (bottom half = [0, split], top = [split, 1])
split = h_lo / fig_height


def to_half(rect, lo, hi):
    """Map an original full-figure rect [l, b, r, t] into a vertical band [lo, hi]."""
    l, b, r, t = rect
    return [l, lo + b * (hi - lo), r, lo + t * (hi - lo)]


# --- Top panel: high coverage ---
gs_hi = fig.add_gridspec(nrat_hi, len(covs_hi))
axes_hi = gs_hi.subplots(sharex=True, squeeze=False)
fill_axes(axes_hi, high, groups_hi, covs_hi, x_hi, off_hi)
fig.legend(
    handles=[Patch(facecolor=colors[t], alpha=0.8, label=t) for t in high['tool_names']],
    loc='lower center',
    bbox_to_anchor=(0.5, split + 0.02 * (1 - split)),
    **high['legend_kw']
)

# --- Bottom panel: low coverage ---
gs_lo = fig.add_gridspec(nrat_lo, len(covs_lo))
axes_lo = gs_lo.subplots(sharex=True, squeeze=False)
fill_axes(axes_lo, low, groups_lo, covs_lo, x_lo, off_lo)
fig.legend(
    handles=[Patch(facecolor=colors[t], alpha=0.8, label=t) for t in low['tool_names']],
    loc='lower center',
    bbox_to_anchor=(0.5, 0.02 * split),
    **low['legend_kw']
)


gs_hi.tight_layout(fig, rect=to_half([0, 0.08, 1, 0.95], split, 1.0))
gs_lo.tight_layout(fig, rect=to_half([0.05, 0.085, 0.99, 0.975], 0.0, split))
fig.align_ylabels(axes_hi[:, 0])
fig.align_ylabels(axes_lo[:, 0])


fig.canvas.draw()
renderer = fig.canvas.get_renderer()
inv = fig.transFigure.inverted()
x0s, x1s = [], []
for ax in axes_lo.ravel():
    bb = ax.get_tightbbox(renderer).transformed(inv)
    x0s.append(bb.x0)
    x1s.append(bb.x1)
content_x0, content_x1 = min(x0s), max(x1s)
delta = 0.5 - (content_x0 + content_x1) / 2.0
delta = max(delta, 0.005 - content_x0)
for ax in axes_lo.ravel():
    pos = ax.get_position()
    ax.set_position([pos.x0 + delta, pos.y0, pos.width, pos.height])

# --- Panel labels: "A" on the top panel, "B" on the bottom panel ---

PANEL_LABEL_SIZE = 170
fig.canvas.draw()
renderer = fig.canvas.get_renderer()


def content_bboxes(axes_panel):
    return [ax.get_tightbbox(renderer).transformed(inv) for ax in axes_panel.ravel()]


label_x = min(bb.x0 for bb in content_bboxes(axes_hi))
for label, axes_panel in (('A', axes_hi), ('B', axes_lo)):
    top = max(bb.y1 for bb in content_bboxes(axes_panel))
    fig.text(label_x, top, label, fontsize=PANEL_LABEL_SIZE, fontweight='bold',
             ha='left', va='top')


plt.savefig(OUT_PDF, bbox_inches='tight', dpi=600)
plt.close()
print(f"Saved {OUT_PDF}")