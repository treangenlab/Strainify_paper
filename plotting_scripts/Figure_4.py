#!/usr/bin/env python3
"""
Combine Strainify, ChronoStrain and StrainGE strain-abundance results into ONE
large multi-panel figure, in the style of the ChronoStrain paper benchmark
figure.

Layout (species are COLUMNS, metric/tool are ROWS):

    row a : L1            -- all three tools in ONE boxplot panel
    row b : RMSE-log      -- all three tools in ONE boxplot panel
    row c : r^2 scatter   -- Strainify        (truth vs. prediction, log-log)
    row d : r^2 scatter   -- ChronoStrain
    row e : r^2 scatter   -- StrainGE

Strainify and StrainGE metrics are computed here directly from their raw output
(same "ground-truth-only" evaluation as plot_strainify_strainge.py). ChronoStrain
is heavy to parse (needs the chronostrain package, marker-ratio pickles, an .ini,
etc.), so we DO NOT re-run it here. Instead we read the CSVs that
`evaluate_chronostrain_cami.py` already writes:

    <chronostrain_eval_dir>/<species_short>/per_sample_metrics.csv
    <chronostrain_eval_dir>/<species_short>/per_genome_truth_vs_pred.csv

Run `evaluate_chronostrain_cami.py ... --out-dir <chronostrain_eval_dir>` first,
then point `chronostrain_eval_dir` (below) at that directory.

The metric calculations are IDENTICAL to those in the ChronoStrain paper (see the evaluate.ipynb in their CAMI analysis repo):

    * calc_errors (L1, RMSE-log) is ported verbatim, with the notebook's globals
      EPSILON_FOR_LOG_PADDING = 1e-7 and filter_lb = 0. Here RENORMALIZE_PER_SPECIES
      is set to True, so truth and every tool's prediction are renormalised to sum
      to 1 over the within-species strains before comparison (common scale). Set
      it back to False to reproduce the notebook's community-scale default.
    * r^2 (log-abundance) = sklearn.metrics.r2_score(log10(truth+eps), log10(pred+eps))
      pooled over all (sample, strain) pairs, exactly as plot_truth_vs_pred does.
    * "# of zeros ~ N per sample" = (count of zero predictions) / 100.

All three tools are routed through the SAME calc_errors, so the comparison is
fair. For ChronoStrain we read the per-genome CSV and recompute its per-sample
errors here, so the RENORMALIZE_PER_SPECIES setting is applied uniformly --
therefore run evaluate_chronostrain_cami.py with `--renorm none` so that CSV
holds raw (un-renormalised) Truth/Pred.

Run `evaluate_chronostrain_cami.py ... --out-dir <chronostrain_eval_dir> --renorm none`
first, then point `chronostrain_eval_dir` (below) at that directory.
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

try:
    import sklearn.metrics
    _HAVE_SKLEARN = True
except Exception:
    _HAVE_SKLEARN = False

# =============================================================================
# Paths  (edit these for your machine)
# =============================================================================
strainify_pred_base = "/dodo/rl152/ChronoStrain/chronostrain_cami/cami_dataset/strainify_cami/force_include_all"

strainge_pred_base = "/dodo/rl152/ChronoStrain/chronostrain_cami/core_genome_only/strainge/results"

gt_dir = "/dodo/rl152/ChronoStrain/chronostrain_cami/cami_dataset/short_read/ground_truth"

# Directory passed as --out-dir to evaluate_chronostrain_cami.py
# (must contain <species_short>/per_sample_metrics.csv and
#  <species_short>/per_genome_truth_vs_pred.csv)
chronostrain_eval_dir = "/dodo/rl152/ChronoStrain/chronostrain_cami/cami_dataset/inference/analyze_output/chronostrain_eval_results"

out_dir = "/dodo/rl152/ChronoStrain/chronostrain_cami/cami_dataset/strainify_cami/plots/combined_collapsed_strains_binned_renormalized"

os.makedirs(out_dir, exist_ok=True)

# =============================================================================
# Species: order = COLUMNS of the figure (matches the paper figure ordering)
# =============================================================================
SPECIES_ORDER = ["ecoli", "saureus", "spneumoniae", "kpneumoniae", "efaecium"]

species_map = {
    "ecoli": "Escherichia coli",
    "saureus": "Staphylococcus aureus",
    "spneumoniae": "Streptococcus pneumoniae",
    "kpneumoniae": "Klebsiella pneumoniae",
    "efaecium": "Enterococcus faecium",
}

# Short italic display names used in column titles (genus abbreviated)
species_display = {
    "ecoli": "E. coli",
    "saureus": "S. aureus",
    "spneumoniae": "S. pneumoniae",
    "kpneumoniae": "K. pneumoniae",
    "efaecium": "E. faecium",
}

# =============================================================================
# Tools.  ROW order for the scatter rows is: Strainify, Strainify (collapsed),
# ChronoStrain, StrainGE. Box order inside the L1 / RMSE panels follows the same
# order. "Strainify (collapsed)" is Strainify re-scored after merging strains
# whose core genome is identical (see load_core_genome_groups); ground truth is
# collapsed the same way, so it is the strain-cluster-level view of Strainify
# (analogous to ChronoStrain's internal clustering).
# Colours echo the paper palette (ChronoStrain navy, StrainGE green, Strainify
# orange); the collapsed variant gets a distinct purple.
# =============================================================================
TOOLS = ["Strainify", "Strainify (collapsed)", "ChronoStrain", "StrainGE"]

TOOL_COLORS = {
    "Strainify": "#E8820E",              # orange
    "Strainify (collapsed)": "#9467BD",  # purple
    "ChronoStrain": "#1f1f6e",           # navy
    "StrainGE": "#2CA02C",               # green
}

sample_ids = range(100)

# =============================================================================
# Metric settings -- copied VERBATIM from the lab notebook (evaluate.ipynb) so
# the numbers are computed identically for all three tools.
#
#   EPSILON_FOR_LOG_PADDING = 1e-7          (notebook cell: "Global settings")
#   RENORMALIZE_PER_SPECIES = False         (notebook cell: "Global settings")
#   FILTER_LB               = 0.0           (calc_errors default; keeps truth > 0)
#
# RENORMALIZE_PER_SPECIES controls within-species composition normalization.
# When True, each tool's FULL predicted vector (over every strain it reports) and
# the truth are divided by their own per-sample totals so they sum to 1 BEFORE any
# strain subsetting -- so all tools are compared on one common relative-abundance
# scale (this also fixes tools whose raw output does not sum to 1, e.g. StrainGE).
# The metrics are then evaluated on the present (ground-truth) strains using those
# already-normalized values (calc_errors/_scatter_pair are called with
# renormalize=False). Consequence: mass a tool places on ABSENT strains stays in
# the denominator, so the predictions over present strains sum to <= 1 and that
# off-target mass counts against the tool (it is not silently dropped). corr^2 is
# unaffected by this choice (shift-invariant); r^2 / L1 / RMSE-log reflect it.
# Set to False to use raw community-scale values (the notebook default).
# =============================================================================
EPSILON_FOR_LOG_PADDING = 1e-7
RENORMALIZE_PER_SPECIES = True
FILTER_LB = 0.0

# Ground-truth column used as "truth" for Strainify / StrainGE evaluation
# (the notebook's `restrict_profile`; "strain_relative_abundance" sums to 1
# within a species, "original_percentage" is the community-scale value).
TRUTH_COLUMN = "strain_relative_abundance"

# The notebook reports zeros as (number of zero predictions) / 100, i.e. per
# the fixed 100-sample CAMI design. Kept identical here.
N_SAMPLES_FOR_ZEROS = 100

# Back-compat alias used throughout the file.
eps = EPSILON_FOR_LOG_PADDING


def safe_name(name):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)


# =============================================================================
# Loaders for Strainify / StrainGE / ground truth  
# =============================================================================
def load_strainify_predictions(species_short):
    path = os.path.join(
        strainify_pred_base, species_short, "abundance_estimates_combined.csv"
    )
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    df = pd.read_csv(path)
    if "strain name" not in df.columns:
        raise ValueError(f"'strain name' column not found in {path}")

    df = df.rename(columns={"strain name": "genome_id"})
    df["genome_id"] = df["genome_id"].astype(str).str.strip()

    sample_cols = [c for c in df.columns if c.startswith("sample_")]
    if not sample_cols:
        raise ValueError(f"No sample columns found in {path}")

    df[sample_cols] = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    # Strainify output is percent scale -> convert to 0-1 fraction
    df[sample_cols] = df[sample_cols] / 100.0
    return df.set_index("genome_id")[sample_cols]


# VCF metadata columns in filtered_variant_matrix.csv (everything else is a strain)
_VCF_META_COLS = {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"}


def load_core_genome_groups(species_short):
    """Group strains whose CORE GENOME is identical, from filtered_variant_matrix.csv.

    The matrix is VCF-style: rows are variant sites, columns are VCF metadata then
    one 0/1 genotype column per strain. Two strains with an identical genotype
    vector across all rows have no distinguishing core-genome variants -> identical
    core genome. Returns dict {strain_name -> representative_name}; strains absent
    from the matrix are not in the dict (treated as their own singleton group).
    """
    path = os.path.join(strainify_pred_base, species_short, "filtered_variant_matrix.csv")
    if not os.path.exists(path):
        return None

    vm = pd.read_csv(path)
    strain_cols = [c for c in vm.columns if c not in _VCF_META_COLS]
    if not strain_cols:
        return None

    # Missing calls (NaN) are treated as a distinct value (-1) so they only match
    # other missing calls -- conservative, never merges on the basis of absent data.
    geno = (
        vm[strain_cols].apply(pd.to_numeric, errors="coerce")
        .fillna(-1).round().astype(np.int16).to_numpy().T  # strains x variants
    )
    _, inv = np.unique(geno, axis=0, return_inverse=True)

    members = {}
    for s, lab in zip(strain_cols, inv):
        members.setdefault(int(lab), []).append(str(s).strip())

    group_map = {}
    for strs in members.values():
        rep = sorted(strs)[0]                 # deterministic representative
        for s in strs:
            group_map[s] = rep

    n_collapsible = sum(1 for v in members.values() if len(v) > 1)
    print(f"  core-genome groups[{species_short}]: {len(strain_cols)} strains "
          f"-> {len(members)} groups ({n_collapsible} multi-strain)")
    return group_map


def collapse_series(series, group_map):
    """Sum abundances of strains that share a core-genome group.

    series: index = strain/genome ids, values = abundance. Strains not in
    group_map keep their own id (singleton group).
    """
    idx = series.index.astype(str).str.strip()
    grp = [group_map.get(s, s) for s in idx]
    return series.groupby(grp).sum()


def load_strainge_predictions(species_short):
    pred_by_sample = {}
    for i in sample_ids:
        sample = f"sample_{i}"
        path = os.path.join(
            strainge_pred_base, sample, species_short, "result.strains.tsv"
        )
        if not os.path.exists(path):
            continue
        try:
            df = pd.read_csv(path, sep="\t")
        except pd.errors.EmptyDataError:
            continue
        if df.empty:
            continue
        if "strain" not in df.columns or "rapct" not in df.columns:
            print(f"  Warning: missing strain or rapct column in {path}")
            continue

        df["strain"] = df["strain"].astype(str).str.strip()
        df["rapct"] = pd.to_numeric(df["rapct"], errors="coerce").fillna(0)
        s = df.groupby("strain")["rapct"].sum()
        total = s.sum()
        s = s / total if total > 0 else s * 0.0
        pred_by_sample[sample] = s

    if not pred_by_sample:
        raise ValueError(f"No StrainGE predictions found for {species_short}")

    pred_df = pd.DataFrame(pred_by_sample).fillna(0)
    pred_df.index = pred_df.index.astype(str).str.strip()
    return pred_df


def load_ground_truth(species_name):
    folder = os.path.join(gt_dir, safe_name(species_name))
    if not os.path.exists(folder):
        raise FileNotFoundError(folder)

    gt_by_sample = {}
    for i in sample_ids:
        sample = f"sample_{i}"
        path = os.path.join(folder, f"{sample}_ground_truth.csv")
        if not os.path.exists(path):
            continue

        df = pd.read_csv(path)
        if "genome_id" not in df.columns or TRUTH_COLUMN not in df.columns:
            print(f"  Warning: missing 'genome_id' or '{TRUTH_COLUMN}' in {path}")
            continue

        df["genome_id"] = df["genome_id"].astype(str).str.strip()
        df[TRUTH_COLUMN] = pd.to_numeric(df[TRUTH_COLUMN], errors="coerce").fillna(0)
        gt_by_sample[sample] = df.groupby("genome_id")[TRUTH_COLUMN].sum()
    return gt_by_sample


# =============================================================================
# Metric core -- PORTED VERBATIM from evaluate.ipynb :: calc_errors
# (operates on already-aligned 1-D truth/pred arrays for the evaluated strains)
# =============================================================================
def normalize_to_sum1(v):
    """Divide a Series/array by its total so the relative abundances sum to 1.

    No-op if the total is <= 0. Applied to each tool's FULL predicted composition
    (and to the truth) once per sample, BEFORE any strain subsetting, so every tool
    is placed on the same relative-abundance scale. Mass a tool puts on absent
    strains stays in the denominator, so it counts against the tool rather than
    being silently dropped.
    """
    arr = np.asarray(v, dtype=float)
    s = float(np.nansum(arr))
    if s <= 0:
        return v
    return v / s


def calc_errors(_truth, _pred, eps=EPSILON_FOR_LOG_PADDING,
                renormalize=RENORMALIZE_PER_SPECIES, filter_lb=FILTER_LB):
    """Returns (l1_error, rms_log_error), identical to the notebook.

    _pred may be 1-D, or 2-D (ChronoStrain posterior samples) -> median trajectory.
    """
    _truth = np.asarray(_truth, dtype=float)
    _pred = np.asarray(_pred, dtype=float)

    if _pred.ndim == 2:
        _pred = np.median(_pred, axis=0)          # ChronoStrain: median trajectory
    elif _pred.ndim != 1:
        raise ValueError(f"bad prediction shape {_pred.shape}")
    assert _pred.shape == _truth.shape, (_pred.shape, _truth.shape)

    if renormalize:
        sp, st = _pred.sum(), _truth.sum()
        _pred = _pred / sp if sp > 0 else _pred
        _truth = _truth / st if st > 0 else _truth

    # Filter by lower bound. (filter_lb = 0.0 keeps strains with truth > 0.)
    (_locs,) = np.where(_truth > filter_lb)
    _truth = _truth[_locs]
    _pred = _pred[_locs]

    l1_error = np.abs(_pred - _truth).sum()
    rms_log_error = np.sqrt(np.mean(np.square(
        np.log10(_pred + eps) - np.log10(_truth + eps)
    )))
    return l1_error, rms_log_error


def _scatter_pair(_truth, _pred, renormalize=RENORMALIZE_PER_SPECIES):
    """Truth/pred pair as recorded for the pooled r^2 scatter.

    Mirrors the notebook's evaluate_each_genome: with renormalize=True each is
    divided by its own sum; otherwise raw values are used. NO truth>0 filter is
    applied here (the scatter includes every evaluated strain, matching the
    notebook). _pred may be 2-D -> median is taken.
    """
    _truth = np.asarray(_truth, dtype=float)
    _pred = np.asarray(_pred, dtype=float)
    if _pred.ndim == 2:
        _pred = np.median(_pred, axis=0)
    if renormalize:
        sp, st = _pred.sum(), _truth.sum()
        _pred = _pred / sp if sp > 0 else _pred
        _truth = _truth / st if st > 0 else _truth
    return _truth, _pred


# =============================================================================
# Per-tool metric computation for Strainify / StrainGE
# (ground-truth strain set per sample; same math as the notebook's calc_errors)
# =============================================================================
def compute_metrics(pred_df, gt_dict, tool_name, group_map=None):
    """Returns (metrics_df, pooled_true, pooled_pred, zeros_per_sample).

    metrics_df has columns: sample, tool, L1, RMSE_log, n_strains_evaluated.
    pooled_true / pooled_pred are 1-D arrays over all (sample, evaluated-strain)
    pairs, recorded the same way the notebook records them for the r^2 scatter.

    If group_map is given, BOTH the prediction and the ground truth are collapsed
    by summing strains that share a core-genome group, before any metric is
    computed -- i.e. evaluation happens at the identical-core-genome cluster level.
    """
    rows = []
    all_y_true, all_y_pred = [], []
    n_zeros_total = 0

    for sample, gt_vec in gt_dict.items():
        if sample not in pred_df.columns:
            continue
        pred_vec = pred_df[sample]

        # Normalize each tool's FULL predicted composition (and the truth) to sum
        # to 1 BEFORE any collapsing/subsetting, so all tools share one scale.
        if RENORMALIZE_PER_SPECIES:
            pred_vec = normalize_to_sum1(pred_vec)
            gt_vec = normalize_to_sum1(gt_vec)

        # collapse strains with identical core genome (truth AND prediction)
        if group_map is not None:
            gt_vec = collapse_series(gt_vec, group_map)
            pred_vec = collapse_series(pred_vec, group_map)

        gt_genomes = sorted(list(gt_vec.index))
        y_true = gt_vec.reindex(gt_genomes).fillna(0).to_numpy(dtype=float)
        y_pred = pred_vec.reindex(gt_genomes).fillna(0).to_numpy(dtype=float)

        if y_true.sum() == 0:
            continue

        # ---- per-sample L1 / RMSE-log via the ported notebook metric ----
        # (vectors are already summed-to-1 above, so renormalize=False here)
        l1, rmse_log = calc_errors(
            y_true, y_pred,
            eps=EPSILON_FOR_LOG_PADDING,
            renormalize=False,
            filter_lb=FILTER_LB,
        )
        rows.append({
            "sample": sample,
            "tool": tool_name,
            "L1": l1,
            "RMSE_log": rmse_log,
            "n_strains_evaluated": len(gt_genomes),
        })

        # ---- pooled truth/pred for the r^2 scatter (notebook convention) ----
        pt, pp = _scatter_pair(y_true, y_pred, renormalize=False)
        all_y_true.extend(pt)
        all_y_pred.extend(pp)
        n_zeros_total += int(np.sum(pp == 0))

    metrics_df = pd.DataFrame(rows)
    zeros_per_sample = n_zeros_total / N_SAMPLES_FOR_ZEROS
    return metrics_df, np.array(all_y_true), np.array(all_y_pred), zeros_per_sample


# =============================================================================
# ChronoStrain loader: read the CSVs written by evaluate_chronostrain_cami.py.
#
# We recompute the per-sample L1 / RMSE-log and the pooled scatter HERE, from the
# raw per-genome Truth/Pred, using the SAME calc_errors / _scatter_pair as the
# other two tools. That way RENORMALIZE_PER_SPECIES / eps / filter_lb are applied
# uniformly to all three methods, independent of how the CSV was generated.
#
# IMPORTANT: run evaluate_chronostrain_cami.py with `--renorm none` so the CSV
# holds RAW (un-renormalised) Truth/Pred; this script then applies the notebook's
# RENORMALIZE_PER_SPECIES setting itself.
# =============================================================================
def load_chronostrain_results(species_short):
    """Returns (metrics_df, pooled_true, pooled_pred, zeros_per_sample) or None."""
    spec_dir = os.path.join(chronostrain_eval_dir, species_short)
    per_sample = os.path.join(spec_dir, "per_sample_metrics.csv")
    per_genome = os.path.join(spec_dir, "per_genome_truth_vs_pred.csv")

    # ---- preferred path: recompute from raw per-genome Truth/Pred ----
    if os.path.exists(per_genome):
        g_df = pd.read_csv(per_genome)
        g_df["Truth"] = pd.to_numeric(g_df["Truth"], errors="coerce").fillna(0)
        g_df["Pred"] = pd.to_numeric(g_df["Pred"], errors="coerce").fillna(0)

        rows, all_t, all_p, n_zeros_total = [], [], [], 0
        for sample_id, sub in g_df.groupby("Sample"):
            y_true = sub["Truth"].to_numpy(dtype=float)
            y_pred = sub["Pred"].to_numpy(dtype=float)
            # normalize this tool's full composition (and truth) to sum to 1
            if RENORMALIZE_PER_SPECIES:
                y_true = normalize_to_sum1(y_true)
                y_pred = normalize_to_sum1(y_pred)
            if y_true.sum() == 0:
                continue
            l1, rmse_log = calc_errors(
                y_true, y_pred,
                eps=EPSILON_FOR_LOG_PADDING,
                renormalize=False,
                filter_lb=FILTER_LB,
            )
            rows.append({
                "sample": str(sample_id),
                "tool": "ChronoStrain",
                "L1": l1,
                "RMSE_log": rmse_log,
                "n_strains_evaluated": int(len(y_true)),
            })
            pt, pp = _scatter_pair(y_true, y_pred, renormalize=False)
            all_t.extend(pt)
            all_p.extend(pp)
            n_zeros_total += int(np.sum(pp == 0))

        metrics_df = pd.DataFrame(rows)
        zeros_per_sample = n_zeros_total / N_SAMPLES_FOR_ZEROS
        return metrics_df, np.array(all_t), np.array(all_p), zeros_per_sample

    # ---- fallback: only per-sample metrics available (no scatter possible) ----
    if os.path.exists(per_sample):
        s_df = pd.read_csv(per_sample)
        metrics_df = pd.DataFrame({
            "sample": s_df["Sample"].astype(str),
            "tool": "ChronoStrain",
            "L1": pd.to_numeric(s_df["L1"], errors="coerce"),
            "RMSE_log": pd.to_numeric(s_df["RMSLE"], errors="coerce"),
        })
        print(f"  ChronoStrain[{species_short}]: per_genome CSV missing -- using "
              f"precomputed per-sample L1/RMSLE (note: its renorm setting may "
              f"differ from RENORMALIZE_PER_SPECIES={RENORMALIZE_PER_SPECIES}).")
        return metrics_df, np.array([]), np.array([]), np.nan

    return None


# =============================================================================
# Pooled r^2 metrics on log10 abundance:
#   pooled_r2    = coefficient of determination (sklearn r2_score) -- paper-style,
#                  distance from the y=x line; can be negative.
#   pooled_corr2 = squared Pearson correlation (corrcoef**2) -- as in the older
#                  plot_strainify_strainge.py; fit of the BEST line, always 0..1.
# =============================================================================
def pooled_r2(true_vals, pred_vals):
    if len(true_vals) < 2:
        return np.nan
    lt = np.log10(true_vals + eps)
    lp = np.log10(pred_vals + eps)
    if _HAVE_SKLEARN:
        return float(sklearn.metrics.r2_score(y_true=lt, y_pred=lp))
    # manual coefficient of determination fallback
    ss_res = np.sum((lt - lp) ** 2)
    ss_tot = np.sum((lt - np.mean(lt)) ** 2)
    return float(1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan


def pooled_corr2(true_vals, pred_vals):
    if len(true_vals) < 2:
        return np.nan
    lt = np.log10(true_vals + eps)
    lp = np.log10(pred_vals + eps)
    corr = np.corrcoef(lt, lp)[0, 1]
    return float(corr ** 2) if np.isfinite(corr) else np.nan



# =============================================================================
# Figure assembly
# =============================================================================
def build_figure(results):
    """results[species_short][tool] = dict(metrics, true, pred, zeros, r2) or None."""
    n_cols = len(SPECIES_ORDER)
    n_scatter = len(TOOLS)
    n_rows = 2 + n_scatter + 1  # a:L1, b:RMSE, scatter per tool, + binned RMSE-log row

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(3.7 * n_cols, 3.5 * n_rows),
        squeeze=False,
    )

    row_letters = [chr(ord("a") + i) for i in range(n_rows)]
    # scatter rows map to tools in TOOLS order; the final row is the binned panel
    scatter_rows = range(2, 2 + n_scatter)
    scatter_row_tool = {2 + i: tool for i, tool in enumerate(TOOLS)}
    scatter_ylabel = {tool: f"{tool} Estimate" for tool in TOOLS}
    binned_row = 2 + n_scatter

    box_positions = list(range(1, len(TOOLS) + 1))

    for j, sp in enumerate(SPECIES_ORDER):
        res = results.get(sp, {})

        # number of GT genomes (union/max evaluated) and samples, for the title
        n_genomes, n_samples = 0, 0
        for tool in TOOLS:
            r = res.get(tool)
            if r is not None and not r["metrics"].empty:
                if "n_strains_evaluated" in r["metrics"].columns:
                    n_genomes = max(n_genomes, int(r["metrics"]["n_strains_evaluated"].max()))
                n_samples = max(n_samples, r["metrics"].shape[0])

        # ---------------- row a: L1 ----------------
        ax = axes[0][j]
        _boxplot_panel(ax, res, "L1", box_positions, y_from_zero=True)
        if j == 0:
            ax.set_ylabel("L1")
        gtxt = f"{n_genomes} GT genomes\n" if n_genomes else ""
        ax.set_title(
            f"$\\it{{{species_display[sp]}}}$: {gtxt}{n_samples} samples",
            fontsize=11,
        )

        # ---------------- row b: RMSE-log ----------------
        ax = axes[1][j]
        _boxplot_panel(ax, res, "RMSE_log", box_positions, y_from_zero=True)
        if j == 0:
            ax.set_ylabel(f"RMSE-log (eps={EPSILON_FOR_LOG_PADDING})")

        # ---------------- scatter rows: one per tool ----------------
        for row in scatter_rows:
            ax = axes[row][j]
            tool = scatter_row_tool[row]
            r = res.get(tool)
            _scatter_panel(
                ax, r, TOOL_COLORS[tool],
                ylabel=scatter_ylabel[tool] if j == 0 else None,
            )

        # ---------------- final row: binned RMSE-log error vs abundance quantiles --
        ax = axes[binned_row][j]
        _binned_panel(ax, res)
        if j == 0:
            ax.set_ylabel("RMSE-log Error")

        # panel labels a1..e5 (top-left, bold)
        for row in range(n_rows):
            axes[row][j].text(
                -0.02, 1.06, f"{row_letters[row]}{j + 1}",
                transform=axes[row][j].transAxes,
                fontsize=13, fontweight="bold", va="bottom", ha="right",
            )

    # ---------------- shared legend ----------------
    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="",
               markerfacecolor=TOOL_COLORS[t], markeredgecolor="black",
               markersize=10, label=t)
        for t in TOOLS
    ]
    fig.legend(
        handles=legend_handles, title="Methods",
        loc="lower center", ncol=len(TOOLS),
        frameon=True, bbox_to_anchor=(0.5, -0.005), fontsize=12, title_fontsize=12,
    )

    fig.tight_layout(rect=[0, 0.035, 1, 1])
    return fig


def _boxplot_panel(ax, res, metric_key, positions, y_from_zero=True):
    """Box + jittered points for the three tools, coloured per tool.

    The y-range is data-driven (no fixed cap): it spans the actual values with a
    little headroom, anchored at 0 when y_from_zero (L1 / RMSE-log are >= 0).
    """
    rng = np.random.default_rng(1)
    data, present_positions, present_colors = [], [], []

    for pos, tool in zip(positions, TOOLS):
        r = res.get(tool)
        if r is None or r["metrics"].empty or metric_key not in r["metrics"].columns:
            continue
        vals = r["metrics"][metric_key].dropna().to_numpy()
        if vals.size == 0:
            continue
        data.append(vals)
        present_positions.append(pos)
        present_colors.append(TOOL_COLORS[tool])

    if not data:
        ax.set_xticks([])
        ax.text(0.5, 0.5, "no data", transform=ax.transAxes,
                ha="center", va="center", color="gray")
        return

    bp = ax.boxplot(
        data, positions=present_positions, widths=0.5,
        patch_artist=True, showfliers=False,
        whis=(2.5, 97.5),  # notebook: whis=(0.025, 0.975)
        medianprops=dict(color="gold", linewidth=2.0),  # notebook gold medians
        whiskerprops=dict(color="black"), capprops=dict(color="black"),
        boxprops=dict(color="black"),
    )
    for patch, color in zip(bp["boxes"], present_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)

    for pos, vals, color in zip(present_positions, data, present_colors):
        x = rng.normal(pos, 0.06, size=len(vals))
        ax.scatter(x, vals, s=14, alpha=0.75, color=color,
                   edgecolor="black", linewidth=0.3, zorder=10)

    ax.set_xticks([])  # tools identified by colour + shared legend
    ax.set_xlim(0.4, len(TOOLS) + 0.6)

    # data-driven y-range (no hard cap): span the real values with headroom
    all_vals = np.concatenate(data)
    dmin, dmax = float(np.min(all_vals)), float(np.max(all_vals))
    if dmax <= dmin:
        dmax = dmin + 1.0
    pad = 0.05 * (dmax - dmin)
    if y_from_zero:
        ax.set_ylim(0.0, dmax + pad)
    else:
        ax.set_ylim(dmin - pad, dmax + pad)


def _scatter_panel(ax, r, color, ylabel=None):
    """truth-vs-prediction log-log scatter with r^2 and zeros annotations."""
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-6, 1.3)
    ax.set_ylim(5e-8, 1.3)
    ax.plot([1e-6, 1.0], [1e-6, 1.0], linestyle="--", linewidth=1.3,
            color="red", alpha=0.7, zorder=1)
    ax.set_xlabel("True Relative Abundance")
    if ylabel:
        ax.set_ylabel(ylabel)

    if r is None or len(r["true"]) == 0:
        ax.text(0.5, 0.5, "no data", transform=ax.transAxes,
                ha="center", va="center", color="gray")
        return

    ax.scatter(r["true"] + eps, r["pred"] + eps, s=7, alpha=0.5,
               linewidths=0, color=color, rasterized=True, zorder=5)

    # Solid (semi-opaque) white box keeps labels readable over the dense cloud.
    box = dict(boxstyle="round,pad=0.3", facecolor="white",
               edgecolor="0.6", linewidth=0.6, alpha=0.92)
    r2, c2 = r.get("r2"), r.get("corr2")
    ax.text(0.04, 0.96, f"$r^2$ (log-abundance) = {r2:.2f}",
            transform=ax.transAxes, fontsize=11, va="top", ha="left",
            zorder=20, color="black", bbox=box)
    if c2 is not None and np.isfinite(c2):
        ax.text(0.04, 0.85, f"corr$^2$ (log-abundance) = {c2:.2f}",
                transform=ax.transAxes, fontsize=11, va="top", ha="left",
                zorder=20, color="#1a1a1a", bbox=box)

    # "# of zeros ~ N per sample" -- only annotated when the tool emits zeros
    z = r.get("zeros")
    if z is not None and np.isfinite(z) and z > 0:
        z_txt = f"{z:.2f}".rstrip("0").rstrip(".")
        ax.text(0.96, 0.07, f"# of zeros ~ {z_txt} per sample",
                transform=ax.transAxes, fontsize=11, va="bottom", ha="right",
                zorder=20, color="black", bbox=box)


def _binned_panel(ax, res, rng=np.random.default_rng(0)):
    """RMSE-log error vs log10 true-abundance quantiles, per the notebook's
    plot_binned_differences (cell 17):

        Bin  = pd.qcut(log10(Truth + eps), q=10)          # 10 equal-count bins
        Diff = |log10(Pred + eps) - log10(Truth + eps)|   # 'RMSE-log Error'

    Grouped boxplots: one box per (bin x tool), gold medians, whiskers at the
    2.5/97.5 percentiles, with a light subsampled strip overlay -- matching the
    paper's last row. Bins are defined on the pooled truth across all tools so the
    quantile edges are shared.
    """
    # gather per-tool pooled (truth, pred) -> long frame
    frames = []
    for tool in TOOLS:
        r = res.get(tool)
        if r is None or len(r.get("true", [])) == 0:
            continue
        t = np.asarray(r["true"], dtype=float)
        p = np.asarray(r["pred"], dtype=float)
        logT = np.log10(t + eps)
        diff = np.abs(np.log10(p + eps) - logT)
        frames.append(pd.DataFrame({"logT": logT, "Diff": diff, "Method": tool}))

    if not frames:
        ax.text(0.5, 0.5, "no data", transform=ax.transAxes,
                ha="center", va="center", color="gray")
        ax.set_xlabel(r"$\log_{10}$ True Abundance Quantiles")
        return

    df = pd.concat(frames, ignore_index=True)
    try:
        df["Bin"] = pd.qcut(df["logT"], q=10, duplicates="drop")
    except (ValueError, IndexError):
        ax.text(0.5, 0.5, "no data", transform=ax.transAxes,
                ha="center", va="center", color="gray")
        return
    bins = list(df["Bin"].cat.categories)

    present = [t for t in TOOLS if (df["Method"] == t).any()]
    n_m = len(present)
    width = 0.8 / max(n_m, 1)
    offsets = (np.arange(n_m) - (n_m - 1) / 2.0) * width

    for gi, b in enumerate(bins):
        for mi, tool in enumerate(present):
            vals = df.loc[(df["Bin"] == b) & (df["Method"] == tool), "Diff"].to_numpy()
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                continue
            pos = gi + offsets[mi]
            color = TOOL_COLORS[tool]
            bp = ax.boxplot(
                [vals], positions=[pos], widths=width * 0.9,
                patch_artist=True, showfliers=False, whis=(2.5, 97.5),
                medianprops=dict(color="gold", linewidth=1.5),
                whiskerprops=dict(color="black", linewidth=0.6),
                capprops=dict(color="black", linewidth=0.6),
                boxprops=dict(color="black", linewidth=0.5),
            )
            bp["boxes"][0].set_facecolor(color)
            bp["boxes"][0].set_alpha(0.55)
            # light, subsampled strip overlay (rasterized) -- paper aesthetic
            if vals.size > 150:
                vals_s = rng.choice(vals, 150, replace=False)
            else:
                vals_s = vals
            x = rng.normal(pos, width * 0.12, size=vals_s.size)
            ax.scatter(x, vals_s, s=2, alpha=0.2, color=color,
                       edgecolor="black", linewidth=0.3, rasterized=True, zorder=3)

    ax.set_xticks(range(len(bins)))
    ax.set_xticklabels([f"({b.left:.2f}, {b.right:.2f}]" for b in bins],
                       rotation=20, ha="right", fontsize=6.5)
    ax.set_xlim(-0.6, len(bins) - 0.4)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(r"$\log_{10}$ True Abundance Quantiles")


# =============================================================================
# Main
# =============================================================================
def main():
    results = {}
    all_metrics_rows = []
    all_r2_rows = []

    for sp in SPECIES_ORDER:
        species_name = species_map[sp]
        print(f"\nProcessing {species_name} ({sp})...")
        results[sp] = {}

        # ground truth (shared by Strainify & StrainGE evaluation)
        try:
            gt_by_sample = load_ground_truth(species_name)
        except Exception as e:
            print(f"  Failed to load ground truth: {e}")
            gt_by_sample = {}

        # ----- Strainify -----
        strainify_pred = None
        try:
            strainify_pred = load_strainify_predictions(sp)
            m, t, p, z = compute_metrics(strainify_pred, gt_by_sample, "Strainify")
            r2, c2 = pooled_r2(t, p), pooled_corr2(t, p)
            results[sp]["Strainify"] = dict(metrics=m, true=t, pred=p, zeros=z, r2=r2, corr2=c2)
            print(f"  Strainify : {m.shape[0]} samples, r2={r2:.3f}, corr2={c2:.3f}, zeros/sample={z:.2f}")
        except Exception as e:
            print(f"  Strainify failed: {e}")
            results[sp]["Strainify"] = None

        # ----- Strainify (collapsed): merge strains with identical core genome -----
        try:
            group_map = load_core_genome_groups(sp)
            if group_map is None:
                print(f"  Strainify (collapsed): no filtered_variant_matrix.csv for {sp} (skipping)")
                results[sp]["Strainify (collapsed)"] = None
            elif strainify_pred is None:
                results[sp]["Strainify (collapsed)"] = None
            else:
                m, t, p, z = compute_metrics(
                    strainify_pred, gt_by_sample, "Strainify (collapsed)",
                    group_map=group_map,
                )
                r2, c2 = pooled_r2(t, p), pooled_corr2(t, p)
                results[sp]["Strainify (collapsed)"] = dict(
                    metrics=m, true=t, pred=p, zeros=z, r2=r2, corr2=c2,)
                print(f"  Strainify (collapsed): {m.shape[0]} samples, "
                      f"r2={r2:.3f}, corr2={c2:.3f}, zeros/sample={z:.2f}")
        except Exception as e:
            print(f"  Strainify (collapsed) failed: {e}")
            results[sp]["Strainify (collapsed)"] = None

        # ----- StrainGE -----
        try:
            pred = load_strainge_predictions(sp)
            m, t, p, z = compute_metrics(pred, gt_by_sample, "StrainGE")
            r2, c2 = pooled_r2(t, p), pooled_corr2(t, p)
            results[sp]["StrainGE"] = dict(metrics=m, true=t, pred=p, zeros=z, r2=r2, corr2=c2)
            print(f"  StrainGE  : {m.shape[0]} samples, r2={r2:.3f}, corr2={c2:.3f}, zeros/sample={z:.2f}")
        except Exception as e:
            print(f"  StrainGE failed: {e}")
            results[sp]["StrainGE"] = None

        # ----- ChronoStrain (read precomputed CSVs) -----
        try:
            cs = load_chronostrain_results(sp)
            if cs is None:
                print(f"  ChronoStrain: no CSVs under {chronostrain_eval_dir}/{sp}/ (skipping)")
                results[sp]["ChronoStrain"] = None
            else:
                m, t, p, z = cs
                r2, c2 = pooled_r2(t, p), pooled_corr2(t, p)
                results[sp]["ChronoStrain"] = dict(metrics=m, true=t, pred=p, zeros=z, r2=r2, corr2=c2)
                print(f"  ChronoStrain: {m.shape[0]} samples, r2={r2:.3f}, corr2={c2:.3f}, zeros/sample={z:.2f}")
        except Exception as e:
            print(f"  ChronoStrain failed: {e}")
            results[sp]["ChronoStrain"] = None

        # collect tidy tables
        for tool in TOOLS:
            rr = results[sp].get(tool)
            if rr is not None and not rr["metrics"].empty:
                mm = rr["metrics"].copy()
                mm["species"] = species_name
                mm["species_short"] = sp
                all_metrics_rows.append(mm)
                all_r2_rows.append({
                    "species": species_name, "species_short": sp, "tool": tool,
                    "r2_log_abundance": rr["r2"],
                    "corr2_log_abundance": rr.get("corr2"),
                    "zeros_per_sample": rr["zeros"],
                    "n_points": len(rr["true"]),
                })

    # ----- build & save the big figure -----
    fig = build_figure(results)
    #out_png = os.path.join(out_dir, "combined_three_tools_benchmark.png")
    out_pdf = os.path.join(out_dir, "combined_three_tools_benchmark.pdf")
    #fig.savefig(out_png, dpi=600, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=600, bbox_inches="tight")
    plt.close(fig)
    #print(f"\nSaved figure: {out_png}")
    print(f"Saved figure: {out_pdf}")

    # ----- save tidy metric tables -----
    if all_metrics_rows:
        combined = pd.concat(all_metrics_rows, ignore_index=True)
        csv_path = os.path.join(out_dir, "combined_three_tools_metrics.csv")
        combined.to_csv(csv_path, index=False)
        print(f"Saved metrics table: {csv_path}")

    if all_r2_rows:
        r2_df = pd.DataFrame(all_r2_rows)
        r2_path = os.path.join(out_dir, "combined_three_tools_r2.csv")
        r2_df.to_csv(r2_path, index=False)
        print(f"Saved r2 table: {r2_path}")

    print(f"\nDone. Results saved to: {out_dir}")


if __name__ == "__main__":
    main()