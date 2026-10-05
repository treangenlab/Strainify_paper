import pandas as pd
import numpy as np
import re
import os
from scipy.spatial.distance import jensenshannon

# ---- Config ----
species_list = ["ecoli", "cacnes", "cdiff", "mtuberculosis", "sepidermidis"]

# StrainGE per-species combined rapct CSVs (from combine_rapct.py):
#   <STRAINGE_DIR>/30_strain_simulated__<species>.csv
STRAINGE_DIR = "/home/Users/rl152/Strainify_dev/StrainGE_benchmarks/combined_rapct"
STRAINGE_TEMPLATE = os.path.join(STRAINGE_DIR, "30_strain_simulated__{species}.csv")

# Ground-truth ratios per species
TRUTH_TEMPLATE = "/home/Users/rl152/Strainify/30_strains/{species}/{species}_30_ratios.csv"

# Output: one JSD CSV per species
OUT_DIR = "/home/Users/rl152/Strainify_dev/StrainGE_benchmarks/30_strain_simulated/jsds"
os.makedirs(OUT_DIR, exist_ok=True)


# Function to normalize index by keeping only the first two parts of genome name
def normalize_index(df):
    def clean(name):
        name = re.sub(r'\.fna$|\.fa$|\.fasta$|\.gz$', '', name)
        name = name.split('.')[0]
        #name = re.sub(r'\.\d+$', '', name)  # Remove version suffix like .1
        parts = name.split('_')
        return '_'.join(parts[:2]) if len(parts) >= 2 else name
    df.index = df.index.to_series().apply(clean)
    return df


for species in species_list:
    est_path = STRAINGE_TEMPLATE.format(species=species)
    truth_path = TRUTH_TEMPLATE.format(species=species)

    if not os.path.exists(est_path):
        print(f"[! Warning] Missing StrainGE CSV for {species}: {est_path}")
        continue
    if not os.path.exists(truth_path):
        print(f"[! Warning] Missing truth CSV for {species}: {truth_path}")
        continue

    # Load data
    estimates_df = pd.read_csv(est_path, index_col=0)
    truth_df = pd.read_csv(truth_path, index_col=0)

    # Make sure estimate values are numeric
    estimates_df = estimates_df.apply(pd.to_numeric, errors="coerce").fillna(0)

    # Normalize indices
    estimates_df = normalize_index(estimates_df)
    truth_df = normalize_index(truth_df)

    # Collapse any strains that normalize to the same name (sum their abundances).
    # Prevents a duplicate-index error at reindex if two genomes share a prefix.
    estimates_df = estimates_df.groupby(level=0).sum()
    truth_df = truth_df.groupby(level=0).sum()

    # Build a mapping from estimate sample name to corresponding ground truth column
    sample_to_ratio = {
        col: re.search(r'(ratio_\d+)', col).group(1)
        for col in estimates_df.columns
        if re.search(r'(ratio_\d+)', col)
    }

    if not sample_to_ratio:
        print(f"[! Warning] No sample columns matched 'ratio_\\d+' for {species}. "
              f"Columns: {list(estimates_df.columns)[:5]} ...")
        continue

    # Align strain indices (take union and fill missing with 0s)
    all_strains = sorted(set(estimates_df.index) | set(truth_df.index))
    estimates_df = estimates_df.reindex(all_strains, fill_value=0)
    truth_df = truth_df.reindex(all_strains, fill_value=0)

    # Compute JSD per sample (matched to ground truth by ratio)
    jsd_results = {}
    for sample, ratio in sample_to_ratio.items():
        if ratio in truth_df.columns:
            # Normalize to probability distributions
            est = estimates_df[sample].clip(lower=0)
            est = est / est.sum() if est.sum() > 0 else est
            true = truth_df[ratio].clip(lower=0)
            true = true / true.sum() if true.sum() > 0 else true

            # Compute JSD (square of Jensen-Shannon distance)
            jsd = jensenshannon(est, true, base=2) ** 2
            jsd_results[sample] = jsd

    # Convert to DataFrame and save
    jsd_df = pd.DataFrame.from_dict(jsd_results, orient="index", columns=["JSD"])
    jsd_df.index.name = "sample"
    out_path = os.path.join(OUT_DIR, f"{species}.csv")
    jsd_df.to_csv(out_path)
    print(f"[\u2713] {species}: {len(jsd_df)} samples -> {out_path}")

print("\n[\u2713] Done.")