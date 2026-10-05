import os
import re
import pandas as pd
import numpy as np
from scipy.spatial.distance import jensenshannon

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

species_list = [
    "cacnes",
    "cdiff",
    "ecoli",
    "mtuberculosis",
    "sepidermidis",
]

pantax_dir = "/dodo/rl152/pantax/30_strains/pantax_abundance_matrices"
truth_base_dir = "/home/Users/rl152/Strainify/30_strains"

out_dir = "/dodo/rl152/pantax/30_strains/jsd"
os.makedirs(out_dir, exist_ok=True)


# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

def normalize_index(df):
    """
    Normalize genome IDs so PanTax and truth tables can be aligned.

    Example:
        GCA_000145195.1_ASM14519v1_genomic -> GCA_000145195
        GCA_000145195.1_ASM14519v1_genomic.fna -> GCA_000145195
    """
    def clean(name):
        name = str(name)

        # Remove common genome file extensions
        name = re.sub(r'\.fna$|\.fa$|\.fasta$|\.gz$', '', name)

        # Keep accession before first dot:
        # GCA_000145195.1_ASM14519v1_genomic -> GCA_000145195
        name = name.split('.')[0]

        # Keep first two underscore-separated parts:
        # GCA_000145195 -> GCA_000145195
        parts = name.split('_')
        return '_'.join(parts[:2]) if len(parts) >= 2 else name

    df = df.copy()
    df.index = df.index.to_series().apply(clean)
    return df


def load_pantax_matrix(path):
    """
    Load PanTax relative abundance matrix.

    PanTax uses 'undetected' for missing strains.
    For JSD, individual 'undetected' values are converted to zero.
    A separate mask is kept so samples with all strains undetected can be labeled
    as 'no_strains_detected' instead of calculating JSD.
    """
    raw_df = pd.read_csv(path, index_col=0)

    undetected_mask = raw_df.astype(str).apply(
        lambda col: col.str.lower().eq("undetected")
    )

    numeric_df = raw_df.replace("undetected", 0)
    numeric_df = numeric_df.apply(pd.to_numeric, errors="coerce").fillna(0)

    return numeric_df, undetected_mask


def get_sample_to_ratio(columns):
    """
    Build mapping from PanTax sample names to truth ratio columns.

    Examples:
        100x_ratio_1 -> ratio_1
        10x_ratio_2  -> ratio_2
        200x_ratio_3 -> ratio_3
    """
    sample_to_ratio = {}

    for col in columns:
        match = re.search(r'(ratio_\d+)', col)
        if match:
            sample_to_ratio[col] = match.group(1)

    return sample_to_ratio


def sample_sort_key(sample_name):
    """
    Sort samples by coverage and ratio number.

    Example order:
        10x_ratio_1, 10x_ratio_2, 10x_ratio_3,
        20x_ratio_1, ...
    """
    match = re.match(r'(\d+)x_ratio_(\d+)', sample_name)
    if match:
        cov, ratio_num = match.groups()
        return int(cov), int(ratio_num)

    return 999999, 999999


def compute_jsd_for_species(species):
    print()
    print("============================================================")
    print(f"[INFO] Processing species: {species}")
    print("============================================================")

    # PanTax estimate file for this species
    pantax_file = os.path.join(
        pantax_dir,
        f"{species}_pantax_relative_abundances.csv"
    )

    # Ground-truth ratio file for this species
    truth_file = os.path.join(
        truth_base_dir,
        species,
        f"{species}_30_ratios.csv"
    )

    if not os.path.exists(pantax_file):
        raise FileNotFoundError(f"Missing PanTax file: {pantax_file}")

    if not os.path.exists(truth_file):
        raise FileNotFoundError(f"Missing truth file: {truth_file}")

    # --------------------------------------------------------
    # Load data
    # --------------------------------------------------------

    estimates_df, undetected_mask = load_pantax_matrix(pantax_file)
    truth_df = pd.read_csv(truth_file, index_col=0)

    # --------------------------------------------------------
    # Normalize genome IDs
    # --------------------------------------------------------

    estimates_df = normalize_index(estimates_df)
    undetected_mask = normalize_index(undetected_mask)
    truth_df = normalize_index(truth_df)

    # If normalization creates duplicate IDs, collapse them.
    estimates_df = estimates_df.groupby(estimates_df.index).sum()

    # For the undetected mask, a normalized strain is considered undetected
    # only if all duplicate rows were undetected.
    undetected_mask = undetected_mask.groupby(undetected_mask.index).all()

    truth_df = truth_df.groupby(truth_df.index).sum()

    # --------------------------------------------------------
    # Match PanTax samples to truth ratio columns
    # --------------------------------------------------------

    sample_to_ratio = get_sample_to_ratio(estimates_df.columns)

    if not sample_to_ratio:
        raise ValueError(f"No ratio columns found in PanTax file: {pantax_file}")

    # --------------------------------------------------------
    # Align strains by union
    # --------------------------------------------------------
    # This guarantees estimate and truth distributions have the same shape.
    # Strains that are missing from either table are filled with zero.

    all_strains = sorted(set(estimates_df.index) | set(truth_df.index))

    estimates_df = estimates_df.reindex(all_strains, fill_value=0)
    truth_df = truth_df.reindex(all_strains, fill_value=0)

    # If a strain is not present in the original PanTax table, treat it as undetected.
    undetected_mask = undetected_mask.reindex(all_strains, fill_value=True)

    print(f"[INFO] Total aligned strains: {len(all_strains)}")
    print(f"[INFO] Samples found: {len(sample_to_ratio)}")

    # --------------------------------------------------------
    # Compute JSD per sample
    # --------------------------------------------------------

    jsd_results = []

    for sample, ratio in sample_to_ratio.items():
        if ratio not in truth_df.columns:
            print(
                f"[WARNING] Skipping {sample}: "
                f"truth column {ratio} not found in {truth_file}"
            )
            continue

        est_raw = estimates_df[sample].clip(lower=0)
        true_raw = truth_df[ratio].clip(lower=0)

        est_sum = est_raw.sum()
        true_sum = true_raw.sum()

        num_estimated_nonzero = int((est_raw > 0).sum())
        num_truth_nonzero = int((true_raw > 0).sum())

        # True when PanTax reported "undetected" for every strain in this sample.
        all_pantax_undetected = bool(undetected_mask[sample].all())

        if all_pantax_undetected or est_sum == 0:
            jsd_value = "no_strains_detected"
        elif true_sum == 0:
            jsd_value = np.nan
        else:
            # Normalize to probability distributions after shape alignment.
            est = est_raw / est_sum
            true = true_raw / true_sum

            # scipy.spatial.distance.jensenshannon returns JS distance.
            # Squaring gives Jensen-Shannon divergence.
            jsd_value = jensenshannon(est, true, base=2) ** 2

        jsd_results.append({
            "sample": sample,
            "truth_ratio": ratio,
            "JSD": jsd_value,
            "num_estimated_nonzero": num_estimated_nonzero,
            "num_truth_nonzero": num_truth_nonzero,
            "all_pantax_undetected": all_pantax_undetected,
        })

    jsd_df = pd.DataFrame(jsd_results)

    jsd_df = jsd_df.sort_values(
        by="sample",
        key=lambda s: s.map(sample_sort_key)
    )

    species_outfile = os.path.join(
        out_dir,
        f"{species}_pantax.csv"
    )

    jsd_df.to_csv(species_outfile, index=False)

    print(f"[INFO] Saved: {species_outfile}")


# ------------------------------------------------------------
# Run all species
# ------------------------------------------------------------

for species in species_list:
    compute_jsd_for_species(species)

print()
print("============================================================")
print("[INFO] Finished all species")
print(f"[INFO] Output directory: {out_dir}")
print("============================================================")