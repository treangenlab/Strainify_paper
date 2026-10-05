#!/usr/bin/env python3

import os
import glob
import pandas as pd

# -----------------------------
# Paths
# -----------------------------
base_dir = "/dodo/rl152/pantax/30_strains"

genome_info_dir = os.path.join(base_dir, "genome_info")
pantax_results_dir = os.path.join(base_dir, "pantax_results")

out_dir = os.path.join(base_dir, "pantax_abundance_matrices")
os.makedirs(out_dir, exist_ok=True)

species_list = [
    "cacnes",
    "cdiff",
    "ecoli",
    "mtuberculosis",
    "sepidermidis",
]

# -----------------------------
# Helper: read expected strains
# -----------------------------
def read_expected_strains(genome_info_file):
    """
    genome_info file:
    first column = genome_ID / strain name
    Handles files with or without a header.
    """
    df = pd.read_csv(genome_info_file, sep="\t", header=None, comment="#")

    if df.empty:
        raise ValueError(f"Genome info file is empty: {genome_info_file}")

    strains = df.iloc[:, 0].astype(str).tolist()

    # Remove header if present
    strains = [s for s in strains if s != "genome_ID"]

    return strains


# -----------------------------
# Main
# -----------------------------
for species in species_list:
    print()
    print("============================================================")
    print(f"[INFO] Processing species: {species}")
    print("============================================================")

    genome_info_file = os.path.join(genome_info_dir, f"{species}_genomes_info.txt")
    species_results_dir = os.path.join(pantax_results_dir, species)
    out_csv = os.path.join(out_dir, f"{species}_pantax_relative_abundances.csv")

    if not os.path.isfile(genome_info_file):
        print(f"[WARN] Missing genome info file, skipping: {genome_info_file}")
        continue

    if not os.path.isdir(species_results_dir):
        print(f"[WARN] Missing Pantax results directory, skipping: {species_results_dir}")
        continue

    expected_strains = read_expected_strains(genome_info_file)

    print(f"[INFO] Expected strains: {len(expected_strains)}")

    files = sorted(
        glob.glob(os.path.join(species_results_dir, "*_strains_abundance.txt"))
    )

    if not files:
        print(f"[WARN] No *_strains_abundance.txt files found for {species}")
        continue

    combined = pd.DataFrame(index=expected_strains)

    for f in files:
        fname = os.path.basename(f)

        # Example:
        # 10x_ratio_1_interleaved_strains_abundance.txt
        # -> 10x_ratio_1
        sample_name = fname.replace("_interleaved_strains_abundance.txt", "")

        print(f"[INFO] Reading {species} / {sample_name}")

        # Default: all expected strains are undetected
        sample_abund = {strain: "undetected" for strain in expected_strains}

        df = pd.read_csv(f, sep="\t")

        if df.empty:
            print(f"[WARN] Empty strain abundance file: {f}")
            combined[sample_name] = [sample_abund[s] for s in expected_strains]
            continue

        required_cols = {"genome_ID", "predicted_abundance"}
        missing_cols = required_cols - set(df.columns)

        if missing_cols:
            raise ValueError(f"{f} is missing required columns: {missing_cols}")

        df["genome_ID"] = df["genome_ID"].astype(str)
        df["predicted_abundance"] = pd.to_numeric(
            df["predicted_abundance"], errors="coerce"
        )

        df = df.dropna(subset=["genome_ID", "predicted_abundance"])

        # Keep only strains expected for this species
        df = df[df["genome_ID"].isin(expected_strains)]

        if df.empty:
            print(f"[WARN] No expected strains detected in {f}")
            combined[sample_name] = [sample_abund[s] for s in expected_strains]
            continue

        # If Pantax somehow reports the same genome more than once, combine rows
        df = (
            df.groupby("genome_ID", as_index=False)["predicted_abundance"]
            .sum()
        )

        # Normalize detected strains only
        total = df["predicted_abundance"].sum()

        if total > 0:
            df["normalized_abundance"] = df["predicted_abundance"] / total
        else:
            print(f"[WARN] Detected abundances sum to 0 in {f}")
            df["normalized_abundance"] = 0.0

        # Fill detected strains; missing expected strains remain "undetected"
        for _, row in df.iterrows():
            sample_abund[row["genome_ID"]] = row["normalized_abundance"]

        combined[sample_name] = [sample_abund[s] for s in expected_strains]

    combined.index.name = "strain"
    combined.to_csv(out_csv)

    print(f"[DONE] Wrote: {out_csv}")