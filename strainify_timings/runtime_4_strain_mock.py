# This script is for the 4-strain mock community dataset.

import pandas as pd
import glob
import os
import re

# === CONFIGURATION ===
species_benchdir = "/home/Users/rl152/Strainify_benchmark_V120/not_parallelized/4_strain_mock/benchmarks"  # change path to your benchmarks folder

# Build-stage rules (top-level files + subfolder for rename_fasta_headers)
build_variant_rules = [
    "bwa_index",
    "faidx_ref",
    "filter_variants",
    "finalized_read_counts",
    "maf2vcf",
    "run_parsnp"
]

# Sample-level rules (subfolders)
sample_rules = [
    "count_reads",
    "filter_and_index_sam",
    "map_reads"
]

# Shared rules (top-level files)
shared_rules = ["get_ref", "compute_abundances"]

# === COLUMN CONFIG ===
num_cols = ["s", "max_rss", "max_vms", "max_uss", "max_pss",
            "io_in", "io_out", "cpu_time"]

sum_cols = ["s", "io_in", "io_out", "cpu_time"]
max_cols = ["max_rss", "max_vms", "max_uss", "max_pss"]

# === LOAD SAMPLE RULE FOLDERS ===
bench_records = []
for rule_dir in sample_rules:
    rule_path = os.path.join(species_benchdir, rule_dir)
    if not os.path.isdir(rule_path):
        continue

    tsv_files = glob.glob(os.path.join(rule_path, "*.tsv"))
    for f in tsv_files:
        try:
            df = pd.read_csv(f, sep="\t")

            # Extract sample name
            if "sample" in df.columns:
                sample_name = str(df["sample"].iloc[0])
            else:
                sample_name = os.path.splitext(os.path.basename(f))[0]

            df["rule"] = rule_dir
            df["sample"] = sample_name
            bench_records.append(df)
        except Exception as e:
            print(f"Skipping {f}: {e}")

if not bench_records:
    raise FileNotFoundError(f"No sample benchmark files found in {species_benchdir}")

bench = pd.concat(bench_records, ignore_index=True)
for col in num_cols:
    bench[col] = pd.to_numeric(bench.get(col, 0), errors="coerce").fillna(0)

# ===  Aggregate BUILD VARIANT MATRIX ===
build_records = []

# 1a. rename_fasta_headers (subfolder)
rename_dir = os.path.join(species_benchdir, "rename_fasta_headers")
if os.path.isdir(rename_dir):
    rename_files = glob.glob(os.path.join(rename_dir, "*.tsv"))
    for f in rename_files:
        try:
            df = pd.read_csv(f, sep="\t")
            df["rule"] = "rename_fasta_headers"
            build_records.append(df)
        except Exception as e:
            print(f"Skipping {f}: {e}")

# 1b. other build rules (top-level .tsv files)
for rule in build_variant_rules:
    fpath = os.path.join(species_benchdir, f"{rule}.tsv")
    if os.path.exists(fpath):
        try:
            df = pd.read_csv(fpath, sep="\t")
            df["rule"] = rule
            build_records.append(df)
        except Exception as e:
            print(f"Skipping {fpath}: {e}")

# Combine and aggregate
if build_records:
    build_df = pd.concat(build_records, ignore_index=True)
    for col in num_cols:
        build_df[col] = pd.to_numeric(build_df.get(col, 0), errors="coerce").fillna(0)
    build_sums = {col: build_df[col].sum() for col in sum_cols}
    for col in max_cols:
        build_sums[col] = build_df[col].max()
else:
    build_sums = {col: 0 for col in num_cols}

build_sums["group"] = "build_variant_matrix"

# ===  Aggregate SAMPLE RULES per sample ===
sample_groups = []
for sample, sdf in bench.groupby("sample"):
    sums = {col: sdf[col].sum() for col in sum_cols}
    for col in max_cols:
        sums[col] = sdf[col].max()
    sums["group"] = sample
    sample_groups.append(sums)

sample_sums = pd.DataFrame(sample_groups)

# ===  Add SHARED RULES ===
for shared_rule in shared_rules:
    fpath = os.path.join(species_benchdir, f"{shared_rule}.tsv")
    if not os.path.exists(fpath) or sample_sums.empty:
        continue

    df_shared = pd.read_csv(fpath, sep="\t")
    for col in num_cols:
        df_shared[col] = pd.to_numeric(df_shared.get(col, 0), errors="coerce").fillna(0)

    n_samples = len(sample_sums)
    for col in sum_cols:
        sample_sums[col] += df_shared[col].sum() / n_samples
    for col in max_cols:
        sample_sums[col] = sample_sums[col].combine(df_shared[col].max(), max)

# ===  Combine and export ===
final = pd.concat([pd.DataFrame([build_sums]), sample_sums], ignore_index=True)
final = final[["group"] + num_cols]

outpath = os.path.join(species_benchdir, "pipeline_runtime_summary.csv")
final.to_csv(outpath, index=False)

print(f"\n Summary written to: {outpath}\n")
print(final)
