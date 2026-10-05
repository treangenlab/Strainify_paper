#!/usr/bin/env python3
"""
Combine StrainGST relative-abundance (rapct) results into per-species CSVs.

Walks the results tree produced by run_straingst_inference.sh:

    <bench_root>/<dataset>/results/<species>/<sample>/result_strains.tsv

For every (dataset, species) it builds a strain x sample matrix from the
`rapct` column (RAW values, no normalization) and writes ONE CSV per
(dataset, species):

    <out_dir>/<dataset>__<species>.csv

Normalization is deliberately left to the downstream plotting step so that all
tools can be put on a common scale together under one uniform rule.

So 4_strain_simulated (single species) yields one CSV, and
30_strain_simulated yields one CSV per species (ecoli, cacnes, cdiff,
mtuberculosis, sepidermidis). The mock dataset is excluded by default.

Rows = strains, columns = samples (read sets). Column headers are the sample
folder names (unique within a species).

Usage:
    python3 extract_results.py
    python3 extract_results.py --bench-root /path/to/strainge_benchmark \
                             --datasets 4_strain_simulated 30_strain_simulated \
                             --out-dir combined_rapct
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

DEFAULT_BENCH_ROOT = "/home/Users/rl152/Strainify_paper/StrainGE" # Please change this to your local path
DEFAULT_DATASETS = ["4_strain_simulated", "30_strain_simulated", "4_strain_1-2-5x"]  # mock excluded


def find_strains_tsv(sample_dir: Path):
    """Return the TSV in sample_dir that holds per-strain results (has rapct)."""
    candidates = sorted(sample_dir.glob("*.tsv"))
    candidates.sort(key=lambda p: (0 if "strain" in p.name.lower() else 1, p.name))
    for p in candidates:
        try:
            header = pd.read_csv(p, sep="\t", nrows=0).columns
        except Exception:
            continue
        if "rapct" in header and "strain" in header:
            return p
    return None


def load_sample_rapct(tsv: Path) -> pd.Series:
    """strain -> rapct series for one sample (duplicate strains summed)."""
    df = pd.read_csv(tsv, sep="\t", usecols=lambda c: c in ("strain", "rapct"))
    df = df.dropna(subset=["strain"])
    return df.groupby("strain")["rapct"].sum()


def write_species_csv(dataset: str, species: str, columns: dict, out_dir: Path):
    """Build and write one CSV of RAW rapct (no normalization).

    Normalization is intentionally left to the downstream plotting step so all
    tools can be harmonized together under one uniform rule.
    """
    matrix = pd.DataFrame(columns).sort_index().fillna(0.0)
    matrix.index.name = "strain"

    out_path = out_dir / f"{dataset}__{species}.csv"
    matrix.to_csv(out_path)

    col_sums = matrix.sum(axis=0)
    print(f"[\u2713] {out_path.name}: {matrix.shape[0]} strains x {matrix.shape[1]} samples "
          f"(raw rapct; per-sample sums {col_sums.min():.3f}-{col_sums.max():.3f})")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bench-root", default=DEFAULT_BENCH_ROOT,
                    help="Root containing <dataset>/results/... (default: %(default)s)")
    ap.add_argument("--datasets", nargs="+", default=DEFAULT_DATASETS,
                    help="Datasets to include (default: %(default)s)")
    ap.add_argument("--out-dir", default="combined_rapct",
                    help="Directory for the per-species CSVs (default: %(default)s)")
    args = ap.parse_args()

    bench_root = Path(args.bench_root)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    wrote_any = False

    for dataset in args.datasets:
        results_root = bench_root / dataset / "results"
        if not results_root.is_dir():
            print(f"[! Warning] Missing results dir: {results_root}", file=sys.stderr)
            continue

        for species_dir in sorted(p for p in results_root.iterdir() if p.is_dir()):
            columns: dict[str, pd.Series] = {}
            for sample_dir in sorted(p for p in species_dir.iterdir() if p.is_dir()):
                tsv = find_strains_tsv(sample_dir)
                if tsv is None:
                    print(f"[! Warning] No strains TSV in {sample_dir}", file=sys.stderr)
                    continue
                columns[sample_dir.name] = load_sample_rapct(tsv)

            if not columns:
                print(f"[! Warning] No samples for {dataset}/{species_dir.name}", file=sys.stderr)
                continue

            write_species_csv(dataset, species_dir.name, columns, out_dir)
            wrote_any = True

    if not wrote_any:
        print("[! Error] Nothing written. Check --bench-root / --datasets.", file=sys.stderr)
        return 1

    print(f"\n[\u2713] CSVs written to: {out_dir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())