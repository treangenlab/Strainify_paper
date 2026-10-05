#!/usr/bin/env python3
"""
Combine StrainGE /usr/bin/time -v logs into one timing CSV per species.

Reads the benchmark tree written by run_straingst_inference.sh:

    <bench_root>/<dataset>/results/<species>/<sample>/timing/*.time

For each sample it sums the wall-clock time across that sample's steps
(read kmerize + straingst run) and takes the peak RSS (max across steps),
then writes ONE CSV per (dataset, species):

    <out_dir>/<dataset>__<species>_strainge_timing.csv

Rows are the two metrics (wall clock time, peak RAM); columns are sample names
(e.g. "10x_ratio_1"). Single-species datasets yield one CSV; the 30-strain
dataset yields one CSV per species.

It also reads the DB-build logs:

    <bench_root>/<dataset>/db/<species>/timing/*.time   (genome kmerize + createdb)

and writes one build CSV per dataset (columns = species):

    <out_dir>/<dataset>_strainge_build_timing.csv

Build wall clock is the SUM of all build steps; since genome kmerize may have
run in parallel, treat it as total compute time (an upper bound on real
elapsed). Build peak RAM is the max across steps.

Usage:
    python3 combine_strainge_timing.py
    python3 combine_strainge_timing.py --bench-root /path/to/StrainGE_benchmarks \
                                       --out-dir /path/to/StrainGE_benchmarks/timing_summary \
                                       --ram-unit gb
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

DEFAULT_BENCH_ROOT = "/home/Users/rl152/Strainify_paper/StrainGE" # Please change this to your local path

_ELAPSED_RE = re.compile(r"Elapsed \(wall clock\) time \([^)]*\):\s*([\d:.]+)")
_MAXRSS_RE = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")

RAM_DIVISOR = {"kb": 1, "mb": 1024, "gb": 1024 * 1024}


def parse_elapsed(val: str) -> float:
    """'h:mm:ss', 'm:ss', 'm:ss.ss' or plain seconds -> float seconds."""
    parts = [float(p) for p in val.split(":")]
    if len(parts) == 3:
        h, m, s = parts
        return h * 3600 + m * 60 + s
    if len(parts) == 2:
        m, s = parts
        return m * 60 + s
    return parts[0]


def parse_time_file(path: Path):
    """Return (wall_seconds, max_rss_kb) from one /usr/bin/time -v log."""
    text = path.read_text(errors="ignore")
    wall = _ELAPSED_RE.search(text)
    rss = _MAXRSS_RE.search(text)
    wall_s = parse_elapsed(wall.group(1)) if wall else None
    rss_kb = int(rss.group(1)) if rss else None
    return wall_s, rss_kb


def sample_sort_key(sample: str):
    """Sort by numeric coverage, then numeric ratio, then name."""
    cov = re.search(r"(\d+)x", sample)
    ratio = re.search(r"ratio_(\d+)", sample)
    return (int(cov.group(1)) if cov else 1_000_000,
            int(ratio.group(1)) if ratio else 1_000_000,
            sample)


def collect_species(species_dir: Path) -> dict:
    """{ '<sample>': (wall_sec, peak_rss_kb) } for one species."""
    columns = {}
    for sample_dir in sorted(p for p in species_dir.iterdir() if p.is_dir()):
        timing_dir = sample_dir / "timing"
        if not timing_dir.is_dir():
            continue
        time_files = sorted(timing_dir.glob("*.time"))
        if not time_files:
            continue

        wall_total = 0.0
        peak_kb = 0
        found = False
        for tf in time_files:
            wall_s, rss_kb = parse_time_file(tf)
            if wall_s is not None:
                wall_total += wall_s
                found = True
            if rss_kb is not None:
                peak_kb = max(peak_kb, rss_kb)
                found = True
        if not found:
            print(f"[! Warning] No parsable metrics in {timing_dir}", file=sys.stderr)
            continue

        columns[sample_dir.name] = (wall_total, peak_kb)
    return columns


def collect_build(db_root: Path) -> dict:
    """{ '<species>': (wall_sum_sec, peak_rss_kb) } for a dataset's DB build.

    Wall clock is the SUM of all build steps (genome kmerize + createdb). Note
    the genome-kmerize steps may have been run in parallel, so this sum is total
    compute time and an upper bound on real elapsed. Peak RAM is the max across
    steps.
    """
    columns = {}
    for species_dir in sorted(p for p in db_root.iterdir() if p.is_dir()):
        timing_dir = species_dir / "timing"
        if not timing_dir.is_dir():
            continue
        time_files = sorted(timing_dir.glob("*.time"))
        if not time_files:
            continue

        wall_total = 0.0
        peak_kb = 0
        found = False
        for tf in time_files:
            wall_s, rss_kb = parse_time_file(tf)
            if wall_s is not None:
                wall_total += wall_s
                found = True
            if rss_kb is not None:
                peak_kb = max(peak_kb, rss_kb)
                found = True
        if not found:
            print(f"[! Warning] No parsable metrics in {timing_dir}", file=sys.stderr)
            continue

        columns[species_dir.name] = (wall_total, peak_kb)
    return columns


def make_metrics_df(cols: dict, ram_row: str, divisor: int, sort_key=None):
    """Build a 2-row (wall clock, peak RAM) DataFrame from {name: (wall, rss_kb)}."""
    data = {
        name: {
            "wall_clock_time_sec": wall,
            ram_row: round(rss_kb / divisor, 4),
        }
        for name, (wall, rss_kb) in cols.items()
    }
    df = pd.DataFrame(data)
    ordered = sorted(df.columns, key=sort_key) if sort_key else sorted(df.columns)
    df = df[ordered]
    df = df.reindex(["wall_clock_time_sec", ram_row])
    df.index.name = "metric"
    return df


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bench-root", default=DEFAULT_BENCH_ROOT,
                    help="Root containing <dataset>/results/... (default: %(default)s)")
    ap.add_argument("--datasets", nargs="+", default=None,
                    help="Datasets to process (default: auto-detect all with a results/ dir)")
    ap.add_argument("--out-dir", default=None,
                    help="Where to write the CSVs (default: <bench_root>/timing_summary)")
    ap.add_argument("--ram-unit", choices=["kb", "mb", "gb"], default="gb",
                    help="Unit for peak RAM row (default: gb)")
    args = ap.parse_args()

    bench_root = Path(args.bench_root)
    out_dir = Path(args.out_dir) if args.out_dir else bench_root / "timing_summary"
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.datasets:
        datasets = args.datasets
    else:
        datasets = sorted(p.name for p in bench_root.iterdir()
                          if p.is_dir() and (p / "results").is_dir())
    if not datasets:
        print(f"[! Error] No datasets with a results/ dir under {bench_root}", file=sys.stderr)
        return 1

    divisor = RAM_DIVISOR[args.ram_unit]
    ram_row = f"peak_ram_{args.ram_unit}"
    wrote_any = False

    for dataset in datasets:
        # --- Inference: one CSV per species (columns = sample names) ---
        results_root = bench_root / dataset / "results"
        if results_root.is_dir():
            for species_dir in sorted(p for p in results_root.iterdir() if p.is_dir()):
                species = species_dir.name
                cols = collect_species(species_dir)
                if not cols:
                    print(f"[! Warning] No inference timing for {dataset}/{species}.",
                          file=sys.stderr)
                    continue
                df = make_metrics_df(cols, ram_row, divisor, sample_sort_key)
                out_path = out_dir / f"{dataset}__{species}_strainge_timing.csv"
                df.to_csv(out_path)
                print(f"[\u2713] inference {dataset}/{species}: {df.shape[1]} samples -> {out_path}")
                wrote_any = True
        else:
            print(f"[! Warning] No results/ for dataset {dataset}.", file=sys.stderr)

        # --- Build: one CSV per dataset (columns = species) ---
        db_root = bench_root / dataset / "db"
        if db_root.is_dir():
            build_cols = collect_build(db_root)
            if build_cols:
                df = make_metrics_df(build_cols, ram_row, divisor)  # species sorted alphabetically
                out_path = out_dir / f"{dataset}_strainge_build_timing.csv"
                df.to_csv(out_path)
                print(f"[\u2713] build     {dataset}: {df.shape[1]} species -> {out_path}")
                wrote_any = True
            else:
                print(f"[! Warning] No build timing for {dataset}.", file=sys.stderr)
        else:
            print(f"[! Warning] No db/ for dataset {dataset}.", file=sys.stderr)

    if not wrote_any:
        print("[! Error] Nothing written.", file=sys.stderr)
        return 1

    print(f"\n[\u2713] CSVs written to: {out_dir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())