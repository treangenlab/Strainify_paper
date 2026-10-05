#!/usr/bin/env bash
set -euo pipefail

# Please change the paths to your local setup
FASTQ_DIR="/dodo/rl152/UMB/PRJNA400628_UMB18/fastq"
CSV="/home/Users/rl152/Strainify/UMB/UMB18_stool_SRRs.csv"
OUTDIR="/dodo/rl152/UMB/PRJNA400628_UMB18/UMB18_selected_fastqs"

mkdir -p "$OUTDIR"

COL=$(head -n 1 "$CSV" | tr ',' '\n' | nl -v1 | grep -w "sample_name" | awk '{print $1}')

awk -F, -v col="$COL" 'NR>1{
  gsub(/"/,"",$col); gsub(/\r/,"",$col)
  if ($col!="") print $col
}' "$CSV" | sort -u > sample_names_unique.txt

# Build a regex pattern from all sample names
PATTERN=$(paste -sd'|' sample_names_unique.txt)

find "$FASTQ_DIR" -type f | grep -E "$PATTERN" | while read -r f; do
  cp -n "$f" "$OUTDIR/"
done

echo "Done."
