#!/usr/bin/env bash
set -euo pipefail

# Please change the paths to your local setup
INPUT="/home/Users/rl152/Strainify/UMB/acc_list_unique.txt"
OUTDIR="/home/Users/rl152/Strainify/UMB/UMB18_genomes"
THREADS=32

mkdir -p "$OUTDIR"

echo "Extracting GCF accessions from column 2 (whitespace-delimited)..."

awk '{
  gsub(/\r/,"",$2)
  if ($2 ~ /^GCF_[0-9]+\.[0-9]+$/) print $2
}' "$INPUT" | sort -u > gcf_accessions.txt

echo "Unique GCF assemblies found: $(wc -l < gcf_accessions.txt)"
head -n 5 gcf_accessions.txt || true

if [[ ! -s gcf_accessions.txt ]]; then
  echo "ERROR: No valid GCF accessions found in column 2."
  exit 1
fi

echo "Downloading RefSeq genomes..."

ncbi-genome-download bacteria \
  --section refseq \
  --assembly-accessions gcf_accessions.txt \
  --formats fasta \
  --output-folder "$OUTDIR" \
  --parallel "$THREADS" \
  --flat-output \
  -v

echo "Done. Genomes saved to: $OUTDIR"

