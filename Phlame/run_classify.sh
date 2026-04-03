#!/usr/bin/env bash
set -euo pipefail

# -------- Paths / settings --------
FASTQ_DIR="/home/Users/rl152/Strainify/4-strain-mock-exp/fastq"
OUT_DIR="/home/Users/rl152/Strainify/phlame/4_strain_mock_community/results_default" 
IDX="/home/Users/rl152/Strainify/phlame/reference_genome_idx/Sakai_idx"
CLASSIFIER="/home/Users/rl152/Strainify/phlame/makedb/ecoli_db.classifier"
REF_FASTA="/home/Users/rl152/Strainify/phlame/genomes/Sakai.fna"

THREADS="${THREADS:-12}"   # override like: THREADS=24 ./run_all.sh

mkdir -p "$OUT_DIR"

shopt -s nullglob

# Loop over *_r1.fq and find the matching *_r2.fq
for r1 in "$FASTQ_DIR"/*_r1.fq; do
  base="$(basename "$r1" _r1.fq)"
  r2="$FASTQ_DIR/${base}_r2.fq"

  if [[ ! -f "$r2" ]]; then
    echo "[WARN] Missing R2 for: $r1  (expected: $r2). Skipping."
    continue
  fi

  sam="$OUT_DIR/${base}.sam"
  bam="$OUT_DIR/${base}.bam"
  freq_csv="$OUT_DIR/${base}_frequencies.csv"
  fitinfo="$OUT_DIR/${base}_fitinfo.data"

  echo "==== Processing: $base ===="

  # 1) Align
  bowtie2 -p "$THREADS" -X 2000 --no-mixed --dovetail \
    -x "$IDX" -1 "$r1" -2 "$r2" -S "$sam"

  # 2) SAM -> sorted BAM
  samtools view -@ "$THREADS" -bS "$sam" \
    | samtools sort -@ "$THREADS" -o "$bam" -

  # 3) Index BAM
  samtools index -@ "$THREADS" "$bam"

  # (optional) remove SAM to save space
  rm -f "$sam"

  # 4) PHLAME classify
  phlame classify -i "$bam" -c "$CLASSIFIER" -r "$REF_FASTA" -m bayesian \
    -o "$freq_csv" -p "$fitinfo"

  echo "==== Done: $base ===="
done

echo "All done."
