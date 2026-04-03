#!/usr/bin/env bash
set -euo pipefail

TIMEBIN="/usr/bin/time"

# Please change all paths in this script to match your local setup
THREADS="${THREADS:-12}"
export THREADS

# =========================
# Part 1: build DB
# =========================
bash <<'PART1'
set -euo pipefail

TIMEBIN="/usr/bin/time"

mkdir -p makedb/db_fastq
mkdir -p makedb/db_sam
mkdir -p makedb/db_bam
mkdir -p makedb/counts
mkdir -p makedb/reference_genome_idx
mkdir -p makedb/time_logs

# ---- simulate reads (serial, timed individually) ----
"$TIMEBIN" -v -o makedb/time_logs/wgsim_E24377A.time.log \
wgsim -e 0.0 -d 500 -N 1000000 -1 150 -2 150 -r 0.0 -R 0.0 -X 0.0 \
  genomes/E24377A.fna \
  makedb/db_fastq/E24377A_r1.fastq makedb/db_fastq/E24377A_r2.fastq

"$TIMEBIN" -v -o makedb/time_logs/wgsim_H10407.time.log \
wgsim -e 0.0 -d 500 -N 1000000 -1 150 -2 150 -r 0.0 -R 0.0 -X 0.0 \
  genomes/H10407.fna \
  makedb/db_fastq/H10407_r1.fastq makedb/db_fastq/H10407_r2.fastq

"$TIMEBIN" -v -o makedb/time_logs/wgsim_Sakai.time.log \
wgsim -e 0.0 -d 500 -N 1000000 -1 150 -2 150 -r 0.0 -R 0.0 -X 0.0 \
  genomes/Sakai.fna \
  makedb/db_fastq/Sakai_r1.fastq makedb/db_fastq/Sakai_r2.fastq

"$TIMEBIN" -v -o makedb/time_logs/wgsim_UTI89.time.log \
wgsim -e 0.0 -d 500 -N 1000000 -1 150 -2 150 -r 0.0 -R 0.0 -X 0.0 \
  genomes/UTI89.fna \
  makedb/db_fastq/UTI89_r1.fastq makedb/db_fastq/UTI89_r2.fastq

# ---- compress (serial, timed individually) ----
"$TIMEBIN" -v -o makedb/time_logs/gzip_E24377A.time.log \
gzip -f makedb/db_fastq/E24377A_r1.fastq makedb/db_fastq/E24377A_r2.fastq

"$TIMEBIN" -v -o makedb/time_logs/gzip_H10407.time.log \
gzip -f makedb/db_fastq/H10407_r1.fastq makedb/db_fastq/H10407_r2.fastq

"$TIMEBIN" -v -o makedb/time_logs/gzip_Sakai.time.log \
gzip -f makedb/db_fastq/Sakai_r1.fastq makedb/db_fastq/Sakai_r2.fastq

"$TIMEBIN" -v -o makedb/time_logs/gzip_UTI89.time.log \
gzip -f makedb/db_fastq/UTI89_r1.fastq makedb/db_fastq/UTI89_r2.fastq

# ---- build reference index ----
"$TIMEBIN" -v -o makedb/time_logs/bowtie2_build_Sakai.time.log \
bowtie2-build -q genomes/Sakai.fna makedb/reference_genome_idx/Sakai_idx

# ---- align (serial, timed individually) ----
"$TIMEBIN" -v -o makedb/time_logs/bowtie2_E24377A.time.log \
bowtie2 -X 2000 --no-mixed --dovetail \
  -x makedb/reference_genome_idx/Sakai_idx \
  -1 makedb/db_fastq/E24377A_r1.fastq.gz \
  -2 makedb/db_fastq/E24377A_r2.fastq.gz \
  -S makedb/db_sam/E24377A.sam

"$TIMEBIN" -v -o makedb/time_logs/bowtie2_UTI89.time.log \
bowtie2 -X 2000 --no-mixed --dovetail \
  -x makedb/reference_genome_idx/Sakai_idx \
  -1 makedb/db_fastq/UTI89_r1.fastq.gz \
  -2 makedb/db_fastq/UTI89_r2.fastq.gz \
  -S makedb/db_sam/UTI89.sam

"$TIMEBIN" -v -o makedb/time_logs/bowtie2_H10407.time.log \
bowtie2 -X 2000 --no-mixed --dovetail \
  -x makedb/reference_genome_idx/Sakai_idx \
  -1 makedb/db_fastq/H10407_r1.fastq.gz \
  -2 makedb/db_fastq/H10407_r2.fastq.gz \
  -S makedb/db_sam/H10407.sam

"$TIMEBIN" -v -o makedb/time_logs/bowtie2_Sakai.time.log \
bowtie2 -X 2000 --no-mixed --dovetail \
  -x makedb/reference_genome_idx/Sakai_idx \
  -1 makedb/db_fastq/Sakai_r1.fastq.gz \
  -2 makedb/db_fastq/Sakai_r2.fastq.gz \
  -S makedb/db_sam/Sakai.sam

# ---- SAM -> sorted BAM (serial, timed individually) ----
"$TIMEBIN" -v -o makedb/time_logs/sam2bam_E24377A.time.log bash -c \
'samtools view -bS makedb/db_sam/E24377A.sam | samtools sort - -o makedb/db_bam/E24377A.bam'

"$TIMEBIN" -v -o makedb/time_logs/sam2bam_H10407.time.log bash -c \
'samtools view -bS makedb/db_sam/H10407.sam | samtools sort - -o makedb/db_bam/H10407.bam'

"$TIMEBIN" -v -o makedb/time_logs/sam2bam_UTI89.time.log bash -c \
'samtools view -bS makedb/db_sam/UTI89.sam | samtools sort - -o makedb/db_bam/UTI89.bam'

"$TIMEBIN" -v -o makedb/time_logs/sam2bam_Sakai.time.log bash -c \
'samtools view -bS makedb/db_sam/Sakai.sam | samtools sort - -o makedb/db_bam/Sakai.bam'

# ---- index BAMs (serial, timed individually) ----
"$TIMEBIN" -v -o makedb/time_logs/samtools_index_E24377A.time.log \
samtools index makedb/db_bam/E24377A.bam

"$TIMEBIN" -v -o makedb/time_logs/samtools_index_H10407.time.log \
samtools index makedb/db_bam/H10407.bam

"$TIMEBIN" -v -o makedb/time_logs/samtools_index_UTI89.time.log \
samtools index makedb/db_bam/UTI89.bam

"$TIMEBIN" -v -o makedb/time_logs/samtools_index_Sakai.time.log \
samtools index makedb/db_bam/Sakai.bam

# ---- counts (serial, timed individually) ----
"$TIMEBIN" -v -o makedb/time_logs/phlame_counts_E24377A.time.log \
phlame counts -i makedb/db_bam/E24377A.bam -r genomes/Sakai.fna -o makedb/counts/E24377A.counts

"$TIMEBIN" -v -o makedb/time_logs/phlame_counts_H10407.time.log \
phlame counts -i makedb/db_bam/H10407.bam -r genomes/Sakai.fna -o makedb/counts/H10407.counts

"$TIMEBIN" -v -o makedb/time_logs/phlame_counts_Sakai.time.log \
phlame counts -i makedb/db_bam/Sakai.bam -r genomes/Sakai.fna -o makedb/counts/Sakai.counts

"$TIMEBIN" -v -o makedb/time_logs/phlame_counts_UTI89.time.log \
phlame counts -i makedb/db_bam/UTI89.bam -r genomes/Sakai.fna -o makedb/counts/UTI89.counts

# ---- downstream DB creation steps, serial and timed ----
"$TIMEBIN" -v -o makedb/time_logs/phlame_cmt.time.log \
phlame cmt -i required_files/counts_files.txt -s required_files/sample_names.txt -r genomes/Sakai.fna -o makedb/ecoli.pickle.gz

"$TIMEBIN" -v -o makedb/time_logs/phlame_tree.time.log \
phlame tree -i makedb/ecoli.pickle.gz -p makedb/ecoli.phylip -r makedb/ecoli_phylip2names.txt -o makedb/ecoli.tre --rescale

"$TIMEBIN" -v -o makedb/time_logs/phlame_makedb.time.log \
phlame makedb -i makedb/ecoli.pickle.gz -t makedb/rescaled_ecoli.tre -o makedb/ecoli_db.classifier -p makedb/ecoli_cladeIDs.txt --min_branchlen 0 --min_leaves 1 --min_snps 1 -c required_files/ecoli_manual_clade_ID.txt
PART1

# =========================
# Part 2: classify
# =========================
"$TIMEBIN" -v -o part2_classify.time.log bash <<'PART2'
set -euo pipefail
shopt -s nullglob

FASTQ_DIR="/home/Users/rl152/Strainify/4-strain-mock-exp/fastq"
OUT_DIR="/home/Users/rl152/Strainify/phlame/4_strain_mock_community/benchmarks/results"
IDX="/home/Users/rl152/Strainify/phlame/4_strain_mock_community/benchmarks/makedb/reference_genome_idx/Sakai_idx"
CLASSIFIER="/home/Users/rl152/Strainify/phlame/4_strain_mock_community/benchmarks/makedb/ecoli_db.classifier"
REF_FASTA="/home/Users/rl152/Strainify/phlame/4_strain_mock_community/benchmarks/genomes/Sakai.fna"
THREADS="${THREADS:-12}"

mkdir -p "$OUT_DIR"

for r1 in "$FASTQ_DIR"/*_r1.fq; do
  base="$(basename "$r1" _r1.fq)"
  r2="$FASTQ_DIR/${base}_r2.fq"

  if [[ ! -f "$r2" ]]; then
    echo "[WARN] Missing R2 for: $r1 (expected: $r2). Skipping."
    continue
  fi

  sam="$OUT_DIR/${base}.sam"
  bam="$OUT_DIR/${base}.bam"
  freq_csv="$OUT_DIR/${base}_frequencies.csv"
  fitinfo="$OUT_DIR/${base}_fitinfo.data"

  echo "==== Processing: $base ===="

  bowtie2 -p "$THREADS" -X 2000 --no-mixed --dovetail \
    -x "$IDX" -1 "$r1" -2 "$r2" -S "$sam"

  samtools view -@ "$THREADS" -bS "$sam" \
    | samtools sort -@ "$THREADS" -o "$bam" -

  samtools index -@ "$THREADS" "$bam"

  rm -f "$sam"

  phlame classify -i "$bam" -c "$CLASSIFIER" -r "$REF_FASTA" -m bayesian \
    -o "$freq_csv" -p "$fitinfo"

  echo "==== Done: $base ===="
done
PART2

echo "Finished."
echo "Part 1 per-step timing logs: makedb/time_logs/"
echo "Part 2 total timing: part2_classify.time.log"