#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------

# Please change these paths to your local setup.

DATA_DIR="/dodo/rl152/ChronoStrain/chronostrain_cami/cami_dataset"

# Folder containing one subfolder per species, each with <species>.hdf5
PER_SPECIES_DB_DIR="/dodo/rl152/ChronoStrain/chronostrain_cami/core_genome_only/strainge/db_per_species"

# Output directory
OUT_BASE="/dodo/rl152/ChronoStrain/chronostrain_cami/core_genome_only/strainge/results"

# Samples to run
sample_start=0
sample_end=99

# Number of samples to run in parallel
# Use one CPU/core per sample.
max_parallel_samples=24

# StrainGST parameters
num_iters=408
min_score=0

# ------------------------------------------------------------
# Run StrainGST for one sample against one species DB
# ------------------------------------------------------------

run_straingst_species() {
  local sample_id=$1
  local species_name=$2
  local db_file=$3

  local sample_outdir="${OUT_BASE}/sample_${sample_id}"
  local outdir="${sample_outdir}/${species_name}"

  local run_breadcrumb="${outdir}/strainGST.DONE"

  if [ -f "${run_breadcrumb}" ]; then
    echo "[! straingst_per_species] ${species_name}, sample_${sample_id} already done."
    return 0
  fi

  local fq1_paired="${DATA_DIR}/reads/extracted/sample_${sample_id}/1_paired.fq.gz"
  local fq2_paired="${DATA_DIR}/reads/extracted/sample_${sample_id}/2_paired.fq.gz"

  if [ ! -f "${fq1_paired}" ] || [ ! -f "${fq2_paired}" ]; then
    echo "[! Error] Read input files not found for sample_${sample_id}. Skipping."
    return 0
  fi

  mkdir -p "${outdir}"

  # Kmerize reads once per sample and share across species DB runs.
  local read_kmer_dir="${sample_outdir}/read_kmers"
  local read_kmers="${read_kmer_dir}/reads.hdf5"

  mkdir -p "${read_kmer_dir}"

  if [ -f "${read_kmers}" ]; then
    echo "[! straingst_per_species] Read kmerization already done for sample_${sample_id}."
  else
    echo "[! straingst_per_species] Kmerizing reads for sample_${sample_id}."
    echo "[! straingst_per_species] Read_1: ${fq1_paired}"
    echo "[! straingst_per_species] Read_2: ${fq2_paired}"

    straingst kmerize \
      -k 23 \
      -o "${read_kmers}" \
      "${fq1_paired}" "${fq2_paired}"
  fi

  local results_tsv="${outdir}/result.tsv"

  echo "[! straingst_per_species] Running StrainGST."
  echo "[! straingst_per_species] Sample: sample_${sample_id}"
  echo "[! straingst_per_species] Species: ${species_name}"
  echo "[! straingst_per_species] DB: ${db_file}"

  straingst run \
    -o "${results_tsv}" \
    -i "${num_iters}" \
    "${db_file}" "${read_kmers}" \
    -s "${min_score}" \
    --separate-output

  touch "${run_breadcrumb}"
}

# ------------------------------------------------------------
# Run all species DBs for one sample
# ------------------------------------------------------------

run_one_sample() {
  local sample_id=$1

  echo
  echo "============================================================"
  echo "[* straingst_per_species] Handling sample_${sample_id}"
  echo "============================================================"

  local fq1_paired="${DATA_DIR}/reads/extracted/sample_${sample_id}/1_paired.fq.gz"
  local fq2_paired="${DATA_DIR}/reads/extracted/sample_${sample_id}/2_paired.fq.gz"

  if [ ! -f "${fq1_paired}" ] || [ ! -f "${fq2_paired}" ]; then
    echo "[! Error] Read input files not found for sample_${sample_id}. Skipping whole sample."
    return 0
  fi

  for species_dir in "${PER_SPECIES_DB_DIR}"/*; do
    [ -d "${species_dir}" ] || continue

    local species_name
    species_name=$(basename "${species_dir}")

    local db_file="${species_dir}/${species_name}.hdf5"

    if [ ! -f "${db_file}" ]; then
      echo "[! Warning] DB file not found for ${species_name}: ${db_file}. Skipping."
      continue
    fi

    run_straingst_species "${sample_id}" "${species_name}" "${db_file}"
  done

  echo "[* straingst_per_species] Finished sample_${sample_id}"
}

# ------------------------------------------------------------
# Main loop: parallelize by sample
# ------------------------------------------------------------

mkdir -p "${OUT_BASE}"

echo "[*] Running samples ${sample_start} to ${sample_end}"
echo "[*] Max parallel samples: ${max_parallel_samples}"

for i in $(seq "${sample_start}" "${sample_end}"); do
  run_one_sample "${i}" &

  # Limit number of concurrently running sample jobs.
  while [ "$(jobs -rp | wc -l)" -ge "${max_parallel_samples}" ]; do
    wait -n
  done
done

# Wait for remaining sample jobs.
wait

echo
echo "[!] All per-species StrainGST runs finished."