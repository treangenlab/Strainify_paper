#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Build one StrainGST DB per species, for a given dataset.
#
# Usage:
#   ./build_straingst_db.sh <dataset>
#     dataset: 4_strain_simulated | 30_strain_simulated | 4_strain_mock
#
# Every straingst call is wrapped with /usr/bin/time -v; resource logs are
# written next to the outputs under <species>/timing/.
# ============================================================================

SCRIPT_DIR="/home/Users/rl152/Strainify_paper/StrainGE"  # Please change this to your local path
source "${SCRIPT_DIR}/strainge_datasets.sh"

DATASET="${1:-}"
if [ -z "${DATASET}" ]; then
  echo "Usage: $0 <dataset>"
  echo "  dataset: 4_strain_simulated | 30_strain_simulated | 4_strain_mock"
  exit 1
fi

load_dataset_config "${DATASET}"

mkdir -p "${DB_BASE}"

echo "============================================================"
echo "[!] Building StrainGST DBs for dataset: ${DATASET_NAME}"
echo "[!] DB output base: ${DB_BASE}"
echo "============================================================"

for species_name in "${SPECIES_LIST[@]}"; do
  genome_dir="${GENOME_DIR[$species_name]}"

  echo
  echo "------------------------------------------------------------"
  echo "[!] Species: ${species_name}"
  echo "[!] Genomes: ${genome_dir}"
  echo "------------------------------------------------------------"

  if [ ! -d "${genome_dir}" ]; then
    echo "[! Warning] Genome dir not found: ${genome_dir}. Skipping."
    continue
  fi

  species_outdir="${DB_BASE}/${species_name}"
  kmer_dir="${species_outdir}/kmers"
  time_dir="${species_outdir}/timing"
  listing="${species_outdir}/${species_name}_references.txt"

  # straingst createdb -o expects a FILE, not a directory.
  db_file="${species_outdir}/${species_name}.hdf5"

  mkdir -p "${species_outdir}" "${kmer_dir}" "${time_dir}"
  : > "${listing}"

  # Collect genome files across all configured extensions.
  shopt -s nullglob
  genome_files=()
  for ext in ${GENOME_EXTS}; do
    genome_files+=( "${genome_dir}"/*."${ext}" )
  done
  shopt -u nullglob

  if [ "${#genome_files[@]}" -eq 0 ]; then
    echo "[! Warning] No genome files (${GENOME_EXTS}) in ${genome_dir}. Skipping ${species_name}."
    continue
  fi

  for fasta_path in "${genome_files[@]}"; do
    [ -e "${fasta_path}" ] || continue

    fasta_file=$(basename "${fasta_path}")
    asm_name="${fasta_file%.*}"
    kmer_file="${kmer_dir}/${asm_name}.hdf5"

    if [ ! -f "${kmer_file}" ]; then
      echo "[!] Kmerizing genome: ${species_name}/${fasta_file}"
      /usr/bin/time -v -o "${time_dir}/${asm_name}.kmerize.time" \
        straingst kmerize -o "${kmer_file}" "${fasta_path}"
    else
      echo "[!] Kmer file already exists: ${kmer_file}"
    fi

    echo "${kmer_file}" >> "${listing}"
  done

  if [ -f "${db_file}" ]; then
    echo "[!] DB already exists: ${db_file}. Skipping createdb for ${species_name}."
  else
    echo "[!] Creating StrainGST database: ${db_file}"
    /usr/bin/time -v -o "${time_dir}/${species_name}.createdb.time" \
      straingst createdb -o "${db_file}" -f "${listing}"
  fi

  echo "[!] Finished DB for ${species_name}: ${db_file}"
done

echo
echo "[!] All per-species StrainGST databases built for dataset ${DATASET_NAME}."