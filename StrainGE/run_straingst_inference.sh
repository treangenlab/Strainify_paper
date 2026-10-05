#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Run StrainGST inference for a given dataset.
#
# Usage:
#   ./run_straingst_inference.sh <dataset>
#     dataset: 4_strain_simulated | 4_strain_1-2-5x | 30_strain_simulated | 4_strain_mock
#
# For each species in the dataset, every read set in that species' fastq
# folder is kmerized and run against that species' DB. Every straingst call is
# wrapped with /usr/bin/time -v; resource logs go under
# <species>/<set>/timing/.
# ============================================================================

SCRIPT_DIR="/home/Users/rl152/Strainify_paper/StrainGE" # Please change this to your local path
source "${SCRIPT_DIR}/strainge_datasets.sh"

DATASET="${1:-}"
if [ -z "${DATASET}" ]; then
  echo "Usage: $0 <dataset>"
  echo "  dataset: 4_strain_simulated | 4_strain_1-2-5x | 30_strain_simulated | 4_strain_mock"
  exit 1
fi

load_dataset_config "${DATASET}"

# StrainGST parameters (matching your originals).
num_iters=408
min_score=0
read_k=23

mkdir -p "${RESULTS_BASE}"

for species_name in "${SPECIES_LIST[@]}"; do
  fastq_dir="${FASTQ_DIR[$species_name]}"
  db_file="${DB_BASE}/${species_name}/${species_name}.hdf5"

  echo
  echo "============================================================"
  echo "[*] Dataset ${DATASET_NAME} | species ${species_name}"
  echo "[*] Fastq dir: ${fastq_dir}"
  echo "[*] DB: ${db_file}"
  echo "============================================================"

  if [ ! -f "${db_file}" ]; then
    echo "[! Warning] DB not found: ${db_file}. Run build_straingst_db.sh first. Skipping."
    continue
  fi
  if [ ! -d "${fastq_dir}" ]; then
    echo "[! Warning] Fastq dir not found: ${fastq_dir}. Skipping."
    continue
  fi

  # Discover the read sets (handles flat files and one-subfolder-per-set).
  mapfile -t read_sets < <(find_read_sets "${fastq_dir}")
  if [ "${#read_sets[@]}" -eq 0 ]; then
    echo "[! Warning] No paired read sets found in ${fastq_dir}. Skipping ${species_name}."
    continue
  fi

  echo "[*] Found ${#read_sets[@]} read set(s)."

  for line in "${read_sets[@]}"; do
    IFS=$'\t' read -r set_name fq1 fq2 <<< "${line}"

    outdir="${RESULTS_BASE}/${species_name}/${set_name}"
    time_dir="${outdir}/timing"
    breadcrumb="${outdir}/strainGST.DONE"

    if [ -f "${breadcrumb}" ]; then
      echo "[!] ${species_name}/${set_name} already done. Skipping."
      continue
    fi

    mkdir -p "${outdir}" "${time_dir}"

    # Kmerize this read set (one set == one sample here).
    read_kmers="${outdir}/reads.hdf5"
    if [ ! -f "${read_kmers}" ]; then
      echo "[!] Kmerizing reads: ${species_name}/${set_name}"
      echo "    R1: ${fq1}"
      echo "    R2: ${fq2}"
      /usr/bin/time -v -o "${time_dir}/reads.kmerize.time" \
        straingst kmerize -k "${read_k}" -o "${read_kmers}" "${fq1}" "${fq2}"
    else
      echo "[!] Read kmers already exist: ${read_kmers}"
    fi

    results_tsv="${outdir}/result.tsv"

    echo "[!] Running StrainGST: ${species_name}/${set_name}"
    /usr/bin/time -v -o "${time_dir}/straingst_run.time" \
      straingst run \
        -o "${results_tsv}" \
        -i "${num_iters}" \
        -s "${min_score}" \
        --separate-output \
        "${db_file}" "${read_kmers}"

    touch "${breadcrumb}"
  done
done

echo
echo "[!] All StrainGST inference runs finished for dataset ${DATASET_NAME}."