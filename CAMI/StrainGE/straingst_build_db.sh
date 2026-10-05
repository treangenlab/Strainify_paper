#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------
# Please change these paths to your local setup.

# Folder containing one subfolder per species.
# Each species folder contains .fna strain genomes.
genomes_base="/dodo/rl152/ChronoStrain/chronostrain_cami/core_genome_only/core_genomes"

# Put output OUTSIDE genomes_base so the script does not accidentally
# treat the output folder as another species folder.
out_base="/dodo/rl152/ChronoStrain/chronostrain_cami/core_genome_only/strainge/db_per_species"

mkdir -p "${out_base}"

# ------------------------------------------------------------
# Build one StrainGST DB per species
# ------------------------------------------------------------

for species_dir in "${genomes_base}"/*; do
  [ -d "${species_dir}" ] || continue

  species_name=$(basename "${species_dir}")

  echo
  echo "============================================================"
  echo "[!] Processing species: ${species_name}"
  echo "============================================================"

  species_outdir="${out_base}/${species_name}"
  kmer_dir="${species_outdir}/kmers"
  listing="${species_outdir}/${species_name}_references.txt"

  # IMPORTANT:
  # straingst createdb -o expects a FILE, not a directory.
  db_file="${species_outdir}/${species_name}.hdf5"

  mkdir -p "${species_outdir}" "${kmer_dir}"

  # Generated automatically from all .fna genomes in this species folder
  > "${listing}"

  found_any=0

  for fasta_path in "${species_dir}"/*.fna; do
    [ -e "${fasta_path}" ] || continue

    found_any=1

    fasta_file=$(basename "${fasta_path}")
    asm_name="${fasta_file%.fna}"

    kmer_file="${kmer_dir}/${asm_name}.hdf5"

    if [ ! -f "${kmer_file}" ]; then
      echo "[!] Kmerizing: ${species_name}/${fasta_file}"
      straingst kmerize -o "${kmer_file}" "${fasta_path}"
    else
      echo "[!] Kmer file already exists: ${kmer_file}"
    fi

    echo "${kmer_file}" >> "${listing}"
  done

  if [ "${found_any}" -eq 0 ]; then
    echo "[! Warning] No .fna files found in ${species_dir}. Skipping ${species_name}."
    continue
  fi

  if [ -f "${db_file}" ]; then
    echo "[!] DB already exists: ${db_file}"
    echo "[!] Skipping createdb for ${species_name}."
  else
    echo "[!] Creating StrainGST database for species: ${species_name}"
    straingst createdb -o "${db_file}" -f "${listing}"
  fi

  echo "[!] Finished DB for ${species_name}: ${db_file}"
done

echo
echo "[!] All per-species StrainGST databases created."