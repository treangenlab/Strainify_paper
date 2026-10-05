#!/usr/bin/env bash
set -euo pipefail
# Please change the paths below to match your system before running this script.
base_dir="/home/Users/rl152/Strainify/30_strains"
out_base="/dodo/rl152/pantax/30_strains/interleaved_fastqs"

species_list=(
    cacnes
    cdiff
    ecoli
    mtuberculosis
    sepidermidis
)

mkdir -p "${out_base}"

for species in "${species_list[@]}"; do
    fastq_dir="${base_dir}/${species}/fastq"
    out_dir="${out_base}/${species}"

    mkdir -p "${out_dir}"

    if [ ! -d "${fastq_dir}" ]; then
        echo "[WARN] Missing fastq directory for ${species}: ${fastq_dir}"
        continue
    fi

    echo
    echo "============================================================"
    echo "[INFO] Processing species: ${species}"
    echo "[INFO] FASTQ dir: ${fastq_dir}"
    echo "[INFO] Output dir: ${out_dir}"
    echo "============================================================"

    for r1 in "${fastq_dir}"/*_r1.fq; do
        [ -e "${r1}" ] || {
            echo "[WARN] No R1 FASTQs found for ${species} in ${fastq_dir}"
            continue
        }

        r2="${r1/_r1.fq/_r2.fq}"
        base=$(basename "${r1}" _r1.fq)

        if [ ! -f "${r2}" ]; then
            echo "[WARN] Missing pair for ${r1}, skipping"
            continue
        fi

        out="${out_dir}/${base}_interleaved.fq"

        echo "[INFO] Interleaving:"
        echo "       R1:  ${r1}"
        echo "       R2:  ${r2}"
        echo "       OUT: ${out}"

        seqtk mergepe "${r1}" "${r2}" > "${out}"
    done
done

echo
echo "[INFO] Done creating interleaved FASTQs."
echo "[INFO] Outputs are in: ${out_base}"