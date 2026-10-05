#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Paths/settings
# -----------------------------
base_dir="/dodo/rl152/pantax/4_strain_simulated"
interleaved_dir="${base_dir}/interleaved_fastqs"

genomes_info="${base_dir}/genomes_info.txt"
db_dir="/dodo/rl152/pantax/4_strain_simulated/4_strain_ecoli_db"

out_base="${base_dir}/pantax_results"
log_dir="${base_dir}/pantax_time_logs"

threads=12

# Gurobi license
export GRB_LICENSE_FILE="/dodo/rl152/pantax/gurobi1103/gurobi.lic"

mkdir -p "${out_base}" "${log_dir}"

# Use /user/bin/time
TIME="/usr/bin/time"

# -----------------------------
# Check inputs
# -----------------------------
if [ ! -f "${genomes_info}" ]; then
    echo "[ERROR] Missing genomes_info file: ${genomes_info}" >&2
    exit 1
fi

if [ ! -d "${interleaved_dir}" ]; then
    echo "[ERROR] Missing interleaved FASTQ directory: ${interleaved_dir}" >&2
    exit 1
fi

if [ ! -f "${GRB_LICENSE_FILE}" ]; then
    echo "[ERROR] Missing Gurobi license: ${GRB_LICENSE_FILE}" >&2
    exit 1
fi

# -----------------------------
# Build DB once
# -----------------------------
echo "[INFO] Building PanTax DB: ${db_dir}"

"${TIME}" -v -o "${log_dir}/build_db.time.log" \
    pantax \
        -f "${genomes_info}" \
        --create \
        -t "${threads}" \
        -d "${db_dir}"

echo "[INFO] DB build finished."

# -----------------------------
# Run all interleaved FASTQs
# -----------------------------
for fq in "${interleaved_dir}"/*_interleaved.fq; do
    [ -e "${fq}" ] || {
        echo "[ERROR] No interleaved FASTQs found in ${interleaved_dir}" >&2
        exit 1
    }

    sample_name=$(basename "${fq}" .fq)
    out_dir="${out_base}/${sample_name}"
    time_log="${log_dir}/${sample_name}.time.log"

    mkdir -p "${out_dir}"

    echo
    echo "============================================================"
    echo "[INFO] Running sample: ${sample_name}"
    echo "[INFO] FASTQ: ${fq}"
    echo "[INFO] Output: ${out_dir}"
    echo "============================================================"

    export GRB_LICENSE_FILE="/dodo/rl152/pantax/gurobi1103/gurobi.lic"

    "${TIME}" -v -o "${time_log}" \
        pantax \
            -d "${db_dir}" \
            -s \
            -p \
            -r "${fq}" \
            --species \
            --strain \
            -t "${threads}" \
            -o "${out_dir}"

    echo "[INFO] Finished sample: ${sample_name}"
    echo "[INFO] Timing log: ${time_log}"
done

echo
echo "[INFO] All PanTax runs finished."
echo "[INFO] Outputs: ${out_base}"
echo "[INFO] Timing logs: ${log_dir}"