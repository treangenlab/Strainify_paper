#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Paths/settings
# -----------------------------
# Please change the paths below to match your system before running this script.
base_dir="/dodo/rl152/pantax/30_strains"

genome_info_dir="${base_dir}/genome_info"
interleaved_base="${base_dir}/interleaved_fastqs"

db_base="${base_dir}/pantax_dbs"
out_base="${base_dir}/pantax_results"
log_base="${base_dir}/pantax_time_logs"

threads=12

species_list=(
    cacnes
    cdiff
    ecoli
    mtuberculosis
    sepidermidis
)

# Gurobi license
export GRB_LICENSE_FILE="/dodo/rl152/pantax/gurobi1103/gurobi.lic"

# Use /usr/bin/time
TIME="/usr/bin/time"

mkdir -p "${db_base}" "${out_base}" "${log_base}"

failed_samples_log="${log_base}/failed_samples.tsv"

# Only write header if file does not exist yet
if [ ! -f "${failed_samples_log}" ]; then
    echo -e "species\tsample\tfastq\tout_dir\ttime_log\tstderr_log" > "${failed_samples_log}"
fi

# -----------------------------
# Check global inputs
# -----------------------------
if [ ! -f "${GRB_LICENSE_FILE}" ]; then
    echo "[ERROR] Missing Gurobi license: ${GRB_LICENSE_FILE}" >&2
    exit 1
fi

if [ ! -x "${TIME}" ]; then
    echo "[ERROR] Cannot find executable time command: ${TIME}" >&2
    exit 1
fi

if [ ! -d "${genome_info_dir}" ]; then
    echo "[ERROR] Missing genome_info directory: ${genome_info_dir}" >&2
    exit 1
fi

if [ ! -d "${interleaved_base}" ]; then
    echo "[ERROR] Missing interleaved FASTQ base directory: ${interleaved_base}" >&2
    exit 1
fi

# -----------------------------
# Build DB and run samples per species
# -----------------------------
for species in "${species_list[@]}"; do
    genomes_info="${genome_info_dir}/${species}_genomes_info.txt"
    interleaved_dir="${interleaved_base}/${species}"

    db_dir="${db_base}/${species}_db"
    species_out_base="${out_base}/${species}"
    species_log_dir="${log_base}/${species}"

    mkdir -p "${species_out_base}" "${species_log_dir}"

    echo
    echo "============================================================"
    echo "[INFO] Processing species: ${species}"
    echo "[INFO] genomes_info: ${genomes_info}"
    echo "[INFO] interleaved_dir: ${interleaved_dir}"
    echo "[INFO] db_dir: ${db_dir}"
    echo "============================================================"

    # -----------------------------
    # Check species inputs
    # -----------------------------
    if [ ! -f "${genomes_info}" ]; then
        echo "[ERROR] Missing genomes_info file for ${species}: ${genomes_info}" >&2
        exit 1
    fi

    if [ ! -d "${interleaved_dir}" ]; then
        echo "[ERROR] Missing interleaved FASTQ directory for ${species}: ${interleaved_dir}" >&2
        exit 1
    fi

    # -----------------------------
    # Build DB for this species
    # Skip if DB already exists and is non-empty
    # -----------------------------
    if [ -d "${db_dir}" ] && [ "$(find "${db_dir}" -mindepth 1 -print -quit)" ]; then
        echo "[INFO] Existing DB found for ${species}, skipping DB build: ${db_dir}"
    else
        echo "[INFO] Building PanTax DB for ${species}: ${db_dir}"

        export GRB_LICENSE_FILE="/dodo/rl152/pantax/gurobi1103/gurobi.lic"

        "${TIME}" -v -o "${species_log_dir}/build_db.time.log" \
            pantax \
                -f "${genomes_info}" \
                --create \
                -t "${threads}" \
                -d "${db_dir}"

        echo "[INFO] DB build finished for ${species}"
    fi

    # -----------------------------
    # Run all interleaved FASTQs for this species
    # -----------------------------
    found_fastq=0

    for fq in "${interleaved_dir}"/*_interleaved.fq; do
        [ -e "${fq}" ] || continue
        found_fastq=1

        sample_name=$(basename "${fq}" .fq)

        out_dir="${species_out_base}/${sample_name}"
        time_log="${species_log_dir}/${sample_name}.time.log"
        stdout_log="${species_log_dir}/${sample_name}.stdout.log"
        stderr_log="${species_log_dir}/${sample_name}.stderr.log"

        done_file="${out_dir}/PANTAX.DONE"
        failed_file="${out_dir}/PANTAX.FAILED"

        # -----------------------------
        # Skip samples already done
        # -----------------------------
        if [ -f "${done_file}" ]; then
            echo "[INFO] Already done, skipping: ${species}/${sample_name}"
            continue
        fi

        mkdir -p "${out_dir}"

        echo
        echo "------------------------------------------------------------"
        echo "[INFO] Running species: ${species}"
        echo "[INFO] Sample: ${sample_name}"
        echo "[INFO] FASTQ: ${fq}"
        echo "[INFO] Output: ${out_dir}"
        echo "------------------------------------------------------------"

        export GRB_LICENSE_FILE="/dodo/rl152/pantax/gurobi1103/gurobi.lic"

        # Remove old failed marker before retrying
        rm -f "${failed_file}"

        if "${TIME}" -v -o "${time_log}" \
            pantax \
                -d "${db_dir}" \
                -s \
                -p \
                -r "${fq}" \
                --species \
                --strain \
                -t "${threads}" \
                -o "${out_dir}" \
                > "${stdout_log}" \
                2> "${stderr_log}"
        then
            touch "${done_file}"
            echo "[INFO] Finished sample: ${species}/${sample_name}"
            echo "[INFO] Timing log: ${time_log}"
        else
            touch "${failed_file}"
            echo "[WARN] Sample crashed, skipping: ${species}/${sample_name}" >&2
            echo "[WARN] stdout log: ${stdout_log}" >&2
            echo "[WARN] stderr log: ${stderr_log}" >&2

            echo -e "${species}\t${sample_name}\t${fq}\t${out_dir}\t${time_log}\t${stderr_log}" >> "${failed_samples_log}"

            continue
        fi
    done

    if [ "${found_fastq}" -eq 0 ]; then
        echo "[ERROR] No interleaved FASTQs found for ${species} in ${interleaved_dir}" >&2
        exit 1
    fi

    echo
    echo "[INFO] Finished all available samples for species: ${species}"
done

echo
echo "[INFO] All PanTax 30-strain runs finished."
echo "[INFO] DBs: ${db_base}"
echo "[INFO] Outputs: ${out_base}"
echo "[INFO] Timing logs: ${log_base}"
echo "[INFO] Failed samples log: ${failed_samples_log}"