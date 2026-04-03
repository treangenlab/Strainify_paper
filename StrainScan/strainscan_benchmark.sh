#!/bin/bash
set -euo pipefail

THREADS=12

# Please change all paths in this script to match your local setup
# ======== Construct Databases ========

DB_BASE="/home/Users/rl152/Strainify_dev/StrainScan_benchmark/db"
LOG_BASE="/home/Users/rl152/Strainify_dev/StrainScan_benchmark/logs"
mkdir -p "$DB_BASE" "$LOG_BASE"

declare -A DB_TO_GENOMES=(
  ["4_strain_ecoli_simulated"]="/home/Users/rl152/Strainify/4_strain_ecoli_simulated/genomes"
  ["4_strain_mock"]="/home/Users/rl152/Strainify/4-strain-mock-exp/genomes"
  ["30_strain_cacnes"]="/home/Users/rl152/Strainify/30_strains/cacnes/cacnes_downloads"
  ["30_strain_cdiff"]="/home/Users/rl152/Strainify/30_strains/cdiff/cdiff_downloads"
  ["30_strain_ecoli"]="/home/Users/rl152/Strainify/30_strains/ecoli/ecoli_downloads"
  ["30_strain_mtuberculosis"]="/home/Users/rl152/Strainify/30_strains/mtuberculosis/mtuberculosis_downloads"
  ["30_strain_sepidermidis"]="/home/Users/rl152/Strainify/30_strains/sepidermidis/sepidermidis_downloads"
)

echo "=== STEP 1: Constructing StrainScan databases ==="

for DB_NAME in "${!DB_TO_GENOMES[@]}"; do
    GENOME_DIR="${DB_TO_GENOMES[$DB_NAME]}"
    DB_OUT="${DB_BASE}/${DB_NAME}"
    DB_LOG_DIR="${LOG_BASE}/${DB_NAME}"
    mkdir -p "$DB_OUT" "$DB_LOG_DIR"

    LOG_FILE="${DB_LOG_DIR}/build.log"

    echo "Building database: ${DB_NAME}"
    echo "Input genomes: ${GENOME_DIR}"
    echo "Output path: ${DB_OUT}"
    echo "Logging to: ${LOG_FILE}"

    /usr/bin/time -v strainscan_build -t $THREADS \
        -i "$GENOME_DIR" \
        -o "$DB_OUT" \
        2> "$LOG_FILE"

    echo "Completed DB: ${DB_NAME} (log: $LOG_FILE)"
    echo "------------------------------------------------------"
done


# ======== Run StrainScan using Built DBs ========

OUTPUT_DIR="/home/Users/rl152/Strainify_dev/StrainScan_benchmark/results"
mkdir -p "$OUTPUT_DIR"

declare -A DB_TO_INPUT_DIR=(
  ["4_strain_ecoli_simulated"]="/home/Users/rl152/Strainify/4_strain_ecoli_simulated/fastq"
  ["4_strain_mock"]="/home/Users/rl152/Strainify/4-strain-mock-exp/fastq"
  ["30_strain_cacnes"]="/home/Users/rl152/Strainify/30_strains/cacnes/fastq"
  ["30_strain_cdiff"]="/home/Users/rl152/Strainify/30_strains/cdiff/fastq"
  ["30_strain_ecoli"]="/home/Users/rl152/Strainify/30_strains/ecoli/fastq"
  ["30_strain_mtuberculosis"]="/home/Users/rl152/Strainify/30_strains/mtuberculosis/fastq"
  ["30_strain_sepidermidis"]="/home/Users/rl152/Strainify/30_strains/sepidermidis/fastq"
)

echo "=== STEP 2: Running StrainScan on samples ==="

for DB_NAME in "${!DB_TO_INPUT_DIR[@]}"; do
    DB_PATH="${DB_BASE}/${DB_NAME}"
    INPUT_DIR="${DB_TO_INPUT_DIR[$DB_NAME]}"
    SPECIES_OUT="${OUTPUT_DIR}/${DB_NAME}"
    SPECIES_LOG_DIR="${LOG_BASE}/${DB_NAME}"
    mkdir -p "$SPECIES_OUT" "$SPECIES_LOG_DIR"

    echo "Processing species: ${DB_NAME}"
    echo "Database path: ${DB_PATH}"
    echo "FASTQ input dir: ${INPUT_DIR}"
    echo "Output dir: ${SPECIES_OUT}"
    echo "Log dir: ${SPECIES_LOG_DIR}"
    echo "------------------------------------------------------"

    if [[ ! -d "$DB_PATH" ]]; then
        echo "Skipping $DB_NAME: database not found."
        continue
    fi
    if [[ ! -d "$INPUT_DIR" ]]; then
        echo "Skipping $DB_NAME: FASTQ dir not found."
        continue
    fi

    for R1_FILE in "$INPUT_DIR"/*_r1.fq; do
        SAMPLE_NAME=$(basename "$R1_FILE" _r1.fq)
        R2_FILE="$INPUT_DIR/${SAMPLE_NAME}_r2.fq"

        if [[ -f "$R2_FILE" ]]; then
            echo "Running sample: $SAMPLE_NAME"

            # --- Default mode (no -l 2) ---
            OUT_DEFAULT="${SPECIES_OUT}/default/${SAMPLE_NAME}"
            mkdir -p "$OUT_DEFAULT"
            LOG_DEFAULT="${SPECIES_LOG_DIR}/${SAMPLE_NAME}_default.log"

            echo "→ Running default mode for ${SAMPLE_NAME}"
            /usr/bin/time -v strainscan \
                -i "$R1_FILE" -j "$R2_FILE" \
                -d "$DB_PATH" \
                -o "$OUT_DEFAULT" \
                2> "$LOG_DEFAULT"

            FINAL_REPORT_DEFAULT="${OUT_DEFAULT}/final_report.txt"
            if [[ -f "$FINAL_REPORT_DEFAULT" ]]; then
                mv "$FINAL_REPORT_DEFAULT" "${SPECIES_OUT}/default/${SAMPLE_NAME}.txt"
                echo "Saved default report: ${SAMPLE_NAME}.txt"
            else
                echo "No default report for ${SAMPLE_NAME}"
            fi

            # --- Low coverage mode (-l 2) ---
            OUT_LOW="${SPECIES_OUT}/low_coverage/${SAMPLE_NAME}"
            mkdir -p "$OUT_LOW"
            LOG_LOW="${SPECIES_LOG_DIR}/${SAMPLE_NAME}_low_coverage.log"

            echo "→ Running low coverage mode for ${SAMPLE_NAME}"
            /usr/bin/time -v strainscan \
                -i "$R1_FILE" -j "$R2_FILE" \
                -d "$DB_PATH" \
                -o "$OUT_LOW" \
                -l 2 \
                2> "$LOG_LOW"

            FINAL_REPORT_LOW="${OUT_LOW}/final_report.txt"
            if [[ -f "$FINAL_REPORT_LOW" ]]; then
                mv "$FINAL_REPORT_LOW" "${SPECIES_OUT}/low_coverage/${SAMPLE_NAME}.txt"
                echo "Saved low coverage report: ${SAMPLE_NAME}.txt"
            else
                echo "No low coverage report for ${SAMPLE_NAME}"
            fi

        else
            echo "Missing _r2.fq for $SAMPLE_NAME"
        fi
    done
done

echo "=== All database builds and dual-mode analyses complete! ==="



