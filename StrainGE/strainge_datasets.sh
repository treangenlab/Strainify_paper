#!/usr/bin/env bash
# ============================================================================
# strainge_datasets.sh
#
# Shared configuration + helpers, sourced by:
#   - build_straingst_db.sh        (kmerize genomes + createdb, per species)
#   - run_straingst_inference.sh   (kmerize reads + run, per read set)
#
# Both scripts take a DATASET name as their first argument:
#   4_strain_simulated | 4_strain_1-2-5x | 30_strain_simulated | 4_strain_mock
#
# Keep this file in the SAME directory as the two scripts above.
# ============================================================================

# ---------------------------------------------------------------------------
# >>> ADJUST HERE if needed <<<
# ---------------------------------------------------------------------------

# Where all StrainGST benchmark output (DBs, results, timing) is written.
BENCH_ROOT="/home/Users/rl152/Strainify_paper/StrainGE"

# Genome FASTA extensions to look for (space separated, uncompressed).
GENOME_EXTS="fna fasta fa"

# Paired-end mate tokens tried IN ORDER to match R1<->R2 files.
# First token whose R1 glob matches wins. Add your convention if it's missing.
_READ_TOKEN_PAIRS=(
  "_r1.:_r2."
  "_r1_:_r2_"
  "_r1:_r2"
  "_R1_:_R2_"
  "_R1.:_R2."
  "_R1:_R2"
  ".R1.:.R2."
  "_1.:_2."
  "_1_:_2_"
  "_1_paired:_2_paired"
  "1_paired:2_paired"
)

# ---------------------------------------------------------------------------
# Dataset configuration
# ---------------------------------------------------------------------------
# Each dataset is described as a set of "species", where every species has:
#   GENOME_DIR[species]  -> folder of reference genomes  (-> one DB)
#   FASTQ_DIR[species]   -> folder of read sets to run against that species DB
# Single-species datasets simply have one entry.
# ---------------------------------------------------------------------------

load_dataset_config() {
  local dataset="$1"

  declare -gA GENOME_DIR
  declare -gA FASTQ_DIR
  GENOME_DIR=()
  FASTQ_DIR=()
  SPECIES_LIST=()

  # Please adjust the paths below to your local setup. 

  case "${dataset}" in
    4_strain_simulated)
      DATASET_NAME="4_strain_simulated"
      SPECIES_LIST=( ecoli )
      GENOME_DIR[ecoli]="/home/Users/rl152/Strainify_paper/simulated_exp/4_strain_simulated/genomes"
      FASTQ_DIR[ecoli]="/home/Users/rl152/Strainify_paper/simulated_exp/4_strain_simulated/fastq"
      ;;

    4_strain_1-2-5x)
      # Same 4-strain E. coli genomes as 4_strain_simulated, but the 1x/2x/5x
      # coverage read sets. Builds its own (identical) DB under this dataset.
      DATASET_NAME="4_strain_1-2-5x"
      SPECIES_LIST=( ecoli )
      GENOME_DIR[ecoli]="/home/Users/rl152/Strainify_paper/simulated_exp/4_strain_simulated/genomes"
      FASTQ_DIR[ecoli]="/home/Users/rl152/Strainify_paper/simulated_exp/4_strain_simulated/fastq_1-2-5x"
      ;;

    30_strain_simulated)
      DATASET_NAME="30_strain_simulated"
      SPECIES_LIST=( ecoli cacnes cdiff mtuberculosis sepidermidis )
      local sp
      for sp in "${SPECIES_LIST[@]}"; do
        GENOME_DIR[$sp]="/home/Users/rl152/Strainify_paper/simulated_exp/30_strain_simulated/${sp}/${sp}_downloads"
        FASTQ_DIR[$sp]="/home/Users/rl152/Strainify_paper/simulated_exp/30_strain_simulated/${sp}/fastq"
      done
      ;;

    4_strain_mock)
      DATASET_NAME="4_strain_mock"
      SPECIES_LIST=( mock )
      GENOME_DIR[mock]="/home/Users/rl152/Strainify_paper/4_strain_mock_community/genomes"
      FASTQ_DIR[mock]="/home/Users/rl152/Strainify_paper/4_strain_mock_community/fastq"
      ;;

    *)
      echo "[! Error] Unknown dataset: '${dataset}'" >&2
      echo "          Valid: 4_strain_simulated | 4_strain_1-2-5x | 30_strain_simulated | 4_strain_mock" >&2
      return 1
      ;;
  esac

  OUT_BASE="${BENCH_ROOT}/${DATASET_NAME}"
  DB_BASE="${OUT_BASE}/db"
  RESULTS_BASE="${OUT_BASE}/results"
}

# ---------------------------------------------------------------------------
# Read-set discovery
# ---------------------------------------------------------------------------

# Emit "set_name<TAB>R1<TAB>R2" for every mate pair directly inside a directory.
#   $1 = directory
#   $2 = optional set-name override (used for the "one subfolder per set" layout)
_emit_pairs_in_dir() {
  local dir="$1"
  local name_override="${2:-}"
  local pair r1tok r2tok
  local -a r1files=()

  local f
  shopt -s nullglob
  for pair in "${_READ_TOKEN_PAIRS[@]}"; do
    r1tok="${pair%%:*}"
    r2tok="${pair##*:}"
    r1files=()
    for f in "${dir}"/*"${r1tok}"*; do
      case "${f}" in
        *.fq|*.fastq|*.fq.gz|*.fastq.gz) r1files+=( "${f}" ) ;;
      esac
    done
    [ "${#r1files[@]}" -gt 0 ] && break
  done
  shopt -u nullglob

  [ "${#r1files[@]}" -eq 0 ] && return 0

  local r1 r2 base setname
  for r1 in "${r1files[@]}"; do
    r2="${r1//${r1tok}/${r2tok}}"
    if [ ! -f "${r2}" ]; then
      echo "[! Warning] No R2 mate for $(basename "${r1}") (expected $(basename "${r2}"))" >&2
      continue
    fi
    if [ -n "${name_override}" ]; then
      setname="${name_override}"
    else
      base=$(basename "${r1}")
      setname="${base%%${r1tok}*}"
      [ -n "${setname}" ] || setname="${base%%.*}"
    fi
    printf '%s\t%s\t%s\n' "${setname}" "${r1}" "${r2}"
  done
}

# Discover read sets in a fastq folder.
#   - If it contains subdirectories, each subdirectory is treated as one set.
#   - Otherwise, flat files in the folder are paired by name.
find_read_sets() {
  local dir="$1"
  local d

  shopt -s nullglob
  local subdirs=( "${dir}"/*/ )
  shopt -u nullglob

  if [ "${#subdirs[@]}" -gt 0 ]; then
    for d in "${subdirs[@]}"; do
      _emit_pairs_in_dir "${d%/}" "$(basename "${d%/}")"
    done
  else
    _emit_pairs_in_dir "${dir}"
  fi
}