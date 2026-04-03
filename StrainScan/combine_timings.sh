#!/bin/bash

set -euo pipefail

# Please change this path to match your local setup
LOG_DIR="/home/Users/rl152/Strainify_dev/StrainScan_benchmark/logs/30_strain_sepidermidis"

OUT_DEFAULT="$LOG_DIR/combined_default.csv"
OUT_LOW="$LOG_DIR/combined_low_coverage.csv"

# Convert h:mm:ss or m:ss(.xx) → seconds (float)
to_seconds() {
  local t="$1"
  IFS=: read -r a b c <<<"$t"
  if [[ -n "${c:-}" ]]; then
    awk -v h="$a" -v m="$b" -v s="$c" 'BEGIN{printf "%.2f", h*3600 + m*60 + s}'
  else
    awk -v m="$a" -v s="$b" 'BEGIN{printf "%.2f", m*60 + s}'
  fi
}

# Parse build.log metrics (if exists)
parse_build_metrics() {
  local file="$1"
  if [[ ! -f "$file" ]]; then
    echo "0,0,0"  # dummy values if missing
    return
  fi

  local wall_raw user sys mem_kb mem_mb wall_s cpu_s
  wall_raw=$(awk -F': ' '/Elapsed \(wall clock\) time/{print $2}' "$file" | tr -d '[:space:]')
  user=$(awk -F': ' '/User time \(seconds\)/{print $2}' "$file")
  sys=$(awk -F': ' '/System time \(seconds\)/{print $2}' "$file")
  mem_kb=$(awk -F': ' '/Maximum resident set size/{print $2}' "$file")

  wall_s=$(to_seconds "$wall_raw")
  cpu_s=$(awk -v u="$user" -v s="$sys" 'BEGIN{printf "%.2f", u+s}')
  mem_mb=$(awk -v kb="$mem_kb" 'BEGIN{printf "%.2f", kb/1024}')

  echo "$wall_s,$cpu_s,$mem_mb"
}

# Combine group (default or low_coverage)
combine_group() {
  local pattern="$1"
  local outfile="$2"

  mapfile -t files < <(ls "$LOG_DIR"/*"${pattern}"*.log 2>/dev/null | sort -V || true)
  if [[ ${#files[@]} -eq 0 ]]; then
    echo "No ${pattern} logs found; skipping $outfile"
    return
  fi

  declare -a samples=()
  declare -a wall=()
  declare -a cpu=()
  declare -a mem=()

  # Include build.log first
  build_file="$LOG_DIR/build.log"
  IFS=, read -r build_wall build_cpu build_mem <<<"$(parse_build_metrics "$build_file")"
  samples+=("build")
  wall+=("$build_wall")
  cpu+=("$build_cpu")
  mem+=("$build_mem")

  # Then all sample logs
  for file in "${files[@]}"; do
    base=$(basename "$file" .log)
    sample="${base/_${pattern}/}"

    wall_raw=$(awk -F': ' '/Elapsed \(wall clock\) time/{print $2}' "$file" | tr -d '[:space:]')
    wall_sec=$(to_seconds "$wall_raw")

    user=$(awk -F': ' '/User time \(seconds\)/{print $2}' "$file")
    sys=$(awk -F': ' '/System time \(seconds\)/{print $2}' "$file")
    cpu_time=$(awk -v u="$user" -v s="$sys" 'BEGIN{printf "%.2f", u+s}')

    mem_kb=$(awk -F': ' '/Maximum resident set size/{print $2}' "$file")
    mem_mb=$(awk -v kb="$mem_kb" 'BEGIN{printf "%.2f", kb/1024}')

    samples+=("$sample")
    wall+=("$wall_sec")
    cpu+=("$cpu_time")
    mem+=("$mem_mb")
  done

  # Write CSV
  {
    printf "metric"
    for s in "${samples[@]}"; do printf ",%s" "$s"; done
    printf "\n"

    printf "Wall clock time (s)"
    for v in "${wall[@]}"; do printf ",%s" "$v"; done
    printf "\n"

    printf "CPU time (s)"
    for v in "${cpu[@]}"; do printf ",%s" "$v"; done
    printf "\n"

    printf "Max RSS (MB)"
    for v in "${mem[@]}"; do printf ",%s" "$v"; done
    printf "\n"
  } > "$outfile"

  echo "Wrote $outfile"
}

combine_group "default" "$OUT_DEFAULT"
combine_group "low_coverage" "$OUT_LOW"
