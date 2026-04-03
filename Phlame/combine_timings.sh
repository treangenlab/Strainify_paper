#!/usr/bin/env bash
set -euo pipefail

# Please change all paths in this script to match your local setup
MAKELOG_DIR="/home/Users/rl152/Strainify/phlame/4_strain_mock_community/benchmarks/makedb/time_logs"
PART2_LOG="/home/Users/rl152/Strainify/phlame/4_strain_mock_community/benchmarks/part2_classify.time.log"
OUT_CSV="/home/Users/rl152/Strainify/phlame/4_strain_mock_community/benchmarks/runtime_ram_summary.csv"

# Convert elapsed time string from /usr/bin/time -v to seconds.
# Supports:
#   h:mm:ss
#   m:ss
#   m:ss.xx
elapsed_to_seconds() {
    local t="$1"
    awk -v t="$t" '
    BEGIN {
        n = split(t, a, ":")
        if (n == 3) {
            print a[1] * 3600 + a[2] * 60 + a[3]
        } else if (n == 2) {
            print a[1] * 60 + a[2]
        } else {
            print t + 0
        }
    }'
}

# Parse one /usr/bin/time -v log file.
# Output:
#   cpu_seconds<TAB>wall_seconds<TAB>max_rss_kb
parse_time_log() {
    local logfile="$1"

    local user_sec sys_sec wall_str max_rss

    user_sec=$(awk -F': *' '/User time \(seconds\)/ {print $2}' "$logfile" | tail -n1)
    sys_sec=$(awk -F': *' '/System time \(seconds\)/ {print $2}' "$logfile" | tail -n1)
    max_rss=$(awk -F': *' '/Maximum resident set size \(kbytes\)/ {print $2}' "$logfile" | tail -n1)

    # Extract only the actual elapsed time value, e.g. 0:07.09
    wall_str=$(sed -n 's/.*):[[:space:]]*//p' "$logfile" | grep -F ':' | tail -n1)

    user_sec="${user_sec:-0}"
    sys_sec="${sys_sec:-0}"
    wall_str="${wall_str:-0}"
    max_rss="${max_rss:-0}"

    local wall_sec
    wall_sec=$(elapsed_to_seconds "$wall_str")

    awk -v u="$user_sec" -v s="$sys_sec" -v w="$wall_sec" -v r="$max_rss" '
    BEGIN {
        printf "%.6f\t%.6f\t%d\n", u+s, w, r
    }'
}

if [[ ! -d "$MAKELOG_DIR" ]]; then
    echo "ERROR: makedb time log directory not found: $MAKELOG_DIR" >&2
    exit 1
fi

if [[ ! -f "$PART2_LOG" ]]; then
    echo "ERROR: part2 log not found: $PART2_LOG" >&2
    exit 1
fi

shopt -s nullglob
make_logs=("$MAKELOG_DIR"/*)
shopt -u nullglob

if [[ ${#make_logs[@]} -eq 0 ]]; then
    echo "ERROR: no makedb log files found in $MAKELOG_DIR" >&2
    exit 1
fi

build_cpu_sum=0
build_wall_sum=0
build_peak_rss=0

for log in "${make_logs[@]}"; do
    parsed=$(parse_time_log "$log")
    cpu_sec=$(echo "$parsed" | cut -f1)
    wall_sec=$(echo "$parsed" | cut -f2)
    rss_kb=$(echo "$parsed" | cut -f3)

    build_cpu_sum=$(awk -v a="$build_cpu_sum" -v b="$cpu_sec" 'BEGIN {printf "%.6f", a+b}')
    build_wall_sum=$(awk -v a="$build_wall_sum" -v b="$wall_sec" 'BEGIN {printf "%.6f", a+b}')

    if (( rss_kb > build_peak_rss )); then
        build_peak_rss="$rss_kb"
    fi
done

part2_parsed=$(parse_time_log "$PART2_LOG")
sample_cpu=$(echo "$part2_parsed" | cut -f1)
sample_wall=$(echo "$part2_parsed" | cut -f2)
sample_peak_rss=$(echo "$part2_parsed" | cut -f3)

# Convert KB to MB
build_peak_rss_mb=$(awk -v r="$build_peak_rss" 'BEGIN {printf "%.2f", r/1024}')
sample_peak_rss_mb=$(awk -v r="$sample_peak_rss" 'BEGIN {printf "%.2f", r/1024}')

cat > "$OUT_CSV" <<EOF
metric,build,SRR13355226
Wall clock time (s),$build_wall_sum,$sample_wall
CPU time (s),$build_cpu_sum,$sample_cpu
Max RSS (MB),$build_peak_rss_mb,$sample_peak_rss_mb
EOF

echo "Wrote $OUT_CSV"
cat "$OUT_CSV"