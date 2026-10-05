# Please change the paths below to match your system before running this script.
mkdir -p /dodo/rl152/pantax/4_strain_simulated/interleaved_fastqs_1-2-5x

fastq_dir="/home/Users/rl152/Strainify/4_strain_ecoli_simulated/fastq_1-2-5x"
out_dir="/dodo/rl152/pantax/4_strain_simulated/interleaved_fastqs_1-2-5x"

for r1 in "${fastq_dir}"/{1..5}x_ratio_*_r1.fq; do
    [ -e "$r1" ] || continue

    r2="${r1/_r1.fq/_r2.fq}"
    base=$(basename "$r1" _r1.fq)

    if [ ! -f "$r2" ]; then
        echo "[WARN] Missing pair for $r1, skipping"
        continue
    fi

    out="${out_dir}/${base}_interleaved.fq"

    echo "[INFO] Interleaving:"
    echo "       R1: $r1"
    echo "       R2: $r2"
    echo "       OUT: $out"

    seqtk mergepe "$r1" "$r2" > "$out"
done