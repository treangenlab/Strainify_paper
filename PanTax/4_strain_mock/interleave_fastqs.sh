# Please change the paths below to match your system before running this script.
mkdir -p /dodo/rl152/pantax/4_strain_mock/interleaved_fastqs

fastq_dir="/home/Users/rl152/Strainify/4-strain-mock-exp/fastq"
out_dir="/dodo/rl152/pantax/4_strain_mock/interleaved_fastqs"

for r1 in "${fastq_dir}"/*_r1.fq; do
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