#!/bin/bash
set -euo pipefail

# -------- Configuration --------
input_csv="mtuberculosis_30_ratios.csv"
read_length=250
fragment_mean=600
fragment_std=150
profile="MSv3"
genome_folder="mtuberculosis_downloads"
output_folder="fastq"
art_bin="path/to/art_illumina"  # Please change this to the actual path of your ART Illumina binary

coverages=(10 20 50 100 200)

# Create output folder once
mkdir -p "$output_folder"

# -------- Parse Header --------
header=$(head -n 1 "$input_csv")
IFS=',' read -ra columns <<< "$header"
samples=("${columns[@]:1}")

# -------- Loop Over Each Coverage --------
for total_coverage in "${coverages[@]}"; do
    echo "===================="
    echo "Simulating ${total_coverage}x"
    echo "===================="

    for ((i=0; i<${#samples[@]}; i++)); do
        sample="${samples[$i]}"
        tag="${total_coverage}x_${sample}"

        echo "Simulating sample: $tag"

        r1_list="${output_folder}/tmp_${tag}_r1.txt"
        r2_list="${output_folder}/tmp_${tag}_r2.txt"
        > "$r1_list"
        > "$r2_list"

        # -------- Process Each Genome --------
        while IFS=',' read -r -a fields; do
            genome="${fields[0]}"
            proportion="${fields[$((i+1))]}"

            [[ -z "$genome" || -z "$proportion" ]] && continue

            base=$(basename "$genome")
            base=${base%.fna}

            indiv_cov=$(echo "$proportion * $total_coverage" | bc -l)

            echo "  - $genome -> ${indiv_cov}x"

            "$art_bin" \
                -ss "$profile" \
                -i "$genome_folder/$genome" \
                -p \
                -l "$read_length" \
                -f "$indiv_cov" \
                -m "$fragment_mean" \
                -s "$fragment_std" \
                -o "$output_folder/${tag}_${base}_sim"

            echo "$output_folder/${tag}_${base}_sim1.fq" >> "$r1_list"
            echo "$output_folder/${tag}_${base}_sim2.fq" >> "$r2_list"

        done < <(tail -n +2 "$input_csv")

        # -------- Merge Reads --------
        echo "Merging reads for $tag"
        xargs cat < "$r1_list" > "$output_folder/${tag}_r1.fq"
        xargs cat < "$r2_list" > "$output_folder/${tag}_r2.fq"

        rm -f "$r1_list" "$r2_list"

        echo "Output: ${tag}_r1.fq, ${tag}_r2.fq"
    done
done

# -------- Cleanup intermediate files --------
echo "Cleaning up intermediate ART files"
rm -f "$output_folder"/*_sim1.fq "$output_folder"/*_sim2.fq "$output_folder"/*.aln