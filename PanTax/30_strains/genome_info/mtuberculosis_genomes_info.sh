#!/usr/bin/env bash
set -euo pipefail

genome_dir="/home/Users/rl152/Strainify/30_strains/mtuberculosis/mtuberculosis_downloads"
out_file="/dodo/rl152/pantax/30_strains/genome_info/mtuberculosis_genomes_info.txt"

echo -e "genome_ID\tstrain_taxid\tspecies_taxid\torganism_name\tid" > "${out_file}"

i=1
for fasta in "${genome_dir}"/*.{fa,fasta,fna}; do
    [ -e "${fasta}" ] || continue

    genome_ID=$(basename "${fasta}")
    genome_ID="${genome_ID%.*}"

    strain_taxid="1773.${i}"
    species_taxid="1773"
    organism_name="Mycobacterium tuberculosis"

    echo -e "${genome_ID}\t${strain_taxid}\t${species_taxid}\t${organism_name}\t${fasta}" >> "${out_file}"

    i=$((i + 1))
done

echo "Wrote: ${out_file}"