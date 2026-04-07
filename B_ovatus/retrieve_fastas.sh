#!/bin/bash

# Step 0: Extract assembly accessions from the AssemblyDetails.txt file
awk -F'\t' '
  $4 ~ /^SAMN[0-9]+$/ {
    samn_num = substr($4, 5)
    if (samn_num >= 11943626 && samn_num <= 11943719) {
      print $1
    }
  }
' PRJNA544527_AssemblyDetails.txt > assemblies.txt


# Step 1: Download only genome files (e.g. .fna)
datasets download genome accession --inputfile assemblies.txt --include genome --filename genomes.zip

# Step 2: Unzip
unzip genomes.zip -d genomes

# Step 3: Create a folder to collect all .fna files
mkdir -p all_fastas

# Step 4: Find and move all .fna files from subfolders into it
find genomes -type f -name "*.fna" -exec mv {} all_fastas/ \;

# Step 5: Clean up
rm -r genomes genomes.zip