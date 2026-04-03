#!/bin/bash
set -euo pipefail

# Please change all paths in this script to match your local setup
DIR="/home/Users/rl152/Strainify/phlame/4_strain_mock_community/results_customed_clades_2"
cd "$DIR"

OUT="combined_relative_abundance_with_novel.csv"

mapfile -t FILES < <(ls *_frequencies.csv 2>/dev/null | sort)
if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "No *_frequencies.csv files found in $DIR" >&2
  exit 1
fi

# Build header
HEADER="strain"
for f in "${FILES[@]}"; do
  name="${f/_frequencies.csv/}"
  HEADER="$HEADER,$name"
done

echo "$HEADER" > combined.tmp

# Merge on strain name using column 2 = Relative abundance
# If novel_strain would be negative, scale all existing strains down so column sum = 1
awk -F',' -v OFS=',' '
FNR==1 { next }  # skip header line of each file

{
  file_index = ARGIND
  strain = $1
  value = $2 + 0

  data[strain, file_index] = value
  strains[strain] = 1
  if (file_index > max_file) max_file = file_index
  colsum[file_index] += value
}

END {
  # Determine scaling factor per file:
  # if sum > 1, scale all strains by 1/sum and set novel_strain=0
  # otherwise keep as is and novel_strain=1-sum
  for (i = 1; i <= max_file; i++) {
    if (colsum[i] > 1) {
      scale[i] = 1 / colsum[i]
      novel[i] = 0
    } else {
      scale[i] = 1
      novel[i] = 1 - colsum[i]
    }
  }

  # Print strain rows
  for (s in strains) {
    printf "%s", s
    for (i = 1; i <= max_file; i++) {
      val = data[s, i]
      if (val == "") val = 0
      val *= scale[i]
      printf OFS "%.16f", val
    }
    printf "\n"
  }

  # Append novel_strain row
  printf "novel_strain"
  for (i = 1; i <= max_file; i++) {
    printf OFS "%.16f", novel[i]
  }
  printf "\n"
}
' "${FILES[@]}" >> combined.tmp

mv combined.tmp "$OUT"
echo "Wrote: $DIR/$OUT"