#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Driver: build DBs + run StrainGST inference for all datasets, end to end.
#
# Usage:
#   ./run_all.sh                       # runs all datasets
#   ./run_all.sh 4_strain_1-2-5x       # runs just one (or several) datasets
#
# All .sh files (this one, build_straingst_db.sh, run_straingst_inference.sh,
# strainge_datasets.sh) must sit in this folder.
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# Datasets to process (override by passing names as arguments).
if [ "$#" -gt 0 ]; then
  DATASETS=( "$@" )
else
  DATASETS=( 4_strain_simulated 4_strain_1-2-5x 30_strain_simulated 4_strain_mock )
fi

for dataset in "${DATASETS[@]}"; do
  echo
  echo "############################################################"
  echo "# DATASET: ${dataset}"
  echo "############################################################"

  echo "[>>] Building databases..."
  bash build_straingst_db.sh "${dataset}"

  echo "[>>] Running inference..."
  bash run_straingst_inference.sh "${dataset}"

  echo "[>>] Done with ${dataset}."
done

echo
echo "############################################################"
echo "# All datasets finished."
echo "############################################################"