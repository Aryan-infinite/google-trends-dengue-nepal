#!/usr/bin/env bash
# Optional and slow (several hours): feeds deliberately falsified data, figures and tables to the checks, which must
# reject them. Run after run_all_figures.sh, from this folder: bash run_mutation_tests.sh
# Each test works in a temporary copy and leaves the files in this folder unchanged.
set -euo pipefail
export PREPRINT_XLSX="$(cd ../data/raw && pwd)/finaldata.xlsx"
for t in 27_mutation_test_figure1 30_mutation_test_figure2 33_mutation_test_images 40_mutation_test_figure3 \
         46_mutation_test_figure4 49_mutation_test_residual 55_mutation_test_figureS1 67_mutation_test_figureS3 \
         76_mutation_test_tables45; do
  echo; echo "=== $t"; python3 -I scripts/$t.py
done
if [ -f data_raw/edcd_district_monthly_2019_2023.xlsx ]; then echo; echo "=== 61_mutation_test_figureS2"; python3 -I scripts/61_mutation_test_figureS2.py; fi
