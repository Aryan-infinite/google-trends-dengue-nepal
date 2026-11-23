#!/usr/bin/env bash
# Rebuilds Figure 1 from the raw files and checks it. Run from this folder: bash run_figure1.sh
set -euo pipefail
export PREPRINT_XLSX=../data/raw/finaldata.xlsx   # the preprint table already in this repository
mkdir -p out manuscript/figures
python3 -I scripts/00_prepare_inputs.py          # preprint table -> out/paper_72month_table.csv
python3 -I scripts/25_figure1.py                 # draws Figure 1, writes the plotted values
Rscript scripts/26_check_figure1.R               # independent rebuild in R; must print 0 failed checks
python3 -I scripts/27_mutation_test_figure1.py   # the R check must catch 5 deliberate errors
