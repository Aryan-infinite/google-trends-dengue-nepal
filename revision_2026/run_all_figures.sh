#!/usr/bin/env bash
# Rebuilds Figures 1-4, S1 and S2 and Tables 1-2 from the raw files and runs every check. Run from this folder:
#   bash run_all_figures.sh
# Stops at the first failed check. Supplementary Figure S3 and the regional results are built only if the EDCD
# district workbook (not public) has been placed at data_raw/edcd_district_monthly_2019_2023.xlsx.
set -euo pipefail
export PREPRINT_XLSX="$(cd ../data/raw && pwd)/finaldata.xlsx"   # the preprint table already in this repository
WB=data_raw/edcd_district_monthly_2019_2023.xlsx
mkdir -p out/final manuscript/figures
step() { echo; echo "=== $*"; }

step "Inputs and analysis"
python3 -I scripts/00_prepare_inputs.py
if [ -f "$WB" ]; then python3 -I scripts/04_dev_regions.py "$WB" > out/04_dev_regions.txt; else echo "(EDCD district workbook not present: regional analysis and Figure S3 skipped)"; fi
python3 -I scripts/06_2024_vintages.py > out/06_2024_vintages.txt
python3 -I scripts/08_external_verification.py > out/08_external_verification.txt
python3 -I scripts/09_final_analysis.py > out/final/run_log.txt
python3 -I scripts/15_nowcast_monthly_refit.py > out/15_nowcast_monthly_refit.txt
python3 -I scripts/19_verification_corrections.py > out/19_verification_corrections.txt
python3 -I scripts/21_reviewer_checks.py > out/21_reviewer_checks.txt
Rscript scripts/22_reviewer_checks.R > out/22_reviewer_checks.txt
python3 -I scripts/23_announcement_timing.py > out/23_announcement_timing.txt
Rscript scripts/24_announcement_timing.R
echo "done"

step "Figure 1"
python3 -I scripts/25_figure1.py
Rscript scripts/26_check_figure1.R

step "Figure 2"
python3 -I scripts/28_figure2.py
Rscript scripts/29_check_figure2.R

step "Figures 1 and 2: reading the TIFF files back (Python, then R)"
python3 -I scripts/31_digitise_figures.py
Rscript scripts/32_digitise_figures.R

step "Figure 3"
python3 -I scripts/35_figure3.py
Rscript scripts/36_check_figure3.R
python3 -I scripts/37_ccf_ar_order_sensitivity.py
Rscript scripts/37_ccf_ar_order_sensitivity.R
python3 -I scripts/38_digitise_figure3.py
Rscript scripts/39_digitise_figure3.R

step "Figure 4"
python3 -I scripts/41_figure4.py
node scripts/43_nowcast_independent.js
Rscript scripts/42_check_figure4.R
python3 -I scripts/44_digitise_figure4.py
Rscript scripts/45_digitise_figure4.R

step "Supplementary Figure S1"
python3 -I scripts/50_figureS1.py
node scripts/52_acf_independent.js
Rscript scripts/51_check_figureS1.R
python3 -I scripts/53_digitise_figureS1.py
Rscript scripts/54_digitise_figureS1.R

step "Supplementary Figure S2"
python3 -I scripts/62_figureS2.py
node scripts/64_s3_independent.js
Rscript scripts/63_check_figureS2.R
python3 -I scripts/65_digitise_figureS2.py
Rscript scripts/66_digitise_figureS2.R

step "Supplementary Figure S3"
if [ -f "$WB" ]; then
  python3 -I scripts/56_figureS3.py
  node scripts/57_regional_independent.js
  Rscript scripts/58_check_figureS3.R
  python3 -I scripts/59_digitise_figureS3.py
  Rscript scripts/60_digitise_figureS3.R
else echo "skipped: needs the EDCD district workbook, which is not public"; fi

step "All figures: text read by OCR, and no stray marks (Python, then R)"
python3 -I scripts/34_ocr_figures.py
python3 -I scripts/47_residual_ink.py
Rscript scripts/48_residual_ink.R

step "Tables 1 and 2 (written to manuscript/Tables_4-5.docx; the file name is from an earlier numbering)"
Rscript scripts/68_table4_runs.R > out/68_table4_runs.txt
python3 -I scripts/69_table4_runs_numpy.py > out/69_table4_runs_numpy.txt
Rscript scripts/70_table5_nowcast.R > out/70_table5_nowcast.txt
python3 -I scripts/71_table5_nowcast_numpy.py > out/71_table5_nowcast_numpy.txt
python3 -I scripts/72_tables45_values.py > out/72_tables45_values.txt
node scripts/73_tables45_docx.js
Rscript scripts/74_check_tables45.R
python3 -I scripts/75_check_tables45_docx.py

echo; echo "ALL STEPS FINISHED; EVERY CHECK PASSED"
