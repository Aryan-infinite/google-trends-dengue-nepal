# Revised analysis (2026)

This folder holds the code and inputs for the revised version of our study of Google search activity and reported dengue cases in Nepal, 2019–2024. The code elsewhere in this repository reproduces the preprint (Zenodo, doi:10.5281/zenodo.21773367) and is left unchanged.

The revised analysis is being added in small, checked steps. This first step contains Figure 1.

## Run it

Requirements: Python 3 with the packages in `requirements.txt`, and R with the `readxl` package. From this folder:

```bash
bash run_figure1.sh
```

This writes `manuscript/figures/Fig1.tif` and must end with `RESULT: 0 failed check(s)` and `caught 5/5`.

## What each script does

| Script | Purpose |
|---|---|
| `scripts/00_prepare_inputs.py` | Reads the preprint's monthly table (path set by `PREPRINT_XLSX` in `run_figure1.sh`) and averages the three Google Trends retrievals of "dengue" ("<1" counted as 0.5). |
| `scripts/fig_common.py` | Shared figure settings: 9 pt Arial (Liberation Sans if Arial is not installed), 170 mm width, 1200 dpi, CMYK TIFF. |
| `scripts/25_figure1.py` | Draws Figure 1. It then reads the plotted values back from the drawn lines and saves them to `out/fig1_plotted_values.csv`. |
| `scripts/26_check_figure1.R` | Separate check in R. Rebuilds the three plotted series from the original files, without any Python output, and compares them with the saved plotted values. |
| `scripts/27_mutation_test_figure1.py` | Feeds five deliberately wrong versions of the plotted values to the R check, which must reject every one. |

## Data in this folder

| File | Source |
|---|---|
| `../data/raw/finaldata.xlsx` (already in this repository) | The preprint's 72-month table: monthly dengue cases reported by the Epidemiology and Disease Control Division (EDCD), plus three Google Trends series. Read unaltered. January–May 2019 were not recorded by EDCD and are treated as missing. |
| `data_raw/gt_dengue_term*.csv` | Google Trends, "dengue", Nepal, January 2019 – December 2024, all categories, web search; retrieved 21, 22 and 23 November 2026. |

`data_raw/SHA256SUMS.txt` lists checksums of the Google Trends files.

The district-by-month EDCD workbook used for the regional analysis is not included, because permission to redistribute it has not been confirmed. It will be described in the paper's data availability statement.

## Colour profile

CMYK conversion uses a FOGRA39 profile. If yours is not at the default TeX Live path, set `CMYK_ICC=/path/to/profile.icc`. Without a profile the script warns and uses an uncalibrated conversion.
