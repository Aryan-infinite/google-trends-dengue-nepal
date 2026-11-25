# Revised analysis (2026)

This folder holds the code, data and analysis plan for the revised study "Google Trends and reported dengue cases in Nepal, 2020–2024: a time-series study" (Neupane R, Neupane B, Bhusal A). The code elsewhere in this repository reproduces the first version, released as a preprint (Zenodo, doi:10.5281/zenodo.21773367), and is left unchanged.

The folder contains everything needed to rebuild every reported result except the regional analyses and the district comparisons in Supplementary Table S1: the main analysis, Figures 1–4, Supplementary Figures S1 and S2, and Tables 1 and 2, together with the programs that check them. Supplementary Figure S3 and the regional results also need a data file that is not public (see "Data that are not included"). `ANALYSIS_PLAN.md` is the analysis plan as fixed before the final analysis (two redactions, marked), and `DATA_DICTIONARY.md` describes every input file and variable.

## Run it

You need:

- Python 3 with the packages in `requirements.txt` (`pip install -r requirements.txt`);
- R with the packages `readxl`, `sandwich` and `MASS`;
- Node.js, with the two packages in `package.json` (run `npm install` in this folder);
- the command-line tools Tesseract (for reading text in the figures) and ImageMagick (`convert`, used by the R image checks);
- a FOGRA39 colour profile (see "Colour profile") and the Arial or Liberation Sans font;
- about 4 GB of free memory: the image checks read 1200-dpi files, and one R check alone uses about 3.5 GB, so do not run two scripts at the same time.

From this folder:

```bash
bash run_all_figures.sh
```

The script stops at the first failed check. A complete run ends with `ALL STEPS FINISHED; EVERY CHECK PASSED`. It writes the figures to `manuscript/figures/` (TIFF, CMYK, 1200 dpi), Tables 1 and 2 to `manuscript/Tables_4-5.docx` (the file name is from an earlier numbering; the tables inside are titled Table 1 and Table 2), and intermediate results to `out/`. These folders are not part of the repository. `run_figure1.sh`, from the first step, still rebuilds Figure 1 on its own.

The versions we used are listed in `requirements.txt`. With them, our own rebuilds gave TIFF files that were byte-identical to the submitted ones. Other versions of the plotting library or fonts can move lines by a pixel or more, and the image checks may then report a failure even though the values are right.

## How the checks work

Every figure is checked in four ways:

1. **The numbers.** The plotting script saves the values it actually drew. A separate program in R recomputes them from the original files and compares the two. For Figures 4, S1, S2 and S3 a third program in JavaScript does the same.
2. **The picture.** Two image readers, one in Python and one in R, open the finished TIFF file, find the axes from the pixels and measure where every line, bar and marker sits. Each must match the recomputed values within a few pixels.
3. **The text.** OCR (Tesseract) reads every label in every figure, and the result must equal the intended text.
4. **Nothing extra.** A residual-ink check in Python and in R confirms that nothing was drawn that should not be there.

Tables 1 and 2 are computed three times: by the main analysis, and separately in R and in NumPy. The three must agree before the Word file is written. The finished Word file is then opened by a check in R and a check in Python, and every cell must equal their own recomputation.

The mutation tests (`bash run_mutation_tests.sh`, optional, several hours) check the checks. They feed deliberately falsified data, figures and tables to the checks and report how many falsifications were rejected. Most work in a temporary copy of this folder; `27` and `30` alter the saved values in place and restore them afterwards.

## What each script does

Script comments refer to the analysis plan (`ANALYSIS_PLAN.md`) and, in places, to the earlier numbering of the display items: Tables 4 and 5 are now Tables 1 and 2, and Methods Tables 1–3 are now Supplementary Tables S1–S3.

**Inputs and analysis**

| Script | Purpose |
|---|---|
| `scripts/00_prepare_inputs.py` | Reads the preprint table and averages the three Google Trends retrievals of "dengue" ("<1" counted as 0.5). |
| `scripts/04_dev_regions.py` | Adds up district cases by the five former development regions. Needs the EDCD district workbook. |
| `scripts/06_2024_vintages.py` | Compares 2024 counts reported in the press during the season with the final counts. |
| `scripts/08_external_verification.py` | Compares the case counts with independently published totals. |
| `scripts/09_final_analysis.py` | Main analysis (Questions 1–6, sensitivity analyses). Writes `out/final/final_results.xlsx`. The regional part runs only when the district workbook is present. |
| `scripts/15_nowcast_monthly_refit.py` | Post hoc: the nowcasts with the models refitted before every test month. |
| `scripts/19_verification_corrections.py` | A correction found during verification (reverse-direction model on 59 months) and checks of alternative inference choices. |
| `scripts/20_peak_margin_full.py` | Post hoc: peak months at margins of 5%, 10% and 20%. Needs the district workbook. |
| `scripts/21_reviewer_checks.py`, `22_reviewer_checks.R` | Additional checks A–D of the second addendum to the plan (autocorrelation; leave-one-year-out, within-year and influential-month analyses; off-season months; robust regression), in Python and separately in R. |
| `scripts/23_announcement_timing.py`, `24_announcement_timing.R` | Check E: the months in which search interest and cases reach half their annual maximum, in Python and in R. |

**Figures** (main figures 1–4; S = supplementary)

| Figure | Draws it | Checks the numbers | Reads the TIFF back | Mutation test |
|---|---|---|---|---|
| 1 | `25_figure1.py` | `26_check_figure1.R` | `31` (Python), `32` (R) | `27`, `33` |
| 2 | `28_figure2.py` | `29_check_figure2.R` | `31`, `32` | `30`, `33` |
| 3 | `35_figure3.py` | `36_check_figure3.R` | `38`, `39` | `40` |
| 4 | `41_figure4.py` | `42_check_figure4.R`, `43_nowcast_independent.js` | `44`, `45` | `46` |
| S1 | `50_figureS1.py` | `51_check_figureS1.R`, `52_acf_independent.js` | `53`, `54` | `55` |
| S2 | `62_figureS2.py` | `63_check_figureS2.R`, `64_s3_independent.js` | `65`, `66` | `67` |
| S3 (needs the district workbook) | `56_figureS3.py` | `58_check_figureS3.R`, `57_regional_independent.js` | `59`, `60` | `61` |

`37_ccf_ar_order_sensitivity.py` and `.R` are a post hoc check of Figure 3: how the cross-correlations change with the order of the autoregressive filter. `34_ocr_figures.py` reads the text of every figure, and `47_residual_ink.py` and `48_residual_ink.R` check every figure for stray marks; `49_mutation_test_residual.py` tests those two. `fig_common.py` holds the shared figure settings (9 pt Arial or Liberation Sans, 170 mm width, 1200 dpi, CMYK TIFF).

**Tables**

| Script | Purpose |
|---|---|
| `68_table4_runs.R`, `69_table4_runs_numpy.py` | Table 1 values for every analysis, in R and separately in NumPy. |
| `70_table5_nowcast.R`, `71_table5_nowcast_numpy.py` | Table 2 values for every schedule and search series, in R and separately in NumPy. |
| `72_tables45_values.py` | Takes the values from the main analysis, stops if they differ from the R and NumPy versions, and formats them. |
| `73_tables45_docx.js` | Writes the Word file. |
| `74_check_tables45.R`, `75_check_tables45_docx.py` | Open the Word file and compare every cell with their own recomputation (R and Python). |
| `76_mutation_test_tables45.py` | Feeds falsified copies of the Word file to both checks. |

## Data in this folder

| File | Source |
|---|---|
| `../data/raw/finaldata.xlsx` (already in this repository) | The preprint's 72-month table: monthly dengue cases reported by the Epidemiology and Disease Control Division (EDCD), plus three Google Trends series. Read unaltered. January–May 2019 were not recorded by EDCD and are treated as missing. |
| `data_raw/gt_dengue_term*.csv` | Google Trends, "dengue", Nepal, January 2019 – December 2024, all categories, web search; retrieved 21, 22 and 23 November 2026. |
| `data_raw/gt_dengue_region_*.csv` | Google Trends, "dengue", for each of Nepal's five former development regions, same settings, one retrieval. |
| `data_raw/gt_dengue_fever_topic.csv`, `gt_dengue_fever_term.csv`, `gt_dengu_term.csv`, `gt_devanagari_dengu_term.csv`, `gt_5term_comparison.csv` | The other queries compared when the primary query was chosen (Methods, Section 2.3): the Dengue fever topic, three single terms, and the five terms retrieved together on one scale. Same settings. |
| `governance/announcement_events.csv` | Dated official statements on dengue found in news and WHO reports, with links (Supplementary Figure S2). The 2024 press totals are typed into `scripts/06_2024_vintages.py` with their sources. |
| `verification/` | The record of the cross-implementation verification. `VERIFICATION_SPEC.md` defines every quantity; four programs computed them: the main analysis (`impl1_pipeline.py` reads its results) and three separate implementations (`impl2_R.R`, `impl3_numpy.py`, `impl4_independent_analysis.py`); `compare.py` compared them; the result is `comparison_table.csv` and is summarised in `VERIFICATION_REPORT.md`. The bootstrap resamples shared by all four are in `boot_*.csv`, with checksums in `INDEX_SHA256.txt`. |

`data_raw/SHA256SUMS.txt` lists checksums of the Google Trends files.

## Data that are not included

The district-by-month EDCD workbook (`edcd_district_monthly_2019_2023.xlsx`) is not included, because permission to redistribute it has not been confirmed. It will be described in the paper's data availability statement. Without it:

- Supplementary Figure S3 and the regional results (Question 5) are not built, and `run_all_figures.sh` says so;
- `scripts/20_peak_margin_full.py` and the four programs in `verification/` cannot be rerun. Their stored results are included (`verification/results_*.json` and `verification/comparison_table.csv`); the table checks read the comparison table.

With the workbook placed at `data_raw/edcd_district_monthly_2019_2023.xlsx`, `run_all_figures.sh` builds and checks Figure S3 as well. The workbook can be requested from the Epidemiology and Disease Control Division, Ministry of Health and Population, Nepal.

## Colour profile

CMYK conversion uses a FOGRA39 profile. If yours is not at the default TeX Live path, set `CMYK_ICC=/path/to/profile.icc`. Without a profile the figure scripts warn and use an uncalibrated conversion. The image checks decode the TIFF files with the same profile, so they need it too.
