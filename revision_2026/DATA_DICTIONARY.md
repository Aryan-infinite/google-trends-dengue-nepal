# Data dictionary

All series are monthly, January 2019 to December 2024 (72 rows), for Nepal unless stated otherwise.

## Case data

**`data/raw/finaldata.xlsx`** (repository root; the table used for the preprint, read unaltered). Header in row 2; 72 data rows.

| Column | Definition | Values |
|---|---|---|
| Months | Month number; 1 = January 2019, 72 = December 2024 | 1–72 |
| Confirmed Cases | Dengue cases reported to the Epidemiology and Disease Control Division (EDCD) in the month. Despite the heading, the counts include suspected as well as confirmed cases (Methods, Section 2.2). January–May 2019 were not recorded by EDCD; they are entered as 0 and are treated as missing in the analysis. 2024 values were taken from the EDCD website in June 2026. | count, 0–27,529 |
| Dengue fever | Google Trends series retrieved for the preprint. The preprint describes it as a search term; it was identified as the Dengue fever topic (Supplementary File S2, Section S2.3). Retrieved together with the two columns below on one scale. | index, 1–100 |
| Dengue Symptoms | Google Trends search term, same retrieval as "Dengue fever" | index, 0–31 |
| Fever | Google Trends search term, same retrieval; used as a non-specific control | index, 9–35 |

The district-by-month workbook supplied by EDCD for 2019–2023 is not included (see the data availability statement).

## Google Trends files (`data_raw/`)

Settings for every file: all categories, web search, monthly, 1 January 2019 to 31 December 2024 (Supplementary Table S2). Each file has two header lines, then a column header `Month` (YYYY-MM) and one column per query. Values are search interest scaled 0–100 to the highest month in that file; "<1" means above zero but below one index point; 0 means too few searches to report.

| File | Query | Retrieved |
|---|---|---|
| `gt_dengue_term.csv` | "dengue", search term (primary series, retrieval 1) | 21 Nov 2026 |
| `gt_dengue_term_day2_2026-11-22.csv` | same (retrieval 2) | 22 Nov 2026 |
| `gt_dengue_term_day3_2026-11-23.csv` | same (retrieval 3) | 23 Nov 2026 |
| `gt_dengue_fever_topic.csv` | Dengue fever, topic | 21–22 Nov 2026 |
| `gt_dengue_fever_term.csv` | "Dengue fever", search term, retrieved alone | 21–22 Nov 2026 |
| `gt_dengu_term.csv` | "dengu" (romanised Nepali), retrieved alone | 21–22 Nov 2026 |
| `gt_devanagari_dengu_term.csv` | "डेंगु" (Devanagari), retrieved alone | 21–22 Nov 2026 |
| `gt_5term_comparison.csv` | dengue, Dengue fever, Dengue Symptoms, dengu and डेंगु, retrieved together on one scale | 21–22 Nov 2026 |
| `gt_dengue_region_{central,eastern,western,midwestern,farwestern}.csv` | "dengue", search term, for each of the five former development regions; each region is scaled to its own maximum | 21–22 Nov 2026 |

`data_raw/SHA256SUMS.txt` gives a checksum for each file.

## Other inputs

| File | Content |
|---|---|
| `governance/announcement_events.csv` | Dated national dengue alerts and official statements, 2019–2024: year, date, event type, issuing body, statement, any figure quoted, source, URL, and whether the page was opened (Supplementary Figure S2) |
| `scripts/06_2024_vintages.py` | Cumulative 2024 case totals reported in the national press, typed in with their sources (Supplementary Table S4) |

## Files created by the code (`out/`)

| File | Content |
|---|---|
| `paper_72month_table.csv` | `finaldata.xlsx` with short column names (m, DF, Fever, DS, cases) and a calendar date |
| `gt_dengue_3day.csv` | The three retrievals of the primary series ("<1" set to 0.5), their mean (the analysed series) and range |
| `final/final_results.xlsx` | Estimates from the main analysis, in sheets by question (Q2–Q6), sensitivity analysis and data check |

## Transformations used in the analysis

Cases: log(cases + 1). Search interest: log(index + 0.5), with "<1" set to 0.5 (0 in a sensitivity analysis). Seasonal adjustment: the mean of each calendar month over the analysis window is subtracted.
