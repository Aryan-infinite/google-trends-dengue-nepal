# Statistical analysis plan (Supplementary File S3 of the paper)

The plan is reproduced as fixed on 23 November 2026 (version 1.0, with addenda v1.1, withdrawn, and v1.2). Two redactions have been made, each marked in square brackets: the name of the EDCD officer who supplied the district data and the path of an internal file. No other text has been changed. References in the plan to sections of the Methods (for example, "Methods 2.3.3") refer to the corresponding sections of Supplementary File S2 (Section S2.3.3).

---

## Statistical analysis plan — Google search activity and reported dengue, Nepal 2019–2024
Version 1.0, fixed 2026-11-23, before the final run (scripts/09_final_analysis.py).

Status note, stated openly in the paper: the data were explored before this plan was written (the original 2026 analysis, the reviewer-led reanalysis of Nov 2026). Every analysis below is therefore exploratory. This plan exists to stop further choices being made after seeing results: every model, window, transformation, threshold and decision rule is fixed here, and any deviation is reported as a deviation.

## 1. Questions
- Q1 Association. Are monthly searches for "dengue" associated with monthly reported dengue cases beyond the shared seasonal pattern?
- Q2 Timing. At monthly resolution, do searches lead or lag cases?
- Q3 Nowcast value. Does adding same-month searches to a case-only model improve estimates of the current month's cases, out of sample?
- Q4 Scale stability. Is the number of reported cases per unit of search interest stable across seasons?
- Q5 Regional. Within Google's five sub-national regions, do regional searches track regional cases in timing and in size?
- Q6 Reporting delay (descriptive). How complete were publicly announced 2024 totals compared with the final monthly counts?

## 2. Data
| Variable | Source | Months |
|---|---|---|
| National reported cases | EDCD district×month workbook [provided by EDCD; officer's name removed] 2019–2023; EDCD website, accessed 26 Jun 2026, for 2024 | 2019-01 to 2024-12 |
| Regional reported cases | EDCD workbook, 77 districts aggregated to the five pre-2015 development regions | 2019-01 to 2023-12 |
| National searches (primary) | Google Trends, "dengue" (search term), Nepal, 2019-01-01 to 2024-12-31, All categories, Web Search; mean of downloads on 21, 22, 23 Nov 2026 | 72 |
| Regional searches | Same settings, one download per region, 22 Nov 2026 | 72 |
| Secondary search series | "Dengue fever" Topic (paper series, Jun–Jul 2026); "Fever" (non-specific control, comparison-scaled); "Dengue Symptoms" (supplement only) | 72 |

Months are Gregorian (verified: Jan–Aug 2022 file total 4,593 vs 4,595 implied by WHO DON412; October 2022 file 17,889 = EDCD's published October 2022 figure).

## 3. Analysis windows
- Primary: 2020-01 to 2024-12 (60 months). Reason: 2019 fails external checks in both directions (annual total −37%; Kathmandu +28%; Sunsari −96%) and Jan–May 2019 were not counted. This is a data-quality reason fixed before the final run.
- Sensitivity A: 2019-06 to 2024-12 (67 months, includes the counted part of 2019).
- Sensitivity B: full 72 months with 2019 Jan–May as recorded (zero), to compare with the original paper.
- Sensitivity C: excluding 2020–2021 (COVID-19 period).
- Sensitivity D: excluding 2022 (the dominant season).
- Regional analyses: 2020-01 to 2023-12 (48 months); 2024 has no district data.

## 4. Transformations and coding
- Cases: y = log(cases + 1). Searches: x = log(search index + 0.5). "<1" in a Google file is coded 0.5 before averaging (sensitivity: 0).
- Seasonal adjustment: subtract each calendar month's mean within the analysis window (equivalent to regressing on 12 month indicators).

## 5. Methods per question
Q1. Pearson r between seasonally adjusted x and y (primary). Also reported for comparison with the original paper: raw-scale Pearson and Spearman. 95% CIs: (a) Fisher z with effective sample size n_eff = n / (1 + 2 Σ_{k=1..12} ρ_x(k) ρ_y(k)); (b) moving-block bootstrap, block length 12 months, 4,000 resamples, seed 20261123. The wider of (a) and (b) is reported as the headline interval.

Q2. (i) Prewhitened cross-correlation: fit an AR(p) model (p ≤ 3, chosen by AIC) to seasonally adjusted x; filter x and y with it; cross-correlate at lags −3…+3 (positive lag = searches earlier). Report r, n and Holm-adjusted p across the 7 lags. (ii) Incremental-prediction test: y_t ~ month indicators + y_{t−1} + x_{t−1}, OLS with Newey–West (3 lags) standard errors; Wald test on x_{t−1}. Reverse direction (x_t ~ … + y_{t−1}) reported too.
Decision rule: "searches lead cases" may be claimed only if the lag +1 prewhitened r is Holm-significant, exceeds the lag −1 value, and the x_{t−1} Wald p < 0.05.

Q3. Models on the log scale, predictions back-transformed with exp(·) − 1:
- M0 naive: ŷ_t = y_{t−1}.
- M1 season + AR: month indicators + y_{t−1}.
- M2 = M1 + x_t (same-month searches).
- M3 = M1 + x_{t−1}.
Primary evaluation: rolling origin, test seasons 2022, 2023, 2024; training = all earlier months from 2019-06. Secondary: leave-one-season-out over 2020–2024. Metrics: MAE and RMSE (cases), MAE (log). Paired difference M2 − M1 in absolute error with a bootstrap 95% CI over test months (4,000 resamples). Outbreak-month detection: a month is an outbreak month if reported cases > 1,000 (sensitivity: > 500); sensitivity and specificity of each model's predictions.
Decision rule: "searches add nowcast value" only if M2 has lower MAE than M1 in at least 2 of 3 rolling-origin seasons AND the paired-difference CI excludes 0.
Caveat built in: all models use final revised counts; a real-time system would see incomplete counts (Q6).

Regression reported for reviewers: negative binomial GLM, cases_t ~ 2 annual harmonics + log(cases_{t−1}+1) + x_t; HAC standard errors; effect expressed as incidence-rate ratio per doubling of search interest; residual checks: Ljung–Box (lag 12) on Pearson residuals, influence (Cook's distance > 4/n flagged).

Q4. Per season: peak-month cases / peak-month search index; season-total cases / season-total search index.

Q5. Regions with < 30 zero months only (Central, Eastern, Western; Mid- and Far-Western descriptive). Seasonally adjusted log correlation with n_eff CI. Peak month per region-season reported only where the top month exceeds the second by ≥ 10%; otherwise "tie". Size: 2023/2022 ratio of peak cases vs peak searches.

Q6. Table of publicly announced cumulative 2024 totals vs final cumulative counts by the same date (within-month linear interpolation).

Multiple testing: Holm within each family (Q2 lags). No other p-value is used to make a claim; estimates and intervals are primary.

## 6. Robustness grid
Every Q1–Q3 headline is re-run under Sensitivities A–D, with the Topic instead of the term, with day-1 only instead of the 3-day mean, and with "<1" = 0. A result is called robust only if its direction holds in every run.

## 7. Software
Python 3.13, pandas 3.0, numpy 2.5, scipy 1.18, statsmodels 0.15, openpyxl 3.1. Seed 20261123.

## 8. Addendum v1.1 (fixed 2026-11-23, before any 2025 data were seen): held-out validation on the 2025 season
Neither the 2025 monthly case counts nor any Google Trends values for 2025 have been downloaded or viewed by the analysts at the time of writing. Disclosure: press reports seen during literature searches gave cumulative 2025 totals of 5,874 cases by 12 October and 8,109 by 3 December 2025 (Kathmandu Post); no monthly breakdown and no search data were seen. This section is therefore the only confirmatory part of the study.

Data to obtain: (a) EDCD national monthly reported cases, Jan–Dec 2025, with access date; (b) Google Trends "dengue", Nepal, 2019-01-01 to 2025-12-31, All categories, Web Search, downloaded on 3 separate days and averaged.
Procedure: the specifications of Section 5 are re-fitted unchanged on 2019-07 to 2024-12 using the new download, then used to predict 2025 month by month (rolling within 2025: each month uses actual earlier 2025 months for the y_{t−1} term). No re-specification is allowed.
Predictions written now:
- P1. The seasonally adjusted association (Q1 method, window 2020–2025) will remain positive.
- P2. The prewhitened lag +1 correlation will not pass Holm's correction (no lead).
- P3. M2 will not reduce 2025 MAE relative to M1 by more than 10%.
- P4. The 2025 peak case month and peak search month will be the same month or adjacent months.
- P5. Cases per peak search point in 2025 will fall outside the 2022–2024 range (275–739), i.e. the scale will not be stable.
Each prediction is reported as confirmed or not confirmed, whatever the outcome.

### Status of §8 (recorded 2026-11-23, same day): WITHDRAWN before any data collection
EDCD has published only part-year 2025 reports, so the monthly 2025 counts the test requires do not exist. No 2025 data were collected or analysed. The paper reports no 2025 analysis and makes no confirmatory claim.

## 9. Addendum v1.2 (fixed 2026-11-23, before any of these analyses were run): checks requested in peer review
These five checks answer the reviews of the earlier submission [internal file path removed]. They are specified after the primary results were known, so they are post hoc relative to Sections 1–7. They will be labelled "added in response to peer review of an earlier version" and used to meet no decision rule.
Already seen when this was written: all results of Sections 5–6; January–May mean search interest by year (Methods 2.3.3); monthly case counts. Not yet computed: anything below.
Common rules: data, transforms (Section 4) and the primary series (3-day mean) are unchanged; window 2020-01..2024-12 unless stated; seed 20261123; 4,000 bootstrap resamples. Every reported number is produced by scripts/21_reviewer_checks.py (Python) and reproduced by scripts/22_reviewer_checks.R. The two must agree within 1e-4 × max(1, |value|), except bootstrap limits, which use shared resamples saved to verification/boot_addendum_*.csv before the run.

### A. Serial dependence (reviewer 1, point 2)
- Series: ly and lx, raw (log scale) and seasonally adjusted (calendar-month means removed), 2020–2024, n = 60.
- ACF at lags 1–20 (standard estimator: overall mean, denominator n), with ±1.96/√n reference bands.
- PACF at lags 1–20 by OLS regression on lags.
- Ljung–Box Q and p at lags 6 and 12 for each of the four series.
- Output: one supplementary figure (4 series × ACF/PACF) and one table. Descriptive; no decision rule.

### B. Influence of single years (reviewer 1, point 8)
- B1. The Q1 primary computation (r, n_eff, Fisher interval, 12-month block bootstrap, union interval), repeated five times, each excluding one of 2020, 2021, 2022, 2023, 2024 (48 months each). Seasonal adjustment is recomputed over the remaining 48 months. As in run (iv), the remaining months are joined across the gap; this is reported as a limitation.
- B2. Within-year Pearson correlation of lx and ly (not seasonally adjusted, since a single year cannot be), for each year 2020–2024, n = 12, with Fisher interval on n − 3. Descriptive; the shared within-year season inflates these values, and this will be stated.
- Pre-stated reading: the Q1 association is called robust to single-year influence if r > 0 in all five B1 runs. The interval limits are reported whatever they are.
- B3. The negative binomial model (Section 5) refitted without the two months flagged by Cook's distance (2020-06, 2022-09), reporting the rate ratio per doubling and its HAC interval.

### C. Off-season resolution (reviewer 4)
- Season = June–November; off-season = December–May (2020–2024: 30 months each).
- C1. For off-season months: the number and share with primary index ≤ 1, and with index ≤ 2; the number of distinct index values; the range of monthly case counts.
- C2. Pearson correlation of the seasonally adjusted lx and ly (adjusted over the full 60 months, as in Q1), computed separately within season and off-season months (n = 30 each). Interval: Fisher on n − 3, labelled as not accounting for serial dependence, because n_eff is undefined for non-contiguous subsets.
- Pre-stated reading: if ≥80% of off-season months have index ≤ 2, the paper states that national Google Trends lacks the resolution to track off-season transmission. Otherwise it states that the off-season index varies, and reports C2.

### D. Robust regression (reviewer 1, point 8)
- Seasonally adjusted ly regressed on seasonally adjusted lx, 2020–2024: (i) OLS slope with HC1 standard error; (ii) Huber M-estimation (tuning constant 1.345, scale by MAD, IRLS), slope with its default (H1) standard error.
- Report both slopes, their 95% intervals and their ratio. Pre-stated reading: robust if both slopes have the same sign and the Huber 95% interval excludes 0.

### E. Public announcement timeline (reviewer 1, point 9)
- Sources: WHO Disease Outbreak News and situation reports; EDCD/MoHP situation reports and press statements; English-language Nepali national news (Kathmandu Post, Himalayan Times, Khabarhub, Republica) reporting official statements.
- Events recorded 2019–2024: first official national alert or "outbreak" statement each year, and dated official case-count announcements. For each: date, source URL, event type and quoted figure.
- Search procedure: web searches per year ("dengue Nepal <year> outbreak", "EDCD dengue <year>", "dengue alert Nepal <year>"). Each query, its date and the events found are logged in governance/announcement_search_log.csv. Only events whose source page is opened and read are kept.
- Analysis, descriptive only: for each year, the month of the first official national statement against the month in which the primary search index first reaches half its annual maximum, and against the peak case month. No test.

### What does not change
Sections 1–7, their results and the decision rules. None of A–E can alter a reported primary estimate.
