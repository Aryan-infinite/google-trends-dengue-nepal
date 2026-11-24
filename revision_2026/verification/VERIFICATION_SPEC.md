# Cross-implementation verification specification

Fixed 2026-11-23 before any implementation was compared. Every implementation computes the keys below from
`data_raw/` (run from the project root) and writes one flat JSON object `{key: value}`, with numbers as floats and lists as JSON lists.
Agreement rules are in Section 9. Shared bootstrap resamples are in this folder (checksums in INDEX_SHA256.txt).

## 1. Data

- **Cases 2019–2023.** Workbook `edcd_district_monthly_2019_2023.xlsx`, one sheet per year (`2019 Monthly` … `2023 Monthly`).
  - Row 1 holds the title and row 2 the header. Column B is the district name ("Total" and "Not specified" also appear), column C the province, and columns D–O January–December.
  - National count = sum of all rows except "Total"; "Not specified" is included. Blank cells = 0.
  - Check that this equals the "Total" row in every month.
- **Cases 2024.** Workbook `paper_finalized_data_RP.xlsx`, first sheet. Row 1 is the header; columns B–F are month index 1–72, Dengue fever Topic, Fever, Dengue Symptoms and cases. Take the 2024 cases from month indices 61–72.
  - Check that indices 1–60 equal the EDCD national counts.
  - The Topic column ("topic", from the preprint) and Fever column ("fever") come from this file.
- **Primary search series x.** The mean of three files: `gt_dengue_term.csv`, `gt_dengue_term_day2_2026-11-22.csv` and `gt_dengue_term_day3_2026-11-23.csv`.
  - Data rows look like `YYYY-MM,value`, January 2019 – December 2024. Treat "<1" as 0.5.
  - `day1` = the first file alone.
- **Regional search.** `gt_dengue_region_{central,eastern,western}.csv`, same format; use 2020–2023.
- **Regional cases.** Assign districts to the pre-2015 development regions by zone (Mechi, Koshi, Sagarmatha = Eastern; Janakpur, Bagmati, Narayani = Central; Gandaki, Dhaulagiri, Lumbini = Western, including both Nawalparasi districts).
  - "Not specified" is excluded.
  - Spellings vary (e.g. CHITAWAN, KAPILBASTU, DHANUSA, KAVREPALANCHOK, TANAHU, both RUKUM and both NAWALPARASI). Verify that each of the 77 districts is assigned to exactly one of the five regions.

## 2. Transforms

- ly = ln(cases + 1); lx = ln(index + 0.5).
- adj(v) = v minus the mean of v over the same calendar month within the analysis window.

## 3. Question 1 (correlation)

For each run R, a = adj(lx(series)), b = adj(ly(cases)) over the run window, and r = Pearson(a, b).

**Effective sample size and interval**
- ρ(k) = Pearson(v[1..n−k], v[1+k..n]): the lag-specific Pearson correlation.
- n_eff = max(4, n / (1 + 2 Σ_{k=1..12} ρa(k) ρb(k))).
- Fisher limits: tanh(atanh(r) ± 1.96 / sqrt(n_eff − 3)).

**Runs**

| Run | Series | Window |
|---|---|---|
| `primary` | x | 2020-01..2024-12 |
| `A` | x | 2019-06..2024-12 |
| `B` | x | 2019-01..2024-12 |
| `C` | x | 2022-01..2024-12 |
| `D` | x | the 48 months of 2020, 2021, 2023 and 2024 (adjustment over these 48) |
| `topic` | topic | 2020–2024 |
| `day1` | day1 | 2020–2024 |
| `fever` | fever | 2020–2024 |

**Keys:** `q1_{run}_r`, `q1_{run}_neff`, `q1_{run}_flo`, `q1_{run}_fhi`.

**Shared block bootstrap** (primary and fever only, n = 60)
- Row j of `boot_q1_block_starts.csv` holds five 0-based starts s.
- idx = concat(s, s+1, …, s+11) over the five starts (length 60).
- r_j = Pearson(a[idx], b[idx]).
- Limits = 2.5th and 97.5th percentiles with linear interpolation (numpy default / R type 7).
- Keys: `q1_primary_blo`, `q1_primary_bhi`, `q1_fever_blo`, `q1_fever_bhi`.

## 4. Question 2 (primary window 2020–2024)

**Autoregressive fit**
- On a = adj(lx(x)), fit AR(p) with a constant by OLS for p = 0..3 on the common sample t = 4..60.
- AIC = n·ln(RSS/n) + 2(p+1) with n = 57. Choose the p with minimum AIC.
- Refit the chosen p with a constant by OLS on t = p+1..60.
- Keys: `q2_ar_order`, `q2_ar_const`, `q2_ar_phi1`, `q2_ar_phi2` (and `q2_ar_phi3` if p = 3).

**Prewhitening and cross-correlations**
- Filter both a and b: e_t = v_t − const − Σ φ_j v_{t−j}, for t = p+1..60.
- For k = −3..3: r_k = Pearson(ea[t], eb[t+k]) over all t where both exist. Positive k means search precedes cases.
- p-value: t-test with n_k − 2 df.
- Holm adjustment across the 7 lags.
- Keys: `q2_ccf_r_{k}`, `q2_ccf_n_{k}`, `q2_ccf_pholm_{k}`, with k written as m3, m2, m1, 0, p1, p2, p3.

**Direction regressions**
- OLS on 2020-02..2024-12 (n = 59) with 12 month dummies (or intercept plus 11):
  - Forward: ly_t ~ months + ly_{t−1} + lx_{t−1}.
  - Reverse: lx_t ~ months + lx_{t−1} + ly_{t−1}.
- Newey–West HAC: Bartlett weights 1 − j/4 for j = 1..3, no small-sample correction, V = (X'X)⁻¹ S (X'X)⁻¹.
- Two-sided normal p-value.
- Keys: `q2_fwd_coef`, `q2_fwd_se`, `q2_fwd_p`, `q2_rev_coef`, `q2_rev_se`, `q2_rev_p`.
- Also: the reverse model on 2020-03..2024-12 (n = 58), keys `q2_rev58_coef`, `q2_rev58_p`.

## 5. Negative binomial model (2020-02..2024-12, n = 59)

- Counts ~ NB2 (variance μ + αμ²) with log link.
- Covariates: intercept, sin(2πm/12), cos(2πm/12), sin(4πm/12), cos(4πm/12), ly_{t−1}, lx_t.
- Fit by maximum likelihood (α estimated).

**Keys**
- `nb_coef_x` and `nb_rr` = 2^coef.
- `nb_alpha`.
- `nb_rr_lo`, `nb_rr_hi`: from HAC standard errors (Bartlett, 3 lags, no correction) of the coefficient of lx. Report the implementation you used.
- `nb_lb_q12`, `nb_lb_p12`: Ljung–Box Q and p at lag 12 (df 12) on Pearson residuals (c − μ)/sqrt(μ + αμ²).
- `nb_cook_months`: list of "YYYY-MM" where Cook's distance > 4/59, from OLS of ly_t on the same covariates.

## 6. Question 3 (nowcast)

**Data and models**
- Data from 2019-07 onward, with lags available: y = ly_t, y1 = ly_{t−1}, x = lx_t, x1 = lx_{t−1}, m = month factor.
- M0: prediction = cases_{t−1}.
- OLS with 12 month dummies (or intercept plus 11):
  - M1: y ~ m + y1.
  - M2: y ~ m + y1 + x.
  - M3: y ~ m + y1 + x1.
- Prediction = max(exp(ŷ) − 1, 0).

**Schedules**

| Schedule | Fitted on | Test months |
|---|---|---|
| `ps` (per season) | for test year Y, all months 2019-07..(Y−1)-12 | each month of 2022, 2023, 2024 (36) |
| `mo` (monthly) | for test month t, all months 2019-07..t−1 | same 36 |
| `loso` | months 2020-02..2024-12 excluding year Y | months of Y = 2020..2024 within 2020-02..2024-12 (59) |

**Keys** (s ∈ {ps, mo, loso}, M ∈ {M0, M1, M2, M3})
- `q3_{s}_mae_{M}`: mean |pred − obs|.
- `q3_{s}_rmse_{M}`.
- `q3_{s}_maelog_{M}`: mean |ln(1+pred) − ln(1+obs)|.
- For ps and mo only:
  - `q3_{s}_season_mae_{M1|M2}_{2022|2023|2024}` and `q3_{s}_seasons_m2_better` (count of seasons where M2 MAE < M1 MAE).
  - `q3_{s}_diff_mean` = mean(|e_M2| − |e_M1|).
  - `q3_{s}_diff_lo`, `q3_{s}_diff_hi`: percentiles of the mean of d[idx] over the rows of `boot_q3_iid_idx.csv` (0-based into the 36 test months in time order).
  - `q3_{s}_diff_b3lo`, `q3_{s}_diff_b3hi`: same, but idx = concat of 3-month blocks from `boot_q3_block3_starts.csv` (12 starts → 36).
  - `q3_{s}_sens_{M}_{1000|500}`, `q3_{s}_spec_{M}_{1000|500}`: an epidemic month is obs > threshold; the model flags pred > threshold.
  - `q3_{s}_rule_met`: 1 if seasons_m2_better ≥ 2 and diff_hi < 0, else 0.

## 7. Questions 4 and 5

**Question 4**, per year Y = 2019..2024 (national):
- `q4_{Y}_peak_ratio` = max monthly cases / max monthly x.
- `q4_{Y}_total_ratio` = sum cases / sum x.

**Question 5**, per region R ∈ {Central, Eastern, Western}, 2020–2023 (as in Q1, with regional search and regional cases):
- `q5_{R}_r`, `q5_{R}_neff`, `q5_{R}_flo`, `q5_{R}_fhi`.
- `q5_{R}_zeros2020_2023`: number of months with index 0.
- `q5_{R}_peakratio_cases` = (max 2023 cases)/(max 2022 cases); `q5_{R}_peakratio_search` likewise for search.

**Peak classification**
- For each series-year (national cases and national x 2020–2024; regional cases and regional search 2020–2023 for the three regions: 34 in total), ratio = largest month / second-largest month.
- Label at margin g ∈ {0.05, 0.10, 0.20}: the month name if ratio ≥ 1+g, else "tie". Second-largest = 0 counts as an infinite ratio.
- Keys: `peak_n_changed` (number of the 34 whose label differs across the three margins), `peak_ties_05`, `peak_ties_10`, `peak_ties_20`.

## 8. Question 6

**Published 2024 cumulative totals** (from news reports, typed in; not recomputed): 2024-07-31: 2,000; 2024-09-03: 6,000; 2024-09-23: 12,644; 2024-10-15: 19,599; 2024-10-25: 27,153.

- Final cumulative count to date d = sum of 2024 monthly cases in earlier months + (cases in d's month) × (day of d) / (days in that month).
- Keys: `q6_final_{date}` and `q6_pct_{date}` = 100 × published / final, with date as YYYYMMDD.

## 9. Agreement rules (fixed before comparison)

1. **Deterministic quantities** must agree within |a − b| ≤ 1e-4 × max(1, |a|). This covers everything above except items 2–4.
2. **Exact match** for counts, AR order, Cook months, peak counts and rule_met.
3. **NB HAC limits** (`nb_rr_lo`/`nb_rr_hi`) are implementation-dependent. They must agree within 0.05 and all must exclude 1.
4. **Shared-index bootstrap limits** must agree within 1e-4 (same resamples). The pipeline's own-RNG bootstrap is compared with the shared-index values within Monte Carlo tolerance: 0.03 for correlation limits and 15% of |limit| for error-difference limits.

Any disagreement outside these rules is reported as a FAIL with both values. Rules are not changed after results are seen.

## Errata

- 2026-11-23, before any comparison: `loso` test months are 2020-02..2024-12 (59), not 60. Sent to the independent implementer.
