# Cross-implementation verification, 2026-11-23

**Specification:** VERIFICATION_SPEC.md, fixed before any comparison. One erratum, for the LOSO test months (59, not 60), was recorded before any comparison.
**Shared bootstrap resamples:** boot_*.csv (SHA256 checksums in INDEX_SHA256.txt).

## Implementations

| # | Program | Written by | Libraries |
|---|---|---|---|
| 1 | Analysis pipeline (scripts/09, 15, 19, 20, 06), read through impl1_pipeline.py | Original analysis | pandas, statsmodels |
| 2 | impl2_R.R | Written from the spec | R stats, MASS::glm.nb, sandwich::NeweyWest |
| 3 | impl3_numpy.py | Written from the spec | NumPy/SciPy only; OLS, HAC, NB2 likelihood, Ljung–Box, Cook's distance and Holm coded by hand |
| 4 | impl4_independent_analysis.py | Written separately from the spec and raw data only, without the scripts, outputs or manuscript | Own NumPy linear algebra, own NB2 likelihood with analytic scores, own HAC |

## Result

- **Agreement:** 218 of 218 keys agree under the pre-set rules; there were 0 disagreements. The median relative spread is 1e-15.
- **Largest spread:** 4e-5, in the Q6 final counts, because the pipeline rounds them to whole cases.
- **The one implementation-dependent value** is the NB rate-ratio interval: 1.199–1.514 in implementations 1, 3 and 4, against 1.172–1.549 in R. R's glm.nb sandwich holds the dispersion fixed. All four intervals exclude 1.
- **Pipeline bootstraps:** its own-RNG bootstraps fall within Monte Carlo tolerance of the shared-resample values (8/8).
- **The comparator itself:** it caught 11 of 12 planted discrepancies. The one it missed (p = 0.0090 changed to 0.00901) is smaller than the pre-set absolute tolerance of 1e-4.

## What agreement does and does not show

- **Shows:** given the specification, the reported numbers are computed correctly. Four programs in two languages and five numerical libraries give the same values, and two of them were written without access to the analysis code.
- **Does not show:** that the specification is the right analysis. All four implementations would reproduce a wrong method identically. The methodological choices were checked separately by two blind reviews of the Methods text against the plan.
- **Does not cover:** the figures from outside publications (Table 1, Q6 press totals), or the Google Trends retrieval settings.
