"""Corrections and inference-choice sensitivity found during verification of Methods v6 (2026-11-23).
scripts/09 is left unchanged; this script recomputes only the affected quantities.
1. Reverse direction regression: scripts/09 lagged x twice through dropna, so the model started in March 2020 (n = 58)
   instead of February 2020 (n = 59). Recomputed on February 2020 – December 2024.
2. Reference distribution / small-sample correction for HAC Wald tests (normal vs t; statsmodels use_correction).
3. n_eff with the autocorrelation estimator of Pyper & Peterman (1998): r(j) = N/(N-j) * sum (x_t - xbar)(x_t+j - xbar) / sum (x_t - xbar)^2,
   compared with the lag-specific Pearson estimator used in scripts/09."""
import numpy as np, pandas as pd, statsmodels.formula.api as smf, scipy.stats as st, warnings
warnings.filterwarnings('ignore')
exec(open('scripts/09_final_analysis.py').read().split('# Q5 regional')[0])
rows = []
x, c = S['dengue_3day']['2020-01':'2024-12'], cases['2020-01':'2024-12']
d = pd.DataFrame({'y': ly(c), 'x': lx(x), 'm': c.index.month}); d['y1'] = d.y.shift(1); d['x1'] = d.x.shift(1); d = d.dropna()
for lab, f, v in [('forward: x(t-1) -> cases', 'y ~ C(m) + y1 + x1', 'x1'), ('reverse: cases(t-1) -> x', 'x ~ C(m) + x1 + y1', 'y1')]:
    for corr in (False, True):
        r = smf.ols(f, d).fit(cov_type='HAC', cov_kwds={'maxlags': 3, 'use_correction': corr}); z = r.params[v] / r.bse[v]
        rows.append(dict(model=lab, n=int(r.nobs), small_sample_correction=corr, coef=r.params[v], z=z,
                         p_normal=2 * st.norm.sf(abs(z)), p_t=2 * st.t.sf(abs(z), r.df_resid)))
reg = pd.DataFrame(rows)

def ac_used(v, k): return np.corrcoef(v[:-k], v[k:])[0, 1]
def ac_pp(v, k): v = v - v.mean(); return len(v) / (len(v) - k) * (v[:-k] * v[k:]).sum() / (v * v).sum()
def ne(a, b, f, wt=False): n = len(a); return max(n / (1 + 2 * sum(((n - k) / n if wt else 1) * f(a, k) * f(b, k) for k in range(1, 13))), 4.0)
def fz(r, n): s = 1 / np.sqrt(n - 3); return np.tanh(np.arctanh(r) - 1.96 * s), np.tanh(np.arctanh(r) + 1.96 * s)
runs = {'primary 2020-24': (S['dengue_3day'], '2020-01', '2024-12'), 'A Jun 2019-2024': (S['dengue_3day'], '2019-06', '2024-12'),
        'B 72 months': (S['dengue_3day'], '2019-01', '2024-12'), 'C 2022-24': (S['dengue_3day'], '2022-01', '2024-12'),
        'Topic': (S['topic'], '2020-01', '2024-12'), 'day-1': (S['dengue_day1'], '2020-01', '2024-12'), 'Fever control': (S['fever'], '2020-01', '2024-12')}
nrows = []
def add(lab, xs, cs):
    a, b = adj(lx(xs)).values, adj(ly(cs)).values; r = np.corrcoef(a, b)[0, 1]; bb = block_boot_ci(a, b)
    for est, f, wt in [('lag-specific Pearson (scripts/09)', ac_used, False), ('overall-mean ACF x n/(n-k)', ac_pp, False), ('lag-specific Pearson, (n-k)/n weight', ac_used, True)]:
        e = ne(a, b, f, wt); lo, hi = fz(r, e)
        nrows.append(dict(run=lab, estimator=est, r=r, n=len(a), n_eff=e, fisher_lo=lo, fisher_hi=hi, union_lo=min(lo, bb[0]), union_hi=max(hi, bb[1])))
for lab, (s, a0, a1) in runs.items(): add(lab, s[a0:a1], cases[a0:a1])
m = (cases.index.year != 2022) & (cases.index.year >= 2020); add('D excl. 2022', S['dengue_3day'][m], cases[m])
for rg, s in (regS.items() if regC is not None else []):       # regional rows need the non-public EDCD district workbook
    if (s == 0).sum() < 30: add(f'Region {rg} 2020-23', s['2020':'2023'], regC[rg]['2020':])
ne_t = pd.DataFrame(nrows)
with pd.ExcelWriter('out/19_verification_corrections.xlsx') as w:
    reg.to_excel(w, 'direction_tests', index=False); ne_t.to_excel(w, 'n_eff_estimators', index=False)
print(reg.round(4).to_string()); print(ne_t.round(3).to_string())
