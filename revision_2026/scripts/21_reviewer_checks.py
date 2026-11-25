"""Checks A-D of ANALYSIS_PLAN.md section 9 (addendum v1.2, fixed 2026-11-23 07:58 UTC before this script was written).
Outputs: out/21_reviewer_checks.json (flat keys, compared with scripts/22 in R) and out/21_reviewer_checks.xlsx."""
import json, numpy as np, pandas as pd, warnings, scipy.stats as st
import statsmodels.api as sm, statsmodels.formula.api as smf
from statsmodels.tsa.stattools import acf
from statsmodels.stats.diagnostic import acorr_ljungbox
warnings.filterwarnings('ignore')
exec(open('scripts/09_final_analysis.py').read().split('# ---------- windows ----------')[0])
R, T = {}, {}
x, c = S['dengue_3day']['2020':'2024'], cases['2020':'2024']
STARTS = np.loadtxt('verification/boot_addendum_n48_block_starts.csv', delimiter=',', dtype=int)

# ---------- A. serial dependence ----------
ser = {'ly_raw': ly(c).values, 'lx_raw': lx(x).values, 'ly_adj': adj(ly(c)).values, 'lx_adj': adj(lx(x)).values}
def pacf_ols(v, K=20):
    out = []
    for k in range(1, K + 1):
        X = np.column_stack([np.ones(len(v) - k)] + [v[k - j:len(v) - j] for j in range(1, k + 1)])
        out.append(np.linalg.lstsq(X, v[k:], rcond=None)[0][k])
    return np.array(out)
rows = []
for nm, v in ser.items():
    a = acf(v, nlags=20, adjusted=False, fft=False)[1:]; p = pacf_ols(v)
    for k in range(20): rows.append(dict(series=nm, lag=k + 1, acf=a[k], pacf=p[k])); R[f'A_{nm}_acf{k+1}'] = a[k]; R[f'A_{nm}_pacf{k+1}'] = p[k]
    lb = acorr_ljungbox(v, lags=[6, 12], return_df=True)
    for L in (6, 12): R[f'A_{nm}_lbq{L}'] = float(lb.loc[L, 'lb_stat']); R[f'A_{nm}_lbp{L}'] = float(lb.loc[L, 'lb_pvalue'])
T['A_acf_pacf'] = pd.DataFrame(rows); R['A_band'] = 1.96 / np.sqrt(60)

# ---------- B. single-year influence ----------
def corr_shared(a, b):
    r = np.corrcoef(a, b)[0, 1]; ne = neff(a, b); f = fisher_ci(r, ne)
    rs = [np.corrcoef(a[i], b[i])[0, 1] for i in (np.concatenate([np.arange(s, s + 12) for s in row])[:len(a)] for row in STARTS)]
    bl, bh = np.percentile(rs, [2.5, 97.5]); return r, ne, f[0], f[1], bl, bh, min(f[0], bl), max(f[1], bh)
b1 = []
for yr in range(2020, 2025):
    m = x.index.year != yr; a, b = adj(lx(x[m])).values, adj(ly(c[m])).values
    vals = corr_shared(a, b); b1.append(dict(excluded=yr, n=len(a), r=vals[0], n_eff=vals[1], fisher_lo=vals[2], fisher_hi=vals[3], boot_lo=vals[4], boot_hi=vals[5], ci_lo=vals[6], ci_hi=vals[7]))
    for k, v in zip(['r', 'neff', 'flo', 'fhi', 'blo', 'bhi', 'lo', 'hi'], vals): R[f'B1_ex{yr}_{k}'] = v
T['B1_leave_one_year_out'] = pd.DataFrame(b1); R['B1_all_positive'] = int(all(d['r'] > 0 for d in b1))
b2 = []
for yr in range(2020, 2025):
    r = np.corrcoef(lx(x[str(yr)]), ly(c[str(yr)]))[0, 1]; lo, hi = fisher_ci(r, 12)
    b2.append(dict(year=yr, n=12, r=r, lo=lo, hi=hi)); R[f'B2_{yr}_r'], R[f'B2_{yr}_lo'], R[f'B2_{yr}_hi'] = r, lo, hi
T['B2_within_year'] = pd.DataFrame(b2)
d = pd.DataFrame({'c': c, 'y1': ly(cases).shift(1)['2020':'2024'], 'x': lx(x)})
d['s1'], d['c1'] = np.sin(2 * np.pi * d.index.month / 12), np.cos(2 * np.pi * d.index.month / 12)
d['s2'], d['c2'] = np.sin(4 * np.pi * d.index.month / 12), np.cos(4 * np.pi * d.index.month / 12)
d = d['2020-02':]; d = d[~d.index.strftime('%Y-%m').isin(['2020-06', '2022-09'])]
nb = smf.negativebinomial('c ~ s1 + c1 + s2 + c2 + y1 + x', d).fit(disp=0, maxiter=500, cov_type='HAC', cov_kwds={'maxlags': 3})
bx, se = nb.params['x'], nb.bse['x']
R['B3_n'], R['B3_coef_x'], R['B3_rr'], R['B3_rr_lo'], R['B3_rr_hi'], R['B3_alpha'] = len(d), bx, 2 ** bx, 2 ** (bx - 1.96 * se), 2 ** (bx + 1.96 * se), nb.params['alpha']
T['B3_nb_without_cook'] = pd.DataFrame([dict(n=len(d), rate_ratio=2 ** bx, lo=R['B3_rr_lo'], hi=R['B3_rr_hi'], alpha=nb.params['alpha'])])

# ---------- C. off-season resolution ----------
off = x.index.month.isin([12, 1, 2, 3, 4, 5]); xo, co = x[off], c[off]
R['C1_n_off'] = int(off.sum()); R['C1_le1'] = int((xo <= 1).sum()); R['C1_le2'] = int((xo <= 2).sum())
R['C1_share_le2'] = R['C1_le2'] / R['C1_n_off']; R['C1_distinct'] = int(xo.round(6).nunique())
R['C1_case_min'], R['C1_case_max'] = float(co.min()), float(co.max())
R['C1_index_min'], R['C1_index_max'] = float(xo.min()), float(xo.max())
a, b = adj(lx(x)), adj(ly(c))
for lab, msk in [('season', ~off), ('off', off)]:
    r = np.corrcoef(a[msk], b[msk])[0, 1]; lo, hi = fisher_ci(r, int(msk.sum()))
    R[f'C2_{lab}_n'], R[f'C2_{lab}_r'], R[f'C2_{lab}_lo'], R[f'C2_{lab}_hi'] = int(msk.sum()), r, lo, hi
R['C_reading'] = 'lacks resolution' if R['C1_share_le2'] >= 0.8 else 'varies'
T['C_offseason'] = pd.DataFrame([{k: v for k, v in R.items() if k.startswith('C')}])

# ---------- D. robust regression ----------
X = sm.add_constant(a.values); y = b.values
o = sm.OLS(y, X).fit(cov_type='HC1'); h = sm.RLM(y, X, M=sm.robust.norms.HuberT(t=1.345)).fit(scale_est='mad', cov='H1', tol=1e-12, maxiter=500)
R['D_ols_slope'], R['D_ols_se'] = o.params[1], o.bse[1]; R['D_ols_lo'], R['D_ols_hi'] = o.params[1] - 1.96 * o.bse[1], o.params[1] + 1.96 * o.bse[1]
R['D_huber_slope'], R['D_huber_se'] = h.params[1], h.bse[1]; R['D_huber_lo'], R['D_huber_hi'] = h.params[1] - 1.96 * h.bse[1], h.params[1] + 1.96 * h.bse[1]
R['D_ratio'] = h.params[1] / o.params[1]; R['D_robust'] = int(np.sign(h.params[1]) == np.sign(o.params[1]) and (R['D_huber_lo'] > 0 or R['D_huber_hi'] < 0))
T['D_robust_regression'] = pd.DataFrame([{k: v for k, v in R.items() if k.startswith('D_')}])

R = {k: (float(v) if isinstance(v, (float, np.floating)) else int(v) if isinstance(v, (int, np.integer, bool, np.bool_)) else v) for k, v in R.items()}
json.dump(R, open('out/21_reviewer_checks.json', 'w'), indent=1)
with pd.ExcelWriter('out/21_reviewer_checks.xlsx') as w:
    for k, t in T.items(): t.to_excel(w, k[:31], index=False)
print(len(R), 'keys')
