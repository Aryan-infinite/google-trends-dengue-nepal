"""POST HOC (not in the analysis plan): AR-order sensitivity of the prewhitened cross-correlations (AR(1), AR(2), AR(3)).
R twin: scripts/37_ccf_ar_order_sensitivity.R. Needs out/paper_72month_table.csv from scripts/00_prepare_inputs.py."""
import numpy as np, pandas as pd, scipy.stats as st
from statsmodels.stats.multitest import multipletests
from statsmodels.tsa.ar_model import AutoReg
p = pd.read_csv('out/paper_72month_table.csv', index_col=0, parse_dates=True); c = p.cases.astype(float).values[12:]
r = lambda f: pd.read_csv(f, skiprows=2).iloc[:, 1].replace('<1', 0.5).astype(float).values
s = np.mean([r('data_raw/gt_dengue_term.csv'), r('data_raw/gt_dengue_term_day2_2026-11-22.csv'), r('data_raw/gt_dengue_term_day3_2026-11-23.csv')], 0)[12:]
mon = np.tile(np.arange(1, 13), 5); dm = lambda v: v - np.array([v[mon == k].mean() for k in mon])
a = pd.Series(dm(np.log(s + .5))); b = pd.Series(dm(np.log(c + 1)))
for P in (1, 2, 3):
    f = AutoReg(a.values, lags=P, trend='c').fit(); cst = f.params[0]; ph = f.params[1:]
    flt = lambda v: (v - cst - sum(ph[j - 1] * v.shift(j) for j in range(1, P + 1))).iloc[P:]
    ex, ey = flt(a), flt(b); rr, pp = [], []
    for k in range(-3, 4):
        q = ex.shift(k); m = q.notna() & ey.notna(); R, pv = st.pearsonr(q[m], ey[m]); rr.append(R); pp.append(pv)
    h = multipletests(pp, method='holm')[1]; nz = np.delete(h, 3); lagnz = np.delete(np.arange(-3, 4), 3)
    print(f'PY AR({P}): r = ' + ' '.join(f'{x:.3f}' for x in rr) +
          f' | Holm p lag0 = {h[3]:.4f}, lag+1 = {h[4]:.3f}, min non-zero = {nz.min():.3f} at lag {lagnz[nz.argmin()]}')
