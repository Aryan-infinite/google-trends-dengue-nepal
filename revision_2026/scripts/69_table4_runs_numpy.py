"""Table 4, third implementation (NumPy from scratch; scipy.special only for the t and normal tail areas).
Same runs and conventions as scripts/68 (R) and the pipeline (scripts/09, statsmodels): OLS by numpy.linalg.lstsq,
Newey-West (Bartlett, 3 lags, no small-sample correction) and Holm coded by hand; csv module for input.
Output: out/table4_runs_np.csv"""
import csv
import numpy as np
from scipy.special import betainc, ndtr

def rd(f, lt=0.5):
    return np.array([lt if r[1].strip() == '<1' else float(r[1]) for r in csv.reader(open('data_raw/' + f, encoding='utf-8-sig'))
                     if len(r) == 2 and len(r[0]) == 7 and r[0][:2] == '20'])
F3 = ['gt_dengue_term.csv', 'gt_dengue_term_day2_2026-11-22.csv', 'gt_dengue_term_day3_2026-11-23.csv']
x = sum(rd(f) for f in F3) / 3; x0 = sum(rd(f, 0.0) for f in F3) / 3; d1 = rd(F3[0])
T = list(csv.DictReader(open('out/paper_72month_table.csv', encoding='utf-8'))); assert len(T) == 72
cases = np.array([float(r['cases']) for r in T]); topic = np.array([float(r['DF']) for r in T]); fever = np.array([float(r['Fever']) for r in T])
YR = np.repeat(np.arange(2019, 2025), 12); MO = np.tile(np.arange(1, 13), 6); YM = YR * 100 + MO
W = lambda lo, hi: np.where((YM >= lo) & (YM <= hi))[0]
P = W(202001, 202412)
RUNS = [('primary', x, P), ('A', x, W(201906, 202412)), ('B', x, W(201901, 202412)), ('C', x, W(202201, 202412)), ('D', x, P[YR[P] != 2022]),
        ('topic', topic, P), ('day1', d1, P), ('lt0', x0, P), ('fever', fever, P)]

def adj(v, m): return v - np.array([v[m == k].mean() for k in m])
def cor(u, v): u = u - u.mean(); v = v - v.mean(); return float(u @ v / np.sqrt((u @ u) * (v @ v)))
def neff(a, b):
    s = sum(cor(a[:-k], a[k:]) * cor(b[:-k], b[k:]) for k in range(1, 13)); return max(len(a) / (1 + 2 * s), 4.0)
def ols(X, y):
    b = np.linalg.lstsq(X, y, rcond=None)[0]; return b, y - X @ b
def t_p(r, n):                                   # two-sided p for a Pearson r: t = r sqrt((n-2)/(1-r^2)), df = n-2
    df = n - 2; t = r * np.sqrt(df / (1 - r * r)); return float(betainc(df / 2, 0.5, df / (df + t * t)))
def holm(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); adjp = np.empty(m); run = 0.0
    for i, j in enumerate(o): run = max(run, min(1.0, (m - i) * p[j])); adjp[j] = run
    return adjp
def nw_se(X, e, L=3):
    n = len(e); XtXi = np.linalg.inv(X.T @ X); u = X * e[:, None]; S = u.T @ u
    for l in range(1, L + 1):
        w = 1 - l / (L + 1); G = u[l:].T @ u[:-l]; S += w * (G + G.T)
    return np.sqrt(np.diag(XtXi @ S @ XtXi))

rows = []
for name, s_all, i in RUNS:
    s, c, mo = s_all[i], cases[i], MO[i]
    a, b = adj(np.log(s + 0.5), mo), adj(np.log(c + 1), mo); n = len(a)
    r = cor(a, b); ne = neff(a, b); z = np.arctanh(r); lo, hi = np.tanh(z - 1.96 / np.sqrt(ne - 3)), np.tanh(z + 1.96 / np.sqrt(ne - 3))
    aic = []
    for p in range(4):                                   # common sample: observations 4..n
        y = a[3:]; X = np.column_stack([np.ones(n - 3)] + [a[3 - j:n - j] for j in range(1, p + 1)])
        _, e = ols(X, y); aic.append((n - 3) * np.log(e @ e / (n - 3)) + 2 * (p + 1))
    p = int(np.argmin(aic))
    if p:
        X = np.column_stack([np.ones(n - p)] + [a[p - j:n - j] for j in range(1, p + 1)]); cf, _ = ols(X, a[p:])
        filt = lambda v: np.array([v[t] - cf[0] - sum(cf[j] * v[t - j] for j in range(1, p + 1)) for t in range(p, n)])
        ea, eb = filt(a), filt(b)
    else: ea, eb = a, b
    m = len(ea); rr, nn, pv = [], [], []
    for L in range(-3, 4):
        u, v = (ea[:m - L], eb[L:]) if L >= 0 else (ea[-L:], eb[:m + L])
        rr.append(cor(u, v)); nn.append(len(u)); pv.append(t_p(rr[-1], len(u)))
    ph = holm(pv)
    LY, LX = np.log(c + 1), np.log(s + 0.5)
    X = np.column_stack([(mo[1:] == k).astype(float) for k in range(1, 13)] + [LY[:-1], LX[:-1]])
    keep = X.any(0); X = X[:, keep]                       # drop calendar months absent from the run (none here)
    bw, e = ols(X, LY[1:]); se = nw_se(X, e)
    rows.append(dict(run=name, n=n, r=r, n_eff=ne, ci_lo=lo, ci_hi=hi, ar_order=p, lag_m1=rr[2], lag_0=rr[3], lag_p1=rr[4],
                     holm_m1=ph[2], holm_0=ph[3], holm_p1=ph[4], n_lag0=nn[3], fwd_coef=bw[-1], fwd_p=2 * ndtr(-abs(bw[-1] / se[-1])), fwd_n=n - 1))
with open('out/table4_runs_np.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader()
    for r_ in rows: w.writerow({k: (repr(float(v)) if isinstance(v, (float, np.floating)) else v) for k, v in r_.items()})
for r_ in rows: print(r_['run'], r_['n'], f"{r_['r']:.4f}", r_['ar_order'], f"{r_['lag_p1']:.4f}", f"{r_['holm_p1']:.4f}", f"{r_['fwd_coef']:+.4f}", f"{r_['fwd_p']:.5f}")
