"""Implementation 3: written from VERIFICATION_SPEC.md only, with NumPy, SciPy and openpyxl.
No pandas or statsmodels. OLS, HAC, NB2 maximum likelihood, Ljung-Box, Cook's distance and Holm are coded by hand."""
import csv, json, math, calendar, numpy as np, openpyxl
from scipy import optimize, stats, special
D = 'data_raw/'; V = 'verification/'   # run from the project root

# ---------- data ----------
wb = openpyxl.load_workbook(D + 'edcd_district_monthly_2019_2023.xlsx', data_only=True)
num = lambda v: 0.0 if v in (None, '') else float(v)
nat, dist = [], {}
for y in range(2019, 2024):
    rows = list(wb[f'{y} Monthly'].iter_rows(min_row=3, values_only=True))
    tot = [num(v) for v in next(r for r in rows if r[1] == 'Total')[3:15]]
    s = [0.0] * 12
    for r in rows:
        if r[1] in (None, 'Total'): continue
        vals = [num(v) for v in r[3:15]]; s = [a + b for a, b in zip(s, vals)]
        if r[1] != 'Not specified': dist.setdefault(str(r[1]).upper().strip(), {})[y] = vals
    assert s == tot, y
    nat += s
pp = [r for r in openpyxl.load_workbook(D + 'paper_finalized_data_RP.xlsx', data_only=True).worksheets[0].iter_rows(min_row=2, values_only=True) if isinstance(r[1], (int, float))]
pp = sorted(pp, key=lambda r: r[1]); assert [int(r[1]) for r in pp] == list(range(1, 73))
assert [float(r[5]) for r in pp[:60]] == nat
cases = np.array(nat + [float(r[5]) for r in pp[60:]]); topic = np.array([float(r[2]) for r in pp]); fever = np.array([float(r[3]) for r in pp])
def gt(f):
    out = []
    for r in csv.reader(open(D + f, encoding='utf-8-sig')):
        if len(r) == 2 and len(r[0]) == 7 and r[0][:2] == '20': out.append(0.5 if r[1].strip() == '<1' else float(r[1]))
    assert len(out) in (72,), f; return np.array(out)
d1, d2, d3 = gt('gt_dengue_term.csv'), gt('gt_dengue_term_day2_2026-11-22.csv'), gt('gt_dengue_term_day3_2026-11-23.csv')
x = (d1 + d2 + d3) / 3
YR = np.repeat(np.arange(2019, 2025), 12); MO = np.tile(np.arange(1, 13), 6)
ly = lambda v: np.log(v + 1.0); lx = lambda v: np.log(v + 0.5)
def adj(v, mo):
    v = np.asarray(v, float); out = np.empty_like(v)
    for m in np.unique(mo): out[mo == m] = v[mo == m] - v[mo == m].mean()
    return out
pear = lambda a, b: float(np.corrcoef(a, b)[0, 1])
def neff(a, b):
    n = len(a); s = sum(pear(a[:-k], a[k:]) * pear(b[:-k], b[k:]) for k in range(1, 13)); return max(4.0, n / (1 + 2 * s))
def fisher(r, ne): z = math.atanh(r); h = 1.96 / math.sqrt(ne - 3); return math.tanh(z - h), math.tanh(z + h)
R = {}
def win(lo, hi): return (YR * 100 + MO >= lo) & (YR * 100 + MO <= hi)

# ---------- Q1 ----------
runs = {'primary': (x, win(202001, 202412)), 'A': (x, win(201906, 202412)), 'B': (x, win(201901, 202412)), 'C': (x, win(202201, 202412)),
        'D': (x, win(202001, 202412) & (YR != 2022)), 'topic': (topic, win(202001, 202412)), 'day1': (d1, win(202001, 202412)), 'fever': (fever, win(202001, 202412))}
AB = {}
for k, (s, w) in runs.items():
    a, b = adj(lx(s[w]), MO[w]), adj(ly(cases[w]), MO[w]); r = pear(a, b); ne = neff(a, b); lo, hi = fisher(r, ne)
    R[f'q1_{k}_r'], R[f'q1_{k}_neff'], R[f'q1_{k}_flo'], R[f'q1_{k}_fhi'] = r, ne, lo, hi; AB[k] = (a, b)
starts = np.loadtxt(V + 'boot_q1_block_starts.csv', delimiter=',', dtype=int)
for k in ('primary', 'fever'):
    a, b = AB[k]; rs = [pear(a[i], b[i]) for i in (np.concatenate([np.arange(s0, s0 + 12) for s0 in row])[:60] for row in starts)]
    R[f'q1_{k}_blo'], R[f'q1_{k}_bhi'] = float(np.percentile(rs, 2.5)), float(np.percentile(rs, 97.5))

# ---------- OLS / HAC helpers ----------
def ols(y, X): XtXi = np.linalg.inv(X.T @ X); b = XtXi @ X.T @ y; return b, y - X @ b, XtXi
def hac(X, e, XtXi, L=3):
    u = X * e[:, None]; S = u.T @ u
    for j in range(1, L + 1): G = u[j:].T @ u[:-j]; S += (1 - j / (L + 1)) * (G + G.T)
    return XtXi @ S @ XtXi
def dummies(mo): return np.column_stack([(mo == m).astype(float) for m in range(1, 13)])

# ---------- Q2 ----------
w = win(202001, 202412); a, b = AB['primary']; n = len(a)
aic = {}
for p in range(4):
    yv = a[3:]; X = np.column_stack([np.ones(n - 3)] + [a[3 - j:n - j] for j in range(1, p + 1)])
    _, e, _ = ols(yv, X); aic[p] = (n - 3) * math.log((e @ e) / (n - 3)) + 2 * (p + 1)
p = min(aic, key=aic.get); R['q2_ar_order'] = p
X = np.column_stack([np.ones(n - p)] + [a[p - j:n - j] for j in range(1, p + 1)]); coef, _, _ = ols(a[p:], X)
R['q2_ar_const'] = coef[0]
for j in range(1, p + 1): R[f'q2_ar_phi{j}'] = coef[j]
filt = lambda v: np.array([v[t] - coef[0] - sum(coef[j] * v[t - j] for j in range(1, p + 1)) for t in range(p, n)])
ea, eb = filt(a), filt(b); m = len(ea); lab = {-3: 'm3', -2: 'm2', -1: 'm1', 0: '0', 1: 'p1', 2: 'p2', 3: 'p3'}
ps_, rs_, ns_ = [], [], []
for k in range(-3, 4):
    u, v = (ea[:m - k], eb[k:]) if k >= 0 else (ea[-k:], eb[:m + k])
    r = pear(u, v); nk = len(u); tt = r * math.sqrt((nk - 2) / (1 - r * r)); ps_.append(2 * stats.t.sf(abs(tt), nk - 2)); rs_.append(r); ns_.append(nk)
order = np.argsort(ps_); adjp = np.empty(7); run = 0.0
for i, j in enumerate(order): run = max(run, min(1.0, (7 - i) * ps_[j])); adjp[j] = run
for i, k in enumerate(range(-3, 4)): R[f'q2_ccf_r_{lab[k]}'], R[f'q2_ccf_n_{lab[k]}'], R[f'q2_ccf_pholm_{lab[k]}'] = rs_[i], ns_[i], float(adjp[i])
LY, LX = ly(cases), lx(x)
def direction(lo, dep, own, other):
    idx = np.where(win(lo, 202412))[0]; X = np.column_stack([dummies(MO[idx]), own[idx - 1], other[idx - 1]])
    bb, e, XtXi = ols(dep[idx], X); se = math.sqrt(hac(X, e, XtXi)[-1, -1]); return bb[-1], se, 2 * stats.norm.sf(abs(bb[-1] / se))
R['q2_fwd_coef'], R['q2_fwd_se'], R['q2_fwd_p'] = direction(202002, LY, LY, LX)
R['q2_rev_coef'], R['q2_rev_se'], R['q2_rev_p'] = direction(202002, LX, LX, LY)
c58, _, p58 = direction(202003, LX, LX, LY); R['q2_rev58_coef'], R['q2_rev58_p'] = c58, p58

# ---------- NB2 ----------
idx = np.where(win(202002, 202412))[0]; yc = cases[idx]
Z = np.column_stack([np.ones(len(idx)), np.sin(2*np.pi*MO[idx]/12), np.cos(2*np.pi*MO[idx]/12), np.sin(4*np.pi*MO[idx]/12), np.cos(4*np.pi*MO[idx]/12), LY[idx - 1], LX[idx]])
def ll_i(th):
    be, la = th[:-1], th[-1]; al = math.exp(la); mu = np.exp(Z @ be); r_ = 1 / al
    return special.gammaln(yc + r_) - special.gammaln(r_) - special.gammaln(yc + 1) + r_ * np.log(r_ / (r_ + mu)) + yc * np.log(mu / (r_ + mu))
b0, _, _ = ols(np.log(yc + 1), Z)
opt = optimize.minimize(lambda th: -ll_i(th).sum(), np.r_[b0, math.log(0.3)], method='BFGS', options={'gtol': 1e-10, 'maxiter': 10000})
opt = optimize.minimize(lambda th: -ll_i(th).sum(), opt.x, method='Nelder-Mead', options={'xatol': 1e-12, 'fatol': 1e-12, 'maxiter': 200000, 'maxfev': 200000}) if not opt.success else opt
th = opt.x; al = math.exp(th[-1]); bx = th[-2]
R['nb_coef_x'], R['nb_rr'], R['nb_alpha'] = bx, 2 ** bx, al
# HAC sandwich on (beta, alpha) with numerical per-observation scores and Hessian (parameterised in alpha, not log alpha)
def ll_a(t2): return ll_i(np.r_[t2[:-1], math.log(t2[-1])])
t2 = np.r_[th[:-1], al]; h = 1e-6
scores = np.column_stack([(ll_a(t2 + h * np.eye(len(t2))[i]) - ll_a(t2 - h * np.eye(len(t2))[i])) / (2 * h) for i in range(len(t2))])
H = np.zeros((len(t2), len(t2)))
def grad(t_): return np.array([(ll_a(t_ + h * np.eye(len(t_))[i]).sum() - ll_a(t_ - h * np.eye(len(t_))[i]).sum()) / (2 * h) for i in range(len(t_))])
for i in range(len(t2)): ei = np.eye(len(t2))[i] * 1e-4; H[:, i] = (grad(t2 + ei) - grad(t2 - ei)) / 2e-4
Hi = np.linalg.inv(-H); S = scores.T @ scores
for j in range(1, 4): G = scores[j:].T @ scores[:-j]; S += (1 - j / 4) * (G + G.T)
se = math.sqrt((Hi @ S @ Hi)[-2, -2]); R['nb_rr_lo'], R['nb_rr_hi'] = 2 ** (bx - 1.96 * se), 2 ** (bx + 1.96 * se)
mu = np.exp(Z @ th[:-1]); pr = (yc - mu) / np.sqrt(mu + al * mu ** 2); prc = pr - pr.mean(); N = len(pr)
acf = [(prc[:-k] @ prc[k:]) / (prc @ prc) for k in range(1, 13)]; Q = N * (N + 2) * sum(r ** 2 / (N - k) for k, r in enumerate(acf, 1))
R['nb_lb_q12'], R['nb_lb_p12'] = Q, float(stats.chi2.sf(Q, 12))
bo, e, XtXi = ols(np.log(yc + 1), Z); hdiag = np.einsum('ij,jk,ik->i', Z, XtXi, Z); s2 = (e @ e) / (N - Z.shape[1])
cook = e ** 2 / (Z.shape[1] * s2) * hdiag / (1 - hdiag) ** 2
R['nb_cook_months'] = [f'{YR[i]}-{MO[i]:02d}' for i, cd in zip(idx, cook) if cd > 4 / N]

# ---------- Q3 ----------
ii = np.where(win(201907, 202412))[0]
def design(rows, M):
    cols = [dummies(MO[rows]), LY[rows - 1][:, None]]
    if M == 'M2': cols.append(LX[rows][:, None])
    if M == 'M3': cols.append(LX[rows - 1][:, None])
    return np.hstack(cols)
def predict(train, test, M):
    bb, _, _ = ols(LY[train], design(train, M)); return np.maximum(np.expm1(design(test, M) @ bb), 0)
iid = np.loadtxt(V + 'boot_q3_iid_idx.csv', delimiter=',', dtype=int); b3 = np.loadtxt(V + 'boot_q3_block3_starts.csv', delimiter=',', dtype=int)
def evaluate(s, tests, trainer):
    P = {M: [] for M in ('M0', 'M1', 'M2', 'M3')}; obs = []
    for grp in tests:
        tr = trainer(grp)
        P['M0'] += list(cases[grp - 1]); obs += list(cases[grp])
        for M in ('M1', 'M2', 'M3'): P[M] += list(predict(tr, grp, M))
    obs = np.array(obs); P = {k: np.array(v) for k, v in P.items()}; yrs = np.concatenate([YR[g] for g in tests])
    for M, pv in P.items():
        e = pv - obs; R[f'q3_{s}_mae_{M}'] = float(np.abs(e).mean()); R[f'q3_{s}_rmse_{M}'] = float(np.sqrt((e ** 2).mean()))
        R[f'q3_{s}_maelog_{M}'] = float(np.abs(np.log1p(pv) - np.log1p(obs)).mean())
        if s != 'loso':
            for th_ in (1000, 500):
                o, f_ = obs > th_, pv > th_; R[f'q3_{s}_sens_{M}_{th_}'] = (o & f_).sum() / o.sum(); R[f'q3_{s}_spec_{M}_{th_}'] = (~o & ~f_).sum() / (~o).sum()
    if s == 'loso': return
    better = 0
    for yv in (2022, 2023, 2024):
        m1, m2 = np.abs(P['M1'] - obs)[yrs == yv].mean(), np.abs(P['M2'] - obs)[yrs == yv].mean()
        R[f'q3_{s}_season_mae_M1_{yv}'], R[f'q3_{s}_season_mae_M2_{yv}'] = m1, m2; better += m2 < m1
    dd = np.abs(P['M2'] - obs) - np.abs(P['M1'] - obs); R[f'q3_{s}_seasons_m2_better'] = int(better); R[f'q3_{s}_diff_mean'] = float(dd.mean())
    bm = [dd[r].mean() for r in iid]; R[f'q3_{s}_diff_lo'], R[f'q3_{s}_diff_hi'] = float(np.percentile(bm, 2.5)), float(np.percentile(bm, 97.5))
    bm3 = [dd[np.concatenate([np.arange(s0, s0 + 3) for s0 in r])].mean() for r in b3]
    R[f'q3_{s}_diff_b3lo'], R[f'q3_{s}_diff_b3hi'] = float(np.percentile(bm3, 2.5)), float(np.percentile(bm3, 97.5))
    R[f'q3_{s}_rule_met'] = int(better >= 2 and R[f'q3_{s}_diff_hi'] < 0)
test_months = np.where(win(202201, 202412))[0]
evaluate('ps', [test_months[YR[test_months] == yv] for yv in (2022, 2023, 2024)], lambda g: ii[ii < g.min()])
evaluate('mo', [np.array([t_]) for t_ in test_months], lambda g: ii[ii < g[0]])
lo_idx = np.where(win(202002, 202412))[0]
evaluate('loso', [lo_idx[YR[lo_idx] == yv] for yv in range(2020, 2025)], lambda g: lo_idx[YR[lo_idx] != YR[g[0]]])

# ---------- Q4 ----------
for yv in range(2019, 2025):
    k = YR == yv; R[f'q4_{yv}_peak_ratio'] = cases[k].max() / x[k].max(); R[f'q4_{yv}_total_ratio'] = cases[k].sum() / x[k].sum()

# ---------- Q5 ----------
REG = {'Eastern': 'TAPLEJUNG PANCHTHAR ILAM JHAPA SANKHUWASABHA TERHATHUM DHANKUTA BHOJPUR MORANG SUNSARI SOLUKHUMBU OKHALDHUNGA KHOTANG UDAYAPUR SAPTARI SIRAHA',
       'Central': 'DOLAKHA RAMECHHAP SINDHULI DHANUSA MAHOTTARI SARLAHI SINDHUPALCHOK KAVRE LALITPUR BHAKTAPUR KATHMANDU NUWAKOT RASUWA DHADING MAKAWANPUR MAKWANPUR RAUTAHAT BARA PARSA CHITAWAN CHITWAN',
       'Western': 'GORKHA LAMJUNG TANAHU SYANGJA KASKI MANANG MUSTANG MYAGDI PARBAT BAGLUNG GULMI PALPA NAWALPARASI RUPANDEHI KAPILBASTU KAPILVASTU ARGHAKHANCHI',
       'MidWestern': 'PYUTHAN ROLPA RUKUM SALYAN DANG DAILEKH JAJARKOT SURKHET BANKE BARDIYA DOLPA JUMLA KALIKOT MUGU HUMLA',
       'FarWestern': 'BAJURA BAJHANG ACHHAM DOTI KAILALI DARCHULA BAITADI DADELDHURA KANCHANPUR'}
def region_of(name):
    tok = name.replace('(', ' ').replace(')', ' ').split()
    hits = {rg for rg, ds in REG.items() for d_ in ds.split() if any(t_.startswith(d_) or d_.startswith(t_) and len(t_) >= 4 for t_ in tok)}
    return hits
assign = {nm: region_of(nm) for nm in dist}
bad = {k: v for k, v in assign.items() if len(v) != 1}; assert not bad, bad
assert len(assign) == 77, len(assign)
mo48 = np.tile(np.arange(1, 13), 4)
for rg in ('Central', 'Eastern', 'Western'):
    rc = np.array([sum(dist[nm][y][mm] for nm in dist if assign[nm] == {rg}) for y in range(2020, 2024) for mm in range(12)])
    sr = []
    for r_ in csv.reader(open(D + f'gt_dengue_region_{rg.lower()}.csv', encoding='utf-8-sig')):
        if len(r_) == 2 and len(r_[0]) == 7 and r_[0][:2] == '20' and 2020 <= int(r_[0][:4]) <= 2023: sr.append(0.5 if r_[1] == '<1' else float(r_[1]))
    sr = np.array(sr); a_, b_ = adj(lx(sr), mo48), adj(ly(rc), mo48); r = pear(a_, b_); ne = neff(a_, b_); lo, hi = fisher(r, ne)
    R[f'q5_{rg}_r'], R[f'q5_{rg}_neff'], R[f'q5_{rg}_flo'], R[f'q5_{rg}_fhi'] = r, ne, lo, hi
    R[f'q5_{rg}_zeros2020_2023'] = int((sr == 0).sum())
    R[f'q5_{rg}_peakratio_cases'] = rc[36:48].max() / rc[24:36].max(); R[f'q5_{rg}_peakratio_search'] = sr[36:48].max() / sr[24:36].max()
    globals()[f'RC_{rg}'], globals()[f'SR_{rg}'] = rc, sr
# peaks
series = [(cases, range(2020, 2025), 2019), (x, range(2020, 2025), 2019)]
series += [(globals()[f'RC_{rg}'], range(2020, 2024), 2020) for rg in ('Central', 'Eastern', 'Western')]
series += [(globals()[f'SR_{rg}'], range(2020, 2024), 2020) for rg in ('Central', 'Eastern', 'Western')]
labels = []
for v, yrs, y0 in series:
    for yv in yrs:
        seg = v[(yv - y0) * 12:(yv - y0) * 12 + 12]; o = np.argsort(-seg, kind='stable'); top, sec = seg[o[0]], seg[o[1]]
        ratio = math.inf if sec == 0 else top / sec
        labels.append([calendar.month_abbr[o[0] + 1] if ratio >= 1 + g else 'tie' for g in (0.05, 0.10, 0.20)])
assert len(labels) == 34
R['peak_n_changed'] = sum(len(set(l)) > 1 for l in labels)
for i, g in enumerate(('05', '10', '20')): R[f'peak_ties_{g}'] = sum(l[i] == 'tie' for l in labels)

# ---------- Q6 ----------
c24 = cases[60:]
for ds, pub in [('20240731', 2000), ('20240903', 6000), ('20240923', 12644), ('20241015', 19599), ('20241025', 27153)]:
    mth, day = int(ds[4:6]), int(ds[6:]); fin = c24[:mth - 1].sum() + c24[mth - 1] * day / calendar.monthrange(2024, mth)[1]
    R[f'q6_final_{ds}'] = float(fin); R[f'q6_pct_{ds}'] = 100 * pub / fin

R = {k: (float(v) if isinstance(v, (float, np.floating)) else int(v) if isinstance(v, (int, np.integer, bool, np.bool_)) else v) for k, v in R.items()}
json.dump(R, open(V + 'results_impl3_numpy.json', 'w'), indent=1); print(len(R), 'keys')
