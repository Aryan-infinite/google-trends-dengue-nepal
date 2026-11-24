"""Final analysis, following ANALYSIS_PLAN.md v1.0 (departures from the plan are listed in the paper's Methods). Outputs in out/final/."""
import os, json, warnings
import numpy as np, pandas as pd, scipy.stats as st
import statsmodels.api as sm, statsmodels.formula.api as smf
from statsmodels.tsa.ar_model import ar_select_order
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.diagnostic import acorr_ljungbox
warnings.filterwarnings('ignore')
SEED = 20261123; B = 4000
OUT = 'out/final'; os.makedirs(OUT, exist_ok=True)
IDX = pd.date_range('2019-01-01', periods=72, freq='MS')

# ---------- data ----------
def gt(path, lt=0.5):
    s = pd.read_csv(path, skiprows=2).iloc[:, 1].replace('<1', lt).astype(float).values
    return pd.Series(s, index=IDX)
paper = pd.read_csv('out/paper_72month_table.csv', index_col=0, parse_dates=True)
cases = paper.cases.astype(float); cases.index = IDX
days = [gt('data_raw/gt_dengue_term.csv'), gt('data_raw/gt_dengue_term_day2_2026-11-22.csv'), gt('data_raw/gt_dengue_term_day3_2026-11-23.csv')]
days0 = [gt(f, 0) for f in ['data_raw/gt_dengue_term.csv', 'data_raw/gt_dengue_term_day2_2026-11-22.csv', 'data_raw/gt_dengue_term_day3_2026-11-23.csv']]
S = {'dengue_3day': sum(days) / 3, 'dengue_day1': days[0], 'dengue_lt1_zero': sum(days0) / 3,
     'topic': pd.Series(paper.DF.values, IDX), 'fever': pd.Series(paper.Fever.values, IDX)}
# Regional case totals come from the EDCD district workbook (scripts/04), which is not public. Without that file the
# regional analysis (Q5) is skipped; nothing else depends on it.
REG_FILE = 'out/edcd_devregion_monthly_2019_2023.csv'
regC = pd.read_csv(REG_FILE, index_col=0, parse_dates=True) if os.path.exists(REG_FILE) else None
regS = {n: gt(f'data_raw/gt_dengue_region_{f}.csv') for f, n in
        [('central', 'Central'), ('eastern', 'Eastern'), ('western', 'Western'), ('midwestern', 'Mid-Western'), ('farwestern', 'Far-Western')]}

ly = lambda c: np.log(c + 1)
lx = lambda s: np.log(s + 0.5)
adj = lambda s: s - s.groupby(s.index.month).transform('mean')

# ---------- helpers ----------
def neff(a, b, K=12):
    a = np.asarray(a) - np.mean(a); b = np.asarray(b) - np.mean(b)
    s = sum(np.corrcoef(a[:-k], a[k:])[0, 1] * np.corrcoef(b[:-k], b[k:])[0, 1] for k in range(1, K + 1))
    return max(len(a) / (1 + 2 * s), 4.0)
def fisher_ci(r, n):
    z, se = np.arctanh(r), 1 / np.sqrt(max(n - 3, 1)); return np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
def block_boot_ci(a, b, L=12, seed=SEED):
    rng = np.random.default_rng(seed); a = np.asarray(a); b = np.asarray(b); n = len(a); rs = []
    starts = np.arange(n - L + 1)
    for _ in range(B):
        idx = np.concatenate([np.arange(s, s + L) for s in rng.choice(starts, int(np.ceil(n / L)))])[:n]
        if np.std(a[idx]) > 0 and np.std(b[idx]) > 0: rs.append(np.corrcoef(a[idx], b[idx])[0, 1])
    return np.percentile(rs, [2.5, 97.5])
def corr_block(x, y):
    r = np.corrcoef(x, y)[0, 1]; ne = neff(x, y); f = fisher_ci(r, ne); bb = block_boot_ci(x, y)
    lo, hi = min(f[0], bb[0]), max(f[1], bb[1])
    return dict(r=r, n=len(x), n_eff=ne, fisher_lo=f[0], fisher_hi=f[1], boot_lo=bb[0], boot_hi=bb[1], ci_lo=lo, ci_hi=hi)

def q1(xs, cs):
    X, Y = lx(xs), ly(cs)
    out = {'primary_adjusted_log': corr_block(adj(X).values, adj(Y).values)}
    out['unadjusted_log'] = corr_block(X.values, Y.values)
    r, p = st.pearsonr(xs, cs); rho, ps = st.spearmanr(xs, cs)
    out['raw_pearson'] = dict(r=r, p=p); out['raw_spearman'] = dict(rho=rho, p=ps)
    return out

def q2(xs, cs):
    X, Y = adj(lx(xs)), adj(ly(cs))
    sel = ar_select_order(X.values, maxlag=3, ic='aic', trend='c'); lags = sel.ar_lags or []
    p = max(lags) if lags else 0
    if p:
        fit = sel.model.fit(); c = fit.params[0]; phi = dict(zip(lags, fit.params[1:]))
        filt = lambda s: (s - c - sum(phi[k] * s.shift(k) for k in lags)).iloc[p:]
        ex, ey = filt(X), filt(Y)
    else:
        ex, ey = X, Y
    rows = []
    for k in range(-3, 4):
        a = ex.shift(k); m = a.notna() & ey.notna(); r, pv = st.pearsonr(a[m], ey[m]); rows.append(dict(lag=k, r=r, n=int(m.sum()), p=pv))
    t = pd.DataFrame(rows); t['p_holm'] = multipletests(t.p, method='holm')[1]; t['ar_order'] = p
    d = pd.DataFrame({'y': ly(cs), 'y1': ly(cs).shift(1), 'x': lx(xs), 'x1': lx(xs).shift(1), 'm': cs.index.month}).dropna()
    w1 = smf.ols('y ~ C(m) + y1 + x1', d).fit(cov_type='HAC', cov_kwds={'maxlags': 3})
    dd = d.assign(xl=d.x.shift(1)).dropna()
    w2 = smf.ols('x ~ C(m) + xl + y1', dd).fit(cov_type='HAC', cov_kwds={'maxlags': 3})
    g = dict(search_to_cases_coef=w1.params['x1'], search_to_cases_p=w1.pvalues['x1'],
             cases_to_search_coef=w2.params['y1'], cases_to_search_p=w2.pvalues['y1'])
    return t, g

def nowcast(xs, cs, scheme, tests, train_start='2019-07-01'):
    d = pd.DataFrame({'y': ly(cs), 'y1': ly(cs).shift(1), 'x': lx(xs), 'x1': lx(xs).shift(1), 'm': cs.index.month, 'c': cs})
    d = d[d.index >= train_start].dropna()
    forms = {'M1': 'y ~ C(m) + y1', 'M2': 'y ~ C(m) + y1 + x', 'M3': 'y ~ C(m) + y1 + x1'}
    preds = []
    for yr in tests:
        te = d.index.year == yr
        tr = (d.index.year < yr) if scheme == 'rolling' else ~te
        P = pd.DataFrame({'year': yr, 'obs': d.c[te], 'M0': np.expm1(d.y1[te])}, index=d.index[te])
        for k, f in forms.items():
            fit = smf.ols(f, d[tr]).fit(); P[k] = np.expm1(fit.predict(d[te]))
        preds.append(P)
    P = pd.concat(preds); P[['M0', 'M1', 'M2', 'M3']] = P[['M0', 'M1', 'M2', 'M3']].clip(lower=0)
    met = []
    for k in ['M0', 'M1', 'M2', 'M3']:
        e = P[k] - P.obs
        for thr in [1000, 500]:
            o, pr = P.obs > thr, P[k] > thr
            sens = (o & pr).sum() / o.sum() if o.sum() else np.nan; spec = (~o & ~pr).sum() / (~o).sum()
            if thr == 1000: s1000, sp1000 = sens, spec
            else: s500, sp500 = sens, spec
        met.append(dict(model=k, MAE=e.abs().mean(), RMSE=np.sqrt((e ** 2).mean()),
                        MAE_log=(np.log1p(P[k]) - np.log1p(P.obs)).abs().mean(),
                        sens_1000=s1000, spec_1000=sp1000, sens_500=s500, spec_500=sp500))
    met = pd.DataFrame(met)
    per = P.groupby('year').apply(lambda g: pd.Series({k: (g[k] - g.obs).abs().mean() for k in ['M0', 'M1', 'M2', 'M3']}))
    diff = (P.M2 - P.obs).abs().values - (P.M1 - P.obs).abs().values
    rng = np.random.default_rng(SEED); bs = [rng.choice(diff, len(diff)).mean() for _ in range(B)]
    pdiff = dict(mean=diff.mean(), lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5),
                 seasons_M2_better=int((per.M2 < per.M1).sum()), n_seasons=len(per))
    return P, met, per, pdiff

def nb_model(xs, cs):
    d = pd.DataFrame({'c': cs, 'y1': ly(cs).shift(1), 'x': lx(xs), 't': np.arange(len(cs))})
    d['s1'], d['c1'] = np.sin(2 * np.pi * d.index.month / 12), np.cos(2 * np.pi * d.index.month / 12)
    d['s2'], d['c2'] = np.sin(4 * np.pi * d.index.month / 12), np.cos(4 * np.pi * d.index.month / 12)
    d = d.dropna()
    m = smf.negativebinomial('c ~ s1 + c1 + s2 + c2 + y1 + x', d).fit(disp=0, maxiter=500)
    mh = smf.negativebinomial('c ~ s1 + c1 + s2 + c2 + y1 + x', d).fit(disp=0, maxiter=500, cov_type='HAC', cov_kwds={'maxlags': 3})
    b, se = mh.params['x'], mh.bse['x']
    mu = m.predict(d); alpha = m.params['alpha']; pr = (d.c - mu) / np.sqrt(mu + alpha * mu ** 2)
    lb = acorr_ljungbox(pr, lags=[12], return_df=True)
    ols = smf.ols('np.log1p(c) ~ s1 + c1 + s2 + c2 + y1 + x', d).fit()
    cook = ols.get_influence().cooks_distance[0]
    return dict(IRR_per_doubling=np.exp(b * np.log(2)), lo=np.exp((b - 1.96 * se) * np.log(2)), hi=np.exp((b + 1.96 * se) * np.log(2)),
                p_HAC=mh.pvalues['x'], alpha=alpha, ljungbox12_p=float(lb.lb_pvalue.iloc[0]),
                n=len(d), n_cook_flagged=int((cook > 4 / len(d)).sum()),
                cook_flagged=[i.strftime('%Y-%m') for i in d.index[cook > 4 / len(d)]])

# ---------- windows ----------
W = {'primary_2020_2024': ('2020-01', '2024-12'), 'A_with_2019JunDec': ('2019-06', '2024-12'), 'B_full72': ('2019-01', '2024-12')}
def sub(s, w): return s[W[w][0]:W[w][1]] if w in W else s
results = {}
x, c = S['dengue_3day'], cases

# Q1 & Q2 primary
P0 = sub(x, 'primary_2020_2024'), sub(c, 'primary_2020_2024')
results['Q1_primary'] = q1(*P0)
q2t, q2g = q2(*P0); results['Q2_wald'] = q2g
# Q3
Pr, met_r, per_r, d_r = nowcast(x, c, 'rolling', [2022, 2023, 2024])
Pl, met_l, per_l, d_l = nowcast(x, c, 'loso', [2020, 2021, 2022, 2023, 2024], train_start='2020-02-01')
results['Q3_rolling_pairdiff_M2_minus_M1'] = d_r; results['Q3_loso_pairdiff_M2_minus_M1'] = d_l
results['NB_primary'] = nb_model(*P0)

# Q4
q4 = []
for yr in range(2019, 2025):
    cy, xy = c[str(yr)], x[str(yr)]
    q4.append(dict(season=yr, cases_total=int(cy.sum()), peak_cases=int(cy.max()), peak_case_month=cy.idxmax().strftime('%b'),
                   peak_search=round(xy.max(), 2), peak_search_month=xy.idxmax().strftime('%b'),
                   cases_per_peak_search_point=round(cy.max() / xy.max()), cases_per_search_unit_total=round(cy.sum() / xy.sum())))
q4 = pd.DataFrame(q4)

# Q5 regional
q5 = []; q5peak = []
for n, s in (regS.items() if regC is not None else []):
    zeros = int((s == 0).sum()); row = dict(region=n, zero_months=zeros, used=zeros < 30)
    if zeros < 30:
        xs, cs = s['2020':'2023'], regC[n]['2020':]
        row.update({f'adj_log_{k}': v for k, v in corr_block(adj(lx(xs)).values, adj(ly(cs)).values).items() if k in ('r', 'n_eff', 'ci_lo', 'ci_hi')})
        row['cases_2023_vs_2022_peak'] = regC[n]['2023'].max() / regC[n]['2022'].max()
        row['search_2023_vs_2022_peak'] = s['2023'].max() / s['2022'].max()
        for yr in [2022, 2023]:
            def pk(z):
                v = z.sort_values(ascending=False); return v.index[0].strftime('%b') if v.iloc[0] >= 1.10 * v.iloc[1] else f'tie ({v.index[0].strftime("%b")}/{v.index[1].strftime("%b")})'
            q5peak.append(dict(region=n, season=yr, case_peak=pk(regC[n][str(yr)]), search_peak=pk(s[str(yr)])))
    q5.append(row)
q5 = pd.DataFrame(q5); q5peak = pd.DataFrame(q5peak)
natpk = []
for yr in range(2019, 2025):
    def pk(z):
        v = z.sort_values(ascending=False); return v.index[0].strftime('%b') if v.iloc[0] >= 1.10 * v.iloc[1] else f'tie ({v.index[0].strftime("%b")}/{v.index[1].strftime("%b")})'
    natpk.append(dict(season=yr, case_peak=pk(c[str(yr)]), search_peak=pk(x[str(yr)])))
natpk = pd.DataFrame(natpk)

# Robustness grid
grid = []
def core(label, xs, cs, rolling=True):
    a = corr_block(adj(lx(xs)).values, adj(ly(cs)).values)
    t, g = q2(xs, cs)
    row = dict(run=label, n=len(cs), adj_r=a['r'], ci_lo=a['ci_lo'], ci_hi=a['ci_hi'],
               lag_minus1=t.loc[t.lag == -1, 'r'].item(), lag0=t.loc[t.lag == 0, 'r'].item(), lag_plus1=t.loc[t.lag == 1, 'r'].item(),
               lag_plus1_holm=t.loc[t.lag == 1, 'p_holm'].item(), wald_x1_p=g['search_to_cases_p'])
    return row
grid.append(core('primary (2020–24, 3-day mean)', *P0))
grid.append(core('A: 2019 Jun–2024', sub(x, 'A_with_2019JunDec'), sub(c, 'A_with_2019JunDec')))
grid.append(core('B: full 72 months', x, c))
grid.append(core('C: excl. 2020–21 (COVID), 2022–24', x['2022':], c['2022':]))
grid.append(core('D: excl. 2022', pd.concat([x['2020':'2021'], x['2023':]]), pd.concat([c['2020':'2021'], c['2023':]])))
grid.append(core('Topic instead of term', sub(S['topic'], 'primary_2020_2024'), P0[1]))
grid.append(core('Day-1 download only', sub(S['dengue_day1'], 'primary_2020_2024'), P0[1]))
grid.append(core('"<1" coded 0', sub(S['dengue_lt1_zero'], 'primary_2020_2024'), P0[1]))
grid.append(core('Control: "Fever"', sub(S['fever'], 'primary_2020_2024'), P0[1]))
grid = pd.DataFrame(grid)
# Q3 robustness: rolling M2 vs M1 with topic and day1
q3rob = []
for lab, s in [('dengue 3-day', x), ('Topic', S['topic']), ('day-1', S['dengue_day1']), ('Fever control', S['fever'])]:
    _, m, per, dd = nowcast(s, c, 'rolling', [2022, 2023, 2024])
    q3rob.append(dict(series=lab, MAE_M1=m.loc[m.model == 'M1', 'MAE'].item(), MAE_M2=m.loc[m.model == 'M2', 'MAE'].item(),
                      MAE_M3=m.loc[m.model == 'M3', 'MAE'].item(), pairdiff=dd['mean'], lo=dd['lo'], hi=dd['hi'], seasons_M2_better=dd['seasons_M2_better']))
q3rob = pd.DataFrame(q3rob)

# Q6
q6 = pd.read_csv('out/06_2024_reporting_lag.csv')

# ---------- decision rules ----------
t1 = q2t.set_index('lag')
lead_claim = bool(t1.loc[1, 'p_holm'] < 0.05 and t1.loc[1, 'r'] > t1.loc[-1, 'r'] and q2g['search_to_cases_p'] < 0.05)
nowcast_claim = bool(d_r['seasons_M2_better'] >= 2 and (d_r['hi'] < 0))
robust_assoc = bool((grid.iloc[:8].adj_r > 0).all())
results['decisions'] = dict(lead_claim_allowed=lead_claim, nowcast_claim_allowed=nowcast_claim, association_robust_direction=robust_assoc)

# ---------- write ----------
def flat(d, pre=''):
    o = {}
    for k, v in d.items():
        if isinstance(v, dict): o.update(flat(v, pre + k + '.'))
        else: o[pre + k] = v
    return o
summ = pd.DataFrame(list(flat(results).items()), columns=['item', 'value'])
with pd.ExcelWriter(f'{OUT}/final_results.xlsx') as w:
    summ.to_excel(w, 'Summary', index=False)
    q2t.to_excel(w, 'Q2_prewhitened_CCF', index=False)
    met_r.to_excel(w, 'Q3_rolling_metrics', index=False); per_r.to_excel(w, 'Q3_rolling_by_season')
    Pr.to_excel(w, 'Q3_rolling_predictions'); met_l.to_excel(w, 'Q3_LOSO_metrics', index=False)
    q3rob.to_excel(w, 'Q3_robustness', index=False)
    q4.to_excel(w, 'Q4_scale_by_season', index=False); natpk.to_excel(w, 'Q4_peak_months_national', index=False)
    if regC is not None: q5.to_excel(w, 'Q5_regional', index=False); q5peak.to_excel(w, 'Q5_regional_peaks', index=False)
    q6.to_excel(w, 'Q6_2024_reporting_lag', index=False)
    grid.to_excel(w, 'Robustness_grid', index=False)
    pd.read_csv('out/08_external_verification.csv').to_excel(w, 'Case_data_verification', index=False)
    pd.read_csv('out/gt_dengue_3day.csv').to_excel(w, 'GT_3day_stability', index=False)
pd.set_option('display.width', 220); pd.set_option('display.max_columns', 30)
print(json.dumps(flat(results), indent=1, default=lambda v: round(float(v), 4) if isinstance(v, (float, np.floating)) else str(v)))
for nm, t in [('Q2', q2t), ('Q3 rolling', met_r), ('Q3 rolling per season', per_r), ('Q3 LOSO', met_l), ('Q3 robustness', q3rob), ('Q4', q4),
              ('national peaks', natpk), ('Q5', q5), ('Q5 peaks', q5peak), ('Robustness grid', grid)]:
    print('\n==', nm); print(t.round(3).to_string())
