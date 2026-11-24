"""Implementation 1 = the analysis pipeline itself (scripts/09 functions, scripts/15, scripts/19, out/06, out/20).
This adapter only calls the pipeline's own functions and reads its outputs; it re-implements nothing.
Pipeline bootstraps use their own RNG, so they are emitted as *_own keys (compared within Monte Carlo tolerance)."""
import json, numpy as np, pandas as pd, warnings, statsmodels.formula.api as smf
from statsmodels.stats.diagnostic import acorr_ljungbox
warnings.filterwarnings('ignore')
import os   # run from the project root
exec(open('scripts/09_final_analysis.py').read().split('# Q5 regional')[0])
exec(open('scripts/15_nowcast_monthly_refit.py').read().split('out = []')[0].split("exec(open('scripts/09_final_analysis.py')")[0])
src15 = open('scripts/15_nowcast_monthly_refit.py').read(); exec(src15[src15.index('def nowcast_monthly'):src15.index('out = []')])
R = {}
x, c = S['dengue_3day'], cases
lab = lambda k: {-3: 'm3', -2: 'm2', -1: 'm1', 0: '0', 1: 'p1', 2: 'p2', 3: 'p3'}[k]

# Q1
m = (c.index.year != 2022) & (c.index.year >= 2020)
runs = {'primary': (x['2020':'2024'], c['2020':'2024']), 'A': (x['2019-06':'2024-12'], c['2019-06':'2024-12']), 'B': (x, c),
        'C': (x['2022':'2024'], c['2022':'2024']), 'D': (x[m], c[m]), 'topic': (S['topic']['2020':'2024'], c['2020':'2024']),
        'day1': (S['dengue_day1']['2020':'2024'], c['2020':'2024']), 'fever': (S['fever']['2020':'2024'], c['2020':'2024'])}
for k, (xs, cs) in runs.items():
    q = q1(xs, cs)['primary_adjusted_log']
    R[f'q1_{k}_r'], R[f'q1_{k}_neff'], R[f'q1_{k}_flo'], R[f'q1_{k}_fhi'] = q['r'], q['n_eff'], q['fisher_lo'], q['fisher_hi']
    if k in ('primary', 'fever'): R[f'q1_{k}_blo_own'], R[f'q1_{k}_bhi_own'] = q['boot_lo'], q['boot_hi']

# Q2
xs, cs = runs['primary']
t, g = q2(xs, cs)
X = adj(lx(xs)).values; sel = ar_select_order(X, maxlag=3, ic='aic', trend='c'); fit = sel.model.fit()
R['q2_ar_order'] = int(max(sel.ar_lags)); R['q2_ar_const'] = fit.params[0]
for i, v in enumerate(fit.params[1:], 1): R[f'q2_ar_phi{i}'] = v
for _, row in t.iterrows():
    k = lab(int(row.lag)); R[f'q2_ccf_r_{k}'], R[f'q2_ccf_n_{k}'], R[f'q2_ccf_pholm_{k}'] = row.r, int(row.n), row.p_holm
d = pd.DataFrame({'y': ly(cs), 'y1': ly(cs).shift(1), 'x': lx(xs), 'x1': lx(xs).shift(1), 'm': cs.index.month}).dropna()
w1 = smf.ols('y ~ C(m) + y1 + x1', d).fit(cov_type='HAC', cov_kwds={'maxlags': 3})
R['q2_fwd_coef'], R['q2_fwd_se'], R['q2_fwd_p'] = w1.params['x1'], w1.bse['x1'], w1.pvalues['x1']
R['q2_rev58_coef'], R['q2_rev58_p'] = g['cases_to_search_coef'], g['cases_to_search_p']          # as run in scripts/09 (n = 58)
c19 = pd.read_excel('out/19_verification_corrections.xlsx', 'direction_tests')                     # corrected n = 59 (scripts/19)
rv = c19[(c19.model.str.startswith('reverse')) & (~c19.small_sample_correction)].iloc[0]
w2 = smf.ols('x ~ C(m) + x1 + y1', d).fit(cov_type='HAC', cov_kwds={'maxlags': 3})
R['q2_rev_coef'], R['q2_rev_se'], R['q2_rev_p'] = rv.coef, w2.bse['y1'], rv.p_normal

# NB
nb = nb_model(*runs['primary'])
dd = pd.DataFrame({'c': cs, 'y1': ly(cs).shift(1), 'x': lx(xs)}); dd['s1'], dd['c1'] = np.sin(2*np.pi*dd.index.month/12), np.cos(2*np.pi*dd.index.month/12)
dd['s2'], dd['c2'] = np.sin(4*np.pi*dd.index.month/12), np.cos(4*np.pi*dd.index.month/12); dd = dd.dropna()
mm = smf.negativebinomial('c ~ s1 + c1 + s2 + c2 + y1 + x', dd).fit(disp=0, maxiter=500)
mu = mm.predict(dd); pr = (dd.c - mu) / np.sqrt(mu + mm.params['alpha'] * mu ** 2)
R['nb_coef_x'] = mm.params['x']; R['nb_rr'] = nb['IRR_per_doubling']; R['nb_alpha'] = nb['alpha']
R['nb_rr_lo'], R['nb_rr_hi'] = nb['lo'], nb['hi']; R['nb_lb_q12'] = float(acorr_ljungbox(pr, lags=[12], return_df=True).lb_stat.iloc[0]); R['nb_lb_p12'] = nb['ljungbox12_p']
R['nb_cook_months'] = nb['cook_flagged']

# Q3
def q3keys(s, P, met, per, dd_):
    for _, r_ in met.iterrows():
        M = r_.model; R[f'q3_{s}_mae_{M}'], R[f'q3_{s}_rmse_{M}'], R[f'q3_{s}_maelog_{M}'] = r_.MAE, r_.RMSE, r_.MAE_log
        if s != 'loso':
            for th in (1000, 500): R[f'q3_{s}_sens_{M}_{th}'], R[f'q3_{s}_spec_{M}_{th}'] = r_[f'sens_{th}'], r_[f'spec_{th}']
    if s != 'loso':
        for yv in (2022, 2023, 2024):
            for M in ('M1', 'M2'): R[f'q3_{s}_season_mae_{M}_{yv}'] = per.loc[yv, M]
        R[f'q3_{s}_seasons_m2_better'] = int(dd_['seasons_M2_better'] if 'seasons_M2_better' in dd_ else dd_['seasons_M2_better'])
        R[f'q3_{s}_diff_mean'] = dd_['mean']; R[f'q3_{s}_diff_lo_own'], R[f'q3_{s}_diff_hi_own'] = dd_['lo'], dd_['hi']
        R[f'q3_{s}_rule_met'] = int(dd_['seasons_M2_better'] >= 2 and dd_['hi'] < 0)
q3keys('ps', *nowcast(x, c, 'rolling', [2022, 2023, 2024]))
P, met, per, dm = nowcast_monthly(x, c); per.index = per.index.astype(int); q3keys('mo', P, met, per, dm)
q3keys('loso', *nowcast(x, c, 'loso', [2020, 2021, 2022, 2023, 2024], train_start='2020-02-01'))

# Q4
for yv in range(2019, 2025):
    cy, xy = c[str(yv)], x[str(yv)]; R[f'q4_{yv}_peak_ratio'] = cy.max() / xy.max(); R[f'q4_{yv}_total_ratio'] = cy.sum() / xy.sum()
# Q5
for rg in ['Central', 'Eastern', 'Western']:
    s_ = regS[rg]['2020':'2023']; cc = regC[rg]['2020':]
    q = corr_block(adj(lx(s_)).values, adj(ly(cc)).values)
    R[f'q5_{rg}_r'], R[f'q5_{rg}_neff'], R[f'q5_{rg}_flo'], R[f'q5_{rg}_fhi'] = q['r'], q['n_eff'], q['fisher_lo'], q['fisher_hi']
    R[f'q5_{rg}_zeros2020_2023'] = int((s_ == 0).sum())
    R[f'q5_{rg}_peakratio_cases'] = cc['2023'].max() / cc['2022'].max(); R[f'q5_{rg}_peakratio_search'] = s_['2023'].max() / s_['2022'].max()
pm = pd.read_csv('out/20_peak_margin_full.csv')
R['peak_n_changed'] = int(pm.changes.sum()); R['peak_ties_05'], R['peak_ties_10'], R['peak_ties_20'] = [int((pm[k] == 'tie').sum()) for k in ('m05', 'm10', 'm20')]
# Q6
q6 = pd.read_csv('out/06_2024_reporting_lag.csv')
for _, r_ in q6.iterrows():
    dt = r_.as_of.replace('-', ''); R[f'q6_final_{dt}'] = float(r_.final_count_by_then); R[f'q6_pct_{dt}'] = 100 * r_.reported_then / r_.final_count_by_then
R = {k: (float(v) if isinstance(v, (np.floating, float)) else int(v) if isinstance(v, (np.integer, int, bool, np.bool_)) else v) for k, v in R.items()}
json.dump(R, open('verification/results_impl1_pipeline.json', 'w'), indent=1); print(len(R), 'keys')
