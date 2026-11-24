"""Rolling-origin nowcast with the models refitted before every test month (POST HOC sensitivity analysis; written
2026-11-23, after the primary results were known). ANALYSIS_PLAN.md Q3 says 'training = all earlier months', which
admits refitting before every month; the primary analysis (scripts/09, implemented before its results were seen)
refits once at the start of each test season. scripts/09 is left unchanged.
Uses scripts/09 data loading and transforms so the only difference is the refit schedule."""
import numpy as np, pandas as pd, statsmodels.formula.api as smf, warnings
warnings.filterwarnings('ignore')
exec(open('scripts/09_final_analysis.py').read().split('# ---------- windows ----------')[0])

def nowcast_monthly(xs, cs, tests=(2022, 2023, 2024), train_start='2019-07-01'):
    d = pd.DataFrame({'y': ly(cs), 'y1': ly(cs).shift(1), 'x': lx(xs), 'x1': lx(xs).shift(1), 'm': cs.index.month, 'c': cs})
    d = d[d.index >= train_start].dropna()
    forms = {'M1': 'y ~ C(m) + y1', 'M2': 'y ~ C(m) + y1 + x', 'M3': 'y ~ C(m) + y1 + x1'}
    rows = []
    for t in d.index[d.index.year.isin(tests)]:
        tr, te = d[d.index < t], d.loc[[t]]
        r = dict(t=t, year=t.year, obs=te.c.item(), M0=np.expm1(te.y1.item()))
        for k, f in forms.items(): r[k] = np.expm1(smf.ols(f, tr).fit().predict(te).item())
        rows.append(r)
    P = pd.DataFrame(rows).set_index('t'); P[['M0', 'M1', 'M2', 'M3']] = P[['M0', 'M1', 'M2', 'M3']].clip(lower=0)
    met = []
    for k in ['M0', 'M1', 'M2', 'M3']:
        e = P[k] - P.obs; row = dict(model=k, MAE=e.abs().mean(), RMSE=np.sqrt((e ** 2).mean()), MAE_log=(np.log1p(P[k]) - np.log1p(P.obs)).abs().mean())
        for thr in [1000, 500]:
            o, pr = P.obs > thr, P[k] > thr
            row[f'sens_{thr}'] = (o & pr).sum() / o.sum(); row[f'spec_{thr}'] = (~o & ~pr).sum() / (~o).sum()
        met.append(row)
    per = P.groupby('year').apply(lambda g: pd.Series({k: (g[k] - g.obs).abs().mean() for k in ['M0', 'M1', 'M2', 'M3']}))
    diff = (P.M2 - P.obs).abs().values - (P.M1 - P.obs).abs().values
    rng = np.random.default_rng(SEED); bs = [rng.choice(diff, len(diff)).mean() for _ in range(B)]
    dd = dict(mean=diff.mean(), lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5), seasons_M2_better=int((per.M2 < per.M1).sum()), n_seasons=len(per))
    return P, pd.DataFrame(met), per, dd

out = []
with pd.ExcelWriter('out/15_nowcast_monthly_refit.xlsx') as xw:
    for lab, s in [('dengue 3-day', S['dengue_3day']), ('Topic', S['topic']), ('day-1', S['dengue_day1']), ('Fever control', S['fever'])]:
        P, met, per, dd = nowcast_monthly(s, cases)
        out.append(dict(series=lab, **{f'MAE_{m}': met.loc[met.model == m, 'MAE'].item() for m in ['M0', 'M1', 'M2', 'M3']},
                        pairdiff=dd['mean'], lo=dd['lo'], hi=dd['hi'], seasons_M2_better=dd['seasons_M2_better'],
                        meets_rule=bool(dd['seasons_M2_better'] >= 2 and dd['hi'] < 0)))
        if lab == 'dengue 3-day':
            met.to_excel(xw, 'metrics', index=False); per.to_excel(xw, 'by_season'); P.to_excel(xw, 'predictions')
    pd.DataFrame(out).to_excel(xw, 'robustness', index=False)
print(pd.DataFrame(out).round(1).to_string())
print(met.round(3).to_string()); print(per.round(1).to_string())
