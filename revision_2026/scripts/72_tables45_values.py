"""Assemble the display values of Tables 4 and 5 (Python). Sources, in order of preference:
  - the analysis pipeline's own outputs (out/final/final_results.xlsx; out/15_nowcast_monthly_refit.xlsx);
  - for values the pipeline did not save (n_eff and lagged-search coefficients of the sensitivity runs): the R and NumPy
    implementations (scripts/68, 69), which must agree with each other to 1e-9 first;
  - 3-month block bootstrap limits (post hoc): the shared-resample values of verification/comparison_table.csv,
    on which three implementations agree.
Every pipeline value used is also compared here with R and NumPy (scripts/68-71); the script stops on any disagreement.
Output: out/tables45_values.json (display strings + raw numbers), read by scripts/73 (docx) and checked by scripts/74 (R)."""
import json, csv, math
import pandas as pd
F = 'out/final/final_results.xlsx'
grid = pd.read_excel(F, 'Robustness_grid'); summ = dict(pd.read_excel(F, 'Summary').values)
R4 = pd.read_csv('out/table4_runs_R.csv'); N4 = pd.read_csv('out/table4_runs_np.csv')
R5 = {r['key']: float(r['value']) for r in csv.DictReader(open('out/table5_R.csv'))}
N5 = {r['key']: float(r['value']) for r in csv.DictReader(open('out/table5_np.csv'))}
CT = {r['key']: r for r in csv.DictReader(open('verification/comparison_table.csv'))}
def same(a, b, tol=1e-9): return abs(a - b) <= tol * max(1.0, abs(a))
for c in R4.columns[1:]:
    assert all(same(float(a), float(b)) for a, b in zip(R4[c], N4[c])), c
MINUS = '−'
def num(v, d): s = f'{v:,.{d}f}'; return s.replace('-', MINUS)
def ci(lo, hi, d): return f'{num(lo, d)} to {num(hi, d)}'
def pv(p): return '<0.001' if p < 0.001 else (f'{p:.3f}' if p < 0.01 else f'{p:.2f}')

# ---------------- Table 4 ----------------
RUNS = [('primary', 'primary (2020–24, 3-day mean)', 'Primary analysis, January 2020–December 2024'),
        ('A', 'A: 2019 Jun–2024', '(i) Including June–December 2019'), ('B', 'B: full 72 months', '(ii) All 72 months, as in the preprint'),
        ('C', 'C: excl. 2020–21 (COVID), 2022–24', '(iii) 2022–2024 only'), ('D', 'D: excl. 2022', '(iv) Excluding 2022'),
        ('topic', 'Topic instead of term', '(v) Preprint Topic series'), ('day1', 'Day-1 download only', '(vi) First retrieval only'),
        ('lt0', '"<1" coded 0', '(vii) Values "<1" coded 0'), ('fever', 'Control: "Fever"', 'Control series "Fever"')]
t4 = []
for key, glab, label in RUNS:
    g = grid[grid.run == glab].iloc[0]; r = R4[R4.run == key].iloc[0]
    for gc, rc in [('n', 'n'), ('adj_r', 'r'), ('lag_minus1', 'lag_m1'), ('lag0', 'lag_0'), ('lag_plus1', 'lag_p1'), ('lag_plus1_holm', 'holm_p1'), ('wald_x1_p', 'fwd_p')]:
        assert same(float(g[gc]), float(r[rc]), 1e-8), (key, gc, g[gc], r[rc])
    assert same(g.ci_lo, r.ci_lo, 1e-8)                                     # lower limit: Fisher in every run
    hi_src = 'fisher'
    if not same(g.ci_hi, r.ci_hi, 1e-8):                                   # union with the bootstrap limit (Fever only)
        assert key == 'fever' and g.ci_hi > r.ci_hi; hi_src = 'bootstrap'
    met = bool(r.holm_p1 < 0.05 and r.lag_p1 > r.lag_m1 and r.fwd_p < 0.05)
    t4.append(dict(key=key, label=label, n=int(g.n), r=float(g.adj_r), lo=float(g.ci_lo), hi=float(g.ci_hi), hi_src=hi_src, n_eff=float(r.n_eff),
                   lag_m1=float(g.lag_minus1), lag_0=float(g.lag0), lag_p1=float(g.lag_plus1), holm_p1=float(g.lag_plus1_holm), holm_0=float(r.holm_0), holm_m1=float(r.holm_m1), fwd_text=f'{num(r.fwd_coef, 2)}; P = {pv(g.wald_x1_p)}',
                   fwd_coef=float(r.fwd_coef), fwd_p=float(g.wald_x1_p), rule_met=met,
                   cells=[label, str(int(g.n)), f'{num(g.adj_r, 2)} ({ci(g.ci_lo, g.ci_hi, 2)})', num(r.n_eff, 1),
                          f'{num(g.lag_minus1, 2)} ({pv(r.holm_m1)})', f'{num(g.lag0, 2)} ({pv(r.holm_0)})', f'{num(g.lag_plus1, 2)} ({pv(g.lag_plus1_holm)})',
                          f'{num(r.fwd_coef, 2)} ({pv(g.wald_x1_p)})', 'Yes' if met else 'No']))
# ---------------- Table 5 ----------------
ps = pd.read_excel(F, 'Q3_rolling_metrics').set_index('model'); lo_ = pd.read_excel(F, 'Q3_LOSO_metrics').set_index('model')
rb = pd.read_excel(F, 'Q3_robustness').set_index('series')
mo = pd.read_excel('out/15_nowcast_monthly_refit.xlsx', 'metrics').set_index('model'); mrb = pd.read_excel('out/15_nowcast_monthly_refit.xlsx', 'robustness').set_index('series')
SERIES = [('dengue', 'dengue 3-day', '"dengue" (primary series)'), ('topic', 'Topic', 'Preprint Topic series'), ('day1', 'day-1', 'First retrieval only'), ('fever', 'Fever control', 'Control series "Fever"')]
t5 = []
def row(group, sched, sk, sl, label, maes, diff, lo, hi, better, nseas, ntest, rule):
    for M, v in zip(('M0', 'M1', 'M2', 'M3'), maes): assert same(v, R5[f'{sk}_{sched}_mae_{M}'], 1e-8) and same(v, N5[f'{sk}_{sched}_mae_{M}'], 1e-8), (sk, sched, M)
    assert same(diff, R5[f'{sk}_{sched}_diff'], 1e-8) and int(better) == int(R5[f'{sk}_{sched}_better']) == int(N5[f'{sk}_{sched}_better'])
    # Monte Carlo check of the pipeline's own-RNG limits against R's own-RNG limits: within 5% of the interval width
    # (a relative-to-limit rule fails for limits near zero, e.g. -65 vs -52 when the interval is 1,850 wide).
    for v, rr in ((lo, R5[f'{sk}_{sched}_bootlo_Rrng']), (hi, R5[f'{sk}_{sched}_boothi_Rrng'])):
        assert abs(v - rr) <= 0.05 * (hi - lo), (sk, sched, v, rr)
    if sk == 'dengue' and sched in ('ps', 'mo'):                                    # and against the shared-resample values (spec rule: 15% of |limit|)
        for v, k in ((lo, 'lo'), (hi, 'hi')): assert abs(v - float(CT[f'q3_{sched}_diff_{k}']['impl2_R'])) <= 0.15 * abs(v), (sched, k)
    t5.append(dict(group=group, sched=sched, series=sk, label=label, maes=list(map(float, maes)), diff=float(diff), lo=float(lo), hi=float(hi),
                   better=int(better), nseasons=nseas, ntest=ntest, rule=rule,
                   cells=[label, str(ntest)] + [num(v, 0) for v in maes] + [f'{num(diff, 0)} ({ci(lo, hi, 0)})', f'{int(better)} of {nseas}', rule]))
for sk, rl, label in SERIES:
    m = rb.loc[rl]; maes = [R5[f'{sk}_ps_mae_M0']] + [m.MAE_M1, m.MAE_M2, m.MAE_M3]
    if sk == 'dengue': assert all(same(a, b, 1e-8) for a, b in zip(maes, ps.MAE.values))
    met = bool(m.seasons_M2_better >= 2 and m.hi < 0)
    row('Refitted once per season (primary schedule)', 'ps', sk, rl, label, maes, m.pairdiff, m.lo, m.hi, m.seasons_M2_better, 3, 36, 'Met' if met else 'Not met')
for sk, rl, label in SERIES:
    m = mrb.loc[rl]; maes = [m.MAE_M0, m.MAE_M1, m.MAE_M2, m.MAE_M3]
    if sk == 'dengue': assert all(same(a, b, 1e-8) for a, b in zip(maes, mo.MAE.values))
    met = bool(m.seasons_M2_better >= 2 and m.hi < 0); assert met == bool(m.meets_rule)
    row('Refitted every month (post hoc)', 'mo', sk, rl, label, maes, m.pairdiff, m.lo, m.hi, m.seasons_M2_better, 3, 36, 'Met' if met else 'Not met')
row('Leave one season out, 2020–2024 (secondary)', 'loso', 'dengue', 'dengue', '"dengue" (primary series)', list(lo_.MAE.values),
    summ['Q3_loso_pairdiff_M2_minus_M1.mean'], summ['Q3_loso_pairdiff_M2_minus_M1.lo'], summ['Q3_loso_pairdiff_M2_minus_M1.hi'],
    summ['Q3_loso_pairdiff_M2_minus_M1.seasons_M2_better'], int(summ['Q3_loso_pairdiff_M2_minus_M1.n_seasons']), int(R5['dengue_loso_ntest']), 'No rule')
blocks = {s: (float(CT[f'q3_{s}_diff_b3lo']['impl2_R']), float(CT[f'q3_{s}_diff_b3hi']['impl2_R'])) for s in ('ps', 'mo')}
for s in ('ps', 'mo'): assert CT[f'q3_{s}_diff_b3lo']['ok'] == 'True' and CT[f'q3_{s}_diff_b3hi']['ok'] == 'True'
json.dump(dict(t4=t4, t5=t5, blocks={k: dict(lo=v[0], hi=v[1], text=ci(v[0], v[1], 0)) for k, v in blocks.items()},
               fever_boot_hi_shared=float(CT['q1_fever_bhi']['impl2_R'])), open('out/tables45_values.json', 'w'), indent=1, ensure_ascii=False)
for r in t4: print(' | '.join(r['cells']))
for r in t5: print(r['group'][:12], '|', ' | '.join(r['cells']))
print('blocks', {k: ci(*v, 0) for k, v in blocks.items()})
