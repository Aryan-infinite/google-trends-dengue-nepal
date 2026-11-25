"""Figure 3. Prewhitened cross-correlation between monthly "dengue" search interest and reported cases,
Nepal, 2020-2024, lags -3 to +3 months (positive lag = search interest precedes cases).

Computation (Methods 2.4; verification/VERIFICATION_SPEC.md "Prewhitening and cross-correlations"):
  log(search + 0.5) and log(cases + 1), January 2020 - December 2024, each minus its calendar-month mean;
  AR order 0-3 chosen by AIC (statsmodels ar_select_order, constant, common sample after 3 months);
  coefficients refitted by OLS on months p+1..60; the SEARCH filter is applied to both series;
  r_k = Pearson(e_search[t], e_cases[t+k]); two-sided t-test with n_k - 2 df; Holm across the 7 lags.
Drawn: bars = r_k; dashed segments = +/- the |r| giving two-sided P = 0.05 (unadjusted) for that lag's n_k;
a label over any lag significant after Holm adjustment; n_k under each lag.
Outputs: manuscript/figures/Fig3.tif, Fig3_preview.png; out/fig3_plotted_values.csv (read back from the drawing)
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd, scipy.stats as st
from statsmodels.tsa.ar_model import ar_select_order
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import fig_common as fc

os.makedirs('manuscript/figures', exist_ok=True)
IDX = pd.date_range('2019-01-01', periods=72, freq='MS')
paper = pd.read_csv('out/paper_72month_table.csv', index_col=0, parse_dates=True)
cases = pd.Series(paper.cases.astype(float).values, IDX)
F = ['data_raw/gt_dengue_term.csv', 'data_raw/gt_dengue_term_day2_2026-11-22.csv', 'data_raw/gt_dengue_term_day3_2026-11-23.csv']
read = lambda f: pd.read_csv(f, skiprows=2).iloc[:, 1].replace('<1', 0.5).astype(float).values
search = pd.Series(np.mean([read(f) for f in F], axis=0), IDX)

c, s = cases['2020-01':'2024-12'], search['2020-01':'2024-12']
adj = lambda v: v - v.groupby(v.index.month).transform('mean')
X, Y = adj(np.log(s + 0.5)), adj(np.log(c + 1))
sel = ar_select_order(X.values, maxlag=3, ic='aic', trend='c'); lags = sel.ar_lags or []
p = max(lags) if lags else 0
if p:
    fit = sel.model.fit(); const = fit.params[0]; phi = dict(zip(lags, fit.params[1:]))
    filt = lambda v: (v - const - sum(phi[k] * v.shift(k) for k in lags)).iloc[p:]
    ex, ey = filt(X), filt(Y)
else:
    ex, ey = X, Y
rows = []
for k in range(-3, 4):
    a = ex.shift(k); m = a.notna() & ey.notna()
    r, pv = st.pearsonr(a[m], ey[m]); n = int(m.sum())
    t = st.t.ppf(0.975, n - 2); rows.append(dict(lag=k, r=r, n=n, p=pv, r_crit=t / np.sqrt(n - 2 + t * t)))
T = pd.DataFrame(rows); T['p_holm'] = multipletests(T.p, method='holm')[1]

font = fc.setup()
fig, ax = plt.subplots(figsize=(120 * fc.MM, 72 * fc.MM))
bars = ax.bar(T.lag, T.r, width=0.6, color=fc.BLUE, linewidth=0, zorder=3)
seg = []
for _, row in T.iterrows():
    for sgn in (1, -1):
        ln, = ax.plot([row.lag - 0.4, row.lag + 0.4], [sgn * row.r_crit] * 2, color=fc.INK2, lw=0.8, ls=(0, (3, 2)), zorder=4)
        seg.append((int(row.lag), sgn, ln))
ax.axhline(0, color=fc.INK2, lw=0.6, zorder=2)
for _, row in T[T.p_holm < 0.05].iterrows():
    txt = 'Holm-adjusted\nP < 0.001' if row.p_holm < 0.001 else f'Holm-adjusted\nP = {row.p_holm:.3f}'
    ax.text(row.lag, row.r + 0.03, txt, ha='center', va='bottom', fontsize=fc.FONT_PT - 1, color=fc.INK)
ax.set_xticks(T.lag)
ax.set_xticklabels([f'{"−" if k < 0 else ("+" if k > 0 else "")}{abs(k)}\nn = {n}' for k, n in zip(T.lag, T.n)])
ax.set_xlim(-3.6, 3.6); ax.set_ylim(-0.4, 0.7); ax.set_yticks(np.arange(-0.4, 0.61, 0.2))
ax.set_yticklabels([f'{v:.1f}'.replace('-', '−') for v in np.arange(-0.4, 0.61, 0.2)])
ax.grid(axis='y', color=fc.GRID, lw=0.4); ax.set_axisbelow(True)
ax.tick_params(axis='x', length=0)
ax.set_xlabel('Lag (months); positive lag = search interest precedes cases')
ax.set_ylabel('Cross-correlation (r)')
ax.legend(handles=[Line2D([], [], color=fc.INK2, lw=0.8, ls=(0, (3, 2)), label='Two-sided P = 0.05 (unadjusted)')],
          loc='lower right', bbox_to_anchor=(1.0, 1.0), frameon=False, borderaxespad=0.2)
fig.subplots_adjust(left=0.13, right=0.98, top=0.9, bottom=0.2)
preview = fc.save_tiff_cmyk(fig, 'manuscript/figures/Fig3.tif')

# Read back what was drawn.
out = pd.DataFrame({'lag': [int(round(b.get_x() + b.get_width() / 2)) for b in bars], 'bar_height': [b.get_height() for b in bars]})
crit = pd.DataFrame([(k, s_, ln.get_ydata()[0]) for k, s_, ln in seg], columns=['lag', 'sign', 'y'])
out['crit_pos'] = out.lag.map(crit[crit.sign == 1].set_index('lag').y); out['crit_neg'] = out.lag.map(crit[crit.sign == -1].set_index('lag').y)
out['tick_label'] = [t.get_text().replace('\n', ' | ') for t in ax.get_xticklabels()]
labels = {round(t.get_position()[0]): t.get_text().replace('\n', ' ') for t in ax.texts}
out['annotation'] = out.lag.map(labels).fillna('')
out['ar_order'] = p
out.to_csv('out/fig3_plotted_values.csv', index=False, float_format='%.17g')
print('font:', font, '| AR order', p, '| AR coef', {'const': round(const, 6), **{k: round(v, 6) for k, v in phi.items()}} if p else {})
print(T.round(6).to_string(index=False))
