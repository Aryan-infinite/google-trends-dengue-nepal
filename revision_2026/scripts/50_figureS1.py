"""Supplementary Figure S1. Autocorrelation (ACF) and partial autocorrelation (PACF), lags 1-20, of the four monthly
series analysed (January 2020 - December 2024, n = 60): log(cases + 1) and log(search + 0.5), each as observed and
after subtracting calendar-month means. Definitions as in scripts/21_reviewer_checks.py (check A):
  ACF_k = sum_t (v_t - mean)(v_{t+k} - mean) / sum_t (v_t - mean)^2          (denominator n, not n - k)
  PACF_k = coefficient on v_{t-k} in an OLS regression of v_t on a constant and v_{t-1}..v_{t-k}
  band = +/- 1.96 / sqrt(60), the approximate 95% limits for white noise.
Ljung-Box Q at lags 6 and 12 is written to out/figS1_values.csv (reported in the caption, not drawn).
Outputs: manuscript/figures/FigS1.tif, FigS1_preview.png; out/figS1_values.csv (computed); out/figS1_plotted_values.csv (read back)
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd, scipy.stats as st
import matplotlib.pyplot as plt
import fig_common as fc

os.makedirs('manuscript/figures', exist_ok=True)
paper = pd.read_csv('out/paper_72month_table.csv', index_col=0, parse_dates=True)
cases = paper.cases.astype(float).values[12:]
F = ['data_raw/gt_dengue_term.csv', 'data_raw/gt_dengue_term_day2_2026-11-22.csv', 'data_raw/gt_dengue_term_day3_2026-11-23.csv']
read = lambda f: pd.read_csv(f, skiprows=2).iloc[:, 1].replace('<1', 0.5).astype(float).values
search = np.mean([read(f) for f in F], axis=0)[12:]
mon = np.tile(np.arange(1, 13), 5)
adj = lambda v: v - np.array([v[mon == m].mean() for m in mon])
ly, lx = np.log(cases + 1), np.log(search + 0.5)
SER = [('ly_raw', 'Log cases', ly), ('ly_adj', 'Log cases, month-adjusted', adj(ly)),
       ('lx_raw', 'Log search interest', lx), ('lx_adj', 'Log search interest, month-adjusted', adj(lx))]
K, n = 20, 60
def acf(v):
    d = v - v.mean(); return np.array([np.sum(d[:-k] * d[k:]) for k in range(1, K + 1)]) / np.sum(d * d)
def pacf(v):
    out = []
    for k in range(1, K + 1):
        X = np.column_stack([np.ones(n - k)] + [v[k - j:n - j] for j in range(1, k + 1)])
        out.append(np.linalg.lstsq(X, v[k:], rcond=None)[0][k])
    return np.array(out)
band = 1.96 / np.sqrt(n)
rows = []
for key, _, v in SER:
    a, p = acf(v), pacf(v)
    for L in (6, 12):
        q = n * (n + 2) * np.sum(a[:L] ** 2 / (n - np.arange(1, L + 1))); rows.append(dict(series=key, stat=f'lbq{L}', value=q))
        rows.append(dict(series=key, stat=f'lbp{L}', value=st.chi2.sf(q, L)))
    rows += [dict(series=key, stat=f'acf{k + 1}', value=a[k]) for k in range(K)] + [dict(series=key, stat=f'pacf{k + 1}', value=p[k]) for k in range(K)]
V = pd.DataFrame(rows); V.to_csv('out/figS1_values.csv', index=False, float_format='%.17g')
get = lambda key, s: V[(V.series == key) & (V.stat == s)].value.item()

font = fc.setup()
fig, axes = plt.subplots(4, 2, figsize=(fc.WIDTH_MM * fc.MM, 190 * fc.MM), sharex=True, sharey=True,
                         gridspec_kw=dict(hspace=0.42, wspace=0.12))
lags = np.arange(1, K + 1); drawn = []
for r, (key, title, v) in enumerate(SER):
    for c, kind in enumerate(('acf', 'pacf')):
        ax = axes[r, c]; vals = np.array([get(key, f'{kind}{k}') for k in lags])
        bars = ax.bar(lags, vals, width=0.6, color=fc.BLUE, linewidth=0, zorder=3)
        lo, = ax.plot([0.4, K + 0.6], [-band, -band], color=fc.INK2, lw=0.8, ls=(0, (3, 2)), zorder=4)
        hi, = ax.plot([0.4, K + 0.6], [band, band], color=fc.INK2, lw=0.8, ls=(0, (3, 2)), zorder=4)
        ax.axhline(0, color=fc.INK2, lw=0.6, zorder=2)
        ax.set_ylim(-1, 1); ax.set_yticks([-1, -0.5, 0, 0.5, 1]); ax.set_xlim(0.4, K + 0.6)
        ax.set_yticklabels(['−1.0', '−0.5', '0.0', '0.5', '1.0'])
        ax.grid(axis='y', color=fc.GRID, lw=0.4); ax.set_axisbelow(True)
        ax.set_title(f'({"abcdefgh"[2 * r + c]}) {title}', fontsize=fc.FONT_PT, loc='left', pad=4)
        for b, k in zip(bars, lags): drawn.append(dict(series=key, kind=kind, lag=k, height=b.get_height()))
        drawn.append(dict(series=key, kind=kind, lag=0, height=hi.get_ydata()[0])); drawn.append(dict(series=key, kind=kind, lag=-1, height=lo.get_ydata()[0]))
for ax in axes[-1]: ax.set_xticks([1, 5, 10, 15, 20]); ax.set_xlabel('Lag (months)')
fig.subplots_adjust(left=0.08, right=0.99, top=0.94, bottom=0.06)
for c, head in enumerate(('Autocorrelation', 'Partial autocorrelation')):
    bb = axes[0, c].get_position(); fig.text((bb.x0 + bb.x1) / 2, 0.985, head, ha='center', va='top', fontweight='bold')
preview = fc.save_tiff_cmyk(fig, 'manuscript/figures/FigS1.tif')
D = pd.DataFrame(drawn)
D['title'] = D.apply(lambda r: axes[[s[0] for s in SER].index(r.series), ('acf', 'pacf').index(r.kind)].get_title(loc='left'), axis=1)
D.to_csv('out/figS1_plotted_values.csv', index=False, float_format='%.17g')
pd.DataFrame({'header': [t.get_text() for t in fig.texts]}).to_csv('out/figS1_headers.csv', index=False)
print('font:', font, '| band', round(band, 6), '| bars drawn', int((D.lag > 0).sum()))
print({key: (round(get(key, 'acf1'), 3), round(get(key, 'lbq12'), 1), f"{get(key, 'lbp12'):.1e}") for key, _, _ in SER})
