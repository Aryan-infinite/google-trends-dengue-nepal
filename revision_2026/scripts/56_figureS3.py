"""Supplementary Figure S3. Regional search interest and regional reported cases, 2020-2023 (Question 5; descriptive).
(a) Seasonally adjusted correlation (log scales, calendar-month means removed within 2020-2023, n = 48) with 95% CI,
    for the three regions with fewer than 30 zero months (Central, Eastern, Western).
    CI = union of the Fisher interval on the effective sample size and the 12-month block-bootstrap interval
    (scripts/09 corr_block); for these regions the union equals the Fisher interval, which is asserted below.
(b) Ratio of the 2023 seasonal peak to the 2022 seasonal peak, for cases and for search interest, per region.
Inputs: out/edcd_devregion_monthly_2019_2023.csv (scripts/04, from the EDCD district workbook, which is not public),
        data_raw/gt_dengue_region_{central,eastern,western}.csv.
Outputs: manuscript/figures/FigS3.tif, FigS3_preview.png; out/figS3_values.csv; out/figS3_plotted_values.csv
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import fig_common as fc

os.makedirs('manuscript/figures', exist_ok=True)
REG = ['Central', 'Eastern', 'Western']
rc = pd.read_csv('out/edcd_devregion_monthly_2019_2023.csv', index_col=0, parse_dates=True)['2020':'2023']
def gt(f):
    d = pd.read_csv(f, skiprows=2); d = d[d.iloc[:, 0].str[:4].isin(['2020', '2021', '2022', '2023'])]
    return d.iloc[:, 1].astype(str).replace('<1', '0.5').astype(float).values
mon = np.tile(np.arange(1, 13), 4)
adj = lambda v: v - np.array([v[mon == m].mean() for m in mon])
def lagr(v, k): return np.corrcoef(v[:-k], v[k:])[0, 1]
def neff(a, b): return max(len(a) / (1 + 2 * sum(lagr(a, k) * lagr(b, k) for k in range(1, 13))), 4.0)
def block_ci(a, b, L=12, seed=20261123, B=4000):                      # as scripts/09 corr_block
    rng = np.random.default_rng(seed); n = len(a); rs = []
    for _ in range(B):
        idx = np.concatenate([np.arange(s, s + L) for s in rng.choice(np.arange(n - L + 1), int(np.ceil(n / L)))])[:n]
        if np.std(a[idx]) > 0 and np.std(b[idx]) > 0: rs.append(np.corrcoef(a[idx], b[idx])[0, 1])
    return np.percentile(rs, [2.5, 97.5])
rows = []
for r_ in REG:
    s = gt(f'data_raw/gt_dengue_region_{r_.lower()}.csv'); c = rc[r_].values.astype(float)
    assert len(s) == 48 and len(c) == 48
    a, b = adj(np.log(s + 0.5)), adj(np.log(c + 1)); r = np.corrcoef(a, b)[0, 1]; ne = neff(a, b)
    z, se = np.arctanh(r), 1 / np.sqrt(ne - 3); flo, fhi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
    blo, bhi = block_ci(a, b)
    assert blo >= flo and bhi <= fhi, f'{r_}: bootstrap interval not inside Fisher interval; union differs from Fisher'
    rows.append(dict(region=r_, r=r, n_eff=ne, ci_lo=min(flo, blo), ci_hi=max(fhi, bhi), boot_lo=blo, boot_hi=bhi, zeros=int((s == 0).sum()),
                     ratio_cases=c[36:].max() / c[24:36].max(), ratio_search=s[36:].max() / s[24:36].max()))
V = pd.DataFrame(rows); V.to_csv('out/figS3_values.csv', index=False, float_format='%.17g')

font = fc.setup()
fig, (a1, a2) = plt.subplots(1, 2, figsize=(fc.WIDTH_MM * fc.MM, 62 * fc.MM), gridspec_kw=dict(wspace=0.55, width_ratios=[1, 1.1]))
ys = np.arange(len(REG))[::-1]
err = a1.errorbar(V.r, ys, xerr=[V.r - V.ci_lo, V.ci_hi - V.r], fmt='o', color=fc.INK, ms=4, mew=0, elinewidth=1.0, capsize=0, zorder=3)
a1.set_yticks(ys); a1.set_yticklabels(REG); a1.set_xlim(0, 1); a1.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
a1.set_xticklabels(['0', '0.25', '0.50', '0.75', '1.00']); a1.set_ylim(-0.6, len(REG) - 0.4)
a1.set_xlabel('Correlation, seasonally adjusted (95% CI)'); a1.grid(axis='x', color=fc.GRID, lw=0.4); a1.set_axisbelow(True)
a1.tick_params(axis='y', length=0)
mk = a2.scatter(V.ratio_cases, ys + 0.12, marker='o', s=22, color=fc.INK, linewidths=0, zorder=3, label='Cases')
ms = a2.scatter(V.ratio_search, ys - 0.12, marker='^', s=26, color=fc.BLUE, linewidths=0, zorder=3, label='Search interest')
a2.axvline(1, color=fc.INK2, lw=0.6, zorder=2)
a2.set_xscale('log'); a2.set_xlim(0.08, 15); a2.set_xticks([0.1, 0.3, 1, 3, 10]); a2.set_xticklabels(['0.1', '0.3', '1', '3', '10'])
a2.minorticks_off(); a2.set_yticks(ys); a2.set_yticklabels(REG); a2.set_ylim(-0.6, len(REG) - 0.4); a2.tick_params(axis='y', length=0)
a2.set_xlabel('2023 peak ÷ 2022 peak (log scale)'); a2.grid(axis='x', color=fc.GRID, lw=0.4); a2.set_axisbelow(True)
a2.legend(handles=[Line2D([], [], ls='', marker='o', color=fc.INK, ms=4.5, label='Cases'), Line2D([], [], ls='', marker='^', color=fc.BLUE, ms=5, label='Search interest')],
          loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, borderaxespad=0.3)
for ax, lab in zip((a1, a2), 'ab'): ax.text(-0.3, 1.02, f'({lab})', transform=ax.transAxes, fontweight='bold', va='bottom')
fig.subplots_adjust(left=0.115, right=0.985, top=0.86, bottom=0.2)
preview = fc.save_tiff_cmyk(fig, 'manuscript/figures/FigS3.tif')
pt = err.lines[0].get_data(); seg = err.lines[2][0].get_segments()
D = pd.DataFrame({'region': REG, 'y': pt[1], 'r_drawn': pt[0], 'lo_drawn': [s_[0][0] for s_ in seg], 'hi_drawn': [s_[1][0] for s_ in seg],
                  'cases_drawn': mk.get_offsets()[:, 0], 'search_drawn': ms.get_offsets()[:, 0],
                  'ytick_a': [t.get_text() for t in a1.get_yticklabels()], 'ytick_b': [t.get_text() for t in a2.get_yticklabels()]})
D['xlabel_a'] = a1.get_xlabel(); D['xlabel_b'] = a2.get_xlabel()
D.to_csv('out/figS3_plotted_values.csv', index=False, float_format='%.17g')
print('font:', font); print(V.round(4).to_string(index=False))
