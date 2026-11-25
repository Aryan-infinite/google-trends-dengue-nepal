"""Figure 2. Monthly "dengue" search interest against reported dengue cases, Nepal, 2020-2024 (n = 60).
(a) log scales; (b) after subtracting each calendar month's 2020-2024 mean (the primary Question 1 analysis).
Months flagged by Cook's distance (> 4/n) in the regression model of scripts/09_final_analysis.py are labelled in (b).

Inputs: out/paper_72month_table.csv, data_raw/gt_dengue_term*.csv (as scripts/09_final_analysis.py)
Outputs: manuscript/figures/Fig2.tif, Fig2_preview.png; out/fig2_plotted_values.csv (read back from the drawn points)
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import fig_common as fc

os.makedirs('manuscript/figures', exist_ok=True)
IDX = pd.date_range('2019-01-01', periods=72, freq='MS')
paper = pd.read_csv('out/paper_72month_table.csv', index_col=0, parse_dates=True)
cases = pd.Series(paper.cases.astype(float).values, IDX)
F = ['data_raw/gt_dengue_term.csv', 'data_raw/gt_dengue_term_day2_2026-11-22.csv',
     'data_raw/gt_dengue_term_day3_2026-11-23.csv']
read = lambda f: pd.read_csv(f, skiprows=2).iloc[:, 1].replace('<1', 0.5).astype(float).values
search = pd.Series(np.mean([read(f) for f in F], axis=0), IDX)

c, s = cases['2020-01':'2024-12'], search['2020-01':'2024-12']           # primary window, n = 60
Y, X = np.log(c + 1), np.log(s + 0.5)
adj = lambda v: v - v.groupby(v.index.month).transform('mean')
Ya, Xa = adj(Y), adj(X)
r_raw, r_adj = np.corrcoef(X, Y)[0, 1], np.corrcoef(Xa, Ya)[0, 1]

# Cook's distance from the same OLS used in scripts/09 (log(cases+1) on harmonics, lagged log cases, log search)
d = pd.DataFrame({'c': c, 'y1': np.log1p(c).shift(1), 'x': X})
mo = d.index.month
d['s1'], d['c1'] = np.sin(2 * np.pi * mo / 12), np.cos(2 * np.pi * mo / 12)
d['s2'], d['c2'] = np.sin(4 * np.pi * mo / 12), np.cos(4 * np.pi * mo / 12)
d = d.dropna()
cook = smf.ols('np.log1p(c) ~ s1 + c1 + s2 + c2 + y1 + x', d).fit().get_influence().cooks_distance[0]
flagged = list(d.index[cook > 4 / len(d)].strftime('%Y-%m'))

late = (c.index.year >= 2022)
font = fc.setup()
fig, axes = plt.subplots(1, 2, figsize=(fc.WIDTH_MM * fc.MM, 78 * fc.MM),
                         gridspec_kw=dict(wspace=0.32))
groups = [(~late, fc.BLUE, 'o', '2020–2021'), (late, fc.ORANGE, '^', '2022–2024')]
drawn = {}
for ax, xs, ys, key in [(axes[0], X, Y, 'raw'), (axes[1], Xa, Ya, 'adj')]:
    for mask, col, mk, lab in groups:
        drawn[(key, lab)] = ax.scatter(xs[mask], ys[mask], s=14, c=col, marker=mk, alpha=0.75,
                                       linewidths=0.5, edgecolors='white', zorder=3)
    ax.grid(color=fc.GRID, lw=0.4); ax.set_axisbelow(True)

a, b = axes
xt = [1, 2, 5, 10, 25, 100]; yt = [1, 10, 100, 1000, 10000]
a.set_xticks(np.log(np.array(xt) + 0.5)); a.set_xticklabels([str(v) for v in xt])
a.set_yticks(np.log(np.array(yt) + 1)); a.set_yticklabels(['1', '10', '100', '1,000', '10,000'])
a.set_xlabel('"dengue" search index (log scale)'); a.set_ylabel('Reported cases (log scale)')
b.set_xlabel('Search, deviation from monthly mean (log)'); b.set_ylabel('Cases, deviation from monthly mean (log)')
b.axhline(0, color=fc.INK2, lw=0.5, zorder=1); b.axvline(0, color=fc.INK2, lw=0.5, zorder=1)
for ax, r in [(a, r_raw), (b, r_adj)]:
    ax.text(0.04, 0.96, f'r = {r:.2f}', transform=ax.transAxes, va='top', color=fc.INK)
for m in flagged:
    t = pd.Timestamp(m)
    right = Xa[t] > 1.5                                   # keep labels inside the panel
    b.annotate(t.strftime('%b %Y'), (Xa[t], Ya[t]), xytext=(4, -50) if right else (14, -16),
               textcoords='offset points', ha='right' if right else 'left', va='center',
               arrowprops=dict(arrowstyle='-', color=fc.INK2, lw=0.5, shrinkA=1, shrinkB=2.5),
               fontsize=fc.FONT_PT - 1, color=fc.INK2)
for ax, lab in zip(axes, 'ab'):
    ax.text(-0.24, 1.04, f'({lab})', transform=ax.transAxes, fontweight='bold', va='bottom')
handles = [Line2D([], [], ls='', marker=mk, color=col, markersize=4.5, markeredgecolor='white', markeredgewidth=0.5, label=lab)
           for _, col, mk, lab in groups]
fig.legend(handles=handles, loc='upper right', ncol=2, frameon=False, bbox_to_anchor=(0.995, 1.0),
           handletextpad=0.3, columnspacing=1.2)
fig.subplots_adjust(left=0.115, right=0.985, top=0.86, bottom=0.15)

preview = fc.save_tiff_cmyk(fig, 'manuscript/figures/Fig2.tif')

# Read back the drawn points (offsets of each scatter collection), not the input variables.
rows = []
for (key, lab), coll in drawn.items():
    mask = ~late if lab == '2020–2021' else late
    for m, (px, py) in zip(c.index[mask].strftime('%Y-%m'), coll.get_offsets()):
        rows.append(dict(panel=key, group=lab, month=m, x=px, y=py))
out = pd.DataFrame(rows).sort_values(['panel', 'month'])
out.to_csv('out/fig2_plotted_values.csv', index=False, float_format='%.17g')
txt = lambda ax: '|'.join(sorted(t.get_text() for t in ax.texts if not t.get_text().startswith('(')))
pd.DataFrame({'quantity': ['r_raw', 'r_adj', 'n', 'flagged', 'text_a', 'text_b'],
              'value': [f'{r_raw:.17g}', f'{r_adj:.17g}', str(len(c)), ';'.join(flagged), txt(a), txt(b)]}) \
  .to_csv('out/fig2_stats.csv', index=False)
print('font:', font, '| r_raw', round(r_raw, 4), '| r_adj', round(r_adj, 4), '| n', len(c), '| Cook flagged', flagged,
      '| points drawn', len(out))
