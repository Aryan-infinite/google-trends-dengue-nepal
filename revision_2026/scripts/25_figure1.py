"""Figure 1. Monthly reported dengue cases and Google search interest, Nepal, January 2019 - December 2024.

Inputs (the same files as scripts/09_final_analysis.py):
  out/paper_72month_table.csv  (reported cases; preprint 'Fever' series)  <- scripts/00_prepare_inputs.py
  data_raw/gt_dengue_term*.csv (three retrievals of "dengue"; "<1" coded 0.5; mean of the three)
Outputs:
  manuscript/figures/Fig1.tif           (CMYK, 1200 dpi, 170 mm wide)
  manuscript/figures/Fig1_preview.png   (RGB preview for checking)
  out/fig1_plotted_values.csv           (values read back from the drawn lines, for scripts/26_check_figure1.R)
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
import matplotlib.pyplot as plt, matplotlib.dates as mdates
from matplotlib.ticker import FixedLocator, NullFormatter
import fig_common as fc

os.makedirs('manuscript/figures', exist_ok=True)
IDX = pd.date_range('2019-01-01', periods=72, freq='MS')

paper = pd.read_csv('out/paper_72month_table.csv', index_col=0, parse_dates=True)
cases = pd.Series(paper.cases.astype(float).values, IDX)
cases[cases.index < '2019-06-01'] = np.nan          # January-May 2019 not counted by EDCD: missing, not zero
fever = pd.Series(paper.Fever.astype(float).values, IDX)
F = ['data_raw/gt_dengue_term.csv', 'data_raw/gt_dengue_term_day2_2026-11-22.csv',
     'data_raw/gt_dengue_term_day3_2026-11-23.csv']
read = lambda f: pd.read_csv(f, skiprows=2).iloc[:, 1].replace('<1', 0.5).astype(float).values
search = pd.Series(np.mean([read(f) for f in F], axis=0), IDX)

font = fc.setup()
fig, axes = plt.subplots(3, 1, figsize=(fc.WIDTH_MM * fc.MM, 150 * fc.MM), sharex=True,
                         gridspec_kw=dict(height_ratios=[1.25, 1, 1], hspace=0.18))
x = IDX + pd.Timedelta(days=14)                      # plot each month at its midpoint

a = axes[0]
la, = a.plot(x, cases.values, color=fc.INK, lw=1.0, marker='o', ms=1.8, mew=0)
a.set_yscale('log'); a.set_ylim(1, 60000)
a.yaxis.set_major_locator(FixedLocator([1, 10, 100, 1000, 10000]))
a.set_yticklabels(['1', '10', '100', '1,000', '10,000'])
a.yaxis.set_minor_formatter(NullFormatter()); a.tick_params(axis='y', which='minor', length=0)
a.axvspan(pd.Timestamp('2019-01-01'), pd.Timestamp('2019-06-01'), color='#ecebe7', lw=0)
a.text(pd.Timestamp('2019-03-16'), 30, 'Not recorded', ha='center', va='center', rotation=90, color=fc.INK2, fontsize=fc.FONT_PT - 1)
a.set_ylabel('Reported cases\n(log scale)')

b = axes[1]
lb, = b.plot(x, search.values, color=fc.BLUE, lw=1.0, marker='o', ms=1.8, mew=0)
b.set_ylim(0, 105); b.set_yticks([0, 25, 50, 75, 100])
b.set_ylabel('"dengue"\nsearch index')

c = axes[2]
lc, = c.plot(x, fever.values, color=fc.ORANGE, lw=1.0, marker='o', ms=1.8, mew=0)
c.set_ylim(0, 40); c.set_yticks([0, 10, 20, 30, 40])
c.set_ylabel('"Fever"\nsearch index')

for ax, lab in zip(axes, 'abc'):
    ax.grid(axis='y', color=fc.GRID, lw=0.4); ax.set_axisbelow(True)
    ax.text(-0.115, 1.0, f'({lab})', transform=ax.transAxes, fontweight='bold', va='top')
    for yr in range(2020, 2025):
        ax.axvline(pd.Timestamp(f'{yr}-01-01'), color=fc.GRID, lw=0.4, zorder=0)
c.set_xlim(pd.Timestamp('2019-01-01'), pd.Timestamp('2025-01-01'))
c.xaxis.set_major_locator(mdates.YearLocator(month=7, day=1))
c.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
c.xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[1]))
c.tick_params(axis='x', which='major', length=0)
c.tick_params(axis='x', which='minor', length=3)
fig.align_ylabels(axes)
fig.subplots_adjust(left=0.12, right=0.985, top=0.985, bottom=0.05)

preview = fc.save_tiff_cmyk(fig, 'manuscript/figures/Fig1.tif')

# Read the values back from the drawn lines, not from the input variables.
def back(line):
    xd, yd = line.get_data()
    return pd.to_datetime(xd).strftime('%Y-%m'), np.asarray(yd, dtype=float)
m, yA = back(la); _, yB = back(lb); _, yC = back(lc)
pd.DataFrame({'month': m, 'cases_plotted': yA, 'dengue_search_plotted': yB, 'fever_search_plotted': yC}) \
  .to_csv('out/fig1_plotted_values.csv', index=False, float_format='%.17g')
print('font:', font, '| preview:', preview)
print('cases plotted (non-missing):', int(np.isfinite(yA).sum()), 'max', np.nanmax(yA))
print('search max', yB.max(), 'at', m[int(np.argmax(yB))], '| fever max', yC.max(), 'at', m[int(np.argmax(yC))])
