"""Figure 4. Observed monthly dengue cases and nowcasts, Nepal, January 2022 - December 2024.
(a) per-season refit (primary, specified before results); (b) monthly refit (post hoc).

Models (Methods 2.4; scripts/09_final_analysis.py nowcast(); scripts/15_nowcast_monthly_refit.py), fitted by OLS to
y = log(cases + 1) on months from July 2019:
  M0: last month's count carried forward
  M1: 12 calendar-month terms + y(t-1)
  M2: M1 + log(search(t) + 0.5)          (search = mean of the three "dengue" retrievals; "<1" = 0.5)
Predictions = exp(fit) - 1, floored at 0. Per-season: one fit per test year on all earlier months.
Monthly: one fit per test month on all earlier months. Fitted here with numpy least squares (no formula interface).
Outputs: manuscript/figures/Fig4.tif, Fig4_preview.png; out/fig4_plotted_values.csv (read back from the drawn lines);
out/fig4_predictions_python.csv (every prediction, for the three-program comparison).
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
import matplotlib.pyplot as plt, matplotlib.dates as mdates
from matplotlib.ticker import FuncFormatter
import fig_common as fc

os.makedirs('manuscript/figures', exist_ok=True)
IDX = pd.date_range('2019-01-01', periods=72, freq='MS')
paper = pd.read_csv('out/paper_72month_table.csv', index_col=0, parse_dates=True)
cases = paper.cases.astype(float).values
F = ['data_raw/gt_dengue_term.csv', 'data_raw/gt_dengue_term_day2_2026-11-22.csv', 'data_raw/gt_dengue_term_day3_2026-11-23.csv']
read = lambda f: pd.read_csv(f, skiprows=2).iloc[:, 1].replace('<1', 0.5).astype(float).values
search = np.mean([read(f) for f in F], axis=0)
LY, LX, MON = np.log(cases + 1), np.log(search + 0.5), IDX.month.values

def design(rows, model):
    X = np.zeros((len(rows), 12)); X[np.arange(len(rows)), MON[rows] - 1] = 1          # month indicators, no intercept
    X = np.c_[X, LY[rows - 1]]
    return np.c_[X, LX[rows]] if model == 'M2' else X
def predict(train, test, model):
    b, *_ = np.linalg.lstsq(design(train, model), LY[train], rcond=None)
    return np.maximum(np.expm1(design(test, model) @ b), 0)

start = list(IDX).index(pd.Timestamp('2019-07-01'))
test = np.arange(list(IDX).index(pd.Timestamp('2022-01-01')), 72)
P = {}
for sched in ('season', 'monthly'):
    out = {'M0': cases[test - 1]}
    for m in ('M1', 'M2'):
        if sched == 'season':
            pr = np.concatenate([predict(np.arange(start, t0), np.array([t for t in test if IDX[t].year == y]), m)
                                 for y, t0 in [(y, list(IDX).index(pd.Timestamp(f'{y}-01-01'))) for y in (2022, 2023, 2024)]])
        else:
            pr = np.array([predict(np.arange(start, t), np.array([t]), m)[0] for t in test])
        out[m] = pr
    P[sched] = out
rows = []
for sched in P:
    for m in ('M0', 'M1', 'M2'):
        for t, v in zip(test, P[sched][m]): rows.append(dict(schedule=sched, model=m, month=IDX[t].strftime('%Y-%m'), pred=v))
pd.DataFrame(rows).to_csv('out/fig4_predictions_python.csv', index=False, float_format='%.17g')
obs = cases[test]
mae = {s: {m: np.mean(np.abs(P[s][m] - obs)) for m in ('M0', 'M1', 'M2')} for s in P}

font = fc.setup()
fig, axes = plt.subplots(2, 1, figsize=(fc.WIDTH_MM * fc.MM, 140 * fc.MM), sharex=True, gridspec_kw=dict(hspace=0.22))
x = IDX[test] + pd.Timedelta(days=14)
GREY = '#9c9b97'
style = [('M0', GREY, 'M0, last month carried forward'), ('M1', fc.BLUE, 'M1, seasonal terms + last month'),
         ('M2', fc.ORANGE, 'M2, M1 + current-month searches')]
drawn = {}
names = {'M0': 'M0', 'M1': 'M1', 'M2': 'M2'}
for ax, sched, title in zip(axes, ('season', 'monthly'), ('Refitted once per season (primary)', 'Refitted every month (post hoc)')):
    for m, col, lab in style:
        drawn[(sched, m)], = ax.plot(x, P[sched][m], color=col, lw=1.0, label=lab, zorder=3)
    drawn[(sched, 'obs')], = ax.plot(x, obs, color=fc.INK, lw=1.2, marker='o', ms=2.2, mew=0, label='Observed', zorder=4)
    ax.set_ylim(0, 40000); ax.set_yticks(range(0, 40001, 10000))      # M1 Sep 2023 reaches 37,044: nothing is clipped
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:,.0f}'))
    ax.grid(axis='y', color=fc.GRID, lw=0.4); ax.set_axisbelow(True)
    for yr in (2023, 2024): ax.axvline(pd.Timestamp(f'{yr}-01-01'), color=fc.GRID, lw=0.4, zorder=0)
    ax.set_ylabel('Cases per month')
    # every series is below 1,000 cases before August 2022 and below 70% of the axis in 2022: the upper-left corner is empty
    ax.text(0.015, 0.97, f'{title}\nMean absolute error: M0 {mae[sched]["M0"]:,.0f}, M1 {mae[sched]["M1"]:,.0f}, M2 {mae[sched]["M2"]:,.0f}',
            transform=ax.transAxes, va='top', ha='left', color=fc.INK, linespacing=1.4)
h, l = axes[0].get_legend_handles_labels(); order = [3, 0, 1, 2]
fig.legend([h[i] for i in order], [l[i] for i in order], loc='upper center', ncol=2, frameon=False,
           bbox_to_anchor=(0.545, 1.0), handlelength=1.8, columnspacing=1.2)
for ax, lab in zip(axes, 'ab'):
    ax.text(-0.095, 1.035, f'({lab})', transform=ax.transAxes, fontweight='bold', va='bottom')
ax = axes[1]
ax.set_xlim(pd.Timestamp('2022-01-01'), pd.Timestamp('2025-01-01'))
ax.set_xticks(pd.date_range('2022-01-01', '2024-10-01', freq='QS'))
ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: mdates.num2date(v).strftime('%b\n%Y') if mdates.num2date(v).month == 1 else mdates.num2date(v).strftime('%b')))
fig.subplots_adjust(left=0.1, right=0.985, top=0.865, bottom=0.08)
preview = fc.save_tiff_cmyk(fig, 'manuscript/figures/Fig4.tif')

# Read back the drawn lines and the legend text.
rec = []
for (sched, m), ln in drawn.items():
    xd, yd = ln.get_data()
    for xx, yy in zip(pd.to_datetime(xd).strftime('%Y-%m'), yd): rec.append(dict(schedule=sched, series=m, month=xx, y=float(yy)))
pd.DataFrame(rec).to_csv('out/fig4_plotted_values.csv', index=False, float_format='%.17g')
key = [t.get_text() for t in fig.legends[0].get_texts()]
boxes = [[t.get_text() for t in a.texts if not t.get_text().startswith('(')][0].replace('\n', ' | ') for a in axes]
pd.DataFrame({'item': ['key', 'panel_a', 'panel_b'], 'text': [' | '.join(key)] + boxes}).to_csv('out/fig4_text.csv', index=False)
print('font:', font, '| max plotted', round(max(max(P[s][m].max() for m in P[s]) for s in P)), '| observed max', obs.max())
for s in P: print(s, {m: round(v, 6) for m, v in mae[s].items()})
