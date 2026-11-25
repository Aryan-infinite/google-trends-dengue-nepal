"""Supplementary Figure S2: timing of search interest, cases and official announcements within each year, 2019-2024.
Each series is divided by its own annual maximum (cases 2019: June-December only; January-May 2019 not recorded).
Events: governance/announcement_events.csv (one row per opened page) and the five dated 2024 counts in
out/06_2024_reporting_lag.csv (the 31 December annual total is not an in-season announcement and is not drawn).
Exports: out/figS2_values.csv (computed), out/figS2_plotted_values.csv (read back from the drawn artists),
out/figS2_events_drawn.csv, out/figS2_text.csv."""
import os, sys, calendar, csv
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd, openpyxl
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import fig_common as fc
fc.setup()

rows = [[c.value for c in r] for r in openpyxl.load_workbook(os.environ.get('PREPRINT_XLSX', 'data_raw/paper_finalized_data_RP.xlsx')).active.iter_rows()]
h = next(i for i, r in enumerate(rows) if 'Confirmed Cases' in r); ci = rows[h].index('Confirmed Cases')
cases = np.array([float(r[ci]) for r in rows[h + 1:] if r[ci] is not None]); assert len(cases) == 72
g = lambda f: pd.read_csv(f, skiprows=2).iloc[:, 1].astype(str).replace('<1', '0.5').astype(float).values
srch = np.mean([g('data_raw/gt_dengue_term.csv'), g('data_raw/gt_dengue_term_day2_2026-11-22.csv'), g('data_raw/gt_dengue_term_day3_2026-11-23.csv')], 0)
assert len(srch) == 72

CAT = {'national response': 'alert', 'national alert': 'alert', 'official count (provincial)': 'count', 'official count': 'count',
       'official count; preparatory actions (no formal alert)': 'count', 'formal outbreak reporting begins': 'report',
       'international report': 'report'}
ev = []
for r in csv.DictReader(open('governance/announcement_events.csv', encoding='utf-8')):
    if r['event_type'] in CAT: ev.append((r['date'], CAT[r['event_type']]))
    else: assert r['event_type'] in ('no national alert found', 'official counts (5 dates)'), r['event_type']
for r in csv.DictReader(open('out/06_2024_reporting_lag.csv', encoding='utf-8')):
    if r['as_of'] != '2024-12-31': ev.append((r['as_of'], 'count'))
assert len(ev) == 15, len(ev)
def xpos(d):
    y, m, dd = map(int, d.split('-')); return y, (m - 1) + (dd - 0.5) / calendar.monthrange(y, m)[1]

V = []
for i, y in enumerate(range(2019, 2025)):
    c = cases[12 * i:12 * i + 12].copy(); s = srch[12 * i:12 * i + 12]
    if y == 2019: c[:5] = np.nan
    for m in range(12):
        V.append(dict(year=y, month=m + 1, x=m + 0.5, cases=c[m], search=s[m], cases_rel=c[m] / np.nanmax(c), search_rel=s[m] / s.max()))
V = pd.DataFrame(V); V.to_csv('out/figS2_values.csv', index=False, float_format='%.17g')

SHAPE = {'alert': dict(marker='v', mfc=fc.INK2, mec=fc.INK2), 'count': dict(marker='o', mfc='white', mec=fc.INK2, mew=0.8),
         'report': dict(marker='s', mfc=fc.INK2, mec=fc.INK2)}
EY, MS = 1.15, 4.2
fig, axs = plt.subplots(3, 2, figsize=(fc.WIDTH_MM * fc.MM, 150 * fc.MM), sharex=True, sharey=True,
                        gridspec_kw=dict(hspace=0.32, wspace=0.12))
P, EVD, TXT = [], [], []
for k, (ax, y) in enumerate(zip(axs.flat, range(2019, 2025))):
    d = V[V.year == y]
    ax.set_xlim(0, 12); ax.set_ylim(0, 1.25); ax.set_yticks([0, 0.5, 1]); ax.set_yticklabels(['0', '0.5', '1'])
    ax.set_xticks(np.arange(12) + 0.5); ax.set_xticklabels(list('JFMAMJJASOND'))
    ax.grid(axis='y', color=fc.GRID, lw=0.4); ax.set_axisbelow(True)
    ls, = ax.plot(d.x, d.search_rel, color=fc.BLUE, lw=1.0, zorder=3)
    lc, = ax.plot(d.x, d.cases_rel, color=fc.INK, lw=1.2, zorder=4)
    for e, (nm, ln) in enumerate((('search', ls), ('cases', lc))):
        for xx, yy in zip(*ln.get_data()):
            P.append(dict(year=y, series=nm, x=xx, y=yy))
    for dt, cat in ev:
        ey, ex = xpos(dt)
        if ey == y:
            mk, = ax.plot([ex], [EY], ls='none', ms=MS, zorder=5, clip_on=False, **SHAPE[cat])
            EVD.append(dict(date=dt, year=y, category=cat, x=mk.get_xdata()[0], y=mk.get_ydata()[0], marker=mk.get_marker()))
    t = ax.set_title(f'({"abcdef"[k]}) {y}', loc='left', fontweight='bold', color=fc.INK, fontsize=fc.FONT_PT, pad=3); TXT.append(('panel', y, t.get_text()))
    notes = {2019: 'Cases not recorded Jan–May', 2020: 'Search index ≤ 3 all year', 2021: 'Search index ≤ 2 all year'}
    if y in notes:
        t = ax.set_title(notes[y], loc='right', color=fc.INK2, fontsize=fc.FONT_PT - 1, pad=3); TXT.append(('note', y, t.get_text()))
    if y == 2021:
        t = ax.text(6.0, EY, 'No national alert found', ha='center', va='center', color=fc.INK2, fontsize=fc.FONT_PT - 1); TXT.append(('strip', y, t.get_text()))
for ax in axs[:, 0]: ax.set_ylabel('Share of annual maximum')
for ax in axs[-1, :]: ax.set_xlabel('Month')
H = [Line2D([], [], color=fc.INK, lw=1.2, label='Cases'), Line2D([], [], color=fc.BLUE, lw=1.0, label='Search interest'),
     Line2D([], [], ls='none', ms=MS, label='National alert or EDCD statement', **SHAPE['alert']),
     Line2D([], [], ls='none', ms=MS, label='Count announced', **SHAPE['count']),
     Line2D([], [], ls='none', ms=MS, label='Daily reports began / WHO report', **SHAPE['report'])]
lg = fig.legend(handles=H, loc='upper center', ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.0), handlelength=1.8, columnspacing=1.6)
for t in lg.get_texts(): TXT.append(('legend', 0, t.get_text()))
fig.subplots_adjust(left=0.085, right=0.985, top=0.885, bottom=0.075)

pd.DataFrame(P).to_csv('out/figS2_plotted_values.csv', index=False, float_format='%.17g')
pd.DataFrame(EVD).to_csv('out/figS2_events_drawn.csv', index=False, float_format='%.17g')
pd.DataFrame(TXT, columns=['kind', 'year', 'text']).to_csv('out/figS2_text.csv', index=False)
os.makedirs('manuscript/figures', exist_ok=True)
fc.save_tiff_cmyk(fig, 'manuscript/figures/FigS2.tif')
print(V.groupby('year').apply(lambda d: pd.Series(dict(search_max=d.search.max(), case_max=np.nanmax(d.cases))), include_groups=False))
print(len(EVD), 'events drawn')
