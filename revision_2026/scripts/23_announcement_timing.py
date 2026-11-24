"""Check E (ANALYSIS_PLAN.md section 9): per year, month the primary search index first reaches half its annual maximum,
search peak month, and case peak month; joined with the announcement events logged in governance/announcement_events.csv."""
import csv, json, numpy as np
def gt(f): return [0.5 if r[1] == '<1' else float(r[1]) for r in csv.reader(open('data_raw/' + f, encoding='utf-8-sig')) if len(r) == 2 and len(r[0]) == 7 and r[0][:2] == '20']
x = np.mean([gt('gt_dengue_term.csv'), gt('gt_dengue_term_day2_2026-11-22.csv'), gt('gt_dengue_term_day3_2026-11-23.csv')], axis=0)
import os, openpyxl
ws = openpyxl.load_workbook(os.environ.get('PREPRINT_XLSX', 'data_raw/paper_finalized_data_RP.xlsx'), data_only=True).worksheets[0]
rows = [list(r) for r in ws.iter_rows(values_only=True)]                       # preprint table; columns found by header name
h = next(i for i, r in enumerate(rows) if 'Confirmed Cases' in r); jm, jc = rows[h].index('Months'), rows[h].index('Confirmed Cases')
c = np.array([float(r[jc]) for r in sorted((r for r in rows[h + 1:] if isinstance(r[jm], (int, float))), key=lambda r: r[jm])])
assert len(x) == 72 and len(c) == 72
R = {}
for i, y in enumerate(range(2019, 2025)):
    xs, cs = x[i*12:(i+1)*12], c[i*12:(i+1)*12]
    if y == 2019: cs = cs.copy(); cs[:5] = np.nan          # January-May 2019 not recorded
    R[f'E_{y}_search_half_month'] = int(np.argmax(xs >= xs.max() / 2) + 1)
    R[f'E_{y}_search_peak_month'] = int(np.argmax(xs) + 1)
    R[f'E_{y}_case_half_month'] = int(np.nanargmax(np.where(np.isnan(cs), -1, cs) >= np.nanmax(cs) / 2) + 1)
    R[f'E_{y}_case_peak_month'] = int(np.nanargmax(np.where(np.isnan(cs), -1, cs)) + 1)
json.dump(R, open('out/23_announcement_timing.json', 'w'), indent=1)
for y in range(2019, 2025): print(y, *(R[f'E_{y}_{k}'] for k in ['search_half_month', 'search_peak_month', 'case_half_month', 'case_peak_month']))
