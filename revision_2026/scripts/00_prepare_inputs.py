"""Build derived inputs: the preprint's 72-row table and the 3-retrieval mean of the primary Google series."""
import os, openpyxl, pandas as pd, numpy as np
# The preprint's 72-month table. Same 72 rows as data/raw/finaldata.xlsx in the public repository;
# set PREPRINT_XLSX to use that copy. Columns are found by their header names.
XLSX = os.environ.get('PREPRINT_XLSX', 'data_raw/paper_finalized_data_RP.xlsx')
ws = openpyxl.load_workbook(XLSX).active
rows = [[c.value for c in r] for r in ws.iter_rows(min_row=1)]
h = next(i for i, r in enumerate(rows) if 'Confirmed Cases' in r)
col = {name: rows[h].index(name) for name in ['Months', 'Dengue fever', 'Fever', 'Dengue Symptoms', 'Confirmed Cases']}
data = [r for r in rows[h + 1:] if r[col['Confirmed Cases']] is not None]
d = pd.DataFrame({k: [r[col[n]] for r in data] for k, n in
                  [('m', 'Months'), ('DF', 'Dengue fever'), ('Fever', 'Fever'), ('DS', 'Dengue Symptoms'), ('cases', 'Confirmed Cases')]})
assert len(d) == 72, len(d)
d['date'] = pd.date_range('2019-01-01', periods=len(d), freq='MS'); d = d.set_index('date')
d.to_csv('out/paper_72month_table.csv')
r = lambda f, lt: pd.read_csv(f, skiprows=2).iloc[:, 1].replace('<1', lt).astype(float).values
F = ['data_raw/gt_dengue_term.csv', 'data_raw/gt_dengue_term_day2_2026-11-22.csv', 'data_raw/gt_dengue_term_day3_2026-11-23.csv']
D = np.vstack([r(f, 0.5) for f in F])
m = pd.date_range('2019-01', periods=72, freq='MS').strftime('%Y-%m')
pd.DataFrame({'month': m, 'day1': D[0], 'day2': D[1], 'day3': D[2], 'mean': D.mean(0).round(3), 'range': D.max(0) - D.min(0)}).to_csv('out/gt_dengue_3day.csv', index=False)
