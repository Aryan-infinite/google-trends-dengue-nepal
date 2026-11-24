"""All national and regional peak-month classifications at margins of 5%, 10% and 20% (post hoc, reported in full)."""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
exec(open('scripts/09_final_analysis.py').read().split('# Q5 regional')[0])
def pk(s, y, mg):
    v = s[str(y)].sort_values(ascending=False); r = v.iloc[0] / v.iloc[1] if v.iloc[1] > 0 else np.inf
    return (v.index[0].strftime('%b') if r >= 1 + mg else 'tie'), r
rows = []
sets = [('national cases', cases, range(2020, 2025)), ('national search', S['dengue_3day'], range(2020, 2025))]
sets += [(f'{rg} cases', regC[rg], range(2020, 2024)) for rg in ['Central', 'Eastern', 'Western']]
sets += [(f'{rg} search', regS[rg], range(2020, 2024)) for rg in ['Central', 'Eastern', 'Western']]
for lab, s, yrs in sets:
    for y in yrs:
        (a, r), (b, _), (c, _) = pk(s, y, .05), pk(s, y, .10), pk(s, y, .20)
        rows.append(dict(series=lab, year=y, top_to_second=r, m05=a, m10=b, m20=c, changes=len({a, b, c}) > 1))
t = pd.DataFrame(rows); t.to_csv('out/20_peak_margin_full.csv', index=False); print(t.round(3).to_string()); print('changed:', int(t.changes.sum()), 'of', len(t))
