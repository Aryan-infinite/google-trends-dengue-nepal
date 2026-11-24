"""Mutation test for scripts/29_check_figure2.R: each deliberately wrong version of the Figure 2 outputs
must make the R check fail. The real files are restored afterwards."""
import shutil, subprocess, numpy as np, pandas as pd
P, S = 'out/fig2_plotted_values.csv', 'out/fig2_stats.csv'
for f in (P, S): shutil.copy(f, f + '.bak')
gp, gs = pd.read_csv(P), pd.read_csv(S, dtype=str)

def mp(f):                                   # mutate plotted points
    def run(): f(gp.copy()).to_csv(P, index=False, float_format='%.17g')
    return run
def ms(q, v):                                # mutate one stats value
    def run(): d = gs.copy(); d.loc[d.quantity == q, 'value'] = v; d.to_csv(S, index=False)
    return run
def shift_point(d): i = d.index[(d.panel == 'adj') & (d.month == '2023-08')][0]; d.loc[i, 'y'] += 0.01; return d
def wrong_group(d): d.loc[(d.month == '2022-01'), 'group'] = '2020–2021'; return d
def drop_point(d): return d.drop(d.index[(d.panel == 'raw') & (d.month == '2024-12')])
def plus_one(d): r = d.panel == 'raw'; d.loc[r, 'x'] = np.log(np.exp(d.loc[r, 'x']) - 0.5 + 1); return d   # log(s+1) instead of log(s+0.5)
mutants = {'adj_point_moved': mp(shift_point), 'group_mislabelled': mp(wrong_group), 'point_dropped': mp(drop_point),
           'log_offset_1_not_0.5': mp(plus_one), 'r_label_wrong': ms('r_adj', '0.8'), 'wrong_month_flagged': ms('flagged', '2020-06;2022-10'),
           'drawn_text_wrong': ms('text_b', 'Jun 2020|Sep 2022|r = 0.82')}
caught = 0
try:
    for name, mutate in mutants.items():
        shutil.copy(P + '.bak', P); shutil.copy(S + '.bak', S); mutate()
        r = subprocess.run(['Rscript', 'scripts/29_check_figure2.R'], capture_output=True, text=True)
        ok = r.returncode != 0; caught += ok
        print(f"{name:22s} {'CAUGHT' if ok else 'MISSED'}  {[l for l in r.stdout.splitlines() if 'FAIL' in l]}")
finally:
    for f in (P, S): shutil.move(f + '.bak', f)
print(f'caught {caught}/{len(mutants)}')
r = subprocess.run(['Rscript', 'scripts/29_check_figure2.R'], capture_output=True, text=True)
print('restored files pass:', r.returncode == 0)
