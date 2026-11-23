"""Mutation test for scripts/26_check_figure1.R: each deliberately wrong version of the
plotted values must make the R check fail. The real file is restored afterwards."""
import shutil, subprocess, pandas as pd, numpy as np
F = 'out/fig1_plotted_values.csv'
shutil.copy(F, F + '.bak')
good = pd.read_csv(F); good[good.columns[1:]] = good[good.columns[1:]].astype(float)

def m1(d): d.loc[44, 'cases_plotted'] += 1; return d                       # one case count off by 1
def m2(d): d['dengue_search_plotted'] = d['dengue_search_plotted'].shift(1).bfill(); return d  # search shifted one month
def m3(d): d.loc[0:4, 'cases_plotted'] = 0; return d                         # Jan-May 2019 drawn as zero
def m4(d): d.loc[10, 'fever_search_plotted'] += 0.001; return d             # tiny change in the control series
def m5(d): d['dengue_search_plotted'] = np.round(d['dengue_search_plotted']); return d          # day-1 style rounding
mutants = dict(case_plus_1=m1, search_shifted=m2, missing_as_zero=m3, fever_tiny=m4, search_rounded=m5)
caught = 0
try:
    for name, f in mutants.items():
        f(good.copy()).to_csv(F, index=False, float_format='%.17g')
        r = subprocess.run(['Rscript', 'scripts/26_check_figure1.R'], capture_output=True, text=True)
        ok = r.returncode != 0; caught += ok
        print(f"{name:16s} {'CAUGHT' if ok else 'MISSED'}  ({[l for l in r.stdout.splitlines() if 'FAIL' in l]})")
finally:
    shutil.move(F + '.bak', F)
print(f'caught {caught}/{len(mutants)}')
r = subprocess.run(['Rscript', 'scripts/26_check_figure1.R'], capture_output=True, text=True)
print('restored file passes:', r.returncode == 0)
