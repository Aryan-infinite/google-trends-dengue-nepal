"""Mutation tests for every Figure 3 check.
Level 1: falsified out/fig3_plotted_values.csv must fail scripts/36_check_figure3.R.
Level 2: falsified Fig3.tif images (re-rendered from a patched copy of scripts/35_figure3.py, or edited directly)
must fail at least the intended image check(s): scripts/38 (Python digitiser), 39 (R digitiser), 34 (OCR).
Runs in a temporary copy; the project is not modified.
"""
import os, shutil, subprocess, tempfile
import pandas as pd
from PIL import Image, ImageDraw
Image.MAX_IMAGE_PIXELS = None
SRC = os.getcwd(); T = tempfile.mkdtemp(prefix='f3mut_')
for d in ('data_raw', 'scripts'): shutil.copytree(os.path.join(SRC, d), os.path.join(T, d))
os.makedirs(os.path.join(T, 'out/final'), exist_ok=True); os.makedirs(os.path.join(T, 'manuscript/figures'), exist_ok=True)
for f in sorted(os.listdir(os.path.join(SRC, 'manuscript/figures'))):        # every figure built (the OCR and residual-ink checks read all)
    if f.endswith('.tif'): shutil.copy(os.path.join(SRC, 'manuscript/figures', f), os.path.join(T, 'manuscript/figures', f))
shutil.copy(os.path.join(SRC, 'out/final/final_results.xlsx'), os.path.join(T, 'out/final/final_results.xlsx'))
run = lambda c: subprocess.run(c, cwd=T, capture_output=True, text=True)
assert run(['python3', '-I', 'scripts/00_prepare_inputs.py']).returncode == 0
assert run(['python3', '-I', 'scripts/35_figure3.py']).returncode == 0
CSV = os.path.join(T, 'out/fig3_plotted_values.csv'); good = pd.read_csv(CSV, dtype={'annotation': str}); good['annotation'] = good['annotation'].fillna('')
TIF = os.path.join(T, 'manuscript/figures/Fig3.tif'); shutil.copy(TIF, TIF + '.good')
SCRIPT = os.path.join(T, 'scripts/35_figure3.py'); ORIG = open(SCRIPT).read()

print('Level 1: plotted-values file vs scripts/36_check_figure3.R')
def c1(name, f):
    d = f(good.copy()); d.to_csv(CSV, index=False, float_format='%.17g')
    r = run(['Rscript', 'scripts/36_check_figure3.R']); print(f"  {name:38s} {'CAUGHT' if r.returncode else 'MISSED'}"); return r.returncode != 0
def setv(col, lag, val):
    def f(d): d.loc[d.lag == lag, col] = val; return d
    return f
L1 = [c1('bar +0.001 at lag +2', setv('bar_height', 2, good.loc[good.lag == 2, 'bar_height'].item() + 0.001)),
      c1('threshold for n-1 at lag 0', setv('crit_pos', 0, 0.2608689)),
      c1('Holm label on lag +1 instead of 0', lambda d: (d.assign(annotation=['', '', '', '', 'Holm-adjusted P < 0.001', '', '']))),
      c1('tick label says n = 57 at lag 0', setv('tick_label', 0, '0 | n = 57')),
      c1('AR order recorded as 3', lambda d: d.assign(ar_order=3))]
good.to_csv(CSV, index=False, float_format='%.17g')

print('Level 2: falsified images vs scripts/38 (Py), 39 (R), 34 (OCR)')
def image_checks():
    out = {}
    for k, cmd in [('Py', ['python3', '-I', 'scripts/38_digitise_figure3.py']), ('R', ['Rscript', 'scripts/39_digitise_figure3.R']),
                   ('OCR', ['python3', '-I', 'scripts/34_ocr_figures.py'])]:
        r = run(cmd); out[k] = [l.split('] ')[1].split('  (')[0] for l in r.stdout.splitlines() if l.startswith('[FAIL]') and (k != 'OCR' or 'F3' in l)]
    return out
def patched(after_T):
    def make():
        open(SCRIPT, 'w').write(ORIG.replace("T = pd.DataFrame(rows); T['p_holm'] = multipletests(T.p, method='holm')[1]",
                                            "T = pd.DataFrame(rows); T['p_holm'] = multipletests(T.p, method='holm')[1]\n" + after_T))
        assert run(['python3', '-I', 'scripts/35_figure3.py']).returncode == 0
        open(SCRIPT, 'w').write(ORIG)
    return make
def stray():
    im = Image.open(TIF); im.load(); fill = im.getpixel((im.width // 2, im.height // 2))
    import numpy as np
    A = np.asarray(im); yy, xx = np.where((A[..., 0] > 150) & (A[..., 1] < 150) & (A[..., 3] < 60))   # bar ink: high cyan, low black
    blue = tuple(int(v) for v in A[yy[len(yy) // 2], xx[len(xx) // 2]])
    ImageDraw.Draw(im).ellipse([int(im.width * 0.8), int(im.height * 0.62), int(im.width * 0.8) + 20, int(im.height * 0.62) + 20], fill=blue)
    im.save(TIF, compression='tiff_lzw', dpi=(1200, 1200))
M2 = {'bar +0.003 at lag +1 (about 6 px)': (patched("T.loc[T.lag == 1, 'r'] += 0.003"), ['Py', 'R']),
      'bar +0.01 at lag +2': (patched("T.loc[T.lag == 2, 'r'] += 0.01"), ['Py', 'R']),
      'thresholds computed with n - 3': (patched("T['r_crit'] = [st.t.ppf(.975, n - 5) / np.sqrt(n - 5 + st.t.ppf(.975, n - 5) ** 2) for n in T.n]"), ['Py', 'R']),
      'Holm label moved to lag +1': (patched("T.loc[T.lag == 1, 'p_holm'] = 0.0001; T.loc[T.lag == 0, 'p_holm'] = 0.5"), ['OCR']),
      'n under lag 0 printed as 57': (patched("T.loc[T.lag == 0, 'n'] = 57"), ['OCR']),
      'stray blue dot (r = 10 px)': (stray, ['Py', 'R']),
      'bar sign flipped at lag +3 (0.009 -> -0.009)': (patched("T.loc[T.lag == 3, 'r'] *= -1"), ['Py', 'R'])}
L2 = []
for name, (make, need) in M2.items():
    shutil.copy(TIF + '.good', TIF); make(); res = image_checks()
    ok = all(res[k] for k in need); L2.append(ok)
    print(f"  {name:44s} {'CAUGHT' if ok else 'MISSED'}  " + '; '.join(f'{k}: {len(v)} fail' for k, v in res.items()))
shutil.copy(TIF + '.good', TIF); res = image_checks()
print(f'\nLevel 1 caught {sum(L1)}/{len(L1)}; level 2 caught {sum(L2)}/{len(L2)} (each by every check it targets); '
      f'true image passes all image checks: {not any(res.values())}')
shutil.rmtree(T)
