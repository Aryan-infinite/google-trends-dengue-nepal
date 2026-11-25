"""Mutation tests for every Figure S1 check (runs in a temporary copy; the project is not modified).
Level 1: falsified value / drawn / title files must fail scripts/51_check_figureS1.R (three-program check).
Level 2: falsified FigS1.tif images must fail the image checks they target: 53 (Py digitiser), 54 (R digitiser),
34 (OCR, S1 lines), 47 / 48 (residual ink, S1 lines)."""
import os, shutil, subprocess, tempfile

def _guarded_run(*a, **k):
    """A check killed by the system (e.g. out of memory) prints nothing; stop instead of counting it as caught or missed."""
    r = subprocess.run(*a, **k)
    if r.returncode < 0 or r.returncode == 137:
        raise SystemExit(f'ABORTED: a check was killed (exit {r.returncode}): {a[0] if a else k.get("args")}; rerun with nothing else running')
    return r
import numpy as np, pandas as pd
from PIL import Image, ImageDraw
Image.MAX_IMAGE_PIXELS = None
SRC = os.getcwd(); T = tempfile.mkdtemp(prefix='s1mut_')
for d in ('data_raw', 'scripts'): shutil.copytree(os.path.join(SRC, d), os.path.join(T, d))
os.makedirs(os.path.join(T, 'out'), exist_ok=True); os.makedirs(os.path.join(T, 'manuscript/figures'), exist_ok=True)
NODE_MODULES = os.path.join(SRC, 'node_modules') if os.path.isdir(os.path.join(SRC, 'node_modules')) else os.environ.get('NODE_PATH', '')  # npm install
for f in sorted(os.listdir(os.path.join(SRC, 'manuscript/figures'))):        # every figure built (the OCR and residual-ink checks read all)
    if f.endswith('.tif'): shutil.copy(os.path.join(SRC, 'manuscript/figures', f), os.path.join(T, 'manuscript/figures', f))
for f in ('21_reviewer_checks.json', '22_reviewer_checks.json'): shutil.copy(os.path.join(SRC, 'out', f), os.path.join(T, 'out'))
run = lambda c, env=None: _guarded_run(c, cwd=T, capture_output=True, text=True, env={**os.environ, **(env or {})})
assert run(['python3', '-I', 'scripts/00_prepare_inputs.py']).returncode == 0
assert run(['python3', '-I', 'scripts/50_figureS1.py']).returncode == 0
assert run(['node', 'scripts/52_acf_independent.js'], {'NODE_PATH': NODE_MODULES}).returncode == 0
assert run(['Rscript', 'scripts/51_check_figureS1.R']).returncode == 0, 'true files must pass first'
O = {f: os.path.join(T, 'out', f) for f in ('figS1_values_js.csv', 'figS1_plotted_values.csv', 'figS1_values.csv')}
for p in O.values(): shutil.copy(p, p + '.good')
print('Level 1: files vs scripts/51_check_figureS1.R')
L1 = []
def c1(name, f, edit):
    for p in O.values(): shutil.copy(p + '.good', p)
    d = pd.read_csv(O[f]); edit(d); d.to_csv(O[f], index=False, float_format='%.17g')
    r = run(['Rscript', 'scripts/51_check_figureS1.R']); L1.append(r.returncode != 0)
    print(f"  {name:50s} {'CAUGHT' if L1[-1] else 'MISSED'}")
c1('JavaScript pacf7 (log search) off by 1e-9', 'figS1_values_js.csv', lambda d: d.__setitem__('value', np.where((d.series == 'lx_raw') & (d.stat == 'pacf7'), d.value + 1e-9, d.value)))
c1('drawn bar acf12 (log cases) +0.001', 'figS1_plotted_values.csv', lambda d: d.__setitem__('height', np.where((d.series == 'ly_raw') & (d.kind == 'acf') & (d.lag == 12), d.height + 0.001, d.height)))
c1('drawn band at 0.25 instead of 0.253', 'figS1_plotted_values.csv', lambda d: d.__setitem__('height', np.where(d.lag == 0, 0.25, d.height)))
c1('Python Ljung-Box Q12 (adjusted cases) +0.1', 'figS1_values.csv', lambda d: d.__setitem__('value', np.where((d.series == 'ly_adj') & (d.stat == 'lbq12'), d.value + 0.1, d.value)))
for p in O.values(): shutil.copy(p + '.good', p)

print('Level 2: falsified images')
S = os.path.join(T, 'scripts/50_figureS1.py'); ORIG = open(S).read(); TIF = os.path.join(T, 'manuscript/figures/FigS1.tif'); shutil.copy(TIF, TIF + '.good')
def patch(a, b):
    def make():
        assert a in ORIG; open(S, 'w').write(ORIG.replace(a, b)); r = run(['python3', '-I', 'scripts/50_figureS1.py']); assert r.returncode == 0, r.stderr; open(S, 'w').write(ORIG)
    return make
def dot(xy, r, fill):
    def make():
        im = Image.open(TIF); im.load(); x, y = xy; ImageDraw.Draw(im).ellipse([x - r, y - r, x + r, y + r], fill=fill); im.save(TIF, compression='tiff_lzw', dpi=(1200, 1200))
    return make
CH = {'Py53': ['python3', '-I', 'scripts/53_digitise_figureS1.py'], 'R54': ['Rscript', 'scripts/54_digitise_figureS1.R'],
      'OCR': ['python3', '-I', 'scripts/34_ocr_figures.py'], 'Py47': ['python3', '-I', 'scripts/47_residual_ink.py'], 'R48': ['Rscript', 'scripts/48_residual_ink.R']}
def failed(k, r): return any(l.startswith('[FAIL]') and (k not in ('OCR', 'Py47', 'R48') or 'FigS1' in l or 'FS1' in l) for l in r.stdout.splitlines())
blue = None
M = [('bar acf3 (log cases) drawn 0.01 too high', patch("vals = np.array([get(key, f'{kind}{k}') for k in lags])",
          "vals = np.array([get(key, f'{kind}{k}') for k in lags]); vals = vals + (0.01 * (lags == 3) if (key, kind) == ('ly_raw', 'acf') else 0)"), ['Py53', 'R54']),
     ('bands drawn for n = 55 instead of 60', patch("band = 1.96 / np.sqrt(n)\nrows = []", "band = 1.96 / np.sqrt(55)\nrows = []"), ['Py53', 'R54']),
     ('panel titles (c) and (e) swapped', patch("SER = [('ly_raw', 'Log cases', ly), ('ly_adj', 'Log cases, month-adjusted', adj(ly)),\n       ('lx_raw', 'Log search interest', lx),",
          "SER = [('ly_raw', 'Log cases', ly), ('ly_adj', 'Log search interest', adj(ly)),\n       ('lx_raw', 'Log cases, month-adjusted', lx),"), ['OCR']),
     ('bar pacf9 (log search) missing', patch("        bars = ax.bar(lags, vals, width=0.6, color=fc.BLUE, linewidth=0, zorder=3)",
          "        bars = ax.bar(lags, np.where((lags == 9) & ((key, kind) == ('lx_raw', 'pacf')), 0, vals), width=0.6, color=fc.BLUE, linewidth=0, zorder=3)"), ['Py53', 'R54']),
     ('stray blue dot (r = 8 px) in panel (a)', None, ['Py53', 'R54', 'Py47', 'R48']),
     ('stray grey dot (r = 8 px) in panel (h)', dot((7600, 7700), 8, (0, 0, 0, 100)), ['Py47', 'R48'])]
def blue_dot():
    A = np.asarray(Image.open(TIF)); yy, xx = np.where((A[..., 0] > 150) & (A[..., 1] < 150) & (A[..., 3] < 60))
    dot((2600, 1650), 8, tuple(int(v) for v in A[yy[len(yy) // 2], xx[len(xx) // 2]]))()
M[4] = (M[4][0], blue_dot, M[4][2])
L2 = []
for name, make, need in M:
    shutil.copy(TIF + '.good', TIF); make()
    res = {k: run(CH[k]) for k in need}; ok = all(failed(k, r) for k, r in res.items()); L2.append(ok)
    print(f"  {name:46s} {'CAUGHT' if ok else 'MISSED'}  " + '; '.join(f"{k}: {'fail' if failed(k, r) else 'pass'}" for k, r in res.items()), flush=True)
shutil.copy(TIF + '.good', TIF); clean = {k: not failed(k, run(c)) for k, c in CH.items()}
print(f'\nLevel 1 caught {sum(L1)}/{len(L1)}; level 2 caught {sum(L2)}/{len(L2)} (each by every check it targets); true image passes all: {all(clean.values())}')
shutil.rmtree(T)
