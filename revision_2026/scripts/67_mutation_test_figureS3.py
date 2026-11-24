"""Mutation tests for every Figure S3 check (runs in a temporary copy; the project is not modified).
Level 1: falsified value / drawn files must fail scripts/63_check_figureS3.R (three-program check).
Level 2: falsified FigS3.tif images must fail the image checks they target: 65 (Py digitiser), 66 (R digitiser),
34 (OCR, S3 lines), 47 / 48 (residual ink, S3 lines)."""
import os, shutil, subprocess, tempfile
import numpy as np, pandas as pd
from PIL import Image, ImageDraw
Image.MAX_IMAGE_PIXELS = None
SRC = os.getcwd(); T = tempfile.mkdtemp(prefix='s3mut_')
for d in ('data_raw', 'scripts', 'governance'): shutil.copytree(os.path.join(SRC, d), os.path.join(T, d))
os.makedirs(os.path.join(T, 'out'), exist_ok=True); os.makedirs(os.path.join(T, 'manuscript/figures'), exist_ok=True)
NODE_MODULES = os.path.join(SRC, 'node_modules') if os.path.isdir(os.path.join(SRC, 'node_modules')) else os.environ.get('NODE_PATH', '')  # npm install
for f in sorted(os.listdir(os.path.join(SRC, 'manuscript/figures'))):        # every figure built (the OCR and residual-ink checks read all)
    if f.endswith('.tif'): shutil.copy(os.path.join(SRC, 'manuscript/figures', f), os.path.join(T, 'manuscript/figures', f))
for f in ('06_2024_reporting_lag.csv', '23_announcement_timing.json', '24_announcement_timing.json'): shutil.copy(os.path.join(SRC, 'out', f), os.path.join(T, 'out'))
run = lambda c, env=None: subprocess.run(c, cwd=T, capture_output=True, text=True, env={**os.environ, **(env or {})})
assert run(['python3', '-I', 'scripts/62_figureS3.py']).returncode == 0
assert run(['node', 'scripts/64_s3_independent.js'], {'NODE_PATH': NODE_MODULES}).returncode == 0
assert run(['Rscript', 'scripts/63_check_figureS3.R']).returncode == 0, 'true files must pass first'
O = {f: os.path.join(T, 'out', f) for f in ('figS3_values.csv', 'figS3_values_js.csv', 'figS3_plotted_values.csv', 'figS3_events_drawn.csv', 'figS3_events_js.csv', 'figS3_text.csv')}
for p in O.values(): shutil.copy(p, p + '.good')
print('Level 1: files vs scripts/63_check_figureS3.R')
L1 = []
def c1(name, f, edit):
    for p in O.values(): shutil.copy(p + '.good', p)
    d = pd.read_csv(O[f]); edit(d); d.to_csv(O[f], index=False, float_format='%.17g')
    r = run(['Rscript', 'scripts/63_check_figureS3.R']); L1.append(r.returncode != 0)
    print(f"  {name:50s} {'CAUGHT' if L1[-1] else 'MISSED'}")
c1('Python search 2023-05 relative value + 1e-9', 'figS3_values.csv', lambda d: d.__setitem__('search_rel', np.where((d.year == 2023) & (d.month == 5), d.search_rel + 1e-9, d.search_rel)))
c1('JavaScript event 2024-10-15 x + 0.01 month', 'figS3_events_js.csv', lambda d: d.__setitem__('x', np.where(d.date == '2024-10-15', d.x + 0.01, d.x)))
c1('drawn case line 2022 Sep - 0.002', 'figS3_plotted_values.csv', lambda d: d.__setitem__('y', np.where((d.series == 'cases') & (d.year == 2022) & (d.x == 8.5), d.y - 0.002, d.y)))
c1('2022-05-06 drawn as a national alert', 'figS3_events_drawn.csv', lambda d: (d.__setitem__('category', np.where(d.date == '2022-05-06', 'alert', d.category)), d.__setitem__('marker', np.where(d.date == '2022-05-06', 'v', d.marker))))
c1('note "\u2264 3" printed as "< 3"', 'figS3_text.csv', lambda d: d.__setitem__('text', d.text.str.replace('\u2264 3', '< 3', regex=False)))
for p in O.values(): shutil.copy(p + '.good', p)

print('Level 2: falsified images')
S = os.path.join(T, 'scripts/62_figureS3.py'); ORIG = open(S).read(); TIF = os.path.join(T, 'manuscript/figures/FigS3.tif'); shutil.copy(TIF, TIF + '.good')
def patch(a, b):
    def make():
        assert a in ORIG; open(S, 'w').write(ORIG.replace(a, b)); r = run(['python3', '-I', 'scripts/62_figureS3.py']); assert r.returncode == 0, r.stderr; open(S, 'w').write(ORIG)
    return make
def dot(xy, r, fill):
    def make():
        im = Image.open(TIF); im.load(); x, y = xy; ImageDraw.Draw(im).ellipse([x - r, y - r, x + r, y + r], fill=fill); im.save(TIF, compression='tiff_lzw', dpi=(1200, 1200))
    return make
CH = {'Py65': ['python3', '-I', 'scripts/65_digitise_figureS3.py'], 'R66': ['Rscript', 'scripts/66_digitise_figureS3.R'],
      'OCR': ['python3', '-I', 'scripts/34_ocr_figures.py'], 'Py47': ['python3', '-I', 'scripts/47_residual_ink.py'], 'R48': ['Rscript', 'scripts/48_residual_ink.R']}
def failed(k, r): return any(l.startswith('[FAIL]') and (k not in ('OCR', 'Py47', 'R48') or 'S3' in l) for l in r.stdout.splitlines())
V0 = "V = pd.DataFrame(V); V.to_csv('out/figS3_values.csv', index=False, float_format='%.17g')"
E0 = "assert len(ev) == 15, len(ev)"
M = [('2023 August cases drawn 0.03 too low', patch(V0, V0 + "; V.loc[(V.year == 2023) & (V.month == 8), 'cases_rel'] -= 0.03"), ['Py65', 'R66']),
     ('2023-05-10 alert drawn one week late', patch(E0, E0 + "; ev = [(('2023-05-17' if d == '2023-05-10' else d), c) for d, c in ev]"), ['Py65', 'R66']),
     ('2022-05-06 count drawn as a triangle (alert)', patch(E0, E0 + "; ev = [(d, ('alert' if d == '2022-05-06' else c)) for d, c in ev]"), ['Py65', 'R66']),
     ('extra marker in 2021 (no alert was found)', patch(E0, E0 + "; ev = ev + [('2021-03-15', 'alert')]"), ['Py65', 'R66']),
     ('month labels M and A swapped', patch("ax.set_xticklabels(list('JFMAMJJASOND'))", "ax.set_xticklabels(list('JFAMMJJASOND'))"), ['OCR']),
     ('stray blue dot (r = 8 px) in 2024 panel', None, ['Py65', 'R66']),
     ('stray grey dot (r = 8 px) in 2021 panel', dot((2000, 3900), 8, (0, 0, 0, 100)), ['Py47', 'R48'])]
def blue_dot():
    A = np.asarray(Image.open(TIF)); yy, xx = np.where((A[..., 0] > 150) & (A[..., 1] < 150) & (A[..., 3] < 60))
    dot((6500, 5500), 8, tuple(int(v) for v in A[yy[len(yy) // 2], xx[len(xx) // 2]]))()
M[5] = (M[5][0], blue_dot, M[5][2])
L2 = []
for name, make, need in M:
    shutil.copy(TIF + '.good', TIF); make()
    res = {k: run(CH[k]) for k in need}; ok = all(failed(k, r) for k, r in res.items()); L2.append(ok)
    print(f"  {name:46s} {'CAUGHT' if ok else 'MISSED'}  " + '; '.join(f"{k}: {'fail' if failed(k, r) else 'pass'}" for k, r in res.items()), flush=True)
shutil.copy(TIF + '.good', TIF); clean = {k: not failed(k, run(c)) for k, c in CH.items()}
print(f'\nLevel 1 caught {sum(L1)}/{len(L1)}; level 2 caught {sum(L2)}/{len(L2)} (each by every check it targets); true image passes all: {all(clean.values())}')
shutil.rmtree(T)
