"""Mutation tests for every Figure S2 check (runs in a temporary copy; the project is not modified).
Level 1: falsified value / drawn files must fail scripts/58_check_figureS2.R (three-program check).
Level 2: falsified FigS2.tif images must fail the image checks they target: 59 (Py digitiser), 60 (R digitiser),
34 (OCR, S2 lines), 47 / 48 (residual ink, S2 lines)."""
import os, shutil, subprocess, tempfile
import numpy as np, pandas as pd
from PIL import Image, ImageDraw
Image.MAX_IMAGE_PIXELS = None
SRC = os.getcwd(); T = tempfile.mkdtemp(prefix='s2mut_')
for d in ('data_raw', 'scripts'): shutil.copytree(os.path.join(SRC, d), os.path.join(T, d))
os.makedirs(os.path.join(T, 'out'), exist_ok=True); os.makedirs(os.path.join(T, 'manuscript/figures'), exist_ok=True)
NODE_MODULES = os.path.join(SRC, 'node_modules') if os.path.isdir(os.path.join(SRC, 'node_modules')) else os.environ.get('NODE_PATH', '')  # npm install
for f in sorted(os.listdir(os.path.join(SRC, 'manuscript/figures'))):        # every figure built (the OCR and residual-ink checks read all)
    if f.endswith('.tif'): shutil.copy(os.path.join(SRC, 'manuscript/figures', f), os.path.join(T, 'manuscript/figures', f))
os.makedirs(os.path.join(T, 'out/final'), exist_ok=True); shutil.copy(os.path.join(SRC, 'out/final/final_results.xlsx'), os.path.join(T, 'out/final'))
os.makedirs(os.path.join(T, 'verification'), exist_ok=True); shutil.copy(os.path.join(SRC, 'verification/comparison_table.csv'), os.path.join(T, 'verification'))
run = lambda c, env=None: subprocess.run(c, cwd=T, capture_output=True, text=True, env={**os.environ, **(env or {})})
for sc in ('00_prepare_inputs.py', '04_dev_regions.py'):
    r = run(['python3', '-I', f'scripts/{sc}'] + (['data_raw/edcd_district_monthly_2019_2023.xlsx'] if sc.startswith('04') else [])); assert r.returncode == 0, r.stderr
assert run(['python3', '-I', 'scripts/56_figureS2.py']).returncode == 0
assert run(['node', 'scripts/57_regional_independent.js'], {'NODE_PATH': NODE_MODULES}).returncode == 0
assert run(['Rscript', 'scripts/58_check_figureS2.R']).returncode == 0, 'true files must pass first'
O = {f: os.path.join(T, 'out', f) for f in ('figS2_values_js.csv', 'figS2_plotted_values.csv', 'figS2_values.csv')}
for p in O.values(): shutil.copy(p, p + '.good')
print('Level 1: files vs scripts/58_check_figureS2.R')
L1 = []
def c1(name, f, edit):
    for p in O.values(): shutil.copy(p + '.good', p)
    d = pd.read_csv(O[f]); edit(d); d.to_csv(O[f], index=False, float_format='%.17g')
    r = run(['Rscript', 'scripts/58_check_figureS2.R']); L1.append(r.returncode != 0)
    print(f"  {name:50s} {'CAUGHT' if L1[-1] else 'MISSED'}")
c1('Python Eastern r off by 1e-9', 'figS2_values.csv', lambda d: d.__setitem__('r', np.where(d.region == 'Eastern', d.r + 1e-9, d.r)))
c1('JavaScript Western search ratio 0.5625 -> 0.57', 'figS2_values_js.csv', lambda d: d.__setitem__('ratio_search', np.where(d.region == 'Western', 0.57, d.ratio_search)))
c1('drawn Central CI upper end +0.005', 'figS2_plotted_values.csv', lambda d: d.__setitem__('hi_drawn', np.where(d.region == 'Central', d.hi_drawn + 0.005, d.hi_drawn)))
c1('row labels Eastern/Western swapped', 'figS2_plotted_values.csv', lambda d: d.__setitem__('ytick_a', ['Central', 'Western', 'Eastern']))
for p in O.values(): shutil.copy(p + '.good', p)

print('Level 2: falsified images')
S = os.path.join(T, 'scripts/56_figureS2.py'); ORIG = open(S).read(); TIF = os.path.join(T, 'manuscript/figures/FigS2.tif'); shutil.copy(TIF, TIF + '.good')
def patch(a, b):
    def make():
        assert a in ORIG; open(S, 'w').write(ORIG.replace(a, b)); r = run(['python3', '-I', 'scripts/56_figureS2.py']); assert r.returncode == 0, r.stderr; open(S, 'w').write(ORIG)
    return make
def dot(xy, r, fill):
    def make():
        im = Image.open(TIF); im.load(); x, y = xy; ImageDraw.Draw(im).ellipse([x - r, y - r, x + r, y + r], fill=fill); im.save(TIF, compression='tiff_lzw', dpi=(1200, 1200))
    return make
CH = {'Py59': ['python3', '-I', 'scripts/59_digitise_figureS2.py'], 'R60': ['Rscript', 'scripts/60_digitise_figureS2.R'],
      'OCR': ['python3', '-I', 'scripts/34_ocr_figures.py'], 'Py47': ['python3', '-I', 'scripts/47_residual_ink.py'], 'R48': ['Rscript', 'scripts/48_residual_ink.R']}
def failed(k, r): return any(l.startswith('[FAIL]') and (k not in ('OCR', 'Py47', 'R48') or 'S2' in l) for l in r.stdout.splitlines())
blue = None
V0 = "V = pd.DataFrame(rows); V.to_csv('out/figS2_values.csv', index=False, float_format='%.17g')"
M = [('Eastern r drawn 0.01 too high', patch(V0, "V = pd.DataFrame(rows); V.loc[V.region == 'Eastern', 'r'] += 0.01; V.to_csv('out/figS2_values.csv', index=False, float_format='%.17g')"), ['Py59', 'R60']),
     ('bootstrap CI drawn instead of the reported CI', patch(V0, "V = pd.DataFrame(rows); V['ci_lo'] = V.boot_lo; V['ci_hi'] = V.boot_hi; V.to_csv('out/figS2_values.csv', index=False, float_format='%.17g')"), ['Py59', 'R60']),
     ('Western cases and search ratios swapped', patch(V0, "V = pd.DataFrame(rows); w = V.region == 'Western'; V.loc[w, ['ratio_cases', 'ratio_search']] = V.loc[w, ['ratio_search', 'ratio_cases']].values; V.to_csv('out/figS2_values.csv', index=False, float_format='%.17g')"), ['Py59', 'R60']),
     ('region labels Eastern/Western swapped', patch("a1.set_yticks(ys); a1.set_yticklabels(REG)", "a1.set_yticks(ys); a1.set_yticklabels(['Central', 'Western', 'Eastern'])"), ['OCR']),
     ('stray blue dot (r = 8 px) in panel (b)', None, ['Py59', 'R60']),
     ('stray grey dot (r = 8 px) in panel (a)', dot((2500, 2000), 8, (0, 0, 0, 100)), ['Py47', 'R48'])]
def blue_dot():
    A = np.asarray(Image.open(TIF)); yy, xx = np.where((A[..., 0] > 150) & (A[..., 1] < 150) & (A[..., 3] < 60))
    dot((7300, 900), 8, tuple(int(v) for v in A[yy[len(yy) // 2], xx[len(xx) // 2]]))()
M[4] = (M[4][0], blue_dot, M[4][2])
L2 = []
for name, make, need in M:
    shutil.copy(TIF + '.good', TIF); make()
    res = {k: run(CH[k]) for k in need}; ok = all(failed(k, r) for k, r in res.items()); L2.append(ok)
    print(f"  {name:46s} {'CAUGHT' if ok else 'MISSED'}  " + '; '.join(f"{k}: {'fail' if failed(k, r) else 'pass'}" for k, r in res.items()), flush=True)
shutil.copy(TIF + '.good', TIF); clean = {k: not failed(k, run(c)) for k, c in CH.items()}
print(f'\nLevel 1 caught {sum(L1)}/{len(L1)}; level 2 caught {sum(L2)}/{len(L2)} (each by every check it targets); true image passes all: {all(clean.values())}')
shutil.rmtree(T)
