"""Mutation tests for every Figure 4 check (runs in a temporary copy; the project is not modified).
Level 1: falsified prediction / drawn-value / text files must fail scripts/42_check_figure4.R (three-program check).
Level 2: falsified Fig4.tif images, re-rendered from a patched copy of scripts/41_figure4.py or edited directly,
must fail the image checks they target: scripts/44 (Python digitiser), 45 (R digitiser), 34 (OCR, Figure 4 lines).
"""
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
SRC = os.getcwd(); T = tempfile.mkdtemp(prefix='f4mut_')
for d in ('data_raw', 'scripts'): shutil.copytree(os.path.join(SRC, d), os.path.join(T, d))
os.makedirs(os.path.join(T, 'out/final'), exist_ok=True); os.makedirs(os.path.join(T, 'manuscript/figures'), exist_ok=True)
for f in sorted(os.listdir(os.path.join(SRC, 'manuscript/figures'))):        # every figure built (the OCR and residual-ink checks read all)
    if f.endswith('.tif'): shutil.copy(os.path.join(SRC, 'manuscript/figures', f), os.path.join(T, 'manuscript/figures', f))
shutil.copy(os.path.join(SRC, 'out/final/final_results.xlsx'), os.path.join(T, 'out/final/'))
shutil.copy(os.path.join(SRC, 'out/15_nowcast_monthly_refit.xlsx'), os.path.join(T, 'out/'))
run = lambda c, env=None: _guarded_run(c, cwd=T, capture_output=True, text=True, env={**os.environ, **(env or {})})
NODE_MODULES = os.path.join(SRC, 'node_modules') if os.path.isdir(os.path.join(SRC, 'node_modules')) else os.environ.get('NODE_PATH', '')  # npm install
JS_ENV = {'NODE_PATH': NODE_MODULES}
assert run(['python3', '-I', 'scripts/00_prepare_inputs.py']).returncode == 0
assert run(['python3', '-I', 'scripts/41_figure4.py']).returncode == 0
assert run(['node', 'scripts/43_nowcast_independent.js'], JS_ENV).returncode == 0
assert run(['Rscript', 'scripts/42_check_figure4.R']).returncode == 0, 'true outputs must pass first'
OUT = {f: os.path.join(T, 'out', f) for f in ('fig4_predictions_python.csv', 'fig4_predictions_js.csv', 'fig4_plotted_values.csv', 'fig4_text.csv')}
for p in OUT.values(): shutil.copy(p, p + '.good')
def restore():
    for p in OUT.values(): shutil.copy(p + '.good', p)

print('Level 1: files vs scripts/42_check_figure4.R')
def c1(name, fname, edit):
    restore(); d = pd.read_csv(OUT[fname]); edit(d); d.to_csv(OUT[fname], index=False, float_format='%.17g')
    r = run(['Rscript', 'scripts/42_check_figure4.R']); ok = r.returncode != 0
    print(f"  {name:52s} {'CAUGHT' if ok else 'MISSED'}  {[l[7:60] for l in r.stdout.splitlines() if l.startswith('[FAIL]')][:2]}"); return ok
def at(d, **kw):
    m = np.ones(len(d), bool)
    for k, v in kw.items(): m &= (d[k] == v).values
    return m
L1 = [c1('JavaScript: one prediction off by 1e-6 (relative)', 'fig4_predictions_js.csv', lambda d: d.__setitem__('pred', np.where(at(d, schedule='monthly', model='M2', month='2023-09'), d.pred * (1 + 1e-6), d.pred))),
      c1('Python: Sep 2023 M2 (season) +0.5%', 'fig4_predictions_python.csv', lambda d: d.__setitem__('pred', np.where(at(d, schedule='season', model='M2', month='2023-09'), d.pred * 1.005, d.pred))),
      c1('drawn observed Aug 2024 changed by 1 case', 'fig4_plotted_values.csv', lambda d: d.__setitem__('y', np.where(at(d, schedule='season', series='obs', month='2024-08'), d.y + 1, d.y))),
      c1('panel (a) text says M2 MAE 1,901', 'fig4_text.csv', lambda d: d.__setitem__('text', d.text.str.replace('M2 1,900', 'M2 1,901')))]
restore()

print('Level 2: falsified images vs scripts/44 (Py), 45 (R), 34 (OCR)')
SCRIPT = os.path.join(T, 'scripts/41_figure4.py'); ORIG = open(SCRIPT).read()
TIF = os.path.join(T, 'manuscript/figures/Fig4.tif'); shutil.copy(TIF, TIF + '.good')
ANCHOR = "mae = {s: {m: np.mean(np.abs(P[s][m] - obs)) for m in ('M0', 'M1', 'M2')} for s in P}"
def patched(code, after_mae=False):
    def make():
        new = ORIG.replace(ANCHOR, ANCHOR + '\n' + code) if after_mae else ORIG.replace(ANCHOR, code + '\n' + ANCHOR)
        assert new != ORIG; open(SCRIPT, 'w').write(new)
        r = run(['python3', '-I', 'scripts/41_figure4.py']); assert r.returncode == 0, r.stderr
        open(SCRIPT, 'w').write(ORIG)
    return make
def paint(shape, fill):
    def make():
        im = Image.open(TIF); im.load(); ImageDraw.Draw(im).ellipse(shape, fill=fill); im.save(TIF, compression='tiff_lzw', dpi=(1200, 1200))
    return make
def checks():
    res = {}
    for k, cmd in [('Py', ['python3', '-I', 'scripts/44_digitise_figure4.py']), ('R', ['Rscript', 'scripts/45_digitise_figure4.R']),
                   ('OCR', ['python3', '-I', 'scripts/34_ocr_figures.py'])]:
        r = run(cmd); res[k] = [l for l in r.stdout.splitlines() if l.startswith('[FAIL]') and (k != 'OCR' or 'F4' in l)]
    return res
i_sep23 = 20; i_mar23 = 14; i_aug24 = 31
M2 = {'M2 (season) Sep 2023 +1% (about 6 px)': (patched(f"P['season']['M2'][{i_sep23}] *= 1.01"), ['Py', 'R']),
      'M1 (season) Mar 2023 +100 cases, where lines overlap': (patched(f"P['season']['M1'][{i_mar23}] += 100"), ['Py', 'R']),
      'M2 (monthly) shifted one month later': (patched("P['monthly']['M2'] = np.r_[P['monthly']['M2'][0], P['monthly']['M2'][:-1]]"), ['Py', 'R']),
      'observed Aug 2024 -2%': (patched(f"obs = obs.copy(); obs[{i_aug24}] *= 0.98"), ['Py', 'R']),
      'M1 and M2 colours swapped': (lambda: (open(SCRIPT, 'w').write(ORIG.replace("('M1', fc.BLUE,", "('M1', fc.ORANGE,").replace("('M2', fc.ORANGE,", "('M2', fc.BLUE,")),
                                            run(['python3', '-I', 'scripts/41_figure4.py']), open(SCRIPT, 'w').write(ORIG)), ['Py', 'R']),
      'text box says M2 MAE 1,800': (patched("mae['season']['M2'] = 1800", after_mae=True), ['OCR']),
      'stray grey dot (r = 10 px)': (paint([5200, 1700, 5220, 1720], (0, 0, 0, 100)), ['Py', 'R']),
      'erase 25 px of the M2 line at its Sep 2023 peak (a)': (lambda: None, ['Py', 'R'])}
# locate the M2 Sep 2023 vertex in panel (a) from the Python digitiser's own calibration for the erase mutant
def erase_m2_peak():
    import importlib.util
    src = open(os.path.join(T, 'scripts/44_digitise_figure4.py')).read().split("    inside = np.zeros(A.shape[:2], bool)")[0]
    g = {}; cwd = os.getcwd(); os.chdir(T)
    try: exec(compile(src.replace("for (y0, y1), sch, lab in zip(panels, ('season', 'monthly'), 'ab'):", "for (y0, y1), sch, lab in list(zip(panels, ('season', 'monthly'), 'ab'))[:1]:"), 'x', 'exec'), g)
    finally: os.chdir(cwd)
    x, y = g['P']['M2'][i_sep23]; paint([int(x) - 25, int(y) - 25, int(x) + 25, int(y) + 25], (0, 0, 0, 0))()
M2['erase 25 px of the M2 line at its Sep 2023 peak (a)'] = (erase_m2_peak, ['Py', 'R'])
L2 = []
for name, (make, need) in M2.items():
    shutil.copy(TIF + '.good', TIF); make(); res = checks(); ok = all(res[k] for k in need); L2.append(ok)
    print(f"  {name:52s} {'CAUGHT' if ok else 'MISSED'}  " + '; '.join(f'{k}: {len(v)}' for k, v in res.items()), flush=True)
shutil.copy(TIF + '.good', TIF); res = checks()
print(f'\nLevel 1 caught {sum(L1)}/{len(L1)}; level 2 caught {sum(L2)}/{len(L2)} (each by every check it targets); '
      f'true image passes all image checks: {not any(res.values())}')
shutil.rmtree(T)
