"""Mutation test for the residual-ink layer (scripts/47 Python, 48 R; for Figure 2 also the geometric rule in 31/32).
Stray marks are painted into copies of the delivered TIFFs in empty parts of a panel; each must be caught by the
Python AND the R implementation. Runs in a temporary copy; the project is not modified."""
import os, shutil, subprocess, tempfile

def _guarded_run(*a, **k):
    """A check killed by the system (e.g. out of memory) prints nothing; stop instead of counting it as caught or missed."""
    r = subprocess.run(*a, **k)
    if r.returncode < 0 or r.returncode == 137:
        raise SystemExit(f'ABORTED: a check was killed (exit {r.returncode}): {a[0] if a else k.get("args")}; rerun with nothing else running')
    return r
import numpy as np
from PIL import Image, ImageDraw
Image.MAX_IMAGE_PIXELS = None
SRC = os.getcwd(); T = tempfile.mkdtemp(prefix='resmut_')
for d in ('data_raw', 'scripts'): shutil.copytree(os.path.join(SRC, d), os.path.join(T, d))
os.makedirs(os.path.join(T, 'manuscript/figures'), exist_ok=True); os.makedirs(os.path.join(T, 'out'), exist_ok=True)
for f in sorted(os.listdir(os.path.join(SRC, 'manuscript/figures'))):        # every figure built (the OCR and residual-ink checks read all)
    if f.endswith('.tif'): shutil.copy(os.path.join(SRC, 'manuscript/figures', f), os.path.join(T, 'manuscript/figures', f))
FIG = {n: os.path.join(T, f'manuscript/figures/{n}.tif') for n in ('Fig1', 'Fig2', 'Fig3', 'Fig4')}
def reset():
    for n, p in FIG.items(): shutil.copy(os.path.join(SRC, f'manuscript/figures/{n}.tif'), p)
run = lambda c: _guarded_run(c, cwd=T, capture_output=True, text=True)
def sample(fig, xy):
    im = Image.open(FIG[fig]); return im.getpixel(xy)
reset()
# CMYK colours sampled from the true images (gridline, axis spine, black ink) plus two colours the figures never use
grid1 = sample('Fig1', (3000, 2823)); ink2_1 = sample('Fig1', (963, 3500))
_A = np.asarray(Image.open(FIG['Fig2'])); _y, _x = np.where(_A[..., 3] > 240)        # K channel ~ 100%: genuine black ink
black2 = tuple(int(v) for v in _A[_y[0], _x[0]]); assert black2[3] > 240, black2
grid4 = sample('Fig4', (3000, 4331))
UNKNOWN_GREY, PALE_PINK = (0, 0, 0, 100), (0, 40, 20, 0)
def dot(fig, xy, r, fill):
    def make():
        im = Image.open(FIG[fig]); im.load(); x, y = xy; ImageDraw.Draw(im).ellipse([x - r, y - r, x + r, y + r], fill=fill)
        im.save(FIG[fig], compression='tiff_lzw', dpi=(1200, 1200))
    return make
CHECKS = {'Py47': ['python3', '-I', 'scripts/47_residual_ink.py'], 'R48': ['Rscript', 'scripts/48_residual_ink.R'],
          'Py31': ['python3', '-I', 'scripts/31_digitise_figures.py'], 'R32': ['Rscript', 'scripts/32_digitise_figures.R'],
          'R45': ['Rscript', 'scripts/45_digitise_figure4.R']}
M = [('Fig1(b) grey dot, colour not in the figure (r=8)', dot('Fig1', (3800, 3050), 8, UNKNOWN_GREY), ['Py47', 'R48']),
     ('Fig1(c) gridline-grey dot off any gridline (r=8)', dot('Fig1', (3800, 5100), 8, grid1), ['Py47', 'R48']),
     ('Fig1(b) dark-grey dot (r=6)', dot('Fig1', (5000, 3450), 6, ink2_1), ['Py47', 'R48']),
     ('Fig2(a) grey dot far from markers (r=8)', dot('Fig2', (2600, 2900), 8, UNKNOWN_GREY), ['Py31', 'R32']),
     ('Fig2(b) dark-grey dot: 17th shape (r=6)', dot('Fig2', (7400, 2800), 6, ink2_1), ['Py47', 'R48']),
     ('Fig2(a) black dot outside the label area (r=5)', dot('Fig2', (3300, 2500), 5, black2), ['Py47', 'R48']),
     ('Fig2(a) darker grey dot (K 60%) far from markers (r=8)', dot('Fig2', (2600, 2900), 8, (0, 0, 0, 153)), ['Py31', 'R32']),
     ('Fig2(b) grey dot far from markers (r=8)', dot('Fig2', (7300, 2900), 8, UNKNOWN_GREY), ['Py31', 'R32']),
     ('Fig3 grey dot, colour not in the figure (r=8)', dot('Fig3', (4600, 700), 8, UNKNOWN_GREY), ['Py47', 'R48']),
     ('Fig3 dark-grey dot away from zero line and dashes (r=6)', dot('Fig3', (1500, 1800), 6, ink2_1), ['Py47', 'R48']),
     ('Fig4(a) the grey dot missed earlier by script 44 (r=10)', dot('Fig4', (5210, 1710), 10, UNKNOWN_GREY), ['Py47', ('R48', 'R45')]),
     ('Fig4(b) gridline-grey dot off any gridline (r=8)', dot('Fig4', (2500, 4300), 8, grid4), ['Py47', 'R48']),
     ('Fig4(a) pale pink dot (r=6)', dot('Fig4', (6000, 1300), 6, PALE_PINK), ['Py47', 'R48'])]
import sys
ONLY = sys.argv[1] if len(sys.argv) > 1 else None
caught = []
for name, make, need in M:
    if ONLY and not name.startswith(ONLY): continue
    reset(); make()
    groups = [g if isinstance(g, tuple) else (g,) for g in need]          # a tuple = caught by any one of these (same language)
    res = {k: run(CHECKS[k]) for g in groups for k in g}
    ok = all(any(res[k].returncode != 0 for k in g) for g in groups); caught.append(ok)
    print(f"{name:58s} {'CAUGHT' if ok else 'MISSED'}  " + '; '.join(f"{k}: {sum(l.startswith('[FAIL]') for l in r.stdout.splitlines())} fail" for k, r in res.items()), flush=True)
reset(); clean = {k: run(c).returncode == 0 for k, c in CHECKS.items()}
print(f'\ncaught {sum(caught)}/{len(caught)} (each by both Python and R); true images pass all four checks: {all(clean.values())}')
shutil.rmtree(T)
