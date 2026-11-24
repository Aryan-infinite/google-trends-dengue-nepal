"""Image-level check of the DELIVERED Figure 3 TIFF (Python; independent of the plotting code).
y scale fitted from the gridlines and the zero line found in the pixels; each bar's value read from its top edge;
threshold segments read from the rows of the dashed lines above and below each bar; no blue ink outside the
seven expected bars. Expected values recomputed here from the raw files (statsmodels AutoReg, AR order 2).
"""
import os, sys
import numpy as np, pandas as pd, scipy.stats as st
from PIL import Image, ImageCms
from statsmodels.tsa.ar_model import AutoReg
from statsmodels.stats.multitest import multipletests
Image.MAX_IMAGE_PIXELS = None
ICC = os.environ.get('CMYK_ICC', '/usr/share/texlive/texmf-dist/tex/generic/colorprofiles/FOGRA39L_coated.icc')
fails = []
def check(n, ok, d=''):
    print(f"[{'PASS' if ok else 'FAIL'}] {n}" + (f'  ({d})' if d else '')); ok or fails.append(n)
def runs(idx):
    if len(idx) == 0: return []
    br = np.where(np.diff(idx) > 1)[0]; return list(zip(np.r_[idx[0], idx[br + 1]], np.r_[idx[br], idx[-1]]))

# expected values from raw files
import openpyxl
rows = [[c.value for c in r] for r in openpyxl.load_workbook(os.environ.get('PREPRINT_XLSX', 'data_raw/paper_finalized_data_RP.xlsx')).active.iter_rows()]
h = next(i for i, r in enumerate(rows) if 'Confirmed Cases' in r); ci = rows[h].index('Confirmed Cases')
cases = np.array([float(r[ci]) for r in rows[h + 1:] if r[ci] is not None])[12:]
g = lambda f: pd.read_csv(f, skiprows=2).iloc[:, 1].astype(str).replace('<1', '0.5').astype(float).values
s = np.mean([g('data_raw/gt_dengue_term.csv'), g('data_raw/gt_dengue_term_day2_2026-11-22.csv'), g('data_raw/gt_dengue_term_day3_2026-11-23.csv')], 0)[12:]
mon = np.tile(np.arange(1, 13), 5); dm = lambda v: v - np.array([v[mon == k].mean() for k in mon])
a, b = pd.Series(dm(np.log(s + .5))), pd.Series(dm(np.log(cases + 1)))
f = AutoReg(a.values, lags=2, trend='c').fit(); cst, ph = f.params[0], f.params[1:]
flt = lambda v: (v - cst - ph[0] * v.shift(1) - ph[1] * v.shift(2)).iloc[2:]
ex, ey = flt(a), flt(b); R, N, Pv = [], [], []
for k in range(-3, 4):
    q = ex.shift(k); m = q.notna() & ey.notna(); r_, p_ = st.pearsonr(q[m], ey[m]); R.append(r_); N.append(int(m.sum())); Pv.append(p_)
R, N = np.array(R), np.array(N); crit = np.array([st.t.ppf(.975, n - 2) / np.sqrt(n - 2 + st.t.ppf(.975, n - 2) ** 2) for n in N])

im = Image.open('manuscript/figures/Fig3.tif'); assert im.mode == 'CMYK'
A = np.asarray(ImageCms.applyTransform(im, ImageCms.buildTransform(ICC, ImageCms.createProfile('sRGB'), 'CMYK', 'RGB'))).astype(np.int16)
near = lambda rgb, t: np.abs(A - rgb).sum(2) <= t
ink2, grid, blue = near([76, 75, 73], 30), near([216, 215, 210], 18), near([57, 115, 185], 60)
cs = ink2.sum(0); sp = [r for r in runs(np.where(cs > 0.5 * cs.max())[0]) if r[1] - r[0] + 1 >= 9]
check('one left spine', len(sp) == 1, str(sp))
x0 = sp[0][1] + 1
yr = [r for r in runs(np.where(ink2[:, sp[0][0]:sp[0][1] + 1].all(1))[0]) if r[1] - r[0] > 500][0]; y0, y1 = yr
W = A.shape[1]
hg = [int(round(np.mean(r))) + y0 for r in runs(np.where(grid[y0:y1, x0 + 10:W - 20].sum(1) > 0.5 * (W - 30 - x0))[0])]
zr = [r for r in runs(np.where(ink2[y0:y1, x0 + 10:W - 20].sum(1) > 0.8 * (W - 30 - x0))[0]) if y0 + r[1] < y1 - 30]
check('one zero line', len(zr) == 1, str(zr))
zy = y0 + np.mean(zr[0])
py_ = np.array(sorted(hg + [zy]), float); vals = np.array([0.6, 0.4, 0.2, 0.0, -0.2])
check('4 gridlines + zero line found', len(py_) == 5, str(py_))
coef = np.polyfit(py_, vals, 1); res = np.max(np.abs(np.polyval(coef, py_) - vals)) / abs(coef[0])
check('y scale fits within 2 px', res < 2, f'{res:.2f} px; 1 px = {abs(coef[0]):.5f} r')
v2p = lambda v: (v - coef[1]) / coef[0]
inside = np.zeros(blue.shape, bool); inside[y0:y1, x0 + 5:W - 5] = True; bm = blue & inside
bars = [r for r in runs(np.where(bm.sum(0) > 0)[0]) if r[1] - r[0] > 40]
check('seven bars', len(bars) == 7, f'{len(bars)}')
cx = np.array([(u + v) / 2 for u, v in bars]); sp_ = np.diff(cx)
check('bars equally spaced (lag positions)', np.max(np.abs(sp_ - sp_.mean())) < 2, f'spacing {sp_.mean():.1f} px, max deviation {np.max(np.abs(sp_ - sp_.mean())):.2f} px')
unit = sp_.mean()                                    # px per lag unit
err, got = [], []
for (u, v), r_exp in zip(bars, R):
    c = (u + v) // 2; col = np.where(bm[:, c - 5:c + 6].any(1))[0]
    top = col.min() if r_exp >= 0 else col.max() + 1
    val = np.polyval(coef, top - 0.5 if r_exp >= 0 else top - 0.5); got.append(val); err.append(abs(top - 0.5 - v2p(r_exp)))
check('each bar top edge within 2 px of recomputed r', max(err) <= 2, f'max {max(err):.2f} px = {max(err) * abs(coef[0]):.4f} in r; read: ' + ', '.join(f'{x:.3f}' for x in got))
check('bar widths = 0.6 lag units', all(abs((v - u + 1) / unit - 0.6) < 0.02 for u, v in bars), ', '.join(f'{(v - u + 1) / unit:.3f}' for u, v in bars))
# dashed threshold segments: INK2 rows within each bar's +/-0.4 span, away from the zero line
d_err = []
for c, cr in zip(cx, crit):
    xs = slice(int(c - 0.38 * unit), int(c + 0.38 * unit)); m = ink2[y0:y1, xs].copy()
    m[int(zy - y0) - 12:int(zy - y0) + 13] = False; m[(y1 - y0) - 25:] = False      # zero line and bottom spine
    rr = [y0 + np.mean(r) for r in runs(np.where(m.sum(1) > 0.3 * m.shape[1])[0])]
    if len(rr) != 2: d_err.append(np.inf); continue
    d_err += [abs(rr[0] - v2p(cr)), abs(rr[1] - v2p(-cr))]
check('dashed thresholds at +/- r_crit for each lag within 2 px', max(d_err) <= 2, f'max {max(d_err):.2f} px')
# no extra blue ink outside the seven expected bars (+3 px)
exp = np.zeros(blue.shape, bool)
for c, r_exp in zip(cx, R):
    ya, yb = sorted([v2p(0), v2p(r_exp)]); exp[int(ya) - 3:int(yb) + 4, int(c - 0.3 * unit) - 3:int(c + 0.3 * unit) + 4] = True
check('no blue ink outside the expected bars', not (bm & ~exp).any(), f'{int((bm & ~exp).sum())} stray px')
print(f'RESULT: {len(fails)} failed check(s)'); sys.exit(1 if fails else 0)
