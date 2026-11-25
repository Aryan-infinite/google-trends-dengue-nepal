"""Image-level check of the DELIVERED Supplementary Figure S3 TIFF (Python; independent of the plotting code).
Panel (a): x scale from the gridlines at 0, 0.25, ..., 1.00; each region's dot (centroid of black blob) and CI ends
(horizontal extent of the black CI line in that row). Panel (b): log x scale from gridlines at 0.1, 0.3, 3, 10 and the
reference line at 1; case dots (black) and search triangles (blue) located by centroid. No black or blue ink elsewhere.
Expected values recomputed here from the regional files (cases from out/edcd_devregion_monthly_2019_2023.csv).
"""
import os, sys
import numpy as np, pandas as pd
from PIL import Image, ImageCms
from scipy import ndimage as ndi
Image.MAX_IMAGE_PIXELS = None
ICC = os.environ.get('CMYK_ICC', '/usr/share/texlive/texmf-dist/tex/generic/colorprofiles/FOGRA39L_coated.icc')
fails = []
def check(n, ok, d=''):
    print(f"[{'PASS' if ok else 'FAIL'}] {n}" + (f'  ({d})' if d else '')); ok or fails.append(n)
def runs(idx):
    if len(idx) == 0: return []
    br = np.where(np.diff(idx) > 1)[0]; return list(zip(np.r_[idx[0], idx[br + 1]], np.r_[idx[br], idx[-1]]))
REG = ['Central', 'Eastern', 'Western']
rc = pd.read_csv('out/edcd_devregion_monthly_2019_2023.csv', index_col=0, parse_dates=True)['2020':'2023']
mon = np.tile(np.arange(1, 13), 4); adj = lambda v: v - np.array([v[mon == m].mean() for m in mon])
lagr = lambda v, k: np.corrcoef(v[:-k], v[k:])[0, 1]
E = {}
for r_ in REG:
    d = pd.read_csv(f'data_raw/gt_dengue_region_{r_.lower()}.csv', skiprows=2); d = d[d.iloc[:, 0].str[:4].isin(['2020', '2021', '2022', '2023'])]
    s = d.iloc[:, 1].astype(str).replace('<1', '0.5').astype(float).values; c = rc[r_].values.astype(float)
    a, b = adj(np.log(s + .5)), adj(np.log(c + 1)); r = np.corrcoef(a, b)[0, 1]
    ne = max(48 / (1 + 2 * sum(lagr(a, k) * lagr(b, k) for k in range(1, 13))), 4); z, se = np.arctanh(r), 1 / np.sqrt(ne - 3)
    E[r_] = dict(r=r, lo=np.tanh(z - 1.96 * se), hi=np.tanh(z + 1.96 * se), rc=c[36:].max() / c[24:36].max(), rs=s[36:].max() / s[24:36].max())
im = Image.open('manuscript/figures/FigS3.tif'); assert im.mode == 'CMYK'
A = np.asarray(ImageCms.applyTransform(im, ImageCms.buildTransform(ICC, ImageCms.createProfile('sRGB'), 'CMYK', 'RGB'))).astype(np.int16)
near = lambda rgb, t: np.abs(A - rgb).sum(2) <= t
ink2, grid, blk, blu = near([76, 75, 73], 30), near([216, 215, 210], 18), near([0, 0, 7], 60), near([57, 115, 185], 60)
cs = ink2.sum(0); bands = [r for r in runs(np.where(cs > 0.5 * cs.max())[0]) if r[1] - r[0] + 1 >= 9]
sp = bands[:2]                                           # the third full-height line is panel (b)'s reference line at 1
check('two left spines + one reference line', len(bands) == 3, str(bands))
res = {}
for p_, (s0, s1) in enumerate(sp):
    yr = [r for r in runs(np.where(ink2[:, s0:s1 + 1].all(1))[0]) if r[1] - r[0] > 300][0]; y0, y1 = yr
    bt = runs(np.where(ink2[y0:y1 + 30, s1 + 20:s1 + 1200].sum(1) > 1000)[0]); yb = int(y0 + np.mean(bt[-1]))
    right = sp[1][0] - 800 if p_ == 0 else A.shape[1] - 5
    gx = [s1 + 5 + np.mean(r) for r in runs(np.where(grid[y0:yb - 10, s1 + 5:right].sum(0) > 0.5 * (yb - 10 - y0))[0])]
    inside = np.zeros(A.shape[:2], bool); inside[y0:yb - 8, s1 + 8:right] = True
    if p_ == 0:
        xs = np.array([(s0 + s1) / 2] + gx); xv = np.array([0, 0.25, 0.5, 0.75, 1.0])
    else:
        ref = [s1 + 5 + np.mean(r) for r in runs(np.where(ink2[y0 + 10:yb - 10, s1 + 5:right].sum(0) > 0.8 * (yb - 20 - y0))[0])]
        check('(b) one vertical reference line', len(ref) == 1, str(ref))
        xs = np.array(sorted(gx + ref)); xv = np.log10([0.1, 0.3, 1, 3, 10])
    ok = len(xs) == len(xv); co = np.polyfit(xs, xv, 1) if ok else None
    rsd = np.max(np.abs(np.polyval(co, xs) - xv)) / abs(co[0]) if ok else np.inf
    check(f'({"ab"[p_]}) x scale from {len(xv)} lines within 2 px', rsd < 2, f'{len(xs)} lines; residual {rsd:.2f} px')
    v2x = (lambda v: (v - co[1]) / co[0]) if p_ == 0 else (lambda v: (np.log10(v) - co[1]) / co[0])
    if p_ == 0:
        m = blk & inside; lbl, n = ndi.label(m, structure=np.ones((3, 3)))
        check('(a) three black objects (dot + CI line per region)', n == 3, f'{n}')
        objs = sorted(ndi.find_objects(lbl), key=lambda sl: sl[0].start)        # top to bottom
        err = []
        for sl, r_ in zip(objs, REG):
            sub = m[sl]; cols = sub.sum(0); dotc = np.where(cols >= cols.max() - 2)[0]   # dot = tallest columns
            xd = sl[1].start + dotc.mean(); xl, xr = sl[1].start, sl[1].stop - 1
            err += [abs(xd - v2x(E[r_]['r'])), abs(xl - v2x(E[r_]['lo'])), abs(xr - v2x(E[r_]['hi']))]
            res[r_] = dict(r=float(np.polyval(co, xd)), lo=float(np.polyval(co, xl)), hi=float(np.polyval(co, xr)))
        check('(a) dots and CI ends within 2 px of recomputed r and Fisher limits', max(err) <= 2,
              f'max {max(err):.2f} px; read: ' + '; '.join(f"{k} {v['r']:.3f} ({v['lo']:.3f}-{v['hi']:.3f})" for k, v in res.items()))
    else:
        err = []
        for col_, key, name in ((blk, 'rc', 'cases'), (blu, 'rs', 'search')):
            m = col_ & inside; lbl, n = ndi.label(m, structure=np.ones((3, 3)))
            check(f'(b) three {name} markers', n == 3, f'{n}')
            cen = sorted(ndi.center_of_mass(m, lbl, range(1, n + 1)))                 # top to bottom
            for (cy, cx), r_ in zip(cen, REG): err.append(abs(cx - v2x(E[r_][key]))); res[f'{r_}_{key}'] = float(10 ** np.polyval(co, cx))
        check('(b) all six markers within 2 px of recomputed peak ratios', max(err) <= 2,
              f"max {max(err):.2f} px; read: " + ', '.join(f'{k} {v:.3f}' for k, v in res.items() if '_' in k))
print(f'RESULT: {len(fails)} failed check(s)'); sys.exit(1 if fails else 0)
