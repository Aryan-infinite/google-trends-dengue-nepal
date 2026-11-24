"""Image-level check of the DELIVERED Supplementary Figure S1 TIFF (Python; independent of the plotting code).
Eight panels found from the axis lines; y scale per panel from gridlines and the zero line; each bar's value read from
its outer edge (top for positive, bottom for negative); dashed bands read from their rows; no blue ink outside the
expected bars. Expected values recomputed here from the raw files (numpy). Bars under 3 px tall cannot be read and are
only required to have no ink beyond 3 px from the zero line.
"""
import os, sys
import numpy as np, pandas as pd
from PIL import Image, ImageCms
Image.MAX_IMAGE_PIXELS = None
ICC = os.environ.get('CMYK_ICC', '/usr/share/texlive/texmf-dist/tex/generic/colorprofiles/FOGRA39L_coated.icc')
fails = []
def check(n, ok, d=''):
    print(f"[{'PASS' if ok else 'FAIL'}] {n}" + (f'  ({d})' if d else '')); ok or fails.append(n)
def runs(idx):
    if len(idx) == 0: return []
    br = np.where(np.diff(idx) > 1)[0]; return list(zip(np.r_[idx[0], idx[br + 1]], np.r_[idx[br], idx[-1]]))
import openpyxl
rows = [[c.value for c in r] for r in openpyxl.load_workbook(os.environ.get('PREPRINT_XLSX', 'data_raw/paper_finalized_data_RP.xlsx')).active.iter_rows()]
h = next(i for i, r in enumerate(rows) if 'Confirmed Cases' in r); ci = rows[h].index('Confirmed Cases')
c = np.array([float(r[ci]) for r in rows[h + 1:] if r[ci] is not None])[12:]
g = lambda f: pd.read_csv(f, skiprows=2).iloc[:, 1].astype(str).replace('<1', '0.5').astype(float).values
s = np.mean([g('data_raw/gt_dengue_term.csv'), g('data_raw/gt_dengue_term_day2_2026-11-22.csv'), g('data_raw/gt_dengue_term_day3_2026-11-23.csv')], 0)[12:]
mon = np.tile(np.arange(1, 13), 5); n = 60; adj = lambda v: v - np.array([v[mon == m].mean() for m in mon])
ly, lx = np.log(c + 1), np.log(s + 0.5)
def acf(v): d = v - v.mean(); return np.array([d[:-k] @ d[k:] for k in range(1, 21)]) / (d @ d)
def pacf(v): return np.array([np.linalg.lstsq(np.column_stack([np.ones(n - k)] + [v[k - j:n - j] for j in range(1, k + 1)]), v[k:], rcond=None)[0][k] for k in range(1, 21)])
ser = [ly, adj(ly), lx, adj(lx)]
EXP = [[acf(v), pacf(v)] for v in ser]; band = 1.96 / np.sqrt(n)

im = Image.open('manuscript/figures/FigS1.tif'); assert im.mode == 'CMYK'
A = np.asarray(ImageCms.applyTransform(im, ImageCms.buildTransform(ICC, ImageCms.createProfile('sRGB'), 'CMYK', 'RGB'))).astype(np.int16)
near = lambda rgb, t: np.abs(A - rgb).sum(2) <= t
ink2, grid, blue = near([76, 75, 73], 30), near([216, 215, 210], 18), near([57, 115, 185], 60)
cs = ink2.sum(0); sp = [r for r in runs(np.where(cs > 0.4 * cs.max())[0]) if r[1] - r[0] + 1 >= 9]
check('two columns of panels (two left-spine bands)', len(sp) == 2, str(sp))
tot = {'read': 0, 'tiny': 0, 'maxerr': 0.0, 'band': 0.0}
for ci_, (s0, s1) in enumerate(sp):
    prs = [r for r in runs(np.where(ink2[:, s0:s1 + 1].all(1))[0]) if r[1] - r[0] > 500]
    check(f'column {ci_ + 1}: four panels', len(prs) == 4, str(prs))
    for ri, (y0, y1) in enumerate(prs):
        lab = 'abcdefgh'[2 * ri + ci_]
        bt = runs(np.where(ink2[y0:y1 + 40, s1 + 20:s1 + 1500].sum(1) > 1200)[0]); yb = y0 + np.mean(bt[-1])
        xr = s1 + 1 + np.where(ink2[int(yb) - 3:int(yb) + 4, s1 + 1:].any(0))[0]; xr = xr[np.r_[np.where(np.diff(xr) > 60)[0], len(xr) - 1][0]]
        sub_g = grid[y0 - 12:int(yb), s1 + 20:xr - 10]
        hg = [y0 - 12 + np.mean(r) for r in runs(np.where(sub_g.sum(1) > 0.5 * sub_g.shape[1])[0])]
        # three dark-grey row bands per panel: upper dashed band, zero line (partly under the bars), lower dashed band
        tri = [y0 + np.mean(r) for r in runs(np.where(ink2[y0:int(yb) - 30, s1 + 20:xr - 10].sum(1) > 0.25 * (xr - 30 - s1))[0])]
        zr = [tri[1]] if len(tri) == 3 else []
        check(f'({lab}) three dark-grey row bands (dash, zero line, dash)', len(tri) == 3, f'{len(tri)}')
        pts = np.array(sorted(hg) + zr + [yb]); vals = np.array([1.0, 0.5, -0.5][:len(hg)] + [0.0] * len(zr) + [-1.0])
        if len(hg) == 2: vals = np.array([0.5, -0.5, 0.0, -1.0])
        ok_cal = len(zr) == 1 and len(hg) in (2, 3)
        co = np.polyfit(pts, vals, 1) if ok_cal else (np.nan, np.nan); res = np.max(np.abs(np.polyval(co, pts) - vals)) / abs(co[0]) if ok_cal else np.inf
        check(f'({lab}) y scale from gridlines + zero line + axis within 2 px', res < 2, f'{len(hg)} gridlines; residual {res:.2f} px')
        v2p = lambda v: (v - co[1]) / co[0]; zy = v2p(0.0)
        bm = blue.copy(); bm[:y0 - 5] = False; bm[int(yb) + 5:] = False; bm[:, :s1 + 5] = False; bm[:, xr:] = False
        exp = EXP[ri][ci_]
        bc = [(u + v) / 2 for u, v in runs(np.where(bm.any(0))[0]) if v - u > 20]           # detected bar centres
        big_k = [k for k in range(20) if abs(v2p(exp[k]) - zy) >= 1]                          # bars at least 1 px tall are drawn
        fitok = len(bc) == len(big_k)
        if fitok:
            cf = np.polyfit(np.array(big_k) + 1, bc, 1); cres = np.max(np.abs(np.polyval(cf, np.array(big_k) + 1) - bc))
        check(f'({lab}) visible bars sit at lags 1-20, evenly spaced (centre fit within 2 px)', fitok and cres < 2,
              f'{len(bc)} bars found for {len(big_k)} expected; ' + (f'spacing {cf[0]:.1f} px, max residual {cres:.2f} px' if fitok else ''))
        px_unit = cf[0] if fitok else np.nan; cx = np.polyval(cf, np.arange(1, 21)) if fitok else np.full(20, np.nan)
        allowed = np.zeros_like(bm)
        for k in range(20):
            hpx = abs(v2p(exp[k]) - zy); col = bm[:, int(cx[k]) - 5:int(cx[k]) + 6].any(1)
            ya, yz = sorted([v2p(exp[k]), zy]); allowed[int(ya) - 3:int(yz) + 4, int(cx[k] - 0.3 * px_unit) - 3:int(cx[k] + 0.3 * px_unit) + 4] = True
            if hpx < 3: tot['tiny'] += 1; continue
            ys = np.where(col)[0]
            edge = ys.min() - 0.5 if exp[k] > 0 else ys.max() + 0.5
            err = abs(edge - v2p(exp[k])); tot['maxerr'] = max(tot['maxerr'], err); tot['read'] += 1
        stray = bm & ~allowed
        check(f'({lab}) no blue ink outside the 20 expected bars', not stray.any(), f'{int(stray.sum())} px')
        rr = [tri[0], tri[2]] if len(tri) == 3 else []
        be = max(abs(rr[0] - v2p(band)), abs(rr[1] - v2p(-band))) if len(rr) == 2 else np.inf; tot['band'] = max(tot['band'], be)
        check(f'({lab}) both dashed bands at +/- 1.96/sqrt(60) within 2 px', be <= 2, f'{be:.2f} px')
check('every readable bar edge within 2 px of its recomputed value', tot['maxerr'] <= 2,
      f"{tot['read']} bars read, max error {tot['maxerr']:.2f} px; {tot['tiny']} bars under 3 px tall (|value| < ~0.006) not readable")
print(f'RESULT: {len(fails)} failed check(s)'); sys.exit(1 if fails else 0)
