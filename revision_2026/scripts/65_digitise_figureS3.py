"""Image-level check of the DELIVERED Supplementary Figure S3 TIFF (Python; independent of the plotting code).
Per panel: y scale from the two gridlines (0.5, 1) and the bottom spine (0); x scale from the 12 month tick marks.
Expected values recomputed here from the raw files; event positions from the events file with this script's own
date arithmetic (datetime) and category rule.
  - tube: every blue pixel / every black pixel (after removing 1-2 px fringes) in a data area lies within half the line
    width + 2 px of its expected line;
  - coverage: every point along each expected line is inked in its own colour (blue may be covered by the black line);
  - readback: every strict turning point clear of the other line (> 35 px) is located from the stroke edge, within 2.5 px;
  - events: marker count per panel, shape (ring / square / triangle) and bounding-box centre within 2 px.
"""
import os, csv, datetime
import numpy as np, pandas as pd, openpyxl
from PIL import Image, ImageCms
from scipy import ndimage as ndi
Image.MAX_IMAGE_PIXELS = None
ICC = os.environ.get('CMYK_ICC', '/usr/share/texlive/texmf-dist/tex/generic/colorprofiles/FOGRA39L_coated.icc')
TIF = os.environ.get('FIGS3_TIF', 'manuscript/figures/FigS3.tif')
fails = []
def check(n, ok, d=''):
    print(f"[{'PASS' if ok else 'FAIL'}] {n}" + (f'  ({d})' if d else '')); ok or fails.append(n)
def runs(idx):
    if len(idx) == 0: return []
    br = np.where(np.diff(idx) > 1)[0]; return list(zip(np.r_[idx[0], idx[br + 1]], np.r_[idx[br], idx[-1]]))

rows = [[c.value for c in r] for r in openpyxl.load_workbook(os.environ.get('PREPRINT_XLSX', 'data_raw/paper_finalized_data_RP.xlsx')).active.iter_rows()]
h = next(i for i, r in enumerate(rows) if 'Confirmed Cases' in r); ci = rows[h].index('Confirmed Cases')
cases = np.array([float(r[ci]) for r in rows[h + 1:] if r[ci] is not None]); assert len(cases) == 72
def gt(f): return np.array([0.5 if r[1] == '<1' else float(r[1]) for r in csv.reader(open(f, encoding='utf-8-sig')) if len(r) == 2 and r[0][:2] == '20' and len(r[0]) == 7])
s = (gt('data_raw/gt_dengue_term.csv') + gt('data_raw/gt_dengue_term_day2_2026-11-22.csv') + gt('data_raw/gt_dengue_term_day3_2026-11-23.csv')) / 3
EXP = {}
for i, y in enumerate(range(2019, 2025)):
    c = cases[12 * i:12 * i + 12].copy(); ss = s[12 * i:12 * i + 12]
    if y == 2019: c[:5] = np.nan
    EXP[y] = {'cases': c / np.nanmax(c), 'search': ss / ss.max()}
EV = {}
for r in csv.DictReader(open('governance/announcement_events.csv', encoding='utf-8')):
    t = r['event_type']
    if t.startswith('no national') or '5 dates' in t: continue
    cat = 'count' if 'count' in t else ('alert' if ('alert' in t or 'response' in t) else 'report')
    EV.setdefault(r['date'], cat)
for r in csv.DictReader(open('out/06_2024_reporting_lag.csv', encoding='utf-8')):
    if not r['as_of'].endswith('12-31'): EV[r['as_of']] = 'count'
def xm(d):
    d = datetime.date.fromisoformat(d); n = ((d.replace(day=28) + datetime.timedelta(days=4)).replace(day=1) - datetime.timedelta(days=1)).day
    return (d.month - 1) + (d.day - 0.5) / n

im = Image.open(TIF); assert im.mode == 'CMYK'
A = np.asarray(ImageCms.applyTransform(im, ImageCms.buildTransform(ICC, ImageCms.createProfile('sRGB'), 'CMYK', 'RGB'))).astype(np.int16)
near = lambda rgb, t: np.abs(A - rgb).sum(2) <= t
ink2, grid, blue = near([76, 75, 73], 30), near([216, 215, 210], 18), near([57, 115, 185], 60)
blk_raw = near([0, 0, 7], 60); blk = ndi.binary_opening(blk_raw, structure=np.ones((3, 3)), iterations=2)
HB, HK = 1200 / 72 * 1.0 / 2, 1200 / 72 * 1.2 / 2

cs = ink2.sum(0); sp = [r for r in runs(np.where(cs > 0.4 * cs.max())[0]) if r[1] - r[0] + 1 >= 8]
check('two left spines', len(sp) == 2, str(sp))
panels = []
for (a, b) in sp:
    pr = [r for r in runs(np.where(ink2[:, a:b + 1].all(1))[0]) if r[1] - r[0] > 800]
    check(f'three panels on spine at x={a}', len(pr) == 3, str(pr))
    panels += [(a, b, t, u) for (t, u) in pr]
panels.sort(key=lambda p: (p[2] // 500, p[0]))                      # row-major: 2019 2020 / 2021 2022 / 2023 2024
check('six panels', len(panels) == 6)
W = A.shape[1]; res_all = {'tube': [], 'cov': [], 'turn': [], 'ev': []}
for (a, b, top, bot), y in zip(panels, range(2019, 2025)):
    right = (W - 40) if a > W / 3 else (sp[1][0] - 80)          # left column ends just before the right panels' y ticks
    bs = runs(np.where(ink2[top:bot + 40, b + 30:(W - 40 if a > W / 3 else sp[1][0] - 80)].mean(1) > 0.5)[0]); yb = top + (bs[-1][0] + bs[-1][1]) / 2  # bottom spine centre
    yt = [top + (r[0] + r[1]) / 2 for r in runs(np.where(ink2[top - 5:bot + 10, a - 40:a - 12].any(1))[0])]          # y ticks 0, 0.5, 1
    yt = [t - 5 for t in yt]
    ky = np.polyfit([1, 0.5, 0][:len(yt)], yt, 1) if len(yt) == 3 else (1, 0)
    gr = [top + (r[0] + r[1]) / 2 for r in runs(np.where(grid[top:bot, b + 5:right].sum(1) > 0.1 * (right - b))[0])]
    if os.environ.get('DBG'): print('   dbg', y, [round(v, 1) for v in yt], round(yb, 1), [round(g_, 1) for g_ in gr])
    ry = max(abs(np.polyval(ky, 0) - yb), *(min(abs(np.polyval(ky, v) - g_) for v in (0.5, 1)) for g_ in gr)) if len(yt) == 3 else 99
    check(f'{y}: 3 y ticks agree with bottom spine and visible gridlines (< 1.5 px)', len(yt) == 3 and len(gr) >= 1 and ry < 1.5, f'{len(yt)} ticks, {len(gr)} gridlines, {ry:.2f} px')
    tk = runs(np.where(ink2[int(yb) + 12:int(yb) + 40, b + 2:right].any(0))[0]); tx = np.array([b + 2 + (r[0] + r[1]) / 2 for r in tk])
    ok = len(tx) == 12
    kx = np.polyfit(np.arange(12) + 0.5, tx, 1) if ok else (1, 0)
    rx = np.abs(np.polyval(kx, np.arange(12) + 0.5) - tx).max() if ok else 99
    check(f'{y}: 12 month ticks, x-scale residual < 1.5 px', ok and rx < 1.5, f'{len(tx)} ticks, {rx:.2f} px')
    X = lambda v: np.polyval(kx, v); Y = lambda v: np.polyval(ky, v)
    y_top = int(Y(1.25)) + 2
    # all per-panel arrays are crops (padding 80 px) to keep memory small; (oy, ox) = crop origin
    oy, ox = y_top - 80, b + 7 - 80; cy, cx = slice(oy, int(yb) + 80), slice(ox, right + 80)
    BL, BK, BKr = blue[cy, cx], blk[cy, cx], blk_raw[cy, cx]
    area = np.zeros(BL.shape, bool); area[80:int(yb) - 7 - oy, 80:right - ox] = True
    pts = {k: np.c_[X(np.arange(12) + 0.5), Y(v)][~np.isnan(v)] for k, v in EXP[y].items()}
    def dist_map(P, H=0):
        """Distance of every crop pixel to polyline P (exact segment distance). H > 0: both ends extended by H along
        the end segments, as drawn by matplotlib's default 'projecting' line caps."""
        P = P - [ox, oy]
        if H:
            P = P.copy(); u0 = (P[0] - P[1]) / np.hypot(*(P[0] - P[1])); u1 = (P[-1] - P[-2]) / np.hypot(*(P[-1] - P[-2])); P[0] += H * u0; P[-1] += H * u1
        yy, xx = np.mgrid[0:BL.shape[0], 0:BL.shape[1]].astype(np.float32); d = np.full(BL.shape, np.inf, np.float32)
        for p0, p1 in zip(P[:-1], P[1:]):
            v = (p1 - p0).astype(np.float32); t = np.clip(((xx - p0[0]) * v[0] + (yy - p0[1]) * v[1]) / (v @ v), 0, 1)
            d = np.minimum(d, np.hypot(xx - p0[0] - t * v[0], yy - p0[1] - t * v[1]))
        return d
    Ds, Dc = dist_map(pts['search']), dist_map(pts['cases'])
    Ts, Tc = dist_map(pts['search'], HB), dist_map(pts['cases'], HK)        # for the tube only
    tb = Ts[BL & area].max() - HB; tk_ = Tc[BK & area].max() - HK
    res_all['tube'] += [tb, tk_]
    if os.environ.get('DBG'):
        for nm, m_, D_ in (('blue', BL, Ts), ('black', BK, Tc)):
            yy, xx = np.nonzero(m_ & area & (D_ > (HB if nm == 'blue' else HK) + 2))
            if len(yy): print('   tube', y, nm, len(yy), 'px; worst at', ox + xx[np.argmax(D_[yy, xx])], oy + yy[np.argmax(D_[yy, xx])])
    # coverage: samples every 2 px along each expected line
    cov = []
    edt_b, edt_k = ndi.distance_transform_edt(~BL), ndi.distance_transform_edt(~BKr)
    for k, P in pts.items():
        for p0, p1 in zip(P[:-1], P[1:]):
            n = int(np.hypot(*(p1 - p0)) // 2) + 1
            for t in np.linspace(0, 1, n):
                q = p0 + t * (p1 - p0); iy, ix = int(round(q[1])) - oy, int(round(q[0])) - ox
                if k == 'search':
                    if Dc[iy, ix] <= HK + 3: continue
                    cov.append(edt_b[iy, ix])
                else: cov.append(edt_k[iy, ix])
    res_all['cov'].append(max(cov))
    # turning points
    for k, P in pts.items():
        m = blue if k == 'search' else blk_raw; other = Dc if k == 'search' else Ds; H_ = HB if k == 'search' else HK
        for i in range(1, len(P) - 1):
            up, dn = P[i, 1] < min(P[i - 1, 1], P[i + 1, 1]), P[i, 1] > max(P[i - 1, 1], P[i + 1, 1])
            ix, iy = int(round(P[i, 0])), int(round(P[i, 1]))
            if not (up or dn) or other[iy - oy, ix - ox] <= 35: continue
            if dn and P[i, 1] + H_ + 3 > yb - 7: continue            # stroke reaches the bottom spine: edge not visible
            col = m[:, ix - 1:ix + 2].any(1); col[:y_top] = False; col[int(yb) - 7:] = False
            r_ = [r for r in runs(np.where(col)[0]) if r[0] - 60 <= iy <= r[1] + 60]
            if not r_: res_all['turn'].append(99); continue
            edge = r_[0][0] + H_ if up else r_[0][1] - H_
            res_all['turn'].append(abs(edge - P[i, 1]))
            if os.environ.get('DBG') and abs(edge - P[i, 1]) > 2.5: print('   turn', y, k, i + 1, 'up' if up else 'dn', round(edge - P[i, 1], 2), r_)
    # events
    want = sorted((xm(d), c) for d, c in EV.items() if d.startswith(str(y)))
    strip = np.zeros_like(ink2); strip[y_top:int(Y(1.05)), b + 7:right] = True
    if y == 2021:
        xs_ = np.where((ink2 & strip).any(0))[0]
        check('2021: no markers; strip ink only where the "No national alert found" text sits (x 2.5-9.5 months)',
              len(want) == 0 and len(xs_) > 0 and xs_.min() > X(2.5) and xs_.max() < X(9.5))
        continue
    lab, nl = ndi.label(ndi.binary_closing(ink2 & strip, iterations=2), structure=np.ones((3, 3)))
    got = []
    for sl in ndi.find_objects(lab):
        reg = lab[sl] > 0; raw = reg.sum(); hh, ww = reg.shape
        if hh < 30 or ww < 30: continue
        kind = 'count' if ndi.binary_fill_holes(reg).sum() > 1.8 * raw else ('report' if raw > 0.85 * hh * ww else 'alert')
        got.append(((sl[1].start + sl[1].stop - 1) / 2, (sl[0].start + sl[0].stop - 1) / 2, kind, sl[0].start - 0.5))
    got.sort()
    ok = len(got) == len(want) and all(g[2] == w_[1] for g, w_ in zip(got, want))
    # x: bounding-box centre (all three shapes are left-right symmetric). y: squares and rings by bounding-box centre;
    # the 'v' triangle by its flat top edge, which lies exactly half the marker size + half the 1 pt edge above the
    # marker centre (the mitred apex below does not affect it).
    TOP = (0.5 * 4.2 + 0.5) * 1200 / 72
    e = [max(abs(g[0] - X(w_[0])), abs((g[3] - (Y(1.15) - TOP)) if g[2] == 'alert' else (g[1] - Y(1.15)))) for g, w_ in zip(got, want)] if ok else [99]
    res_all['ev'] += e; None
    check(f'{y}: {len(want)} event markers, shapes = categories', ok, f'got {[g[2] for g in got]}')
check('tube: every blue / black data pixel within half line width + 2 px of its expected line', max(res_all['tube']) <= 2,
      f'max excess {max(res_all["tube"]):.2f} px')
check('coverage: every point of each expected line inked (centreline within 3 px)', max(res_all['cov']) <= 3, f'max {max(res_all["cov"]):.2f} px')
check(f'readback: {len(res_all["turn"])} turning points within 2.5 px (minima touching the bottom spine skipped)', len(res_all['turn']) >= 15 and max(res_all['turn']) <= 2.5,
      f'max {max(res_all["turn"]):.2f} px')
check(f'events: {len(res_all["ev"])} marker centres within 2 px', len(res_all['ev']) == 15 and max(res_all['ev']) <= 2, f'max {max(res_all["ev"]):.2f} px')
print(f'RESULT: {len(fails)} failed check(s)'); raise SystemExit(1 if fails else 0)
