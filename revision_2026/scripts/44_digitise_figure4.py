"""Image-level check of the DELIVERED Figure 4 TIFF (Python; independent of the plotting code).
x scale from the quarterly tick marks under each panel, y scale from the gridlines; four overlapping line series.
  - readback: vertices lying clear of the other series (> 35 px from their lines) are located by a local
    line(+marker) model fit and must be within 2 px: every such vertex of the observed series (it has markers),
    and the turning points (>= 30 degree bend) of the model lines, the only places a marker-less vertex is pinned down;
  - tube: every model-line pixel within half the line width + 2 px of its expected line (bounds the whole stroke);
  - no extra ink: every pixel of a series colour (after removing 1-2 px anti-aliasing fringes by morphological
    opening) lies on that series' expected line;
  - coverage: every point along every expected line is inked in its own colour or by a series drawn on top of it.
Expected values: observed cases from the workbook; predictions recomputed here from the raw files (numpy QR).
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

# ---- expected values (own implementation: QR via numpy.linalg.qr) ----
import openpyxl
rows = [[c.value for c in r] for r in openpyxl.load_workbook(os.environ.get('PREPRINT_XLSX', 'data_raw/paper_finalized_data_RP.xlsx')).active.iter_rows()]
h = next(i for i, r in enumerate(rows) if 'Confirmed Cases' in r); ci = rows[h].index('Confirmed Cases')
cases = np.array([float(r[ci]) for r in rows[h + 1:] if r[ci] is not None]); assert len(cases) == 72
g = lambda f: pd.read_csv(f, skiprows=2).iloc[:, 1].astype(str).replace('<1', '0.5').astype(float).values
s = np.mean([g('data_raw/gt_dengue_term.csv'), g('data_raw/gt_dengue_term_day2_2026-11-22.csv'), g('data_raw/gt_dengue_term_day3_2026-11-23.csv')], 0)
LY, LX, MO = np.log(cases + 1), np.log(s + 0.5), np.tile(np.arange(1, 13), 6)
def D(r, M):
    X = np.c_[(MO[r][:, None] == np.arange(1, 13)).astype(float), LY[r - 1]]
    return np.c_[X, LX[r]] if M == 'M2' else X
def pr(tr, te, M):
    Q, Rm = np.linalg.qr(D(tr, M)); b = np.linalg.solve(Rm, Q.T @ LY[tr]); return np.maximum(np.expm1(D(te, M) @ b), 0)
test = np.arange(36, 72)
E = {}
for sch in ('season', 'monthly'):
    E[(sch, 'M0')] = cases[test - 1]; E[(sch, 'obs')] = cases[test]
    for M in ('M1', 'M2'):
        E[(sch, M)] = (np.concatenate([pr(np.arange(6, 36 + 12 * k), test[12 * k:12 * k + 12], M) for k in range(3)]) if sch == 'season'
                       else np.array([pr(np.arange(6, t), np.array([t]), M)[0] for t in test]))

im = Image.open('manuscript/figures/Fig4.tif'); assert im.mode == 'CMYK'
A = np.asarray(ImageCms.applyTransform(im, ImageCms.buildTransform(ICC, ImageCms.createProfile('sRGB'), 'CMYK', 'RGB'))).astype(np.int16)
near = lambda rgb, t: np.abs(A - rgb).sum(2) <= t
ink2, grid = near([76, 75, 73], 30), near([216, 215, 210], 18)
COL = {'M0': [154, 153, 150], 'M1': [57, 115, 185], 'M2': [232, 99, 44], 'obs': [0, 0, 7]}
ORDER = ['M0', 'M1', 'M2', 'obs']                                          # drawing order (later = on top)
HALF = {'M0': 1200 / 72 / 2, 'M1': 1200 / 72 / 2, 'M2': 1200 / 72 / 2, 'obs': 1.2 * 1200 / 72 / 2}
DISC = 2.2 * 1200 / 72 / 2                                                 # observed-series marker radius
# Opening removes 1-2 px fringes only for grey M0 and black (edges of the black line pass through grey); blue and
# orange have no look-alike fringes, so their raw pixels are tested (this catches small shifts where lines overlap).
masks = {k: (ndi.binary_opening(near(v, 60), structure=np.ones((3, 3)), iterations=2) if k in ('M0', 'obs') else near(v, 60)) for k, v in COL.items()}

cs = ink2.sum(0); sp = [r for r in runs(np.where(cs > 0.5 * cs.max())[0]) if r[1] - r[0] + 1 >= 9]
check('one left spine band', len(sp) == 1, str(sp)); xs0, xs1 = sp[0]
panels = [r for r in runs(np.where(ink2[:, xs0:xs1 + 1].all(1))[0]) if r[1] - r[0] > 500]
check('two panels', len(panels) == 2, str(panels))
W = A.shape[1]
for (y0, y1), sch, lab in zip(panels, ('season', 'monthly'), 'ab'):
    bot = runs(np.where(ink2[y0:y1 + 40, xs1 + 20:xs1 + 2000].sum(1) > 1500)[0]); ybot = y0 + np.mean(bot[-1])
    tk = ink2[int(ybot) + 8:int(ybot) + 50, xs1 + 3:W - 2].all(0)
    tx = [(xs0 + xs1) / 2] + [xs1 + 3 + np.mean(r) for r in runs(np.where(tk)[0])]     # Jan 2022 tick lies on the left spine (x-axis starts there)
    check(f'({lab}) 12 quarterly x calibration points (spine + 11 tick marks)', len(tx) == 12, f'{len(tx)}')
    days = np.array([(d - pd.Timestamp('1970-01-01')).days for d in pd.date_range('2022-01-01', '2024-10-01', freq='QS')], float)
    xc = np.polyfit(np.array(tx), days, 1); xres = np.max(np.abs(np.polyval(xc, tx) - days)) / abs(xc[0])
    hy = [y0 - 12 + np.mean(r) for r in runs(np.where(grid[y0 - 12:y1, xs1 + 20:W - 20].sum(1) > 0.5 * (W - 40 - xs1))[0])]
    yv = np.array([40000, 30000, 20000, 10000, 0], float); yp = np.array(sorted(hy) + [ybot])
    check(f'({lab}) 4 gridlines + axis', len(hy) == 4, str([round(v) for v in hy]))
    yc = np.polyfit(yp, yv, 1); yres = np.max(np.abs(np.polyval(yc, yp) - yv)) / abs(yc[0])
    check(f'({lab}) x and y scales fit within 2 px', xres < 2 and yres < 2, f'x {xres:.2f} px, y {yres:.2f} px; 1 px = {abs(yc[0]):.1f} cases')
    xpx = (np.array([((pd.Timestamp(f'{2022 + i // 12}-{i % 12 + 1:02d}-01') + pd.Timedelta(days=14)) - pd.Timestamp('1970-01-01')).days for i in range(36)], float) - xc[1]) / xc[0]
    P = {k: np.c_[xpx, (E[(sch, k)] - yc[1]) / yc[0]] for k in ORDER}
    inside = np.zeros(A.shape[:2], bool); inside[y0 - 30:int(ybot) + 25, xs1 + 6:W - 5] = True
    def dist_poly(pts, Q):
        dm = np.full(len(pts), np.inf)
        for p, q in zip(Q[:-1], Q[1:]):
            d = q - p; t = np.clip(((pts - p) @ d) / (d @ d), 0, 1); dm = np.minimum(dm, np.hypot(*(pts - (p + t[:, None] * d)).T))
        return dm
    tot_read = tot_skip = 0; worst = 0.0
    for k in ORDER:
        m = masks[k] & inside
        yy, xx = np.where(m); pts = np.c_[xx, yy].astype(float)
        lim = (DISC if k == 'obs' else HALF[k]) + 2
        dm = dist_poly(pts, P[k]) if len(pts) else np.array([0.0])
        if k == 'obs':
            tb = (pts[:, 1] < y0 + 0.2 * (ybot - y0)) & (pts[:, 0] < xs1 + 0.55 * (W - xs1))      # MAE text box region
            off = pts[(dm > lim) & ~tb]
            check(f'({lab}) obs: black ink is on the observed line or inside the text box only', len(off) == 0 and dm[~tb].max() <= lim,
                  f'{len(pts)} px; max distance outside text box {dm[~tb].max():.1f} (limit {lim:.1f}); text box ink {int(tb.sum())} px')
        else:
            check(f'({lab}) {k}: no pixel of its colour more than 2 px outside its expected stroke', dm.max() <= lim,
                  f'{len(pts)} px; max distance {dm.max():.1f} (limit {lim:.1f} = half line width + 2)')
        # coverage: own colour, or hidden under a series drawn later
        later = np.zeros_like(m)
        for k2 in ORDER[ORDER.index(k) + 1:]: later |= near(COL[k2], 60)
        cover = ndi.binary_dilation(near(COL[k], 60) | later, iterations=3)
        Q = P[k]; samp = np.vstack([p + np.outer(np.linspace(0, 1, max(2, int(np.hypot(*(q - p)) / 2))), q - p) for p, q in zip(Q[:-1], Q[1:])])
        ok = cover[np.round(samp[:, 1]).astype(int), np.round(samp[:, 0]).astype(int)]
        check(f'({lab}) {k}: expected line fully inked (own colour or under a later line)', ok.all(), f'{len(ok)} samples; {int((~ok).sum())} without ink')
        # readback of isolated vertices
        others = [P[k2] for k2 in ORDER if k2 != k]
        own = near(COL[k], 60) & inside
        for i in range(36):
            v = P[k][i]
            if min(dist_poly(v[None, :], Q2)[0] for Q2 in others) <= 35: tot_skip += 1; continue
            if k != 'obs':
                if i in (0, 35): tot_skip += 1; continue
                u1, u2 = P[k][i] - P[k][i - 1], P[k][i + 1] - P[k][i]
                bend = np.degrees(np.arccos(np.clip(u1 @ u2 / np.hypot(*u1) / np.hypot(*u2), -1, 1)))
                if bend < 30: tot_skip += 1; continue                 # vertex not identifiable along a nearly straight stroke
            cx, cy = v; WR = 26
            Yw, Xw = np.mgrid[int(cy) - WR:int(cy) + WR + 1, int(cx) - WR:int(cx) + WR + 1]
            win = np.hypot(Xw - cx, Yw - cy) <= WR; obs_ = own[Yw, Xw]; best = None
            nb = [P[k][j] for j in (i - 1, i + 1) if 0 <= j < 36]
            for dy in range(-8, 9):
                for dx in range(-4, 5):
                    c = np.array([cx + dx, cy + dy]); model = np.zeros_like(obs_)
                    if k == 'obs': model |= np.hypot(Xw - c[0], Yw - c[1]) <= DISC
                    for q in nb:
                        d = q - c; t = np.clip(((Xw - c[0]) * d[0] + (Yw - c[1]) * d[1]) / (d @ d), 0, 1)
                        model |= np.hypot(Xw - c[0] - t * d[0], Yw - c[1] - t * d[1]) <= HALF[k]
                    sc = int(((model ^ obs_) & win).sum())
                    if best is None or sc < best[0]: best = (sc, dx, dy)
            worst = max(worst, abs(best[1]), abs(best[2])); tot_read += 1
    check(f'({lab}) every isolated vertex within 2 px of its value', worst <= 2, f'{tot_read} vertices read back (observed: all clear of other lines; model lines: turning points); max shift {worst:.0f} px; {tot_skip} not read singly')
print(f'RESULT: {len(fails)} failed check(s)'); sys.exit(1 if fails else 0)
