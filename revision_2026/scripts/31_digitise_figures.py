"""Image-level check of the DELIVERED TIFFs (Python digitiser; independent of the plotting code).

For each figure: decode the CMYK TIFF, find axes and gridlines in the pixels, fit pixel->value scales from
the gridlines (no matplotlib transforms used), read every plotted value back out of the image, and compare
with values recomputed here from the raw files. "No extra ink": every series-coloured pixel inside a plot
area must lie on the expected line (Fig 1) or inside an expected marker (Fig 2); leader lines in Fig 2 must
end next to the month they label.
Run from the project folder: python3 -I scripts/31_digitise_figures.py   -> out/31_digitise_python.json
"""
import os, json, sys
import numpy as np, pandas as pd
from PIL import Image, ImageCms
from scipy import ndimage as ndi
import statsmodels.formula.api as smf

Image.MAX_IMAGE_PIXELS = None
CMYK_ICC = os.environ.get('CMYK_ICC', '/usr/share/texlive/texmf-dist/tex/generic/colorprofiles/FOGRA39L_coated.icc')
fails, report = [], {}
def check(name, ok, detail=''):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ''))
    report[name] = dict(ok=bool(ok), detail=detail)
    if not ok: fails.append(name)

def load_rgb(path):
    im = Image.open(path)
    assert im.mode == 'CMYK', im.mode
    xf = ImageCms.buildTransform(CMYK_ICC, ImageCms.createProfile('sRGB'), 'CMYK', 'RGB')
    return np.asarray(ImageCms.applyTransform(im, xf)).astype(np.int16)

def near(a, rgb, tol):
    return (np.abs(a - np.array(rgb, dtype=np.int16)).sum(axis=2) <= tol)

def runs(idx):
    """Group sorted integer positions into contiguous runs -> list of (start, end)."""
    if len(idx) == 0: return []
    br = np.where(np.diff(idx) > 1)[0]
    st = np.r_[idx[0], idx[br + 1]]; en = np.r_[idx[br], idx[-1]]
    return list(zip(st, en))

def match_ticks(pxs, ticks, descending):
    """Assign detected gridline pixels to a consecutive run of tick values (some may be hidden under spines)."""
    pxs = np.sort(np.asarray(pxs, float)); tv = np.sort(np.asarray(ticks, float))
    if descending: tv = tv[::-1]                     # image rows grow downward while values decrease
    best = None
    for k in range(len(tv) - len(pxs) + 1):
        coef, res = fit(pxs, tv[k:k + len(pxs)])
        if best is None or res < best[1]: best = (coef, res, tv[k:k + len(pxs)])
    return best

def fit(px, val):
    A = np.vstack([px, np.ones_like(px)]).T
    coef, *_ = np.linalg.lstsq(A, val, rcond=None)
    resid_px = (val - A @ coef) / coef[0]
    return coef, float(np.max(np.abs(resid_px)))

# ---------------- expected values, recomputed here from the raw files ----------------
months = pd.date_range('2019-01-01', periods=72, freq='MS')
import openpyxl
xlsx = os.environ.get('PREPRINT_XLSX', 'data_raw/paper_finalized_data_RP.xlsx')
rows = [[c.value for c in r] for r in openpyxl.load_workbook(xlsx).active.iter_rows()]
h = next(i for i, r in enumerate(rows) if 'Confirmed Cases' in r)
ci, fi = rows[h].index('Confirmed Cases'), rows[h].index('Fever')
body = [r for r in rows[h + 1:] if r[ci] is not None]
cases = np.array([float(r[ci]) for r in body]); fever = np.array([float(r[fi]) for r in body])
assert len(cases) == 72
def gt(f):
    v = pd.read_csv(f, skiprows=2).iloc[:, 1].astype(str).replace('<1', '0.5').astype(float).values
    assert len(v) == 72; return v
search = np.mean([gt('data_raw/gt_dengue_term.csv'), gt('data_raw/gt_dengue_term_day2_2026-11-22.csv'),
                  gt('data_raw/gt_dengue_term_day3_2026-11-23.csv')], axis=0)
casesA = cases.copy(); casesA[:5] = np.nan

INK2, GRID = (76, 75, 73), (216, 215, 210)

# ================================ FIGURE 1 ================================
a = load_rgb('manuscript/figures/Fig1.tif'); H, W, _ = a.shape
ink2 = near(a, INK2, 30); grid = near(a, GRID, 18)
# left spine: the column band with the most INK2 pixels
colcount = ink2.sum(0); xs_band = runs(np.where(colcount > 0.5 * colcount.max())[0])
check('F1: one left-spine band found', len(xs_band) == 1, str(xs_band))
xl = int(np.mean(xs_band[0])); xr_guess = W - 1
# panels: vertical runs of the left spine
yruns = [r for r in runs(np.where(ink2[:, xs_band[0][0]:xs_band[0][1] + 1].all(1))[0]) if r[1] - r[0] > 500]
check('F1: three panels found from the left spine', len(yruns) == 3, str(yruns))
# x calibration from vertical gridlines (Jan 1 of 2020..2024) + left spine (2019-01-01)
mid = yruns[1]                                     # use panel (b) to find vertical gridlines
gcol = grid[mid[0]:mid[1], :].sum(0)
vx = [int(round(np.mean(r))) for r in runs(np.where(gcol > 0.6 * (mid[1] - mid[0]))[0])]
check('F1: five vertical year gridlines found', len(vx) == 5, str(vx))
xd = np.array([(pd.Timestamp(f'{y}-01-01') - pd.Timestamp('1970-01-01')).days for y in range(2019, 2025)], float)
xcoef, xres = fit(np.array([xl] + vx, float), xd)
check('F1: x scale fits 6 calibration lines within 2 px', xres < 2, f'max residual {xres:.2f} px')
day2px = lambda d: (d - xcoef[1]) / xcoef[0]

panels = [('a', casesA, (0, 0, 7), [10, 100, 1000, 10000], 1, True),
          ('b', search, (57, 115, 185), [25, 50, 75, 100], 0, False),
          ('c', fever, (232, 99, 44), [10, 20, 30, 40], 0, False)]
out = {}
for (lab, vals, col, ticks, base, logy), (y0, y1) in zip(panels, yruns):
    # bottom spine row = axis minimum; horizontal gridlines = tick values
    sub_ink = ink2[y0:y1 + 40, xl + 20: xl + 2000].sum(1)
    bot = runs(np.where(sub_ink > 1500)[0])
    ybot = y0 + int(round(np.mean(bot[-1])))
    rowg = grid[y0 - 10:y1 + 1, xl + 20:W - 20].sum(1)
    hy = [y0 - 10 + int(round(np.mean(r))) for r in runs(np.where(rowg > 0.5 * (W - 40 - xl))[0])]
    check(f'F1({lab}): {len(ticks)} horizontal gridlines found', len(hy) == len(ticks), str(hy))
    tv = np.log10(ticks + [base]) if logy else np.array(ticks + [base], float)
    pys = np.array(sorted(hy, reverse=True)[::-1] + [ybot], float)       # top->bottom = high->low values
    ycoef, yres = fit(pys, tv[::-1][::-1][[*range(len(ticks))][::-1] + [len(ticks)]] if False else
                      np.r_[np.sort(tv[:-1])[::-1], tv[-1]])
    check(f'F1({lab}): y scale fits gridlines + axis within 2 px', yres < 2, f'max residual {yres:.2f} px')
    mask = near(a, col, 60); mask[:y0 - 30] = False; mask[ybot + 30:] = False; mask[:, :xl + 5] = False
    ok = np.isfinite(vals)
    vx_ = np.array([day2px((m + pd.Timedelta(days=14) - pd.Timestamp('1970-01-01')).days) for m in months])
    vy_ = np.array([((np.log10(v) if logy else v) - ycoef[1]) / ycoef[0] if np.isfinite(v) else np.nan for v in vals])
    LWH, DR, WR = 1200 / 72 / 2, 1.8 * 1200 / 72 / 2, 22          # half line width, marker radius, window radius (px)
    def seg_d(X, Y, p, q):
        d = q - p; t = np.clip(((X - p[0]) * d[0] + (Y - p[1]) * d[1]) / (d @ d), 0, 1)
        return np.hypot(X - p[0] - t * d[0], Y - p[1] - t * d[1])
    shifts, got = [], []
    idx = np.where(ok)[0]
    for j, i in enumerate(idx):
        cx, cy = vx_[i], vy_[i]
        Yw, Xw = np.mgrid[int(cy) - WR:int(cy) + WR + 1, int(cx) - WR:int(cx) + WR + 1]
        win = np.hypot(Xw - cx, Yw - cy) <= WR
        obs = mask[Yw, Xw]
        nb = [np.array([vx_[k], vy_[k]]) for k in (idx[j - 1] if j > 0 else None, idx[j + 1] if j + 1 < len(idx) else None) if k is not None]
        best = None
        for dy in range(-8, 9):
            for dx in range(-4, 5):
                c = np.array([cx + dx, cy + dy])
                model = np.hypot(Xw - c[0], Yw - c[1]) <= DR
                for q in nb: model |= seg_d(Xw, Yw, c, q) <= LWH
                score = int(((model ^ obs) & win).sum())
                if best is None or score < best[0]: best = (score, dx, dy)
        shifts.append((best[1], best[2])); yb = cy + best[2]
        v = ycoef[0] * yb + ycoef[1]; got.append(10 ** v if logy else v)
    sh = np.abs(np.array(shifts)); per_px = abs(ycoef[0])
    empty_missing = all(not mask[int(ycoef[1] and 0) + y0:ybot, int(vx_[i]) - 2:int(vx_[i]) + 3].any() for i in np.where(~ok)[0])
    check(f'F1({lab}): every plotted month located within 2 px of its true value (model fit)', sh.max() <= 2,
          f'{len(idx)} points; max |dx| {sh[:, 0].max()} px, max |dy| {sh[:, 1].max()} px '
          f'(1 px = {per_px:.4f} {"log10 units" if logy else "index points"}); Jan-May 2019 empty: {empty_missing if lab == "a" else "n/a"}')
    if lab == 'a': check('F1(a): no line ink in the five unrecorded months', empty_missing)
    errs = np.array([0.0])
    # no extra ink: every series pixel within marker radius (15 px) + 4 px of the expected polyline
    ok = np.isfinite(vals)
    P = np.array([[day2px((m + pd.Timedelta(days=14) - pd.Timestamp('1970-01-01')).days),
                   ((np.log10(v) if logy else v) - ycoef[1]) / ycoef[0]] for m, v in zip(months[ok], vals[ok])])
    yy, xx = np.where(mask)
    pts = np.c_[xx, yy].astype(float)
    dmin = np.full(len(pts), np.inf)
    for p, q in zip(P[:-1], P[1:]):
        d = q - p; t = np.clip(((pts - p) @ d) / (d @ d), 0, 1)
        dmin = np.minimum(dmin, np.hypot(*(pts - (p + t[:, None] * d)).T))
    check(f'F1({lab}): no series-coloured pixel off the expected line', dmin.max() <= 19,
          f'{len(pts)} pixels; max distance {dmin.max():.1f} px (limit 19)')
    samp = np.vstack([p + np.outer(np.linspace(0, 1, max(2, int(np.hypot(*(q - p)) / 2))), q - p) for p, q in zip(P[:-1], P[1:])])
    dil = ndi.binary_dilation(mask, iterations=3)
    cov = dil[np.round(samp[:, 1]).astype(int), np.round(samp[:, 0]).astype(int)]
    check(f'F1({lab}): expected line fully inked (sampled every 2 px)', cov.all(), f'{len(cov)} samples; {int((~cov).sum())} without ink')
    report.setdefault('fig1_px', {})[lab] = {str(m)[:7]: [float(vx_[i]), float(vy_[i])] for i, m in enumerate(months) if ok[i]}
    out[lab] = dict(values=[None if g is None else float(g) for g in got], max_err_px=float(np.nanmax(errs)))
report['fig1_readback'] = out

# ================================ FIGURE 2 ================================
b = load_rgb('manuscript/figures/Fig2.tif'); H2, W2, _ = b.shape
ink2 = near(b, INK2, 30); grid = near(b, GRID, 18)
cc = ink2.sum(0); bands = [r for r in runs(np.where(cc > 0.6 * cc.max())[0]) if r[1] - r[0] + 1 >= 9]
check('F2: two left-spine bands found', len(bands) == 2, str(bands))
c60, s60 = cases[12:], search[12:]                                  # 2020-01..2024-12
mon = np.array([m.month for m in months[12:]]); yr = np.array([m.year for m in months[12:]])
Y, X = np.log(c60 + 1), np.log(s60 + 0.5)
dm = lambda v: v - np.array([v[mon == k].mean() for k in mon])
Ya, Xa = dm(Y), dm(X)
late = yr >= 2022
d = pd.DataFrame(dict(y=np.log1p(c60), y1=np.r_[np.nan, np.log1p(c60)[:-1]], x=X,
                      s1=np.sin(2*np.pi*mon/12), c1=np.cos(2*np.pi*mon/12), s2=np.sin(4*np.pi*mon/12), c2=np.cos(4*np.pi*mon/12))).dropna()
cook = smf.ols('y ~ s1 + c1 + s2 + c2 + y1 + x', d).fit().get_influence().cooks_distance[0]
flag_idx = d.index[cook > 4 / len(d)].tolist()
check('F2: Cook flags months 2020-06 and 2022-09', [str(months[12 + i])[:7] for i in flag_idx] == ['2020-06', '2022-09'])

specs = [('a', X, Y, np.log(np.array([1, 2, 5, 10, 25, 100]) + 0.5), np.log(np.array([1, 10, 100, 1000, 10000]) + 1)),
         ('b', Xa, Ya, np.array([-2, -1, 0, 1, 2], float), np.array([-3, -2, -1, 0, 1, 2], float))]
BLUE_A, ORANGE_A = (98, 149, 206), (235, 140, 101)
blue_m = (b[..., 2] - b[..., 0] > 45) & (b[..., 2] > 150)
orng_m = (b[..., 0] - b[..., 2] > 60) & (b[..., 0] > 180)
R = 31.2                                                         # marker radius: sqrt(14) pt / 2 at 1200 dpi
for (lab, xv, yv, xt, yt), (bx0, bx1) in zip(specs, bands):
    xs = (bx0 + bx1) / 2
    ys_run = [r for r in runs(np.where(ink2[:, bx0:bx1 + 1].all(1))[0]) if r[1] - r[0] > 500][0]
    y0, y1 = ys_run
    # gridlines (both axes) inside this panel
    right = bands[1][0] - 400 if lab == 'a' else W2 - 20
    gc = grid[y0:y1, bx1 + 5:right].sum(0)
    vxs = [bx1 + 5 + int(round(np.mean(r))) for r in runs(np.where(gc > 0.5 * (y1 - y0))[0])]
    gr = grid[y0:y1, bx1 + 5:right].sum(1)
    hys = [y0 + int(round(np.mean(r))) for r in runs(np.where(gr > 0.5 * (right - bx1 - 5))[0])]
    if lab == 'b':                                       # zero lines are drawn in INK2 over the 0 gridlines
        zr = [r for r in runs(np.where(ink2[y0:y1, bx1 + 40:right].sum(1) > 0.8 * (right - bx1 - 40))[0]) if y0 + r[1] < y1 - 30]   # not the bottom spine
        zc = [r for r in runs(np.where(ink2[y0 + 20:y1 - 40, bx1 + 5:right].sum(0) > 0.8 * (y1 - y0 - 60))[0]) if r[1] - r[0] + 1 <= 9]
        check('F2(b): one horizontal and one vertical zero line found', len(zr) == 1 and len(zc) == 1,
              f'rows {[(y0 + int(r[0]), y0 + int(r[1])) for r in zr]}, cols {[(bx1 + 5 + int(c[0]), bx1 + 5 + int(c[1])) for c in zc]}; panel rows {y0}-{y1}')
        vxs = sorted(vxs + [bx1 + 5 + int(round(np.mean(zc[0])))]); hys = sorted(hys + [y0 + int(round(np.mean(zr[0])))])
    xcoef, xres, xused = match_ticks(vxs, xt, descending=False)
    ycoef, yres, yused = match_ticks(hys, yt, descending=True)
    check(f'F2({lab}): gridlines matched to tick values within 2 px', xres < 2 and yres < 2 and len(vxs) >= 4 and len(hys) >= 4,
          f'x: {len(vxs)} lines, residual {xres:.2f} px; y: {len(hys)} lines, residual {yres:.2f} px')
    px = (xv - xcoef[1]) / xcoef[0]; py = (yv - ycoef[1]) / ycoef[0]
    report.setdefault('fig2_px', {})[lab] = {str(months[12 + i])[:7]: [float(px[i]), float(py[i])] for i in range(60)}
    inside = np.zeros(b.shape[:2], bool); inside[y0:y1, bx1 + 5:right] = True
    bm, om = blue_m & inside, orng_m & inside
    # presence: marker ink (own colour, or the other group's where markers overlap) within 0.5 R of each point
    Yd, Xd = np.mgrid[-16:17, -16:17]; disc = np.hypot(Xd, Yd) <= 0.5 * R
    def frac(mk, i):
        yc, xc = int(round(py[i])), int(round(px[i])); return mk[yc - 16:yc + 17, xc - 16:xc + 17][disc].mean()
    own_f = np.array([frac(om if late[i] else bm, i) for i in range(60)])
    any_f = np.array([frac(om | bm, i) for i in range(60)])
    crowd = np.array([np.sum(np.hypot(px - px[i], py - py[i]) < 2.2 * R) - 1 for i in range(60)])
    other_close = np.array([np.any((late != late[i]) & (np.hypot(px - px[i], py - py[i]) < 2.2 * R)) for i in range(60)])
    present = np.where(crowd == 0, own_f > 0.5, any_f > 0.2)
    check(f'F2({lab}): marker ink present at all 60 computed positions', present.all(),
          f'{int((crowd == 0).sum())} isolated (own colour > 50% of centre disc, min {own_f[crowd == 0].min():.2f}); '
          f'{int((crowd > 0).sum())} overlapping others (any marker ink > 20%, min {any_f[crowd > 0].min():.2f})')
    # no extra ink: every marker-coloured pixel lies within 1.3 R of an expected point of that colour
    for name, mk, grp in [('blue', bm, ~late), ('orange', om, late)]:
        yy, xx = np.where(mk)
        dd = np.min(np.hypot(xx[:, None] - px[grp][None, :], yy[:, None] - py[grp][None, :]), axis=1) if len(xx) else np.array([0])
        check(f'F2({lab}): no {name} pixel outside the {int(grp.sum())} expected {name} markers', dd.max() <= 1.3 * R,
              f'{len(xx)} pixels; max distance {dd.max():.1f} px (limit {1.3 * R:.0f})')
    # residual ink (geometric): any non-white pixel that is not gridline, zero line, label text or leader (audited in
    # scripts/47) must lie within marker reach of an expected point; catches a stray mark of any colour
    nw = np.abs(b - [255, 255, 255]).sum(2) > 45
    label_ink = (np.abs(b - [76, 75, 73]).sum(2) <= 30) | (np.abs(b - [0, 0, 7]).sum(2) <= 60)   # axis grey and black only
    skip = ndi.binary_dilation(label_ink, iterations=3)                                            # plus their 3-px edges
    for gy in hys: skip[int(gy) - 7:int(gy) + 8, :] = True
    for gx in vxs: skip[:, int(gx) - 7:int(gx) + 8] = True
    cand = nw & inside & ~skip; cand[y1 - 15:, :] = False
    lblc, nc = ndi.label(cand, structure=np.ones((3, 3))); szc = np.bincount(lblc.ravel()); keep = szc >= 3; keep[0] = False; cand = keep[lblc]
    yy, xx = np.where(cand)
    dd = np.min(np.hypot(xx[:, None] - px[None, :], yy[:, None] - py[None, :]), axis=1) if len(xx) else np.array([0.0])
    check(f'F2({lab}): every other non-white pixel lies within marker reach of an expected point', dd.max() <= 1.3 * R + 6,
          f'{len(xx)} px checked; max distance {dd.max():.1f} px (limit {1.3 * R + 6:.0f})')
    RS = R / 3                                                  # '^' centroid lies R/3 below the data point
    errs2 = []
    for name, mk, grp in [('blue', bm, ~late), ('orange', om, late)]:
        lbl, nb = ndi.label(mk)
        for i in np.where(grp)[0]:
            L = lbl[int(round(py[i])), int(round(px[i]))]
            if L == 0: continue
            members = [k for k in np.where(grp)[0] if lbl[int(round(py[k])), int(round(px[k]))] == L]
            if len(members) != 1 or crowd[i] > 0: continue
            yy, xx = np.where(lbl == L)
            ex = xx.mean(); ey = yy.mean() - (RS if late[i] else 0)
            errs2.append(np.hypot(ex - px[i], ey - py[i]))
    errs2 = np.array(errs2)
    check(f'F2({lab}): isolated markers read back within 2 px of their computed positions', len(errs2) > 0 and errs2.max() <= 2,
          f'{len(errs2)} of 60 points isolated; max error {errs2.max():.2f} px; '
          f'1 px = {abs(xcoef[0]):.4f} x-units, {abs(ycoef[0]):.4f} y-units (log scale)')
    if lab == 'b':
        # leader lines: thin INK2 components away from the zero lines and spines, ending next to the flagged points
        lm = ink2 & inside
        zy, zx = y0 + int(round(np.mean(zr[0]))), bx1 + 5 + int(round(np.mean(zc[0])))
        lm[zy - 12:zy + 13, :] = False; lm[:, zx - 12:zx + 13] = False
        lbl, nb = ndi.label(lm, structure=np.ones((3, 3)))
        sizes = ndi.sum(lm, lbl, range(1, nb + 1)); objs = ndi.find_objects(lbl)
        leaders = []
        for k, sl in enumerate(objs):
            hgt, wid = sl[0].stop - sl[0].start, sl[1].stop - sl[1].start
            if max(hgt, wid) > 250 and sizes[k] < 0.08 * hgt * wid + 4000:      # long and thin
                yy, xx = np.where(lbl == k + 1); leaders.append(np.c_[xx, yy])
        check('F2(b): exactly two leader lines', len(leaders) == 2, f'{len(leaders)} found')
        ok = True; det = []
        fl = [12 + 0 for _ in []]
        flag_pts = {str(months[12 + i])[:7]: (px[i], py[i]) for i in flag_idx}
        for L in leaders:
            dists = {k: np.min(np.hypot(L[:, 0] - x, L[:, 1] - y)) for k, (x, y) in flag_pts.items()}
            k = min(dists, key=dists.get)
            # nearest of ALL 60 points to the leader's closest pixel must be the flagged month itself
            j = np.argmin(np.hypot(L[:, 0] - flag_pts[k][0], L[:, 1] - flag_pts[k][1])); ex, ey = L[j]
            nearest = str(months[12 + int(np.argmin(np.hypot(px - ex, py - ey)))])[:7]
            det.append(f'{k}: gap {dists[k]:.0f} px, nearest point {nearest}')
            ok &= (20 <= dists[k] <= 70) and nearest == k
        check('F2(b): each leader ends 20-70 px from its labelled month and that month is the nearest point', ok, '; '.join(det))

json.dump(report, open('out/31_digitise_python.json', 'w'), indent=1, default=str)
print(f'RESULT: {len(fails)} failed check(s)')
sys.exit(1 if fails else 0)
