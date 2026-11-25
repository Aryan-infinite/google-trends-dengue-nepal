"""Residual-ink check for the DELIVERED figures (1-4, S1-S3) (Python). Closes the blind spot of the per-colour checks:
a mark in a colour that is not one of the figure's colours. Inside every plot area:
  (1) every non-white pixel is a figure colour, or an anti-aliasing blend lying within 3 px of a figure colour;
  (2) gridline-grey pixels lie on a detected gridline (+/- 7 px);
  (3) black and dark-grey (axis colour) ink only where the figure draws it:
      Fig 1 - dark grey only inside the 'Not recorded' band (panel a); none in panels b and c;
      Fig 2 - black only in the top-left 'r = ...' label area; dark grey: none in panel a; in panel b, apart from
              the zero lines, exactly 16 connected shapes = 2 leader lines + 14 characters of 'Jun 2020' and 'Sep 2022';
      Fig 3 - dark grey only on the zero line and the dashed threshold rows;
      Fig 4 - dark grey: none inside the panels.
Series-coloured ink is already tied to its expected geometry by scripts 31, 38 and 44 (and 32, 39, 45 in R).
For Fig 2, translucent markers with white outlines create pale tints and outline-over-line greys that no colour rule can
separate from a stray mark, so its unknown-ink test is geometric (scripts 31 and 32: every such pixel near an expected marker).
Groups smaller than 3 px are ignored (single-pixel mixes where two coloured edges cross).
"""
import os, sys
import numpy as np
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
def load(p):
    im = Image.open(p); assert im.mode == 'CMYK'
    return np.asarray(ImageCms.applyTransform(im, ImageCms.buildTransform(ICC, ImageCms.createProfile('sRGB'), 'CMYK', 'RGB'))).astype(np.int16)
WHITE, BLACK, INK2, GRID = (255, 255, 255), (0, 0, 7), (76, 75, 73), (216, 215, 210)
PAL = {'Fig1': [BLACK, INK2, GRID, (236, 234, 231), (57, 115, 185), (232, 99, 44)],
       'Fig2': [BLACK, INK2, GRID],                                          # + translucent marker ink (chromatic), see below
       'Fig3': [BLACK, INK2, GRID, (57, 115, 185)],
       'Fig4': [BLACK, INK2, GRID, (154, 153, 150), (57, 115, 185), (232, 99, 44)],
       'FigS1': [BLACK, INK2, GRID, (57, 115, 185)],
       'FigS3': [BLACK, INK2, GRID, (57, 115, 185)],
       'FigS2': [BLACK, INK2, GRID, (57, 115, 185)]}

def axes_boxes(a):
    """Plot areas from the pixels: (x_left_spine_end, y_top, x_right, y_bottom_spine) per panel."""
    ink2 = np.abs(a - INK2).sum(2) <= 30; grid = np.abs(a - GRID).sum(2) <= 18
    cs = ink2.sum(0); bands = [r for r in runs(np.where(cs > 0.5 * cs.max())[0]) if r[1] - r[0] + 1 >= 9]
    bands = [b for b in bands if not any(o[0] < b[0] < o[2] for o in [])]
    out = []
    for s, e in bands:
        for y0, y1 in [r for r in runs(np.where(ink2[:, s:e + 1].all(1))[0]) if r[1] - r[0] > 500]:
            rows = np.where(ink2[y0:y1 + 1, e + 20:e + 1500].sum(1) > 1200)[0]   # bottom spine = last long INK2 row run
            yb = y0 + runs(rows)[-1][0] if len(rows) else y1
            cols = np.where(ink2[yb:yb + 6, e + 1:].any(0))[0]                       # bottom spine, this panel only
            br = np.where(np.diff(cols) > 60)[0]; xr = e + 1 + (cols[br[0]] if len(br) else cols[-1])
            out.append((e + 1, y0, int(xr), int(yb)))
    return out

def big(m, k=3):
    """Drop connected groups smaller than k px (isolated single-pixel colour mixes where two line edges cross)."""
    lbl, n = ndi.label(m, structure=np.ones((3, 3)))
    if n == 0: return m
    sz = np.bincount(lbl.ravel()); keep = sz >= k; keep[0] = False
    return keep[lbl]

def residual(name, a, box, extra_ok=None, reach=3):
    x0, y0, x1, y1 = box
    sub = a[y0:y1, x0:x1]
    nonwhite = np.abs(sub - WHITE).sum(2) > 45
    pure = np.zeros(nonwhite.shape, bool)
    for c in PAL[name]: pure |= np.abs(sub - c).sum(2) <= 30
    if extra_ok is not None: pure |= extra_ok[y0:y1, x0:x1]
    near_pure = ndi.binary_dilation(pure & nonwhite, iterations=reach)
    bad = big(nonwhite & ~pure & ~near_pure)
    return bad, pure, sub

def black_fringe(sub):
    """Pixels within 3 px of black ink: anti-aliased edges of black lines and text pass through both greys."""
    return ndi.binary_dilation(np.abs(sub - BLACK).sum(2) <= 60, iterations=3)

def dark_neutral_fringe(sub):
    """Pixels within 3 px of any darker neutral ink (black, axis grey, M0 grey): their edges fade through gridline grey."""
    neutral = (sub.max(2) - sub.min(2) <= 12) & (sub.mean(2) < 170)
    return ndi.binary_dilation(neutral, iterations=3)

def grid_rule(sub, frac=0.5):
    """frac: share of the row a gridline must show (FigS2 2021: the flat search line covers most of the 0.5 gridline)."""
    grid = (np.abs(sub - GRID).sum(2) <= 18) & ~dark_neutral_fringe(sub)
    rows = [r for r in runs(np.where(grid.sum(1) > frac * grid.shape[1])[0])]
    cols = [r for r in runs(np.where(grid.sum(0) > 0.5 * grid.shape[0])[0])]
    ok = np.zeros(grid.shape, bool)
    for u, v in rows: ok[max(0, u - 7):v + 8, :] = True
    for u, v in cols: ok[:, max(0, u - 7):v + 8] = True
    return big(grid & ~ok)

for name in ('Fig1', 'Fig2', 'Fig3', 'Fig4', 'FigS1', 'FigS3', 'FigS2'):
    if name == 'FigS3' and not os.path.exists('manuscript/figures/FigS3.tif'):    # S2 needs the non-public EDCD district workbook
        print('[SKIP] FigS3: not built (needs the EDCD district workbook, which is not public)'); continue
    a = load(f'manuscript/figures/{name}.tif'); boxes = axes_boxes(a)
    extra = None
    if name == 'Fig2':                                       # translucent markers and their stacks: any saturated ink
        ch = a.max(2) - a.min(2); extra = ch > 40
    for k, box in enumerate(boxes):
        lab = 'abcdefgh'[k]
        if name == 'FigS2':                                                  # boxes come spine by spine; name them by year
            lab = str(2019 + 2 * sorted({b[1] for b in boxes}).index(box[1]) + (0 if box[0] < a.shape[1] // 2 else 1))
        if name == 'Fig2':
            sub = a[box[1]:box[3], box[0]:box[2]]           # unknown-colour and gridline tests for Fig 2 use marker geometry: scripts 31/32
        else:
            bad, pure, sub = residual(name, a, box, extra)
            check(f'{name}({lab}) no ink of an unknown colour', not bad.any(),
                  f'{int(bad.sum())} px' + (f' at x {box[0] + np.where(bad)[1].min()}-{box[0] + np.where(bad)[1].max()}, y {box[1] + np.where(bad)[0].min()}-{box[1] + np.where(bad)[0].max()}' if bad.any() else ''))
            g = grid_rule(sub, 0.15 if name == 'FigS2' else 0.5)
            check(f'{name}({lab}) gridline grey only on gridlines', not g.any(), f'{int(g.sum())} px off gridlines')
        ink2 = (np.abs(sub - INK2).sum(2) <= 30) & ~black_fringe(sub); blk = np.abs(sub - BLACK).sum(2) <= 30
        H, Wd = ink2.shape
        if name == 'Fig1':
            if lab == 'a':
                band = np.abs(sub - (236, 234, 231)).sum(2) <= 30; cols = np.where(band.sum(0) > 0.3 * H)[0]
                outside = ink2.copy(); outside[:, :cols.max() + 1] = False if len(cols) else outside[:, :0]
                outside[H - 15:, :] = False                                         # bottom spine
                check('Fig1(a) dark grey only inside the Not recorded band', not outside.any(), f'{int(outside.sum())} px outside')
            else:
                o = ink2.copy(); o[H - 15:, :] = False
                check(f'Fig1({lab}) no dark grey inside the panel', not o.any(), f'{int(o.sum())} px');
        if name == 'Fig2':
            tl = np.zeros_like(blk); tl[:int(0.15 * H), :int(0.3 * Wd)] = True
            check(f'Fig2({lab}) black only in the top-left r label area', not (blk & ~tl).any(), f'{int((blk & ~tl).sum())} px elsewhere')
            o = ink2.copy(); o[H - 15:, :] = False
            if lab == 'a':
                check('Fig2(a) no dark grey inside the panel', not o.any(), f'{int(o.sum())} px')
            else:
                zr = [r for r in runs(np.where(o.sum(1) > 0.8 * Wd)[0])]; zc = [r for r in runs(np.where(o[20:H - 40].sum(0) > 0.8 * (H - 60))[0])]
                for u, v in zr: o[u - 12:v + 13, :] = False
                for u, v in zc: o[:, u - 12:v + 13] = False
                o = ndi.binary_opening(o, iterations=1)                              # drop 1-px fringe crumbs left at cut edges
                lbl, n = ndi.label(o, structure=np.ones((3, 3)))
                check('Fig2(b) dark-grey shapes besides the zero lines = 2 leaders + 14 characters', n == 16, f'{n} shapes')
        if name == 'Fig3':
            o = ink2.copy(); o[H - 15:, :] = False
            zr = [r for r in runs(np.where(o.sum(1) > 0.5 * Wd)[0])]                 # zero line (partly under bars)
            for u, v in zr: o[u - 12:v + 13, :] = False
            dr = [r for r in runs(np.where(o.sum(1) > 0.15 * Wd)[0])]                # the two dashed-threshold row bands
            for u, v in dr: o[u - 6:v + 7, :] = False                                # +/- 6 px: anti-aliased edge rows of the dashes
            check('Fig3 dark grey only on the zero line and the dashed threshold rows', not o.any() and len(dr) == 2,
                  f'{len(dr)} dashed row bands; {int(o.sum())} px elsewhere')
        if name == 'FigS1':                                                  # zero line (under the bars) and two dashed bands only
            o = ink2.copy(); o[H - 15:, :] = False
            tri = runs(np.where(o.sum(1) > 0.25 * Wd)[0])
            for u, v in tri: o[u - 8:v + 9, :] = False
            check(f'FigS1({lab}) dark grey only on the zero line and the two dashed bands', len(tri) == 3 and not o.any(), f'{len(tri)} bands; {int(o.sum())} px elsewhere')
        if name == 'FigS3':                                                  # panel (b): the vertical reference line at 1 only
            o = ink2.copy(); o[H - 15:, :] = False
            vl = runs(np.where(o.sum(0) > 0.8 * H)[0])
            for u, v in vl: o[:, u - 8:v + 9] = False
            check(f'FigS3({lab}) dark grey only on the reference line (panel b)', len(vl) == (1 if lab == 'b' else 0) and not o.any(), f'{len(vl)} vertical lines; {int(o.sum())} px elsewhere')
        if name == 'FigS2':                                                  # top strip (above 1.05 on a 0-1.25 axis) holds the event
            st = int(0.17 * H); o = ink2.copy(); o[H - 15:, :] = False      # markers and the 2021 note; the lines never reach it
            o[:st, :] = False; bl = np.abs(sub - (57, 115, 185)).sum(2) <= 60
            check(f'FigS2({lab}) dark grey only in the event strip; no line ink in the strip', not o.any() and not blk[:st].any() and not bl[:st].any(),
                  f'{int(o.sum())} grey px below the strip; {int(blk[:st].sum())} black, {int(bl[:st].sum())} blue px in it')
        if name == 'Fig4':
            o = ink2.copy(); o[H - 15:, :] = False
            check(f'Fig4({lab}) no dark grey inside the panel', not o.any(), f'{int(o.sum())} px')
print(f'RESULT: {len(fails)} failed check(s)'); sys.exit(1 if fails else 0)
