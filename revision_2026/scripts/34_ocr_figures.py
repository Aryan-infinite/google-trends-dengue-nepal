"""OCR check of every tick label and text label printed on the DELIVERED figure TIFFs (Figures 1-4, S1-S3; tesseract 5).
Crops are located from the axis lines found in the pixels; each crop's OCR text must equal the expected text.
Run from the project folder: python3 -I scripts/34_ocr_figures.py
"""
import subprocess, tempfile, os, re, sys
import numpy as np
from PIL import Image, ImageCms, ImageOps
Image.MAX_IMAGE_PIXELS = None
ICC = os.environ.get('CMYK_ICC', '/usr/share/texlive/texmf-dist/tex/generic/colorprofiles/FOGRA39L_coated.icc')
fails = []
def check(name, ok, detail):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}  (OCR: {detail})"); ok or fails.append(name)

def load(p):
    im = Image.open(p); rgb = ImageCms.applyTransform(im, ImageCms.buildTransform(ICC, ImageCms.createProfile('sRGB'), 'CMYK', 'RGB'))
    return rgb, np.asarray(rgb).astype(np.int16)

def ocr(img, box, psm=6, rotate=0, whitelist=None):
    box = (max(0, box[0]), max(0, box[1]), min(img.width, box[2]), min(img.height, box[3]))   # never crop past the edge (padding reads as ink)
    c = ImageOps.grayscale(img.crop(box))
    if rotate: c = c.rotate(rotate, expand=True, fillcolor=255)
    c = c.resize((max(1, c.width // 4), max(1, c.height // 4)), Image.LANCZOS)       # 1200 -> 300 dpi
    f = tempfile.mktemp(suffix='.png'); c.save(f, dpi=(300, 300))
    args = ['tesseract', f, '-', '--psm', str(psm)] + (['-c', f'tessedit_char_whitelist={whitelist}'] if whitelist else [])
    out = subprocess.run(args, capture_output=True, text=True).stdout; os.remove(f)
    return ' '.join(out.split())

def norm(s): return s.replace('−', '-').replace('—', '-').replace('–', '-')

def spines(a):
    ink2 = np.abs(a - [76, 75, 73]).sum(2) <= 30; cs = ink2.sum(0)
    cols = np.where(cs > 0.5 * cs.max())[0]; br = np.where(np.diff(cols) > 1)[0]
    bands = [(s, e) for s, e in zip(np.r_[cols[0], cols[br + 1]], np.r_[cols[br], cols[-1]]) if e - s + 1 >= 9]
    out = []
    for s, e in bands:
        rows = np.where(ink2[:, s:e + 1].all(1))[0]; rb = np.where(np.diff(rows) > 1)[0]
        for y0, y1 in zip(np.r_[rows[0], rows[rb + 1]], np.r_[rows[rb], rows[-1]]):
            if y1 - y0 > 500: out.append((s, e, y0, y1))
    return out

NUM = '0123456789,'
def ticks_y(a, s, y0, y1):
    """Rows of the y tick marks: short INK2 runs just left of the spine at column s."""
    ink2 = np.abs(a[y0 - 20:y1 + 20, s - 50:s - 8] - [76, 75, 73]).sum(2) <= 30
    r = np.where(ink2.all(1))[0]; br = np.where(np.diff(r) > 1)[0]
    return [y0 - 20 + (u + v) // 2 for u, v in zip(np.r_[r[0], r[br + 1]], np.r_[r[br], r[-1]])]
def ticks_x(a, x0, x1, ybot):
    """Columns of the x tick marks: short INK2 runs just below the bottom spine."""
    ink2 = np.abs(a[ybot + 10:ybot + 50, x0:x1] - [76, 75, 73]).sum(2) <= 30
    c = np.where(ink2.all(0))[0]; br = np.where(np.diff(c) > 1)[0]
    return [x0 + (u + v) // 2 for u, v in zip(np.r_[c[0], c[br + 1]], np.r_[c[br], c[-1]])]
def ylabels(img, a, s, y0, y1, wl=NUM, left=700):
    return [norm(ocr(img, (s - left, t - 100, s - 70, t + 100), psm=7, whitelist=wl)) for t in ticks_y(a, s, y0, y1)]
def xlabels(img, a, x0, x1, ybot, wl=None):
    tx = ticks_x(a, x0, x1, ybot); half = min(300, int(np.min(np.diff(tx)) / 2) - 10) if len(tx) > 1 else 300
    return [norm(ocr(img, (t - half, ybot + 65, t + half, ybot + 260), psm=7, whitelist=wl)) for t in tx]
# ---------------- Figure 1 ----------------
img, a = load('manuscript/figures/Fig1.tif'); sp = spines(a)
print('Figure 1 axes (x0, x1, y0, y1):', sp)
exp_y = [['10,000', '1,000', '100', '10', '1'], ['100', '75', '50', '25', '0'], ['40', '30', '20', '10', '0']]
for (s, e, y0, y1), exp, lab in zip(sp, exp_y, 'abc'):
    got = ylabels(img, a, s, y0, y1)
    check(f'F1({lab}) y tick labels = {exp}', got == exp, got)
    got = ocr(img, (40, y0 - 40, 420, y0 + 160), psm=7)
    check(f'F1({lab}) panel letter', got == f'({lab})', got)
s, e, y0, y1 = sp[2]
got = ocr(img, (s, y1 + 15, img.width, y1 + 260), psm=7, whitelist='0123456789')
check('F1 x labels = 2019..2024', got.split() == [str(y) for y in range(2019, 2025)], got)
ylab = {'a': 'Reported cases (log scale)', 'b': '"dengue" search index', 'c': '"Fever" search index'}
for (s, e, y0, y1), lab in zip(sp, 'abc'):
    got = ocr(img, (0, y0 + 200, s - 600, y1), psm=6, rotate=-90)
    check(f'F1({lab}) y axis title', norm(got).replace('”', '"').replace('“', '"') == ylab[lab], got)
s, e, y0, y1 = sp[0]
got = ocr(img, (s + 130, y0 + 1000, s + 320, y1 - 20), psm=7, rotate=-90)
check('F1(a) band label "Not recorded"', got == 'Not recorded', got)

# ---------------- Figure 2 ----------------
img, a = load('manuscript/figures/Fig2.tif'); sp = spines(a)
print('Figure 2 axes (x0, x1, y0, y1):', sp)
(sa, ea, ya0, ya1), (sb, eb, yb0, yb1) = sp
got = ylabels(img, a, sa, ya0, ya1, left=540)
check("F2(a) y tick labels = ['10,000','1,000','100','10','1']", got == ['10,000', '1,000', '100', '10', '1'], got)
got = xlabels(img, a, sa + 12, sb - 900, ya1, wl='0123456789')
check("F2(a) x tick labels = ['1','2','5','10','25','100']", got == ['1', '2', '5', '10', '25', '100'], got)
got = ylabels(img, a, sb, yb0, yb1, wl=None, left=335)
check("F2(b) y tick labels = ['2','1','0','-1','-2','-3']", got == ['2', '1', '0', '-1', '-2', '-3'], got)
got = xlabels(img, a, sb + 12, img.width - 5, yb1)
check("F2(b) x tick labels = ['-2','-1','0','1','2']", got == ['-2', '-1', '0', '1', '2'], got)
got = ocr(img, (sa + 80, ya0 + 20, sa + 900, ya0 + 330), psm=7)
check("F2(a) annotation 'r = 0.87'", got.replace(' ', '') == 'r=0.87', got)
got = ocr(img, (sb + 80, yb0 + 20, sb + 900, yb0 + 330), psm=7)
check("F2(b) annotation 'r = 0.81'", got.replace(' ', '') == 'r=0.81', got)
got = norm(ocr(img, (int(img.width * 0.62), 0, img.width, ya0 - 120), psm=7))
check("F2 legend '2020-2021' '2022-2024' (marker glyphs ignored)", re.findall(r'\d{4}-\d{4}', got) == ['2020-2021', '2022-2024'], got)
got = ocr(img, (sb + 15, yb0, img.width - 10, yb1 - 15), psm=11)
words = got.split()
check("F2(b) month labels 'Jun 2020' and 'Sep 2022' present", 'Jun' in words and 'Sep' in words and '2020' in words and '2022' in words, got[:120])
for (s, y0), lab in [((sa, ya0), 'a'), ((sb, yb0), 'b')]:
    got = ocr(img, (s - 820, y0 - 330, s - 520, y0 - 20), psm=7)
    check(f'F2({lab}) panel letter', got == f'({lab})', got)

# ---------------- Figure 3 ----------------
img, a = load('manuscript/figures/Fig3.tif'); sp = spines(a)
print('Figure 3 axes (x0, x1, y0, y1):', sp)
s3, e3, y30, y31 = sp[0]
got = ylabels(img, a, s3, y30, y31, wl=None, left=420)
check("F3 y tick labels = ['0.6','0.4','0.2','0.0','-0.2','-0.4']", got == ['0.6', '0.4', '0.2', '0.0', '-0.2', '-0.4'], got)
bl = (a[..., 2] - a[..., 0] > 60) & (a[..., 2] > 150); bl[:, :s3 + 20] = False; bl[y31:] = False
cols = np.where(bl.sum(0) > 0)[0]; br = np.where(np.diff(cols) > 1)[0]
bars = [(u, v) for u, v in zip(np.r_[cols[0], cols[br + 1]], np.r_[cols[br], cols[-1]]) if v - u > 40]
cxs = [(u + v) // 2 for u, v in bars]; half = int(np.mean(np.diff(cxs)) / 2) - 10
want = [('-3', 'n=55'), ('-2', 'n=56'), ('-1', 'n=57'), ('0', 'n=58'), ('+1', 'n=57'), ('+2', 'n=56'), ('+3', 'n=55')]
got = []
for c in cxs:                                    # line 1 = lag (signs and digits only), line 2 = n
    lag = norm(ocr(img, (c - half, y31 + 20, c + half, y31 + 215), psm=7, whitelist='+-−0123456789'))
    nn = ocr(img, (c - half, y31 + 200, c + half, y31 + 420), psm=7).replace(' ', '')
    got.append((lag, nn))
check('F3 x tick labels (lag and n under each bar)', got == want, got)
# all black text inside the plot area must be ONE block, centred on the lag 0 bar, reading 'Holm-adjusted P < 0.001'
ink = (np.abs(a - [0, 0, 7]).sum(2) <= 60); ink[:, :s3 + 20] = False; ink[y31 - 20:] = False; ink[:y30] = False
yy, xx = np.where(ink)
box = (xx.min(), yy.min(), xx.max(), yy.max()); cx_txt = (box[0] + box[2]) / 2
txt = ocr(img, (box[0] - 30, box[1] - 30, box[2] + 30, box[3] + 30), psm=6)
check('F3 only text in plot area = "Holm-adjusted P < 0.001", centred on the lag 0 bar',
      txt.replace(' ', '') == 'Holm-adjustedP<0.001' and abs(cx_txt - cxs[3]) < 15 and box[3] < np.where(bl[:, cxs[3]])[0].min(),
      f'OCR "{txt}"; text centre {cx_txt:.0f} px vs lag 0 bar centre {cxs[3]} px; text bottom row {box[3]} above bar top {np.where(bl[:, cxs[3]])[0].min()}')
got = ocr(img, (int(img.width * 0.4), 0, img.width, y30 - 10), psm=7)
check("F3 legend 'Two-sided P = 0.05 (unadjusted)'", got.endswith('Two-sided P = 0.05 (unadjusted)'), got)
got = ocr(img, (s3, y31 + 430, img.width, img.height), psm=7)
check('F3 x axis title', got == 'Lag (months); positive lag = search interest precedes cases', got)
got = ocr(img, (0, y30, s3 - 460, y31), psm=6, rotate=-90)
check('F3 y axis title', got == 'Cross-correlation (r)', got)

# ---------------- Figure 4 ----------------
def bottom_row(a, s, y1):
    """Centre row of the bottom spine: the long horizontal INK2 run just above y1 (y1 may include a tick under the spine)."""
    ink2 = np.abs(a[y1 - 200:y1 + 1, s + 20:s + 2020] - [76, 75, 73]).sum(2) <= 30
    r = np.where(ink2.sum(1) > 1500)[0]; br = np.where(np.diff(r) > 1)[0]; last = (np.r_[r[0], r[br + 1]][-1], r[-1])
    return y1 - 200 + (last[0] + last[1]) // 2
MO_ = lambda t: re.sub(r'\bMO\b', 'M0', t)        # OCR reads the zero in "M0" as a capital O in this typeface
img, a = load('manuscript/figures/Fig4.tif'); sp = spines(a)
print('Figure 4 axes (x0, x1, y0, y1):', sp)
for (s4, e4, y40, y41), lab, title, maes in zip(sp, 'ab', ('Refitted once per season (primary)', 'Refitted every month (post hoc)'),
                                                 ('M0 2,973, M1 2,025, M2 1,900', 'M0 2,973, M1 2,069, M2 1,266')):
    got = ylabels(img, a, s4, y40, y41, left=620)
    check(f"F4({lab}) y tick labels", got == ['40,000', '30,000', '20,000', '10,000', '0'], got)
    got = ocr(img, (s4 + 20, y40 + 20, s4 + 3800, y40 + 480), psm=6)
    check(f'F4({lab}) text box (M0/MO read as the same)', MO_(got) == f'{title} Mean absolute error: {maes}', got)
    got = ocr(img, (s4 - 720, y40 - 300, s4 - 440, y40 - 30), psm=7)
    check(f'F4({lab}) panel letter', got == f'({lab})', got)
    got = ocr(img, (0, y40 + 300, s4 - 640, y41), psm=7, rotate=-90)
    check(f'F4({lab}) y axis title', got == 'Cases per month', got)
s4, e4, y40, y41 = sp[1]; y41 = bottom_row(a, s4, y41)
tx = [s4 + 5] + ticks_x(a, s4 + 12, img.width - 5, y41)
want = ['Jan 2022', 'Apr', 'Jul', 'Oct', 'Jan 2023', 'Apr', 'Jul', 'Oct', 'Jan 2024', 'Apr', 'Jul', 'Oct']
got = [ocr(img, (t - 270, y41 + 60, t + 270, y41 + 460), psm=6) for t in tx]
check('F4 x tick labels (quarters, years under January)', got == want, got)
got = ocr(img, (0, 0, img.width, sp[0][2] - 380), psm=6)
check('F4 key (all four series names; M0/MO read as the same)', all(t in MO_(got) for t in ['Observed', 'M0, last month carried forward', 'M1, seasonal terms + last month', 'M2, M1 + current-month searches']), got)
# ---------------- Figure S1 ----------------
img, a = load('manuscript/figures/FigS1.tif'); sp = spines(a)
titles = ['(a) Log cases', '(b) Log cases', '(c) Log cases, month-adjusted', '(d) Log cases, month-adjusted', '(e) Log search interest',
          '(f) Log search interest', '(g) Log search interest, month-adjusted', '(h) Log search interest, month-adjusted']
order = sorted(sp, key=lambda t: (t[2], t[0]))                                   # row by row, left then right
okt, okyt, bad = 0, 0, []
for k, (s5, e5, y50, y51) in enumerate(order):
    t = ocr(img, (s5 - 40, y50 - 260, s5 + 3700, y50 - 10), psm=7)
    if t == titles[k]: okt += 1
    else: bad.append(t)
    yl = ylabels(img, a, s5, y50, y51, wl=None, left=520) if s5 < img.width // 2 else None
    if yl is not None: okyt += yl == ['1.0', '0.5', '0.0', '-0.5', '-1.0']
check('FS1 eight panel titles', okt == 8, f'{okt}/8 exact' + (f'; differs: {bad}' if bad else ''))
check('FS1 y tick labels on the four left panels', okyt == 4, f'{okyt}/4 exact')
got = ocr(img, (0, 0, img.width, order[0][2] - 270), psm=6)
check('FS1 column headers', got == 'Autocorrelation Partial autocorrelation', got)
for s5, e5, y50, y51 in [t for t in order if t[2] == max(o[2] for o in order)]:
    yb5 = bottom_row(a, s5, y51); tx = ticks_x(a, s5 + 12, s5 + 3800, yb5)
    got = [ocr(img, (t - 150, yb5 + 60, t + 150, yb5 + 250), psm=7, whitelist='0123456789') for t in tx]
    check(f'FS1 x tick labels (panel at x {s5})', got == ['1', '5', '10', '15', '20'], got)
    got = ocr(img, (s5, yb5 + 250, s5 + 3800, yb5 + 520), psm=7)
    check(f'FS1 x axis title (panel at x {s5})', got == 'Lag (months)', got)
# ---------------- Figure S2 ----------------
if not os.path.exists('manuscript/figures/FigS2.tif'):                            # S2 needs the non-public EDCD district workbook
    print('[SKIP] FS2: manuscript/figures/FigS2.tif not built (needs the EDCD district workbook, which is not public)')
else:
    img, a = load('manuscript/figures/FigS2.tif'); sp = spines(a)[:2]                  # third full-height line = reference line
    for (s6, e6, y60, y61), lab, xl, xt in zip(sp, 'ab', ('Correlation, seasonally adjusted (95% CI)', '2023 peak ÷ 2022 peak (log scale)'),
                                               (['0', '0.25', '0.50', '0.75', '1.00'], ['0.1', '0.3', '1', '3', '10'])):
        yb6 = bottom_row(a, s6, y61)
        got = ocr(img, (s6 - 1100, y60, s6 - 40, yb6), psm=6).split()
        check(f'FS2({lab}) row labels', got == ['Central', 'Eastern', 'Western'], got)
        right = sp[1][0] - 800 if lab == 'a' else img.width - 5
        tx = ([s6 + 5] if lab == 'a' else []) + ticks_x(a, s6 + 12, right, yb6)
        def lone_zero(t):
            """A single glyph read as 'O' or nothing is accepted as '0' only if it has a zero's proportions (width/height < 0.8;
            in this typeface a zero is about 0.67 as wide as tall, a capital O about 1.0)."""
            arr = np.asarray(ImageOps.grayscale(img.crop((t - 230, yb6 + 60, t + 230, yb6 + 250)))); yy, xx = np.where(arr < 128)
            return '0' if len(xx) and (xx.max() - xx.min() + 1) / (yy.max() - yy.min() + 1) < 0.8 else '?'
        got = [norm(ocr(img, (t - 230, yb6 + 60, t + 230, yb6 + 250), psm=7, whitelist='0123456789.')) for t in tx]
        got = [g if g not in ('', 'O') else lone_zero(t) for g, t in zip(got, tx)]
        check(f'FS2({lab}) x tick labels', got == xt, got)
        got = ocr(img, (s6 - 500, yb6 + 250, right, yb6 + 520), psm=7).replace('(95% Cl)', '(95% CI)').replace(' peak + 2022', ' peak ÷ 2022')  # OCR: I/l and ÷/+ ambiguous; exact strings checked in scripts/58
        check(f'FS2({lab}) x axis title (CI/Cl read as the same)', got == xl, got)
        lx6 = int(s6 - 0.3 * (right - s6))                                              # panel letters sit at -0.3 axes widths
        got = ocr(img, (lx6 - 150, y60 - 330, lx6 + 520, y60 - 10), psm=7)
        check(f'FS2({lab}) panel letter', got == f'({lab})', got)
    got = ocr(img, (sp[1][0], 0, img.width, sp[1][2] - 30), psm=7)
    check('FS2 legend', all(t in got for t in ('Cases', 'Search interest')), got)
# ---------------- Figure S3 ----------------
img, a = load('manuscript/figures/FigS3.tif'); order = sorted(spines(a), key=lambda t: (t[2], t[0]))       # 2019 2020 / 2021 2022 / 2023 2024
check('FS3 six panels found', len(order) == 6, str(len(order)))
NOTE = {2019: 'Cases not recorded Jan-May', 2020: 'Search index < 3 all year', 2021: 'Search index < 2 all year'}   # OCR reads ≤ as <; exact bytes in scripts/63
okt, okn, bad = 0, 0, []
for k, (s7, e7, y70, y71) in enumerate(order):
    y = 2019 + k; lim = order[1][0] - 150 if s7 < img.width // 2 else img.width
    t = ocr(img, (s7 - 120, y70 - 330, s7 + 1500, y70 - 10), psm=7); okt += t == f'({"abcdef"[k]}) {y}'
    blank = not (np.abs(a[y70 - 330:y70 - 10, s7 + 1500:lim] - 255).sum(2) > 60).any()     # tesseract invents text on blank crops
    n = '' if blank else norm(ocr(img, (s7 + 1500, y70 - 330, lim, y70 - 10), psm=7)).replace('≤', '<'); okn += n == NOTE.get(y, '')
    if t != f'({"abcdef"[k]}) {y}' or n != NOTE.get(y, ''): bad.append((t, n))
check('FS3 six panel titles', okt == 6, f'{okt}/6' + (f'; {bad}' if bad else ''))
check('FS3 notes on 2019-2021 only, text as intended', okn == 6, f'{okn}/6')
oky = 0
for s7, e7, y70, y71 in [t for t in order if t[0] < img.width // 2]:
    tys = ticks_y(a, s7, y70, y71)
    # lone glyphs ('1', '0') are read in single-character mode; '0.5' in line mode (a dropped decimal point, '05', is
    # accepted as '0.5' here; the exact label strings are checked in scripts/63)
    got = [ocr(img, (s7 - 300, t - 110, s7 - 70, t + 110), psm=10, whitelist='0123456789.') if i != 1 else
           ocr(img, (s7 - 520, t - 100, s7 - 70, t + 100), psm=7, whitelist='0123456789.').replace('05', '0.5') for i, t in enumerate(tys)]
    oky += got == ['1', '0.5', '0']
    got_t = ocr(img, (0, y70 - 200, s7 - 330, y71 + 200), psm=7, rotate=-90)
    check(f'FS3 y axis title (panel row at y {y70})', got_t == 'Share of annual maximum', got_t)
check('FS3 y tick labels on the three left panels', oky == 3, f'{oky}/3 exact')
for s7, e7, y70, y71 in order[-2:]:
    yb7 = bottom_row(a, s7, y71); tx = ticks_x(a, s7 + 12, s7 + 3800, yb7)
    got = [ocr(img, (t - 150, yb7 + 60, t + 150, yb7 + 260), psm=10, whitelist='JFMASOND') for t in tx]
    check(f'FS3 month tick labels (panel at x {s7})', got == list('JFMAMJJASOND'), ''.join(got))
    got = ocr(img, (s7, yb7 + 250, s7 + 3750, yb7 + 520), psm=7)
    check(f'FS3 x axis title (panel at x {s7})', got == 'Month', got)
s7, e7, y70, y71 = order[2]
got = ocr(img, (s7 + 300, y70 + 40, s7 + 3600, y70 + 230), psm=7)
check('FS3 2021 strip text', got == 'No national alert found', got)
top = a[:order[0][2] - 380]; r = np.where((np.abs(top - [255, 255, 255]).sum(2) > 60).any(1))[0]; br = np.where(np.diff(r) > 1)[0]
got = ' | '.join(ocr(img, (0, u - 30, img.width, v + 30), psm=7) for u, v in zip(np.r_[r[0], r[br + 1]], np.r_[r[br], r[-1]]))
check('FS3 legend (all five entries)', all(t in got for t in ('Cases', 'Search interest', 'National alert or EDCD statement', 'Count announced',
                                                            'Daily reports began / WHO report')), got)
print(f'RESULT: {len(fails)} failed check(s)')
sys.exit(1 if fails else 0)
