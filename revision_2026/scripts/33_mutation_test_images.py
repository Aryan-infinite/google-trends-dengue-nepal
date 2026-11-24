"""Mutation test for the two image digitisers (scripts/31 Python, scripts/32 R).
Each mutant is a falsified Figure 1 or 2 TIFF: either re-rendered from subtly altered data, or the delivered TIFF
edited directly (stray dot, erased marker, erased stretch of line, leader pointing at the wrong month).
Both digitisers then judge the falsified image against the TRUE raw data; each must report a failure.
Runs in a temporary copy; nothing in the project is modified. Usage: python3 -I scripts/33_mutation_test_images.py
"""
import os, shutil, subprocess, tempfile, json
from PIL import Image, ImageDraw
Image.MAX_IMAGE_PIXELS = None
SRC = os.getcwd()
T = tempfile.mkdtemp(prefix='imgmut_')
for d in ('data_raw', 'scripts'): shutil.copytree(os.path.join(SRC, d), os.path.join(T, d))
os.makedirs(os.path.join(T, 'out'), exist_ok=True); os.makedirs(os.path.join(T, 'manuscript/figures'), exist_ok=True)
TRUE = {f: os.path.join(SRC, 'manuscript/figures', f) for f in ('Fig1.tif', 'Fig2.tif')}
coords = json.load(open(os.path.join(SRC, 'out/31_digitise_python.json')))
# The preprint table (PREPRINT_XLSX, as in the other scripts) is copied into the temporary folder; mutants edit the copy.
PRE = os.path.join(SRC, os.environ.get('PREPRINT_XLSX', 'data_raw/paper_finalized_data_RP.xlsx'))
PRE_T = os.path.join(T, 'preprint_copy.xlsx')
ENV = {**os.environ, 'PREPRINT_XLSX': PRE_T}
run = lambda cmd: subprocess.run(cmd, cwd=T, capture_output=True, text=True, env=ENV)

def reset():
    for f, p in TRUE.items(): shutil.copy(p, os.path.join(T, 'manuscript/figures', f))
    for f in os.listdir(os.path.join(SRC, 'data_raw')): shutil.copy(os.path.join(SRC, 'data_raw', f), os.path.join(T, 'data_raw', f))
    shutil.copy(PRE, PRE_T)
    shutil.copy(os.path.join(SRC, 'scripts/28_figure2.py'), os.path.join(T, 'scripts/28_figure2.py'))

def rerender(fig, edit):
    """Alter inputs, re-render one figure, then restore the TRUE inputs (the digitisers read the raw files)."""
    edit()
    assert run(['python3', '-I', 'scripts/00_prepare_inputs.py']).returncode == 0
    r = run(['python3', '-I', f'scripts/{"25_figure1" if fig == 1 else "28_figure2"}.py']); assert r.returncode == 0, r.stderr
    keep = os.path.join(T, f'manuscript/figures/Fig{fig}.tif'); shutil.copy(keep, keep + '.mut')
    reset(); shutil.move(keep + '.mut', keep)

def gt_edit(fname, month, new):
    def f():
        p = os.path.join(T, 'data_raw', fname); L = open(p).read().splitlines()
        L = [f'{month},{new}' if l.startswith(month + ',') else l for l in L]; open(p, 'w').write('\n'.join(L) + '\n')
    return f

def xlsx_edit(month_index, col_name, new):
    def f():
        import openpyxl
        p = PRE_T; wb = openpyxl.load_workbook(p); ws = wb.active
        rows = list(ws.iter_rows()); h = next(i for i, r in enumerate(rows) if any(c.value == 'Confirmed Cases' for c in r))
        col = [c.value for c in rows[h]].index(col_name)
        data = [r for r in rows[h + 1:] if r[col].value is not None]; data[month_index][col].value = new; wb.save(p)
    return f

def paint(fig, shape, colour_from=None, fill=None):
    p = os.path.join(T, f'manuscript/figures/Fig{fig}.tif'); im = Image.open(p); im.load()
    if colour_from: fill = im.getpixel(colour_from)
    ImageDraw.Draw(im).ellipse(shape, fill=fill); im.save(p, compression='tiff_lzw', dpi=(1200, 1200))

c1, c2 = coords['fig1_px'], coords['fig2_px']
def px1(panel, m): return [int(round(v)) for v in c1[panel][m]]
def px2(panel, m): return [int(round(v)) for v in c2[panel][m]]
def disc(c, r): return [c[0] - r, c[1] - r, c[0] + r, c[1] + r]

mutants = {
 # Figure 1: values re-rendered from altered data
 'F1b: Aug 2023 search +1/3 point (one retrieval 19->20)': lambda: rerender(1, gt_edit('gt_dengue_term.csv', '2023-08', 20)),
 'F1a: Mar 2023 cases +2.6% (195->200)': lambda: rerender(1, xlsx_edit(50, 'Confirmed Cases', 200)),
 'F1c: Jan 2022 Fever +1 (11->12)': lambda: rerender(1, xlsx_edit(36, 'Fever', 12)),
 # Figure 1: image edits
 'F1b: stray blue dot (r=10 px) in empty area': lambda: paint(1, disc([px1('b', '2021-06')[0], px1('b', '2022-09')[1] + 600], 10), colour_from=tuple(px1('b', '2022-09'))),
 'F1a: tiny stray black dot (r=3 px)': lambda: paint(1, disc([px1('a', '2022-04')[0], px1('a', '2022-09')[1] + 120], 3), colour_from=tuple(px1('a', '2022-09'))),
 'F1b: erase marker and line around Aug 2023 (r=25 px)': lambda: paint(1, disc(px1('b', '2023-08'), 25), fill=(0, 0, 0, 0)),
 # Figure 2: values re-rendered from altered data
 'F2: Oct 2024 cases +2.2% (12332->12600)': lambda: rerender(2, xlsx_edit(69, 'Confirmed Cases', 12600)),
 'F2: Aug 2023 search +1/3 point': lambda: rerender(2, gt_edit('gt_dengue_term.csv', '2023-08', 20)),
 # Figure 2: image edits and a wrong leader
 'F2a: stray orange dot (r=12 px) in empty area': lambda: paint(2, disc([px2('a', '2020-01')[0] + 200, px2('a', '2022-09')[1] + 300], 12), colour_from=(int(c2['a']['2022-09'][0]), int(c2['a']['2022-09'][1]))),
 'F2b: erase the Sep 2022 marker': lambda: paint(2, disc(px2('b', '2022-09'), 34), fill=(0, 0, 0, 0)),
 'F2b: erase the Aug 2023 marker': lambda: paint(2, disc(px2('b', '2023-08'), 34), fill=(0, 0, 0, 0)),
 'F2b: leader points at Oct 2022 instead of Sep 2022': lambda: rerender(2, lambda: open(os.path.join(T, 'scripts/28_figure2.py'), 'w').write(
     open(os.path.join(SRC, 'scripts/28_figure2.py')).read().replace("flagged = list(d.index[cook > 4 / len(d)].strftime('%Y-%m'))", "flagged = ['2020-06', '2022-10']"))),
}
rows = []
for name, make in mutants.items():
    reset(); make()
    py = run(['python3', '-I', 'scripts/31_digitise_figures.py']); r = run(['Rscript', 'scripts/32_digitise_figures.R'])
    fp = [l.split('] ')[1].split('  (')[0] for l in py.stdout.splitlines() if l.startswith('[FAIL]')]
    fr = [l.split('] ')[1].split('  (')[0] for l in r.stdout.splitlines() if l.startswith('[FAIL]')]
    rows.append((name, py.returncode != 0, r.returncode != 0))
    print(f"{name}\n   Python: {'CAUGHT' if py.returncode else 'MISSED'} {fp[:3]}\n   R:      {'CAUGHT' if r.returncode else 'MISSED'} {fr[:3]}", flush=True)
reset()
py = run(['python3', '-I', 'scripts/31_digitise_figures.py']); r = run(['Rscript', 'scripts/32_digitise_figures.R'])
print(f"\nPython caught {sum(x[1] for x in rows)}/{len(rows)}; R caught {sum(x[2] for x in rows)}/{len(rows)}; "
      f"either caught {sum(x[1] or x[2] for x in rows)}/{len(rows)}")
print('true images still pass:', py.returncode == 0 and r.returncode == 0)
shutil.rmtree(T)
