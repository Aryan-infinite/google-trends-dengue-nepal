"""Second check of manuscript/Tables_4-5.docx (Python; python-docx reader). Every data cell must equal the value
formatted here from the NumPy implementations (scripts/69, 71), the pipeline's bootstrap limits and the shared-resample
block limits; decision-rule labels recomputed here. Independent of scripts/72-73 (no shared code) and of scripts/74 (R)."""
import csv, os, sys
import openpyxl
from docx import Document
fails = []
def check(n, ok, d=''): print(f"[{'PASS' if ok else 'FAIL'}] {n}" + (f'  ({d})' if d else '')); ok or fails.append(n)
M = '−'
def f(v, d): return format(round(v, d) + 0.0, f',.{d}f').replace('-', M) if d else format(int(round(v)), ',').replace('-', M)
def p(v): return '<0.001' if v < 0.001 else ('%.3f' % v if v < 0.01 else '%.2f' % v)
doc = Document(os.environ.get('TABLES45_DOCX', 'manuscript/Tables_4-5.docx'))
check('two tables', len(doc.tables) == 2)
def rows(t):
    out = []
    for r in t.rows:
        seen, cells = set(), []
        for c in r.cells:                                     # python-docx repeats a merged cell once per grid column
            if id(c._tc) in seen: continue
            seen.add(id(c._tc)); cells.append(c.text.strip())
        out.append(cells)
    return out
N4 = list(csv.DictReader(open('out/table4_runs_np.csv'))); fl = lambda r, k: float(r[k])
wb = openpyxl.load_workbook('out/final/final_results.xlsx', read_only=True, data_only=True)
grid = {r[0]: r for r in wb['Robustness_grid'].iter_rows(min_row=2, values_only=True)}
fever_hi = grid['Control: "Fever"'][4]
L4 = ['Primary analysis, January 2020–December 2024', '(i) Including June–December 2019', '(ii) All 72 months, as in the preprint', '(iii) 2022–2024 only',
      '(iv) Excluding 2022', '(v) Preprint Topic series', '(vi) First retrieval only', '(vii) Values "<1" coded 0', 'Control series "Fever"']
want = []
for lab, r in zip(L4, N4):
    hi = fever_hi if r['run'] == 'fever' else fl(r, 'ci_hi')
    met = fl(r, 'holm_p1') < 0.05 and fl(r, 'lag_p1') > fl(r, 'lag_m1') and fl(r, 'fwd_p') < 0.05
    want.append([lab, r['n'], f'{f(fl(r, "r"), 2)} ({f(fl(r, "ci_lo"), 2)} to {f(hi, 2)})', f(fl(r, 'n_eff'), 1),
                 f'{f(fl(r, "lag_m1"), 2)} ({p(fl(r, "holm_m1"))})', f'{f(fl(r, "lag_0"), 2)} ({p(fl(r, "holm_0"))})', f'{f(fl(r, "lag_p1"), 2)} ({p(fl(r, "holm_p1"))})',
                 f'{f(fl(r, "fwd_coef"), 2)} ({p(fl(r, "fwd_p"))})', 'Yes' if met else 'No'])
got = [r for r in rows(doc.tables[0]) if len(r) == 9][1:]
bad = [i for i, (g, w) in enumerate(zip(got, want)) if g != w]
check('Table 4: every data cell = NumPy value, formatted here', len(got) == 9 and not bad, f'rows differing: {bad}' if bad else '')
N5 = {r['key']: float(r['value']) for r in csv.DictReader(open('out/table5_np.csv'))}
rb = {r[0]: r for r in wb['Q3_robustness'].iter_rows(min_row=2, values_only=True)}
mo = {r[0]: r for r in openpyxl.load_workbook('out/15_nowcast_monthly_refit.xlsx', read_only=True)['robustness'].iter_rows(min_row=2, values_only=True)}
sm = {r[0]: r[1] for r in wb['Summary'].iter_rows(min_row=2, values_only=True)}
SER = [('dengue', 'dengue 3-day', '"dengue" (primary series)'), ('topic', 'Topic', 'Preprint Topic series'), ('day1', 'day-1', 'First retrieval only'), ('fever', 'Fever control', 'Control series "Fever"')]
want5 = []
for sch, src in (('ps', rb), ('mo', mo), ('loso', None)):
    for sk, key, lab in (SER if sch != 'loso' else SER[:1]):
        if src is None: lo, hi = sm['Q3_loso_pairdiff_M2_minus_M1.lo'], sm['Q3_loso_pairdiff_M2_minus_M1.hi']
        else: r = src[key]; lo, hi = (r[5], r[6]) if sch == 'ps' else (r[6], r[7])
        b, ns = int(N5[f'{sk}_{sch}_better']), int(N5[f'{sk}_{sch}_nseasons'])
        rule = 'No rule' if sch == 'loso' else ('Met' if b >= 2 and hi < 0 else 'Not met')
        want5.append([lab, str(int(N5[f'{sk}_{sch}_ntest']))] + [f(N5[f'{sk}_{sch}_mae_{m}'], 0) for m in ('M0', 'M1', 'M2', 'M3')] +
                     [f'{f(N5[f"{sk}_{sch}_diff"], 0)} ({f(lo, 0)} to {f(hi, 0)})', f'{b} of {ns}', rule])
got5 = [r for r in rows(doc.tables[1]) if len(r) == 9][1:]
bad = [i for i, (g, w) in enumerate(zip(got5, want5)) if g != w]
check('Table 5: every data cell = NumPy value, formatted here', len(got5) == 9 and not bad, f'rows differing: {bad}' if bad else '')
CT = {r['key']: r for r in csv.DictReader(open('verification/comparison_table.csv'))}
blk = {s: f"{f(float(CT[f'q3_{s}_diff_b3lo']['impl3_numpy']), 0)} to {f(float(CT[f'q3_{s}_diff_b3hi']['impl3_numpy']), 0)}" for s in ('ps', 'mo')}
notes = '\n'.join(pp.text for pp in doc.paragraphs)
check('Table 5 note: block-bootstrap intervals (NumPy shared-resample values)', f"were {blk['ps']} (once per season) and {blk['mo']} (every month)" in notes, str(blk))
check('Table 4 note: run (iv) coefficient', f"(−0.25; P = 0.04)" in notes and f(float(N4[4]['fwd_coef']), 2) == '−0.25' and p(float(N4[4]['fwd_p'])) == '0.04')
print(f'RESULT: {len(fails)} failed check(s)'); sys.exit(1 if fails else 0)
