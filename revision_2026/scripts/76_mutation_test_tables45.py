"""Mutation test for the Word-table checks: each falsified copy of manuscript/Tables_4-5.docx (one text change in
word/document.xml) must fail BOTH scripts/74 (R) and scripts/75 (Python). Runs on temporary copies."""
import os, re, subprocess, tempfile, zipfile, shutil
SRC = 'manuscript/Tables_4-5.docx'; X = zipfile.ZipFile(SRC).read('word/document.xml').decode('utf-8')
def mutate(a, b, nth=1):
    i = -1
    for _ in range(nth): i = X.index(a, i + 1)
    return X[:i] + b + X[i + len(a):]
M = [('Table 5: M2 MAE 1,900 -> 1,901', mutate('>1,900<', '>1,901<')),
     ('Table 4: primary CI upper 0.96 -> 0.97', mutate('0.81 (0.24 to 0.96)', '0.81 (0.24 to 0.97)')),
     ('Table 5: primary rule Not met -> Met', mutate('>Not met<', '>Met<')),
     ('Table 5 note: block interval 378 -> 387', mutate('−787 to 378', '−787 to 387')),
     ('Table 4: primary lag +1 Holm 0.89 -> 0.05', mutate('0.16 (0.89)', '0.16 (0.05)')),
     ('Table 4: run (iii) months 36 -> 35', mutate('>36<', '>35<')),
     ('Table 5: leave-one-season-out 4 of 5 -> 5 of 5', mutate('>4 of 5<', '>5 of 5<')),
     ('Table 4: run (iv) note coefficient -0.25 -> -0.52', mutate('−0.25; ', '−0.52; '))]
T = tempfile.mkdtemp(prefix='t45mut_'); caught = []
for name, xml in M:
    assert xml != X
    out = os.path.join(T, 'mut.docx')
    with zipfile.ZipFile(SRC) as zi, zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as zo:
        for it in zi.infolist(): zo.writestr(it, xml.encode('utf-8') if it.filename == 'word/document.xml' else zi.read(it.filename))
    env = {**os.environ, 'TABLES45_DOCX': out}
    r1 = subprocess.run(['Rscript', 'scripts/74_check_tables45.R'], env=env, capture_output=True, text=True).returncode != 0
    r2 = subprocess.run(['python3', '-I', 'scripts/75_check_tables45_docx.py'], env=env, capture_output=True, text=True).returncode != 0
    caught.append(r1 and r2); print(f"  {name:52s} {'CAUGHT' if r1 and r2 else 'MISSED'}  R74: {'fail' if r1 else 'pass'}; Py75: {'fail' if r2 else 'pass'}")
ok = [subprocess.run(c, capture_output=True).returncode == 0 for c in (['Rscript', 'scripts/74_check_tables45.R'], ['python3', '-I', 'scripts/75_check_tables45_docx.py'])]
print(f'caught {sum(caught)}/{len(caught)} by both checks; true document passes both: {all(ok)}'); shutil.rmtree(T)
