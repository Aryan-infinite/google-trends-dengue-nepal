"""Compare implementations under the agreement rules fixed in VERIFICATION_SPEC.md Section 9. Usage: python3 compare.py results_*.json"""
import json, sys, math, csv
files = sys.argv[1:]; I = {f.replace('results_', '').replace('.json', '').split('/')[-1]: json.load(open(f)) for f in files}
names = list(I); keys = sorted(set().union(*[set(v) for v in I.values()]))
EXACT = ('_ar_order', 'seasons_m2_better', 'rule_met', 'zeros2020_2023', 'peak_n_changed', 'peak_ties', 'cook_months', '_ccf_n_')
SHARED = ('_blo', '_bhi', '_diff_lo', '_diff_hi', '_diff_b3lo', '_diff_b3hi')
rows, fails, missing = [], [], []
for k in keys:
    if k.endswith('_own'): continue
    have = {n: I[n][k] for n in names if k in I[n]}
    if len(have) < len(names): missing.append((k, sorted(set(names) - set(have))))
    if len(have) < 2: continue
    vals = list(have.values()); ref = vals[0]
    if any(t in k for t in EXACT):
        ok = all(v == ref for v in vals); rule = 'exact'
    elif k in ('nb_rr_lo', 'nb_rr_hi'):
        ok = max(vals) - min(vals) <= 0.05 and (all(v > 1 for v in vals) if k == 'nb_rr_lo' else True); rule = '0.05, excl. 1'
    elif k.endswith(SHARED):
        ok = max(vals) - min(vals) <= 1e-4; rule = 'shared resamples 1e-4'
    else:
        ok = max(vals) - min(vals) <= 1e-4 * max(1, abs(ref)); rule = '1e-4 rel'
    spread = 0 if isinstance(ref, list) else max(vals) - min(vals)
    rows.append(dict(key=k, rule=rule, ok=ok, spread=spread, **{n: have.get(n) for n in names}))
    if not ok: fails.append(rows[-1])
# pipeline own-RNG bootstraps vs shared-resample values (Monte Carlo tolerance)
mc = []
for n in names:
    for k, v in I[n].items():
        if not k.endswith('_own'): continue
        base = k[:-4]; others = [I[m][base] for m in names if base in I[m]]
        if not others: continue
        ref = sum(others) / len(others); tol = 0.03 if base.startswith('q1_') else 0.15 * abs(ref)
        mc.append(dict(key=k, own=v, shared=ref, tol=tol, ok=abs(v - ref) <= tol))
with open('comparison_table.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['key', 'rule', 'ok', 'spread'] + names); w.writeheader()
    for r in rows: w.writerow({k: (json.dumps(v) if isinstance(v, list) else v) for k, v in r.items()})
print(f'implementations: {names}')
print(f'keys compared: {len(rows)}   agree: {sum(r["ok"] for r in rows)}   disagree: {len(fails)}')
print(f'pipeline own-RNG bootstrap vs shared resamples: {sum(m["ok"] for m in mc)}/{len(mc)} within Monte Carlo tolerance')
for m in mc:
    if not m['ok']: print('  MC FAIL', m)
for r in fails: print('  FAIL', r['key'], r['rule'], {n: r[n] for n in names})
if missing: print('  keys missing in some implementations:', missing[:40], '...' if len(missing) > 40 else '')
