"""Table 5, third implementation (NumPy; least squares by numpy.linalg.lstsq; csv module for input).
Same keys as scripts/70 (R). Output: out/table5_np.csv"""
import csv
import numpy as np
def rd(f): return np.array([0.5 if r[1].strip() == '<1' else float(r[1]) for r in csv.reader(open('data_raw/' + f, encoding='utf-8-sig')) if len(r) == 2 and len(r[0]) == 7 and r[0][:2] == '20'])
F3 = ['gt_dengue_term.csv', 'gt_dengue_term_day2_2026-11-22.csv', 'gt_dengue_term_day3_2026-11-23.csv']
T = list(csv.DictReader(open('out/paper_72month_table.csv', encoding='utf-8')))
cases = np.array([float(r['cases']) for r in T])
SER = {'dengue': sum(rd(f) for f in F3) / 3, 'topic': np.array([float(r['DF']) for r in T]), 'day1': rd(F3[0]), 'fever': np.array([float(r['Fever']) for r in T])}
YR = np.repeat(np.arange(2019, 2025), 12); MO = np.tile(np.arange(1, 13), 6); YM = YR * 100 + MO
idx = lambda lo, hi: [i for i in range(72) if lo <= YM[i] <= hi]
LY = np.log(cases + 1); K = {}
for sn, s in SER.items():
    LX = np.log(s + 0.5)
    def des(r, M):
        r = np.asarray(r); cols = [(MO[r] == k).astype(float) for k in range(1, 13)] + [LY[r - 1]]
        if M == 'M2': cols.append(LX[r])
        if M == 'M3': cols.append(LX[r - 1])
        return np.column_stack(cols)
    def pred(tr, te, M):
        b = np.linalg.lstsq(des(tr, M), LY[np.asarray(tr)], rcond=None)[0]; return np.maximum(np.expm1(des(te, M) @ b), 0)
    ii, tm, lo = idx(201907, 202412), idx(202201, 202412), idx(202002, 202412)
    SCH = {'ps': ([[t for t in tm if YR[t] == y] for y in (2022, 2023, 2024)], lambda g: [i for i in ii if i < min(g)]),
           'mo': ([[t] for t in tm], lambda g: [i for i in ii if i < g[0]]),
           'loso': ([[t for t in lo if YR[t] == y] for y in range(2020, 2025)], lambda g: [i for i in lo if YR[i] != YR[g[0]]])}
    for sc, (groups, trainer) in SCH.items():
        P = {m: [] for m in ('M0', 'M1', 'M2', 'M3')}; obs, yrs = [], []
        for g in groups:
            tr = trainer(g); P['M0'] += list(cases[np.asarray(g) - 1]); obs += list(cases[g]); yrs += list(YR[g])
            for M in ('M1', 'M2', 'M3'): P[M] += list(pred(tr, g, M))
        obs, yrs = np.array(obs), np.array(yrs)
        for M, p in P.items():
            p = np.array(p); e = p - obs
            K[f'{sn}_{sc}_mae_{M}'] = np.mean(np.abs(e)); K[f'{sn}_{sc}_rmse_{M}'] = np.sqrt(np.mean(e ** 2))
            K[f'{sn}_{sc}_maelog_{M}'] = np.mean(np.abs(np.log1p(p) - np.log1p(obs)))
            for y in sorted(set(yrs)): K[f'{sn}_{sc}_mae{y}_{M}'] = np.mean(np.abs(e[yrs == y]))
            for thr in (1000, 500):                                                    # 500: sensitivity threshold (added 2026-11-24)
                o, f = obs > thr, p > thr
                K[f'{sn}_{sc}_sens{thr}_{M}'] = (o & f).sum() / o.sum() if o.sum() else np.nan
                K[f'{sn}_{sc}_spec{thr}_{M}'] = (~o & ~f).sum() / (~o).sum(); K[f'{sn}_{sc}_nepi{thr}'] = int(o.sum())
        a1, a2 = np.abs(np.array(P['M1']) - obs), np.abs(np.array(P['M2']) - obs)
        K[f'{sn}_{sc}_diff'] = np.mean(a2 - a1); K[f'{sn}_{sc}_better'] = sum(a2[yrs == y].mean() < a1[yrs == y].mean() for y in set(yrs))
        K[f'{sn}_{sc}_nseasons'] = len(set(yrs)); K[f'{sn}_{sc}_ntest'] = len(obs)
with open('out/table5_np.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['key', 'value'])
    for k, v in K.items(): w.writerow([k, repr(float(v))])
print(len(K), 'keys')
