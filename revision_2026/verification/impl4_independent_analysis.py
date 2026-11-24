"""Independent implementation of VERIFICATION_SPEC.md (written from scratch).
OLS, HAC, NB2 ML, Ljung-Box, Cook's distance and bootstraps implemented with numpy/scipy.
"""
import json, re, calendar, math
import numpy as np
import openpyxl
from scipy import stats, optimize, special

RAW = "data_raw/"   # run from the project root
VER = "verification/"   # run from the project root
OUT = "verification/results_impl4_independent.json"
CHECKS = []

def check(name, cond):
    CHECKS.append((name, bool(cond)))
    if not cond:
        print("CHECK FAILED:", name)

MONTHS = [f"{y}-{m:02d}" for y in range(2019, 2025) for m in range(1, 13)]   # t=0..71
TI = {k: i for i, k in enumerate(MONTHS)}
CAL = np.array([int(k[5:]) for k in MONTHS])
YR = np.array([int(k[:4]) for k in MONTHS])
MNAME = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

def win(a, b):
    return list(range(TI[a], TI[b] + 1))

# ---------------------------------------------------------------- data
ZONE_REGION = {}
EASTERN = "TAPLEJUNG PANCHTHAR ILAM JHAPA SANKHUWASABHA TERHATHUM DHANKUTA BHOJPUR MORANG SUNSARI " \
          "SOLUKHUMBU OKHALDHUNGA KHOTANG UDAYAPUR SAPTARI SIRAHA".split()
CENTRAL = "DOLAKHA RAMECHHAP SINDHULI DHANUSA MAHOTTARI SARLAHI SINDHUPALCHOK RASUWA NUWAKOT DHADING " \
          "KATHMANDU BHAKTAPUR LALITPUR KAVREPALANCHOK MAKWANPUR CHITAWAN RAUTAHAT BARA PARSA".split()
WESTERN = "GORKHA LAMJUNG TANAHU SYANGJA KASKI MANANG MUSTANG MYAGDI PARBAT BAGLUNG GULMI ARGHAKHANCHI " \
          "PALPA RUPANDEHI KAPILBASTU".split() + ["NAWALPARASI EAST", "NAWALPARASI WEST"]
MIDWEST = "ROLPA PYUTHAN DANG SALYAN DAILEKH JAJARKOT SURKHET BANKE BARDIYA DOLPA MUGU HUMLA JUMLA " \
          "KALIKOT".split() + ["RUKUM EAST", "RUKUM WEST"]
FARWEST = "BAJURA BAJHANG DOTI ACHHAM KAILALI DARCHULA BAITADI DADELDHURA KANCHANPUR".split()
REGION_LISTS = {"Eastern": EASTERN, "Central": CENTRAL, "Western": WESTERN, "MidWestern": MIDWEST, "FarWestern": FARWEST}

def load_edcd():
    wb = openpyxl.load_workbook(RAW + "edcd_district_monthly_2019_2023.xlsx", data_only=True)
    nat = np.zeros(60)
    reg = {r: np.zeros(60) for r in REGION_LISTS}
    district_names = None
    for yi, y in enumerate(range(2019, 2024)):
        ws = wb[f"{y} Monthly"]
        rows = list(ws.iter_rows(values_only=True))
        total = None
        s = np.zeros(12)
        names = []
        for r in rows[2:]:
            name = r[1]
            if name is None:
                continue
            vals = np.array([0.0 if v is None else float(v) for v in r[3:15]])
            nm = str(name).strip()
            if nm.upper() == "TOTAL":
                total = vals
                continue
            s += vals
            if nm.upper() == "NOT SPECIFIED":
                continue
            dname = re.sub(r"^\d+\s+", "", nm).upper()
            names.append(dname)
            hits = [rg for rg, lst in REGION_LISTS.items() if dname in lst]
            check(f"{y} district {dname} assigned to exactly one region", len(hits) == 1)
            reg[hits[0]][yi * 12:(yi + 1) * 12] += vals
        check(f"{y} 77 districts", len(names) == 77 and len(set(names)) == 77)
        check(f"{y} district sum incl. Not specified == Total row each month", total is not None and np.allclose(s, total))
        if district_names is None:
            district_names = set(names)
        check(f"{y} same district set", set(names) == district_names)
        nat[yi * 12:(yi + 1) * 12] = s
    allreg = set(sum(REGION_LISTS.values(), []))
    check("region lists cover exactly the 77 districts", allreg == district_names and sum(map(len, REGION_LISTS.values())) == 77)
    return nat, reg

def load_paper():
    ws = openpyxl.load_workbook(RAW + "paper_finalized_data_RP.xlsx", data_only=True).worksheets[0]
    rows = list(ws.iter_rows(values_only=True))
    # locate header row (spec says row 1; actual file has a blank row 1 and header in row 2)
    hdr = [i for i, r in enumerate(rows) if r[1] == "Months"][0]
    data = {}
    for r in rows[hdr + 1:]:
        if isinstance(r[1], (int, float)):
            data[int(r[1])] = [float(r[2]), float(r[3]), float(r[4]), float(r[5])]
    check("paper file has month indices 1..72", sorted(data) == list(range(1, 73)))
    arr = np.array([data[i] for i in range(1, 73)])
    return arr[:, 0], arr[:, 1], arr[:, 3], hdr

def load_gt(fn):
    out = {}
    with open(RAW + fn) as f:
        for line in f:
            m = re.match(r"^(\d{4}-\d{2}),(.*)$", line.strip())
            if m:
                v = m.group(2).strip()
                out[m.group(1)] = 0.5 if v == "<1" else float(v)
    check(f"{fn} covers 2019-01..2024-12", all(k in out for k in MONTHS))
    return np.array([out[k] for k in MONTHS])

nat5, regcases = load_edcd()
topic, fever, paper_cases, paper_hdr = load_paper()
check("paper cases indices 1-60 == EDCD national 2019-2023", np.array_equal(paper_cases[:60], nat5))
cases = np.concatenate([nat5, paper_cases[60:72]])
day1 = load_gt("gt_dengue_term.csv")
day2 = load_gt("gt_dengue_term_day2_2026-11-22.csv")
day3 = load_gt("gt_dengue_term_day3_2026-11-23.csv")
x = (day1 + day2 + day3) / 3.0
regsearch = {r: load_gt(f"gt_dengue_region_{r.lower()}.csv") for r in ["Central", "Eastern", "Western"]}

ly = np.log(cases + 1.0)
def lxf(v):
    return np.log(v + 0.5)
lx = lxf(x)

def adj(v, idx):
    """v[idx] minus calendar-month mean over idx."""
    idx = np.asarray(idx)
    vv = v[idx].astype(float)
    cm = CAL[idx]
    out = vv.copy()
    for m in range(1, 13):
        s = cm == m
        if s.any():
            out[s] = vv[s] - vv[s].mean()
    return out

def pearson(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a - a.mean(); b = b - b.mean()
    return float((a @ b) / math.sqrt((a @ a) * (b @ b)))

def rho(v, k):
    n = len(v)
    return pearson(v[:n - k], v[k:])

def corr_block(a, b):
    n = len(a)
    r = pearson(a, b)
    ssum = sum(rho(a, k) * rho(b, k) for k in range(1, 13))
    neff = max(4.0, n / (1 + 2 * ssum))
    z = math.atanh(r); h = 1.96 / math.sqrt(neff - 3)
    return r, neff, math.tanh(z - h), math.tanh(z + h)

R = {}

# ---------------------------------------------------------------- Q1
runs = {
    "primary": (x, win("2020-01", "2024-12")),
    "A": (x, win("2019-06", "2024-12")),
    "B": (x, win("2019-01", "2024-12")),
    "C": (x, win("2022-01", "2024-12")),
    "D": (x, [i for i in range(72) if YR[i] in (2020, 2021, 2023, 2024)]),
    "topic": (topic, win("2020-01", "2024-12")),
    "day1": (day1, win("2020-01", "2024-12")),
    "fever": (fever, win("2020-01", "2024-12")),
}
check("run D has 48 months", len(runs["D"][1]) == 48)
q1ab = {}
for name, (ser, idx) in runs.items():
    a = adj(lxf(ser), idx); b = adj(ly, idx)
    q1ab[name] = (a, b)
    r, ne, lo, hi = corr_block(a, b)
    R[f"q1_{name}_r"] = r; R[f"q1_{name}_neff"] = ne; R[f"q1_{name}_flo"] = lo; R[f"q1_{name}_fhi"] = hi

starts = np.loadtxt(VER + "boot_q1_block_starts.csv", delimiter=",", dtype=int)
check("q1 block starts within 0..48", starts.min() >= 0 and starts.max() <= 48 and starts.shape[1] == 5)
for name in ["primary", "fever"]:
    a, b = q1ab[name]
    rs = []
    for row in starts:
        idx = np.concatenate([np.arange(s, s + 12) for s in row])
        rs.append(pearson(a[idx], b[idx]))
    R[f"q1_{name}_blo"] = float(np.percentile(rs, 2.5)); R[f"q1_{name}_bhi"] = float(np.percentile(rs, 97.5))

# ---------------------------------------------------------------- OLS helpers
def ols(X, y):
    XtX_inv = np.linalg.inv(X.T @ X)
    beta = XtX_inv @ X.T @ y
    resid = y - X @ beta
    return beta, resid, XtX_inv

def hac_cov(X, u, XtX_inv, L=3):
    g = X * u[:, None]
    S = g.T @ g
    for j in range(1, L + 1):
        w = 1 - j / (L + 1)
        G = g[j:].T @ g[:-j]
        S += w * (G + G.T)
    return XtX_inv @ S @ XtX_inv

def month_dummies(idx):
    return np.array([[1.0 if CAL[i] == m else 0.0 for m in range(1, 13)] for i in idx])

# ---------------------------------------------------------------- Q2
P = win("2020-01", "2024-12")
a2 = adj(lx, P); b2 = adj(ly, P)
n = 60
def ar_design(v, p, t0):
    # t0 = first 0-based t used; rows t = t0..59
    cols = [np.ones(n - t0)] + [v[t0 - j:n - j] for j in range(1, p + 1)]
    return np.column_stack(cols), v[t0:]
aic = {}
for p in range(4):
    X, y = ar_design(a2, p, 3)        # t = 4..60 in 1-based
    beta, res, _ = ols(X, y)
    rss = res @ res
    aic[p] = 57 * math.log(rss / 57) + 2 * (p + 1)
pbest = min(aic, key=aic.get)
X, y = ar_design(a2, pbest, pbest)
arb, _, _ = ols(X, y)
R["q2_ar_order"] = int(pbest); R["q2_ar_const"] = float(arb[0])
for j in range(1, 4):
    if j <= pbest:
        R[f"q2_ar_phi{j}"] = float(arb[j])
AIC_TABLE = aic

def prewhiten(v):
    return np.array([v[t] - arb[0] - sum(arb[j] * v[t - j] for j in range(1, pbest + 1)) for t in range(pbest, n)])
ea = prewhiten(a2); eb = prewhiten(b2)
m = len(ea)
klab = {-3: "m3", -2: "m2", -1: "m1", 0: "0", 1: "p1", 2: "p2", 3: "p3"}
pr = {}
for k in range(-3, 4):
    if k >= 0:
        u, w = ea[:m - k], eb[k:]
    else:
        u, w = ea[-k:], eb[:m + k]
    r = pearson(u, w); nk = len(u)
    tstat = r * math.sqrt((nk - 2) / (1 - r * r))
    pv = 2 * stats.t.sf(abs(tstat), nk - 2)
    R[f"q2_ccf_r_{klab[k]}"] = r; R[f"q2_ccf_n_{klab[k]}"] = int(nk); pr[k] = pv
# Holm
order = sorted(pr, key=pr.get)
padj = {}; running = 0.0
for i, k in enumerate(order):
    running = max(running, min(1.0, (7 - i) * pr[k]))
    padj[k] = running
for k in range(-3, 4):
    R[f"q2_ccf_pholm_{klab[k]}"] = float(padj[k])

def direction(dep, lag1_own, lag1_other, idx):
    idx = np.asarray(idx)
    X = np.column_stack([month_dummies(idx), lag1_own[idx - 1], lag1_other[idx - 1]])
    y = dep[idx]
    beta, u, XtX_inv = ols(X, y)
    V = hac_cov(X, u, XtX_inv, 3)
    c = beta[-1]; se = math.sqrt(V[-1, -1])
    return float(c), float(se), float(2 * stats.norm.sf(abs(c / se)))
I59 = win("2020-02", "2024-12")
R["q2_fwd_coef"], R["q2_fwd_se"], R["q2_fwd_p"] = direction(ly, ly, lx, I59)
R["q2_rev_coef"], R["q2_rev_se"], R["q2_rev_p"] = direction(lx, lx, ly, I59)
c58, _, p58 = direction(lx, lx, ly, win("2020-03", "2024-12"))
R["q2_rev58_coef"], R["q2_rev58_p"] = c58, p58

# ---------------------------------------------------------------- NB2
idx = np.array(I59)
mm = CAL[idx].astype(float)
Xnb = np.column_stack([np.ones(len(idx)), np.sin(2 * np.pi * mm / 12), np.cos(2 * np.pi * mm / 12),
                       np.sin(4 * np.pi * mm / 12), np.cos(4 * np.pi * mm / 12), ly[idx - 1], lx[idx]])
ynb = cases[idx]

def nb_ll_obs(beta, alpha, X=Xnb, y=ynb):
    mu = np.exp(X @ beta); r = 1 / alpha
    return (special.gammaln(y + r) - special.gammaln(r) - special.gammaln(y + 1)
            + r * np.log(r / (r + mu)) + y * np.log(mu / (r + mu)))

def nb_scores(theta, X=Xnb, y=ynb):
    """per-observation scores wrt (beta, alpha)."""
    beta, alpha = theta[:-1], theta[-1]
    mu = np.exp(X @ beta); r = 1 / alpha
    sb = ((y - mu) / (1 + alpha * mu))[:, None] * X
    dr = special.digamma(y + r) - special.digamma(r) + np.log(r) + 1 - np.log(r + mu) - (r + y) / (r + mu)
    sa = dr * (-1 / alpha ** 2)
    return np.column_stack([sb, sa])

def negll(th):
    return -nb_ll_obs(th[:-1], math.exp(th[-1])).sum()
def negll_grad(th):
    s = nb_scores(np.r_[th[:-1], math.exp(th[-1])]).sum(0)
    s[-1] *= math.exp(th[-1])
    return -s
# start: OLS of log(y+1)
b0 = np.linalg.lstsq(Xnb, np.log(ynb + 1), rcond=None)[0]
best = None
for la0 in [-2.0, -1.0, 0.0, 1.0]:
    res = optimize.minimize(negll, np.r_[b0, la0], jac=negll_grad, method="BFGS", options={"gtol": 1e-10, "maxiter": 10000})
    if best is None or res.fun < best.fun:
        best = res
th = best.x
# Newton polishing in (beta, alpha)
theta = np.r_[th[:-1], math.exp(th[-1])]
def total_hess(theta, eps=1e-6):
    k = len(theta); H = np.zeros((k, k))
    for i in range(k):
        e = np.zeros(k); e[i] = eps * max(1.0, abs(theta[i]))
        H[:, i] = (nb_scores(theta + e).sum(0) - nb_scores(theta - e).sum(0)) / (2 * e[i])
    return (H + H.T) / 2
for _ in range(20):
    g = nb_scores(theta).sum(0); H = total_hess(theta)
    step = np.linalg.solve(H, g)
    theta = theta - step
    if np.max(np.abs(step)) < 1e-12:
        break
NB_GRAD_MAX = float(np.max(np.abs(nb_scores(theta).sum(0))))
beta_nb, alpha_nb = theta[:-1], theta[-1]
R["nb_coef_x"] = float(beta_nb[-1]); R["nb_rr"] = float(2 ** beta_nb[-1]); R["nb_alpha"] = float(alpha_nb)

# HAC sandwich, full parameter vector (beta, alpha): bread = inverse observed information
Hn = total_hess(theta)
Ainv = np.linalg.inv(-Hn)
sc = nb_scores(theta)
def hac_meat(g, L=3):
    S = g.T @ g
    for j in range(1, L + 1):
        w = 1 - j / (L + 1); G = g[j:].T @ g[:-j]; S += w * (G + G.T)
    return S
V_full = Ainv @ hac_meat(sc) @ Ainv
se_full = math.sqrt(V_full[6, 6])
R["nb_rr_lo"] = float(2 ** (beta_nb[-1] - 1.96 * se_full)); R["nb_rr_hi"] = float(2 ** (beta_nb[-1] + 1.96 * se_full))
# alternative (reported only): beta-only sandwich with alpha fixed, expected-information bread (glm.nb + sandwich style)
mu_nb = np.exp(Xnb @ beta_nb)
Wt = mu_nb / (1 + alpha_nb * mu_nb)
Binv = np.linalg.inv(Xnb.T @ (Xnb * Wt[:, None]))
V_b = Binv @ hac_meat(sc[:, :-1]) @ Binv
se_b = math.sqrt(V_b[6, 6])
NB_ALT = (2 ** (beta_nb[-1] - 1.96 * se_b), 2 ** (beta_nb[-1] + 1.96 * se_b), se_full, se_b)

# Ljung-Box on Pearson residuals
pres = (ynb - mu_nb) / np.sqrt(mu_nb + alpha_nb * mu_nb ** 2)
def ljung_box(e, h):
    nn = len(e); d = e - e.mean(); den = d @ d
    Q = nn * (nn + 2) * sum(((d[:nn - k] @ d[k:]) / den) ** 2 / (nn - k) for k in range(1, h + 1))
    return float(Q), float(stats.chi2.sf(Q, h))
R["nb_lb_q12"], R["nb_lb_p12"] = ljung_box(pres, 12)

# Cook's distance from OLS of ly_t on same covariates
yc = ly[idx]
bc, uc, XtXi = ols(Xnb, yc)
hlev = np.einsum("ij,jk,ik->i", Xnb, XtXi, Xnb)
pk = Xnb.shape[1]; s2 = (uc @ uc) / (len(yc) - pk)
cook = uc ** 2 / (pk * s2) * hlev / (1 - hlev) ** 2
R["nb_cook_months"] = [MONTHS[i] for i, c in zip(idx, cook) if c > 4 / 59]

# ---------------------------------------------------------------- Q3
def design(model, ids):
    ids = np.asarray(ids)
    cols = [month_dummies(ids), ly[ids - 1][:, None]]
    if model == "M2": cols.append(lx[ids][:, None])
    if model == "M3": cols.append(lx[ids - 1][:, None])
    return np.hstack(cols)

def predict(model, train, test):
    test = np.asarray(test)
    if model == "M0":
        return cases[test - 1]
    beta, _, _ = ols(design(model, train), ly[np.asarray(train)])
    yhat = design(model, test) @ beta
    return np.maximum(np.exp(yhat) - 1, 0)

START = TI["2019-07"]
TEST36 = win("2022-01", "2024-12")
MODELS = ["M0", "M1", "M2", "M3"]
preds = {s: {M: np.zeros(0) for M in MODELS} for s in ["ps", "mo", "loso"]}
tests = {"ps": TEST36, "mo": TEST36, "loso": win("2020-02", "2024-12")}   # loso per coordinator erratum (59)
for M in MODELS:
    pp = []
    for Y in (2022, 2023, 2024):
        tr = win("2019-07", f"{Y-1}-12"); te = win(f"{Y}-01", f"{Y}-12")
        pp.append(predict(M, tr, te))
    preds["ps"][M] = np.concatenate(pp)
    preds["mo"][M] = np.array([predict(M, list(range(START, t)), [t])[0] for t in TEST36])
    pl = []
    allL = win("2020-02", "2024-12")
    for Y in range(2020, 2025):
        te = [i for i in allL if YR[i] == Y]; tr = [i for i in allL if YR[i] != Y]
        pl.append(predict(M, tr, te))
    preds["loso"][M] = np.concatenate(pl)
check("loso has 59 test months", len(preds["loso"]["M1"]) == 59)

iid = np.loadtxt(VER + "boot_q3_iid_idx.csv", delimiter=",", dtype=int)
b3 = np.loadtxt(VER + "boot_q3_block3_starts.csv", delimiter=",", dtype=int)
check("iid idx in 0..35, 36 cols", iid.min() >= 0 and iid.max() <= 35 and iid.shape[1] == 36)
check("block3 starts in 0..33, 12 cols", b3.min() >= 0 and b3.max() <= 33 and b3.shape[1] == 12)
for s in ["ps", "mo", "loso"]:
    obs = cases[np.asarray(tests[s])]
    for M in MODELS:
        e = preds[s][M] - obs
        R[f"q3_{s}_mae_{M}"] = float(np.mean(np.abs(e)))
        R[f"q3_{s}_rmse_{M}"] = float(np.sqrt(np.mean(e ** 2)))
        R[f"q3_{s}_maelog_{M}"] = float(np.mean(np.abs(np.log1p(preds[s][M]) - np.log1p(obs))))
    if s == "loso":
        continue
    yrs = YR[np.asarray(tests[s])]
    better = 0
    for Y in (2022, 2023, 2024):
        mae = {}
        for M in ("M1", "M2"):
            mae[M] = float(np.mean(np.abs(preds[s][M] - obs)[yrs == Y]))
            R[f"q3_{s}_season_mae_{M}_{Y}"] = mae[M]
        better += int(mae["M2"] < mae["M1"])
    R[f"q3_{s}_seasons_m2_better"] = better
    d = np.abs(preds[s]["M2"] - obs) - np.abs(preds[s]["M1"] - obs)
    R[f"q3_{s}_diff_mean"] = float(d.mean())
    bm = np.array([d[r].mean() for r in iid])
    R[f"q3_{s}_diff_lo"] = float(np.percentile(bm, 2.5)); R[f"q3_{s}_diff_hi"] = float(np.percentile(bm, 97.5))
    bm3 = np.array([d[np.concatenate([np.arange(st, st + 3) for st in r])].mean() for r in b3])
    R[f"q3_{s}_diff_b3lo"] = float(np.percentile(bm3, 2.5)); R[f"q3_{s}_diff_b3hi"] = float(np.percentile(bm3, 97.5))
    for M in MODELS:
        for thr in (1000, 500):
            ep = obs > thr; fl = preds[s][M] > thr
            R[f"q3_{s}_sens_{M}_{thr}"] = float((ep & fl).sum() / ep.sum()) if ep.sum() else float("nan")
            R[f"q3_{s}_spec_{M}_{thr}"] = float((~ep & ~fl).sum() / (~ep).sum()) if (~ep).sum() else float("nan")
    R[f"q3_{s}_rule_met"] = int(better >= 2 and R[f"q3_{s}_diff_hi"] < 0)

# ---------------------------------------------------------------- Q4
for Y in range(2019, 2025):
    s = YR == Y
    R[f"q4_{Y}_peak_ratio"] = float(cases[s].max() / x[s].max())
    R[f"q4_{Y}_total_ratio"] = float(cases[s].sum() / x[s].sum())

# ---------------------------------------------------------------- Q5
W5 = win("2020-01", "2023-12")
for Rg in ["Central", "Eastern", "Western"]:
    srch = regsearch[Rg]; cs = regcases[Rg]
    cs72 = np.r_[cs, np.full(12, np.nan)]
    a = adj(lxf(srch), W5); b = adj(np.log(cs72 + 1), W5)
    r, ne, lo, hi = corr_block(a, b)
    R[f"q5_{Rg}_r"] = r; R[f"q5_{Rg}_neff"] = ne; R[f"q5_{Rg}_flo"] = lo; R[f"q5_{Rg}_fhi"] = hi
    R[f"q5_{Rg}_zeros2020_2023"] = int(np.sum(srch[W5] == 0))
    R[f"q5_{Rg}_peakratio_cases"] = float(cs72[YR == 2023].max() / cs72[YR == 2022].max())
    R[f"q5_{Rg}_peakratio_search"] = float(srch[YR == 2023].max() / srch[YR == 2022].max())

# ---------------------------------------------------------------- peaks
series_years = []
for Y in range(2020, 2025):
    series_years.append((f"cases_nat_{Y}", cases[YR == Y]))
    series_years.append((f"x_nat_{Y}", x[YR == Y]))
for Rg in ["Central", "Eastern", "Western"]:
    cs72 = np.r_[regcases[Rg], np.full(12, np.nan)]
    for Y in range(2020, 2024):
        series_years.append((f"cases_{Rg}_{Y}", cs72[YR == Y]))
        series_years.append((f"search_{Rg}_{Y}", regsearch[Rg][YR == Y]))
check("34 peak series-years", len(series_years) == 34)
changed = 0; ties = {0.05: 0, 0.10: 0, 0.20: 0}; PEAK_DETAIL = []; NEAR = []
for nm, v in series_years:
    o = np.argsort(-v, kind="stable")
    top, sec = v[o[0]], v[o[1]]
    ratio = math.inf if sec == 0 else top / sec
    labs = []
    for g in (0.05, 0.10, 0.20):
        lab = MNAME[o[0]] if ratio >= 1 + g else "tie"
        if lab == "tie": ties[g] += 1
        if abs(ratio - (1 + g)) < 1e-9: NEAR.append((nm, ratio, g))
        labs.append(lab)
    if len(set(labs)) > 1: changed += 1
    PEAK_DETAIL.append((nm, float(top), float(sec), ratio, labs))
R["peak_n_changed"] = changed
R["peak_ties_05"] = ties[0.05]; R["peak_ties_10"] = ties[0.10]; R["peak_ties_20"] = ties[0.20]

# ---------------------------------------------------------------- Q6
pub = {"2024-07-31": 2000, "2024-09-03": 6000, "2024-09-23": 12644, "2024-10-15": 19599, "2024-10-25": 27153}
c24 = cases[YR == 2024]
for d, v in pub.items():
    yy, mo_, dd = map(int, d.split("-"))
    final = c24[:mo_ - 1].sum() + c24[mo_ - 1] * dd / calendar.monthrange(yy, mo_)[1]
    key = d.replace("-", "")
    R[f"q6_final_{key}"] = float(final); R[f"q6_pct_{key}"] = float(100 * v / final)

# ---------------------------------------------------------------- output
def clean(v):
    if isinstance(v, (np.integer,)): return int(v)
    if isinstance(v, (np.floating,)): return float(v)
    return v
R = {k: clean(v) for k, v in R.items()}
with open(OUT, "w") as f:
    json.dump(R, f, indent=1, allow_nan=True)

if __name__ == "__main__":
    print("keys:", len(R))
    print("checks failed:", [c for c in CHECKS if not c[1]], "of", len(CHECKS))
    print("paper header row (0-based):", paper_hdr)
    print("AIC:", AIC_TABLE)
    print("NB grad max:", NB_GRAD_MAX, "converged BFGS:", best.success)
    print("NB alt (beta-only sandwich) RR lo/hi, se_full, se_beta:", NB_ALT)
    print("near-boundary peak ratios:", NEAR)
    for p in PEAK_DETAIL: print(p)
    print("nan keys:", [k for k, v in R.items() if isinstance(v, float) and math.isnan(v)])
