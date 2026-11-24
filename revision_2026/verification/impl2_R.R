# Implementation 2: R (stats, MASS, sandwich), written from VERIFICATION_SPEC.md. Outputs results_impl2_R.json.
suppressMessages({library(readxl); library(MASS); library(sandwich)})
D <- "data_raw/"; V <- "verification/"   # run from the project root
R <- list(); put <- function(k, v) R[[k]] <<- v

# ---- data ----
ed <- list(); dist <- list()
for (y in 2019:2023) {
  t <- suppressMessages(as.data.frame(read_excel(paste0(D, "edcd_district_monthly_2019_2023.xlsx"), sheet = paste(y, "Monthly"), col_names = FALSE)))
  t <- t[-(1:2), ]; m <- suppressWarnings(apply(as.matrix(t[, 4:15]), 2, as.numeric)); m[is.na(m)] <- 0
  nm <- toupper(trimws(t[, 2])); keep <- !is.na(nm) & nm != "TOTAL"
  s <- colSums(m[keep, , drop = FALSE]); stopifnot(all(s == m[which(nm == "TOTAL"), ])); ed[[y - 2018]] <- s
  dd <- keep & nm != "NOT SPECIFIED"; dist[[as.character(y)]] <- list(name = nm[dd], m = m[dd, , drop = FALSE])
}
pp <- suppressMessages(as.data.frame(read_excel(paste0(D, "paper_finalized_data_RP.xlsx"), col_names = FALSE)))
pp <- pp[suppressWarnings(!is.na(as.numeric(pp[, 1]))), ]; pp[] <- lapply(pp, as.numeric); pp <- pp[order(pp[, 1]), ]
stopifnot(nrow(pp) == 72, all(pp[1:60, 5] == unlist(ed)))
cases <- c(unlist(ed), pp[61:72, 5]); topic <- pp[, 2]; fever <- pp[, 3]
gt <- function(f) { L <- readLines(paste0(D, f), encoding = "UTF-8"); L <- L[grepl("^20[0-9]{2}-[0-9]{2},", L)]; v <- sub("^[^,]*,", "", L); v[trimws(v) == "<1"] <- "0.5"; as.numeric(v) }
d1 <- gt("gt_dengue_term.csv"); x <- (d1 + gt("gt_dengue_term_day2_2026-11-22.csv") + gt("gt_dengue_term_day3_2026-11-23.csv")) / 3
YR <- rep(2019:2024, each = 12); MO <- rep(1:12, 6); YM <- YR * 100 + MO
ly <- function(v) log(v + 1); lx <- function(v) log(v + 0.5); adj <- function(v, m) v - ave(v, m)
lagcor <- function(v, k) cor(v[1:(length(v) - k)], v[(1 + k):length(v)])
neff <- function(a, b) max(4, length(a) / (1 + 2 * sum(sapply(1:12, function(k) lagcor(a, k) * lagcor(b, k)))))
fisher <- function(r, ne) tanh(atanh(r) + c(-1, 1) * 1.96 / sqrt(ne - 3))
pc <- function(v, p) as.numeric(quantile(v, p, type = 7))

# ---- Q1 ----
W <- function(lo, hi) YM >= lo & YM <= hi
runs <- list(primary = list(x, W(202001, 202412)), A = list(x, W(201906, 202412)), B = list(x, W(201901, 202412)), C = list(x, W(202201, 202412)),
  D = list(x, W(202001, 202412) & YR != 2022), topic = list(topic, W(202001, 202412)), day1 = list(d1, W(202001, 202412)), fever = list(fever, W(202001, 202412)))
AB <- list()
for (k in names(runs)) { s <- runs[[k]][[1]]; w <- runs[[k]][[2]]; a <- adj(lx(s[w]), MO[w]); b <- adj(ly(cases[w]), MO[w]); r <- cor(a, b); ne <- neff(a, b); f <- fisher(r, ne)
  put(paste0("q1_", k, "_r"), r); put(paste0("q1_", k, "_neff"), ne); put(paste0("q1_", k, "_flo"), f[1]); put(paste0("q1_", k, "_fhi"), f[2]); AB[[k]] <- list(a, b) }
st <- as.matrix(read.csv(paste0(V, "boot_q1_block_starts.csv"), header = FALSE))
for (k in c("primary", "fever")) { a <- AB[[k]][[1]]; b <- AB[[k]][[2]]
  rs <- apply(st, 1, function(s0) { i <- as.vector(sapply(s0, function(s) (s + 1):(s + 12)))[1:60]; cor(a[i], b[i]) })
  put(paste0("q1_", k, "_blo"), pc(rs, .025)); put(paste0("q1_", k, "_bhi"), pc(rs, .975)) }

# ---- Q2 ----
a <- AB$primary[[1]]; b <- AB$primary[[2]]; n <- length(a)
aic <- sapply(0:3, function(p) { y <- a[4:n]; X <- if (p == 0) NULL else sapply(1:p, function(j) a[(4 - j):(n - j)])
  e <- if (p == 0) resid(lm(y ~ 1)) else resid(lm(y ~ X)); (n - 3) * log(sum(e^2) / (n - 3)) + 2 * (p + 1) })
p <- which.min(aic) - 1; put("q2_ar_order", p)
Xp <- sapply(1:p, function(j) a[(p + 1 - j):(n - j)]); cf <- coef(lm(a[(p + 1):n] ~ Xp)); put("q2_ar_const", cf[[1]])
for (j in 1:p) put(paste0("q2_ar_phi", j), cf[[j + 1]])
filt <- function(v) sapply((p + 1):n, function(t) v[t] - cf[[1]] - sum(cf[-1] * v[t - (1:p)]))
ea <- filt(a); eb <- filt(b); m <- length(ea); lab <- c("m3", "m2", "m1", "0", "p1", "p2", "p3")
rr <- nn <- pv <- numeric(7)
for (i in 1:7) { k <- i - 4; if (k >= 0) { u <- ea[1:(m - k)]; v <- eb[(1 + k):m] } else { u <- ea[(1 - k):m]; v <- eb[1:(m + k)] }
  rr[i] <- cor(u, v); nn[i] <- length(u); tt <- rr[i] * sqrt((nn[i] - 2) / (1 - rr[i]^2)); pv[i] <- 2 * pt(-abs(tt), nn[i] - 2) }
ph <- p.adjust(pv, "holm")
for (i in 1:7) { put(paste0("q2_ccf_r_", lab[i]), rr[i]); put(paste0("q2_ccf_n_", lab[i]), nn[i]); put(paste0("q2_ccf_pholm_", lab[i]), ph[i]) }
LY <- ly(cases); LX <- lx(x)
direction <- function(lo, dep, own, other) { i <- which(W(lo, 202412)); f <- lm(dep[i] ~ 0 + factor(MO[i]) + own[i - 1] + other[i - 1])
  V <- NeweyWest(f, lag = 3, prewhite = FALSE, adjust = FALSE); k <- length(coef(f)); b <- coef(f)[k]; se <- sqrt(V[k, k]); c(b, se, 2 * pnorm(-abs(b / se))) }
o <- direction(202002, LY, LY, LX); put("q2_fwd_coef", o[1]); put("q2_fwd_se", o[2]); put("q2_fwd_p", o[3])
o <- direction(202002, LX, LX, LY); put("q2_rev_coef", o[1]); put("q2_rev_se", o[2]); put("q2_rev_p", o[3])
o <- direction(202003, LX, LX, LY); put("q2_rev58_coef", o[1]); put("q2_rev58_p", o[3])

# ---- NB2 ----
i <- which(W(202002, 202412)); dn <- data.frame(c = cases[i], s1 = sin(2 * pi * MO[i] / 12), c1 = cos(2 * pi * MO[i] / 12), s2 = sin(4 * pi * MO[i] / 12), c2 = cos(4 * pi * MO[i] / 12), y1 = LY[i - 1], x = LX[i])
nb <- glm.nb(c ~ s1 + c1 + s2 + c2 + y1 + x, data = dn, control = glm.control(epsilon = 1e-12, maxit = 200))
bx <- coef(nb)[["x"]]; al <- 1 / nb$theta; se <- sqrt(NeweyWest(nb, lag = 3, prewhite = FALSE, adjust = FALSE)["x", "x"])
put("nb_coef_x", bx); put("nb_rr", 2^bx); put("nb_alpha", al); put("nb_rr_lo", 2^(bx - 1.96 * se)); put("nb_rr_hi", 2^(bx + 1.96 * se))
mu <- fitted(nb); pr <- (dn$c - mu) / sqrt(mu + al * mu^2); lb <- Box.test(pr, lag = 12, type = "Ljung-Box"); put("nb_lb_q12", unname(lb$statistic)); put("nb_lb_p12", lb$p.value)
ol <- lm(log1p(c) ~ s1 + c1 + s2 + c2 + y1 + x, dn); ck <- cooks.distance(ol); put("nb_cook_months", sprintf("%d-%02d", YR[i], MO[i])[ck > 4 / nrow(dn)])

# ---- Q3 ----
iid <- as.matrix(read.csv(paste0(V, "boot_q3_iid_idx.csv"), header = FALSE)) + 1; b3 <- as.matrix(read.csv(paste0(V, "boot_q3_block3_starts.csv"), header = FALSE))
des <- function(r, M) { X <- cbind(matrix(sapply(1:12, function(k) as.numeric(MO[r] == k)), nrow = length(r)), LY[r - 1]); if (M == "M2") X <- cbind(X, LX[r]); if (M == "M3") X <- cbind(X, LX[r - 1]); X }
pred <- function(tr, te, M) { b <- qr.solve(des(tr, M), LY[tr]); pmax(expm1(des(te, M) %*% b), 0) }
ii <- which(W(201907, 202412))
evaluate <- function(s, groups, trainer) {
  P <- list(M0 = c(), M1 = c(), M2 = c(), M3 = c()); obs <- c(); yrs <- c()
  for (g in groups) { tr <- trainer(g); P$M0 <- c(P$M0, cases[g - 1]); obs <- c(obs, cases[g]); yrs <- c(yrs, YR[g]); for (M in c("M1", "M2", "M3")) P[[M]] <- c(P[[M]], pred(tr, g, M)) }
  for (M in names(P)) { e <- P[[M]] - obs; put(sprintf("q3_%s_mae_%s", s, M), mean(abs(e))); put(sprintf("q3_%s_rmse_%s", s, M), sqrt(mean(e^2)))
    put(sprintf("q3_%s_maelog_%s", s, M), mean(abs(log1p(P[[M]]) - log1p(obs))))
    if (s != "loso") for (th in c(1000, 500)) { o <- obs > th; f <- P[[M]] > th; put(sprintf("q3_%s_sens_%s_%d", s, M, th), sum(o & f) / sum(o)); put(sprintf("q3_%s_spec_%s_%d", s, M, th), sum(!o & !f) / sum(!o)) } }
  if (s == "loso") return(invisible())
  better <- 0
  for (yv in 2022:2024) { m1 <- mean(abs(P$M1 - obs)[yrs == yv]); m2 <- mean(abs(P$M2 - obs)[yrs == yv]); put(sprintf("q3_%s_season_mae_M1_%d", s, yv), m1); put(sprintf("q3_%s_season_mae_M2_%d", s, yv), m2); better <- better + (m2 < m1) }
  d <- abs(P$M2 - obs) - abs(P$M1 - obs); put(sprintf("q3_%s_seasons_m2_better", s), better); put(sprintf("q3_%s_diff_mean", s), mean(d))
  bm <- apply(iid, 1, function(r) mean(d[r])); put(sprintf("q3_%s_diff_lo", s), pc(bm, .025)); put(sprintf("q3_%s_diff_hi", s), pc(bm, .975))
  bm3 <- apply(b3, 1, function(r) mean(d[as.vector(sapply(r, function(s0) (s0 + 1):(s0 + 3)))])); put(sprintf("q3_%s_diff_b3lo", s), pc(bm3, .025)); put(sprintf("q3_%s_diff_b3hi", s), pc(bm3, .975))
  put(sprintf("q3_%s_rule_met", s), as.integer(better >= 2 && pc(bm, .975) < 0)) }
tm <- which(W(202201, 202412))
evaluate("ps", lapply(2022:2024, function(yv) tm[YR[tm] == yv]), function(g) ii[ii < min(g)])
evaluate("mo", as.list(tm), function(g) ii[ii < g])
lo <- which(W(202002, 202412)); evaluate("loso", lapply(2020:2024, function(yv) lo[YR[lo] == yv]), function(g) lo[YR[lo] != YR[g[1]]])

# ---- Q4 ----
for (yv in 2019:2024) { k <- YR == yv; put(sprintf("q4_%d_peak_ratio", yv), max(cases[k]) / max(x[k])); put(sprintf("q4_%d_total_ratio", yv), sum(cases[k]) / sum(x[k])) }

# ---- Q5 ----
REG <- list(Eastern = "TAPLEJUNG|PANCHTHAR|ILAM|JHAPA|SANKHUWASABHA|TERHATHUM|DHANKUTA|BHOJPUR|MORANG|SUNSARI|SOLUKHUMBU|OKHALDHUNGA|KHOTANG|UDAYAPUR|SAPTARI|SIRAHA",
  Central = "DOLAKHA|RAMECHHAP|SINDHULI|DHANUSA|DHANUSHA|MAHOTTARI|SARLAHI|SINDHUPALCHOK|KAVRE|LALITPUR|BHAKTAPUR|KATHMANDU|NUWAKOT|RASUWA|DHADING|MAKWANPUR|MAKAWANPUR|RAUTAHAT|\\bBARA\\b|PARSA|CHITWAN|CHITAWAN",
  Western = "GORKHA|LAMJUNG|TANAHU|SYANGJA|KASKI|MANANG|MUSTANG|MYAGDI|PARBAT|BAGLUNG|GULMI|PALPA|NAWALPARASI|RUPANDEHI|KAPILBASTU|KAPILVASTU|ARGHAKHANCHI",
  MidWestern = "PYUTHAN|ROLPA|RUKUM|SALYAN|\\bDANG\\b|DAILEKH|JAJARKOT|SURKHET|BANKE|BARDIYA|DOLPA|JUMLA|KALIKOT|MUGU|HUMLA",
  FarWestern = "BAJURA|BAJHANG|ACHHAM|DOTI|KAILALI|DARCHULA|BAITADI|DADELDHURA|KANCHANPUR")
nm23 <- dist[["2023"]]$name; hits <- sapply(nm23, function(z) sum(sapply(REG, function(p) grepl(p, z)))); stopifnot(length(nm23) == 77, all(hits == 1))
SR <- list(); RC <- list(); m48 <- rep(1:12, 4)
for (rg in c("Central", "Eastern", "Western")) {
  rc <- unlist(lapply(2020:2023, function(y) { dd <- dist[[as.character(y)]]; colSums(dd$m[grepl(REG[[rg]], dd$name), , drop = FALSE]) }))
  L <- readLines(sprintf("%sgt_dengue_region_%s.csv", D, tolower(rg)), encoding = "UTF-8"); L <- L[grepl("^202[0-3]-[0-9]{2},", L)]; v <- sub("^[^,]*,", "", L); v[v == "<1"] <- "0.5"; sr <- as.numeric(v)
  a <- adj(lx(sr), m48); b <- adj(ly(rc), m48); r <- cor(a, b); ne <- neff(a, b); f <- fisher(r, ne)
  put(sprintf("q5_%s_r", rg), r); put(sprintf("q5_%s_neff", rg), ne); put(sprintf("q5_%s_flo", rg), f[1]); put(sprintf("q5_%s_fhi", rg), f[2])
  put(sprintf("q5_%s_zeros2020_2023", rg), sum(sr == 0)); put(sprintf("q5_%s_peakratio_cases", rg), max(rc[37:48]) / max(rc[25:36])); put(sprintf("q5_%s_peakratio_search", rg), max(sr[37:48]) / max(sr[25:36]))
  SR[[rg]] <- sr; RC[[rg]] <- rc }
labs <- list(); cls <- function(v, yrs, y0) for (yv in yrs) { s <- v[((yv - y0) * 12 + 1):((yv - y0) * 12 + 12)]; o <- order(-s); ra <- if (s[o[2]] == 0) Inf else s[o[1]] / s[o[2]]
  labs[[length(labs) + 1]] <<- sapply(c(.05, .10, .20), function(g) if (ra >= 1 + g) month.abb[o[1]] else "tie") }
cls(cases, 2020:2024, 2019); cls(x, 2020:2024, 2019); for (rg in c("Central", "Eastern", "Western")) cls(RC[[rg]], 2020:2023, 2020); for (rg in c("Central", "Eastern", "Western")) cls(SR[[rg]], 2020:2023, 2020)
stopifnot(length(labs) == 34); put("peak_n_changed", sum(sapply(labs, function(l) length(unique(l)) > 1)))
for (j in 1:3) put(c("peak_ties_05", "peak_ties_10", "peak_ties_20")[j], sum(sapply(labs, function(l) l[j] == "tie")))

# ---- Q6 ----
c24 <- cases[61:72]
for (q in list(c("20240731", 2000), c("20240903", 6000), c("20240923", 12644), c("20241015", 19599), c("20241025", 27153))) {
  mth <- as.integer(substr(q[1], 5, 6)); day <- as.integer(substr(q[1], 7, 8)); nd <- as.integer(format(as.Date(sprintf("2024-%02d-01", mth %% 12 + 1)) - 1, "%d"))
  fin <- sum(c24[seq_len(mth - 1)]) + c24[mth] * day / nd; put(paste0("q6_final_", q[1]), fin); put(paste0("q6_pct_", q[1]), 100 * as.numeric(q[2]) / fin) }

# ---- write JSON (no jsonlite) ----
js <- sapply(names(R), function(k) { v <- R[[k]]; val <- if (is.character(v)) paste0("[", paste0('"', v, '"', collapse = ","), "]") else formatC(v, digits = 15, format = "g"); sprintf('"%s": %s', k, val) })
writeLines(paste0("{\n", paste(js, collapse = ",\n"), "\n}"), paste0(V, "results_impl2_R.json")); cat(length(R), "keys\n")
