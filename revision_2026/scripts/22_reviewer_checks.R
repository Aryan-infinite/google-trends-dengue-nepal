# Checks A-D of ANALYSIS_PLAN.md section 9, reimplemented in R from the addendum text (not from scripts/21).
# Output: out/22_reviewer_checks.json with the same keys as scripts/21.
suppressMessages({library(readxl); library(MASS); library(sandwich)})
D <- "data_raw/"   # run from the project root (as run_all.sh does)
R <- list(); put <- function(k, v) R[[k]] <<- v
# Preprint table (PREPRINT_XLSX), columns found by header name. If the EDCD district workbook is present, its 2019-2023
# national totals must equal the preprint's (checked); the workbook is not public, so without it the preprint values are used.
xl <- Sys.getenv("PREPRINT_XLSX", paste0(D, "paper_finalized_data_RP.xlsx"))
pp <- as.data.frame(read_excel(xl, col_names = FALSE, .name_repair = "minimal"))
hr <- which(apply(pp, 1, function(r) any(r %in% "Confirmed Cases"))); hd <- as.character(unlist(pp[hr, ])); pp <- pp[(hr + 1):nrow(pp), ]
mo_ <- suppressWarnings(as.numeric(pp[[match("Months", hd)]])); cs_ <- suppressWarnings(as.numeric(pp[[match("Confirmed Cases", hd)]]))
cs_ <- cs_[!is.na(mo_)][order(mo_[!is.na(mo_)])]; stopifnot(length(cs_) == 72, !anyNA(cs_))
wbk <- paste0(D, "edcd_district_monthly_2019_2023.xlsx")
if (file.exists(wbk)) {
  ed <- unlist(lapply(2019:2023, function(y) { t <- suppressMessages(as.data.frame(read_excel(wbk, sheet = paste(y, "Monthly"), col_names = FALSE)))[-(1:2), ]
    m <- suppressWarnings(apply(as.matrix(t[, 4:15]), 2, as.numeric)); m[is.na(m)] <- 0; nm <- toupper(trimws(t[, 2])); colSums(m[!is.na(nm) & nm != "TOTAL", , drop = FALSE]) }))
  stopifnot(all(cs_[1:60] == ed))
}
cases_all <- cs_
gt <- function(f) { L <- readLines(paste0(D, f), encoding = "UTF-8"); L <- L[grepl("^20[0-9]{2}-[0-9]{2},", L)]; v <- sub("^[^,]*,", "", L); v[trimws(v) == "<1"] <- "0.5"; as.numeric(v) }
x_all <- (gt("gt_dengue_term.csv") + gt("gt_dengue_term_day2_2026-11-22.csv") + gt("gt_dengue_term_day3_2026-11-23.csv")) / 3
YR <- rep(2019:2024, each = 12); MO <- rep(1:12, 6); w <- YR >= 2020
x <- x_all[w]; cs <- cases_all[w]; yr <- YR[w]; mo <- MO[w]; n <- 60
ly <- function(v) log(v + 1); lx <- function(v) log(v + 0.5); adj <- function(v, m) v - ave(v, m)
lagcor <- function(v, k) cor(v[1:(length(v) - k)], v[(1 + k):length(v)])
neff <- function(a, b) max(4, length(a) / (1 + 2 * sum(sapply(1:12, function(k) lagcor(a, k) * lagcor(b, k)))))
fisher <- function(r, nn) tanh(atanh(r) + c(-1, 1) * 1.96 / sqrt(nn - 3))
pc <- function(v, p) as.numeric(quantile(v, p, type = 7))

# A
S <- list(ly_raw = ly(cs), lx_raw = lx(x), ly_adj = adj(ly(cs), mo), lx_adj = adj(lx(x), mo))
for (nm in names(S)) { v <- S[[nm]]; a <- as.numeric(acf(v, lag.max = 20, plot = FALSE, demean = TRUE)$acf)[-1]
  p <- sapply(1:20, function(k) { X <- sapply(1:k, function(j) v[(k + 1 - j):(n - j)]); coef(lm(v[(k + 1):n] ~ X))[[k + 1]] })
  for (k in 1:20) { put(sprintf("A_%s_acf%d", nm, k), a[k]); put(sprintf("A_%s_pacf%d", nm, k), p[k]) }
  for (L in c(6, 12)) { b <- Box.test(v, lag = L, type = "Ljung-Box"); put(sprintf("A_%s_lbq%d", nm, L), unname(b$statistic)); put(sprintf("A_%s_lbp%d", nm, L), b$p.value) } }
put("A_band", 1.96 / sqrt(60))

# B1, B2
ST <- as.matrix(read.csv("verification/boot_addendum_n48_block_starts.csv", header = FALSE))
pos <- c()
for (ex in 2020:2024) { k <- yr != ex; a <- adj(lx(x[k]), mo[k]); b <- adj(ly(cs[k]), mo[k]); r <- cor(a, b); ne <- neff(a, b); f <- fisher(r, ne)
  rs <- apply(ST, 1, function(s0) { i <- as.vector(sapply(s0, function(s) (s + 1):(s + 12)))[1:length(a)]; cor(a[i], b[i]) }); bl <- pc(rs, .025); bh <- pc(rs, .975)
  v <- c(r = r, neff = ne, flo = f[1], fhi = f[2], blo = bl, bhi = bh, lo = min(f[1], bl), hi = max(f[2], bh)); for (nm in names(v)) put(sprintf("B1_ex%d_%s", ex, nm), v[[nm]]); pos <- c(pos, r > 0) }
put("B1_all_positive", as.integer(all(pos)))
for (y in 2020:2024) { k <- yr == y; r <- cor(lx(x[k]), ly(cs[k])); f <- fisher(r, 12); put(sprintf("B2_%d_r", y), r); put(sprintf("B2_%d_lo", y), f[1]); put(sprintf("B2_%d_hi", y), f[2]) }
# B3
i <- which(YR >= 2020 & !(YR * 100 + MO) %in% c(202001, 202006, 202209))
dn <- data.frame(c = cases_all[i], s1 = sin(2 * pi * MO[i] / 12), c1 = cos(2 * pi * MO[i] / 12), s2 = sin(4 * pi * MO[i] / 12), c2 = cos(4 * pi * MO[i] / 12), y1 = ly(cases_all[i - 1]), x = lx(x_all[i]))
nb <- glm.nb(c ~ s1 + c1 + s2 + c2 + y1 + x, data = dn, control = glm.control(epsilon = 1e-12, maxit = 200)); bx <- coef(nb)[["x"]]
se <- sqrt(NeweyWest(nb, lag = 3, prewhite = FALSE, adjust = FALSE)["x", "x"])
put("B3_n", nrow(dn)); put("B3_coef_x", bx); put("B3_rr", 2^bx); put("B3_rr_lo", 2^(bx - 1.96 * se)); put("B3_rr_hi", 2^(bx + 1.96 * se)); put("B3_alpha", 1 / nb$theta)

# C
off <- mo %in% c(12, 1:5); xo <- x[off]; co <- cs[off]
put("C1_n_off", sum(off)); put("C1_le1", sum(xo <= 1)); put("C1_le2", sum(xo <= 2)); put("C1_share_le2", sum(xo <= 2) / sum(off)); put("C1_distinct", length(unique(round(xo, 6))))
put("C1_case_min", min(co)); put("C1_case_max", max(co)); put("C1_index_min", min(xo)); put("C1_index_max", max(xo))
a <- adj(lx(x), mo); b <- adj(ly(cs), mo)
for (lab in c("season", "off")) { k <- if (lab == "off") off else !off; r <- cor(a[k], b[k]); f <- fisher(r, sum(k))
  put(sprintf("C2_%s_n", lab), sum(k)); put(sprintf("C2_%s_r", lab), r); put(sprintf("C2_%s_lo", lab), f[1]); put(sprintf("C2_%s_hi", lab), f[2]) }
put("C_reading", if (sum(xo <= 2) / sum(off) >= 0.8) "lacks resolution" else "varies")

# D: OLS with HC1; Huber (k = 1.345, MAD scale about 0, IRLS) with the H1 covariance computed explicitly
o <- lm(b ~ a); X <- model.matrix(o); e <- resid(o); XtXi <- solve(crossprod(X))
V1 <- XtXi %*% crossprod(X * e) %*% XtXi * n / (n - 2); put("D_ols_slope", coef(o)[[2]]); put("D_ols_se", sqrt(V1[2, 2]))
put("D_ols_lo", coef(o)[[2]] - 1.96 * sqrt(V1[2, 2])); put("D_ols_hi", coef(o)[[2]] + 1.96 * sqrt(V1[2, 2]))
h <- rlm(b ~ a, psi = psi.huber, k = 1.345, scale.est = "MAD", maxit = 500, acc = 1e-14)
sr <- resid(h) / h$s; psi <- pmax(-1.345, pmin(1.345, sr)); dpsi <- as.numeric(abs(sr) <= 1.345); p <- 2
K <- 1 + p / n * (mean(dpsi^2) - mean(dpsi)^2) / mean(dpsi)^2
VH <- K^2 * (sum(psi^2) / (n - p) * h$s^2) / (mean(dpsi))^2 * XtXi; seh <- sqrt(VH[2, 2]); bh <- coef(h)[[2]]
put("D_huber_slope", bh); put("D_huber_se", seh); put("D_huber_lo", bh - 1.96 * seh); put("D_huber_hi", bh + 1.96 * seh)
put("D_ratio", bh / coef(o)[[2]]); put("D_robust", as.integer(sign(bh) == sign(coef(o)[[2]]) && (bh - 1.96 * seh > 0 || bh + 1.96 * seh < 0)))

js <- sapply(names(R), function(k) { v <- R[[k]]; sprintf('"%s": %s', k, if (is.character(v)) paste0('"', v, '"') else formatC(v, digits = 15, format = "g")) })
writeLines(paste0("{\n", paste(js, collapse = ",\n"), "\n}"), "out/22_reviewer_checks.json"); cat(length(R), "keys\n")
