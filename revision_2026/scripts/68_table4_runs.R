# Table 4, second implementation (R): Questions 1 and 2 in the primary analysis and every sensitivity run.
# Extends verification/impl2_R.R (which covered the primary run only) to all nine runs. Conventions as in the
# pipeline (scripts/09): each run is a series of consecutive positions; run D joins December 2021 to January 2023,
# so its lags span that gap, as stated in Methods 2.5.
# Output: out/table4_runs_R.csv
suppressMessages(library(sandwich))
rd <- function(f, lt = "0.5") { L <- readLines(paste0("data_raw/", f), encoding = "UTF-8"); L <- L[grepl("^20[0-9]{2}-[0-9]{2},", L)]
  v <- trimws(sub("^[^,]*,", "", L)); v[v == "<1"] <- lt; as.numeric(v) }
f3 <- c("gt_dengue_term.csv", "gt_dengue_term_day2_2026-11-22.csv", "gt_dengue_term_day3_2026-11-23.csv")
x <- Reduce(`+`, lapply(f3, rd)) / 3; x0 <- Reduce(`+`, lapply(f3, rd, lt = "0")) / 3; d1 <- rd(f3[1])
pt_ <- read.csv("out/paper_72month_table.csv"); stopifnot(nrow(pt_) == 72)
cases <- pt_$cases; topic <- pt_$DF; fever <- pt_$Fever
YR <- rep(2019:2024, each = 12); MO <- rep(1:12, 6); YM <- YR * 100 + MO
W <- function(lo, hi) which(YM >= lo & YM <= hi)
P <- W(202001, 202412)
runs <- list(primary = list(x, P), A = list(x, W(201906, 202412)), B = list(x, W(201901, 202412)), C = list(x, W(202201, 202412)),
             D = list(x, P[YR[P] != 2022]), topic = list(topic, P), day1 = list(d1, P), lt0 = list(x0, P), fever = list(fever, P))
ly <- function(v) log(v + 1); lx <- function(v) log(v + 0.5); adj <- function(v, m) v - ave(v, m)
lagcor <- function(v, k) cor(v[1:(length(v) - k)], v[(1 + k):length(v)])
neff <- function(a, b) max(4, length(a) / (1 + 2 * sum(sapply(1:12, function(k) lagcor(a, k) * lagcor(b, k)))))
out <- NULL
for (k in names(runs)) {
  i <- runs[[k]][[2]]; s <- runs[[k]][[1]][i]; cc <- cases[i]; mo <- MO[i]
  a <- adj(lx(s), mo); b <- adj(ly(cc), mo); n <- length(a); r <- cor(a, b); ne <- neff(a, b); f <- tanh(atanh(r) + c(-1, 1) * 1.96 / sqrt(ne - 3))
  aic <- sapply(0:3, function(p) { y <- a[4:n]; X <- if (p == 0) NULL else sapply(1:p, function(j) a[(4 - j):(n - j)])
    e <- if (p == 0) resid(lm(y ~ 1)) else resid(lm(y ~ X)); (n - 3) * log(sum(e^2) / (n - 3)) + 2 * (p + 1) })
  p <- which.min(aic) - 1
  if (p > 0) { Xp <- sapply(1:p, function(j) a[(p + 1 - j):(n - j)]); cf <- coef(lm(a[(p + 1):n] ~ Xp))
    filt <- function(v) sapply((p + 1):n, function(t) v[t] - cf[[1]] - sum(cf[-1] * v[t - (1:p)])); ea <- filt(a); eb <- filt(b) } else { ea <- a; eb <- b }
  m <- length(ea); rr <- nn <- pv <- numeric(7)
  for (j in 1:7) { L <- j - 4; if (L >= 0) { u <- ea[1:(m - L)]; v <- eb[(1 + L):m] } else { u <- ea[(1 - L):m]; v <- eb[1:(m + L)] }
    rr[j] <- cor(u, v); nn[j] <- length(u); tt <- rr[j] * sqrt((nn[j] - 2) / (1 - rr[j]^2)); pv[j] <- 2 * pt(-abs(tt), nn[j] - 2) }
  ph <- p.adjust(pv, "holm")
  LY <- ly(cc); LX <- lx(s); t2 <- 2:n
  fw <- lm(LY[t2] ~ 0 + factor(mo[t2]) + LY[t2 - 1] + LX[t2 - 1]); V <- NeweyWest(fw, lag = 3, prewhite = FALSE, adjust = FALSE)
  q <- length(coef(fw)); bw <- coef(fw)[[q]]; sw <- sqrt(V[q, q])
  out <- rbind(out, data.frame(run = k, n = n, r = r, n_eff = ne, ci_lo = f[1], ci_hi = f[2], ar_order = p,
                               lag_m1 = rr[3], lag_0 = rr[4], lag_p1 = rr[5], holm_m1 = ph[3], holm_0 = ph[4], holm_p1 = ph[5], n_lag0 = nn[4],
                               fwd_coef = bw, fwd_p = 2 * pnorm(-abs(bw / sw)), fwd_n = length(t2)))
}
write.csv(out, "out/table4_runs_R.csv", row.names = FALSE)
print(format(out, digits = 4), row.names = FALSE)
