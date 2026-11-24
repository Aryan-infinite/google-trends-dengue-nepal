# Independent check of Figure 3 (R; recomputes the prewhitened cross-correlations from the raw files with base R)
# and compares with the values read back from the drawing (out/fig3_plotted_values.csv) and, where present,
# the logged primary result (out/final/final_results.xlsx, sheet Q2_prewhitened_CCF).
suppressMessages(library(readxl))
fails <- 0
check <- function(label, ok, detail = "") {
  cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", label, if (nzchar(detail)) paste0("  (", detail, ")") else ""))
  if (!isTRUE(ok)) fails <<- fails + 1
}
tol <- 1e-9
mk <- format(seq(as.Date("2019-01-01"), by = "month", length.out = 72), "%Y-%m")
x <- as.data.frame(read_excel(Sys.getenv("PREPRINT_XLSX", "data_raw/paper_finalized_data_RP.xlsx"), col_names = FALSE, .name_repair = "minimal"))
hr <- which(apply(x, 1, function(r) any(r %in% "Confirmed Cases"))); hd <- as.character(unlist(x[hr, ]))
cases <- as.numeric(x[(hr + 1):nrow(x), match("Confirmed Cases", hd)]); cases <- cases[!is.na(cases)]
gt <- function(f) { d <- read.csv(f, skip = 2, stringsAsFactors = FALSE); stopifnot(identical(d[[1]], mk)); v <- d[[2]]; v[v == "<1"] <- "0.5"; as.numeric(v) }
search <- (gt("data_raw/gt_dengue_term.csv") + gt("data_raw/gt_dengue_term_day2_2026-11-22.csv") + gt("data_raw/gt_dengue_term_day3_2026-11-23.csv")) / 3
w <- 13:72; mon <- as.integer(substr(mk[w], 6, 7))
a <- log(search[w] + 0.5); a <- a - ave(a, mon)               # search, month-adjusted
b <- log(cases[w] + 1);   b <- b - ave(b, mon)                # cases, month-adjusted
n <- length(a)

# AR order by AIC on the common sample t = 4..60 (first three months held back), OLS with a constant
aic <- sapply(0:3, function(p) {
  yy <- a[4:n]; if (p == 0) f <- lm(yy ~ 1) else f <- lm(yy ~ sapply(1:p, function(j) a[(4 - j):(n - j)]))
  m <- length(yy); m * log(sum(resid(f)^2) / m) + 2 * (p + 1) })
p <- which.min(aic) - 1
cat(sprintf("AIC (p = 0..3): %s; next best minus best = %.2f\n", paste(sprintf("%.2f", aic), collapse = ", "), sort(aic)[2] - min(aic)))
check("AR order selected = 2", p == 2)
Xp <- sapply(1:p, function(j) a[(p + 1 - j):(n - j)]); cf <- coef(lm(a[(p + 1):n] ~ Xp))
flt <- function(v) { e <- rep(NA, n); for (t in (p + 1):n) e[t] <- v[t] - cf[1] - sum(cf[-1] * v[t - (1:p)]); e }
ea <- flt(a); eb <- flt(b)                                      # the SEARCH filter applied to both series
res <- t(sapply(-3:3, function(k) {
  i <- (p + 1):n; i <- i[i + k >= p + 1 & i + k <= n]               # pairs (search at t, cases at t + k)
  ct <- cor.test(ea[i], eb[i + k]); m <- length(i); tq <- qt(0.975, m - 2)
  c(lag = k, r = unname(ct$estimate), n = m, p = ct$p.value, r_crit = tq / sqrt(m - 2 + tq^2)) }))
res <- as.data.frame(res); res$p_holm <- p.adjust(res$p, "holm")
print(format(res, digits = 6), row.names = FALSE)

d <- read.csv("out/fig3_plotted_values.csv", stringsAsFactors = FALSE); d <- d[order(d$lag), ]
check("seven bars at lags -3..3", identical(as.integer(d$lag), -3:3))
check("bar heights equal recomputed r", max(abs(d$bar_height - res$r)) < tol, sprintf("max diff %.2e", max(abs(d$bar_height - res$r))))
check("dashed segments equal +/- P=0.05 threshold for each n", max(abs(d$crit_pos - res$r_crit)) < tol && max(abs(d$crit_neg + res$r_crit)) < tol)
lab <- ifelse(res$lag < 0, paste0("−", abs(res$lag)), ifelse(res$lag > 0, paste0("+", res$lag), "0"))
check("tick labels show lag and n", identical(d$tick_label, sprintf("%s | n = %d", lab, as.integer(res$n))), paste(d$tick_label, collapse = "; "))
want <- ifelse(res$p_holm < 0.001, "Holm-adjusted P < 0.001", ifelse(res$p_holm < 0.05, sprintf("Holm-adjusted P = %.3f", res$p_holm), ""))
d$annotation[is.na(d$annotation)] <- ""
check("Holm label on exactly the Holm-significant lags", identical(d$annotation, want), sprintf("labelled lags: %s", paste(d$lag[d$annotation != ""], collapse = ",")))
check("AR order recorded in figure output = 2", all(d$ar_order == 2))
if (file.exists("out/final/final_results.xlsx")) {
  L <- as.data.frame(read_excel("out/final/final_results.xlsx", sheet = "Q2_prewhitened_CCF")); L <- L[order(L$lag), ]
  check("r, n, Holm p equal the logged primary analysis",
        max(abs(L$r - res$r)) < tol && all(L$n == res$n) && max(abs(L$p_holm - res$p_holm)) < tol,
        sprintf("max r diff %.1e; max Holm p diff %.1e", max(abs(L$r - res$r)), max(abs(L$p_holm - res$p_holm))))
} else cat("[SKIP] logged results not present (repository copy)\n")
cat(sprintf("facts: lag 0 r = %.3f (n = %d, Holm p = %.4f); lag +1 r = %.3f (Holm p = %.3f); min Holm p at non-zero lags = %.3f\n",
    res$r[4], as.integer(res$n[4]), res$p_holm[4], res$r[5], res$p_holm[5], min(res$p_holm[-4])))
cat(sprintf("RESULT: %d failed check(s)\n", fails))
quit(status = if (fails == 0) 0 else 1)
