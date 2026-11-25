# Independent check of Supplementary Figure S1 (R) and three-program comparison.
# R recomputes ACF (stats::acf), PACF (lm on a constant and lags 1..k) and Ljung-Box (Box.test) from the raw files, then
# compares: Python (out/figS1_values.csv), JavaScript (out/figS1_values_js.csv), the drawn bars and bands
# (out/figS1_plotted_values.csv), panel titles and column headers, and the logged check-A runs (out/21_..., out/22_...).
suppressMessages(library(readxl))
fails <- 0
check <- function(l, ok, d = "") { cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", l, if (nzchar(d)) paste0("  (", d, ")") else "")); if (!isTRUE(ok)) fails <<- fails + 1 }
mk <- format(seq(as.Date("2019-01-01"), by = "month", length.out = 72), "%Y-%m")
x <- as.data.frame(read_excel(Sys.getenv("PREPRINT_XLSX", "data_raw/paper_finalized_data_RP.xlsx"), col_names = FALSE, .name_repair = "minimal"))
hr <- which(apply(x, 1, function(r) any(r %in% "Confirmed Cases"))); hd <- as.character(unlist(x[hr, ]))
cs <- as.numeric(x[(hr + 1):nrow(x), match("Confirmed Cases", hd)]); cs <- cs[!is.na(cs)][13:72]
gt <- function(f) { d <- read.csv(f, skip = 2, stringsAsFactors = FALSE); stopifnot(identical(d[[1]], mk)); v <- d[[2]]; v[v == "<1"] <- "0.5"; as.numeric(v) }
s <- ((gt("data_raw/gt_dengue_term.csv") + gt("data_raw/gt_dengue_term_day2_2026-11-22.csv") + gt("data_raw/gt_dengue_term_day3_2026-11-23.csv")) / 3)[13:72]
mo <- rep(1:12, 5); n <- 60; ly <- log(cs + 1); lx <- log(s + 0.5)
S <- list(ly_raw = ly, ly_adj = ly - ave(ly, mo), lx_raw = lx, lx_adj = lx - ave(lx, mo))
Rv <- do.call(rbind, lapply(names(S), function(k) { v <- S[[k]]
  a <- as.numeric(acf(v, lag.max = 20, plot = FALSE, demean = TRUE)$acf)[-1]
  p <- sapply(1:20, function(L) { X <- sapply(1:L, function(j) v[(L + 1 - j):(n - j)]); coef(lm(v[(L + 1):n] ~ X))[[L + 1]] })
  b6 <- Box.test(v, 6, "Ljung-Box"); b12 <- Box.test(v, 12, "Ljung-Box")
  data.frame(series = k, stat = c(paste0("acf", 1:20), paste0("pacf", 1:20), "lbq6", "lbq12", "lbp6", "lbp12"),
             value = c(a, p, b6$statistic, b12$statistic,                       # Box.test's p (1 - cdf) underflows to 0 here,
                       pchisq(b6$statistic, 6, lower.tail = FALSE), pchisq(b12$statistic, 12, lower.tail = FALSE))) }))  # so take the upper tail directly
key <- function(d) paste(d$series, d$stat)
py <- read.csv("out/figS1_values.csv"); js <- read.csv("out/figS1_values_js.csv")
q <- Rv$stat != "lbp6" & Rv$stat != "lbp12"                       # p-values are compared separately (relative scale)
dp <- max(abs(py$value[match(key(Rv), key(py))][q] - Rv$value[q])); dj <- max(abs(js$value[match(key(Rv[q, ]), key(js))] - Rv$value[q]))
check("160 ACF/PACF + 8 Ljung-Box Q: Python = R", dp < 1e-10, sprintf("max abs difference %.1e", dp))
check("160 ACF/PACF + 8 Ljung-Box Q: JavaScript = R", dj < 1e-10, sprintf("max abs difference %.1e", dj))
pp <- py$value[match(key(Rv), key(py))][!q]; rp <- Rv$value[!q]
check("Ljung-Box p-values: Python = R (relative)", max(abs(pp - rp) / rp) < 1e-6, sprintf("max relative difference %.1e; largest p %.1e", max(abs(pp - rp) / rp), max(rp)))
j21 <- readLines("out/21_reviewer_checks.json", warn = FALSE); j22 <- readLines("out/22_reviewer_checks.json")
grab <- function(js, k) { l <- grep(sprintf('"%s"', k), js, value = TRUE, fixed = FALSE)[1]; as.numeric(sub(",?\\s*$", "", sub('.*:\\s*', "", l))) }
lk <- paste0("A_", Rv$series[q], "_", Rv$stat[q])
d21 <- max(abs(sapply(lk, function(k) grab(j21, k)) - Rv$value[q])); d22 <- max(abs(sapply(lk, function(k) grab(j22, k)) - Rv$value[q]))
check("equal to the logged check A, Python run (script 21)", d21 < 1e-9, sprintf("max diff %.1e", d21))
check("equal to the logged check A, R run (script 22)", d22 < 1e-9, sprintf("max diff %.1e", d22))
d <- read.csv("out/figS1_plotted_values.csv", stringsAsFactors = FALSE)
b <- d[d$lag > 0, ]; ex <- Rv$value[match(paste(b$series, paste0(b$kind, b$lag)), key(Rv))]
check("160 drawn bars equal the R values", nrow(b) == 160 && max(abs(b$height - ex)) < 1e-10, sprintf("max diff %.1e", max(abs(b$height - ex))))
bd <- d[d$lag <= 0, ]
check("16 band lines at +/- 1.96/sqrt(60)", nrow(bd) == 16 && max(abs(abs(bd$height) - 1.96 / sqrt(60))) < 1e-12 && all(sign(bd$height) == ifelse(bd$lag == 0, 1, -1)))
want <- c(ly_raw = "Log cases", ly_adj = "Log cases, month-adjusted", lx_raw = "Log search interest", lx_adj = "Log search interest, month-adjusted")
tt <- unique(d[, c("series", "kind", "title")]); tt$want <- sprintf("(%s) %s", letters[match(paste(tt$series, tt$kind), paste(rep(names(want), each = 2), c("acf", "pacf")))], want[tt$series])
check("8 panel titles: letter order and series names", nrow(tt) == 8 && identical(tt$title, tt$want), paste(tt$title, collapse = "; "))
check("column headers", identical(read.csv("out/figS1_headers.csv")$header, c("Autocorrelation", "Partial autocorrelation")))
cat(sprintf("facts: lag-1 ACF %s; Ljung-Box Q12 %s (all p < %.0e)\n", paste(sprintf("%s %.3f", names(S), sapply(names(S), function(k) Rv$value[Rv$series == k & Rv$stat == "acf1"])), collapse = ", "),
    paste(sprintf("%.1f", Rv$value[Rv$stat == "lbq12"]), collapse = ", "), 10^ceiling(log10(max(Rv$value[Rv$stat %in% c("lbp6", "lbp12")])))))
cat(sprintf("RESULT: %d failed check(s)\n", fails)); quit(status = if (fails == 0) 0 else 1)
