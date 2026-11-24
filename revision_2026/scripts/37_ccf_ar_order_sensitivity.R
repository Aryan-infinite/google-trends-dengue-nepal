# POST HOC (not in the analysis plan): sensitivity of the prewhitened cross-correlations to the AR filter order.
# AIC chose AR(2); AR(3) was 1.61 AIC points behind. Recomputes the CCF with AR(1), AR(2), AR(3). Python twin: scripts/37_ccf_ar_order_sensitivity.py
suppressMessages(library(readxl))
mk <- format(seq(as.Date("2019-01-01"), by = "month", length.out = 72), "%Y-%m")
x <- as.data.frame(read_excel(Sys.getenv("PREPRINT_XLSX", "data_raw/paper_finalized_data_RP.xlsx"), col_names = FALSE, .name_repair = "minimal"))
hr <- which(apply(x, 1, function(r) any(r %in% "Confirmed Cases"))); hd <- as.character(unlist(x[hr, ]))
cases <- as.numeric(x[(hr + 1):nrow(x), match("Confirmed Cases", hd)]); cases <- cases[!is.na(cases)]
gt <- function(f) { d <- read.csv(f, skip = 2); v <- d[[2]]; v[v == "<1"] <- "0.5"; as.numeric(v) }
s <- (gt("data_raw/gt_dengue_term.csv") + gt("data_raw/gt_dengue_term_day2_2026-11-22.csv") + gt("data_raw/gt_dengue_term_day3_2026-11-23.csv")) / 3
w <- 13:72; mon <- as.integer(substr(mk[w], 6, 7)); a <- log(s[w] + .5); a <- a - ave(a, mon); b <- log(cases[w] + 1); b <- b - ave(b, mon); n <- 60
for (p in 1:3) { Xp <- sapply(1:p, function(j) a[(p + 1 - j):(n - j)]); cf <- coef(lm(a[(p + 1):n] ~ Xp))
  f <- function(v) { e <- rep(NA, n); for (t in (p + 1):n) e[t] <- v[t] - cf[1] - sum(cf[-1] * v[t - (1:p)]); e }; ea <- f(a); eb <- f(b)
  r <- sapply(-3:3, function(k) { i <- (p + 1):n; i <- i[i + k >= p + 1 & i + k <= n]; ct <- cor.test(ea[i], eb[i + k]); c(ct$estimate, ct$p.value) })
  ph <- p.adjust(r[2, ], "holm"); cat(sprintf("R  AR(%d): r = %s | Holm p lag0 = %.4f, lag+1 = %.3f, min non-zero = %.3f at lag %d\n", p, paste(sprintf("%.3f", r[1, ]), collapse = " "), ph[4], ph[5], min(ph[-4]), (-3:3)[-4][which.min(ph[-4])])) }
