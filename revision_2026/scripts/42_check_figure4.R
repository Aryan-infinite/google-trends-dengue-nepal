# Independent check of Figure 4 (R) and three-program comparison of every nowcast.
# R fits each model by QR decomposition (qr.solve) from the raw files, then compares:
#   (1) R vs Python (out/fig4_predictions_python.csv, numpy least squares) vs JavaScript (out/fig4_predictions_js.csv,
#       normal equations + Gaussian elimination): all 216 predictions;
#   (2) the values drawn in the figure (out/fig4_plotted_values.csv), the legend and panel text (out/fig4_text.csv);
#   (3) the logged primary pipeline where present (statsmodels formulas: out/final/final_results.xlsx, out/15_nowcast_monthly_refit.xlsx).
suppressMessages(library(readxl))
fails <- 0
check <- function(l, ok, d = "") { cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", l, if (nzchar(d)) paste0("  (", d, ")") else "")); if (!isTRUE(ok)) fails <<- fails + 1 }
mk <- format(seq(as.Date("2019-01-01"), by = "month", length.out = 72), "%Y-%m")
x <- as.data.frame(read_excel(Sys.getenv("PREPRINT_XLSX", "data_raw/paper_finalized_data_RP.xlsx"), col_names = FALSE, .name_repair = "minimal"))
hr <- which(apply(x, 1, function(r) any(r %in% "Confirmed Cases"))); hd <- as.character(unlist(x[hr, ]))
cases <- as.numeric(x[(hr + 1):nrow(x), match("Confirmed Cases", hd)]); cases <- cases[!is.na(cases)]
gt <- function(f) { d <- read.csv(f, skip = 2, stringsAsFactors = FALSE); stopifnot(identical(d[[1]], mk)); v <- d[[2]]; v[v == "<1"] <- "0.5"; as.numeric(v) }
s <- (gt("data_raw/gt_dengue_term.csv") + gt("data_raw/gt_dengue_term_day2_2026-11-22.csv") + gt("data_raw/gt_dengue_term_day3_2026-11-23.csv")) / 3
LY <- log(cases + 1); LX <- log(s + 0.5); MO <- rep(1:12, 6); YR <- rep(2019:2024, each = 12)
des <- function(r, M) { X <- cbind(outer(MO[r], 1:12, "==") * 1, LY[r - 1]); if (M == "M2") X <- cbind(X, LX[r]); X }
prd <- function(tr, te, M) pmax(expm1(des(te, M) %*% qr.solve(des(tr, M), LY[tr])), 0)
start <- 7; test <- 37:72
R <- list()
for (sch in c("season", "monthly")) for (M in c("M0", "M1", "M2")) {
  v <- if (M == "M0") cases[test - 1] else if (sch == "season")
    unlist(lapply(2022:2024, function(y) prd(start:(which(mk == sprintf("%d-01", y)) - 1), test[YR[test] == y], M)))
  else sapply(test, function(t) prd(start:(t - 1), t, M))
  R[[paste(sch, M)]] <- as.numeric(v) }
key <- function(d) paste(d$schedule, d$model, d$month)
rv <- data.frame(schedule = rep(c("season", "monthly"), each = 108), model = rep(rep(c("M0", "M1", "M2"), each = 36), 2),
                 month = rep(mk[test], 6), pred = unlist(R[c(paste("season", c("M0", "M1", "M2")), paste("monthly", c("M0", "M1", "M2")))]))
py <- read.csv("out/fig4_predictions_python.csv", stringsAsFactors = FALSE); js <- read.csv("out/fig4_predictions_js.csv", stringsAsFactors = FALSE)
py <- py[match(key(rv), key(py)), ]; js <- js[match(key(rv), key(js)), ]
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1))
check("216 predictions present in all three programs", nrow(rv) == 216 && !anyNA(py$pred) && !anyNA(js$pred))
check("Python (numpy lstsq) = R (QR)", rel(py$pred, rv$pred) < 1e-8, sprintf("max relative difference %.1e", rel(py$pred, rv$pred)))
check("JavaScript (normal equations) = R (QR)", rel(js$pred, rv$pred) < 1e-8, sprintf("max relative difference %.1e", rel(js$pred, rv$pred)))
check("JavaScript = Python", rel(js$pred, py$pred) < 1e-8, sprintf("max relative difference %.1e", rel(js$pred, py$pred)))

obs <- cases[test]
mae <- sapply(split(rv, paste(rv$schedule, rv$model)), function(d) mean(abs(d$pred - obs)))
cat(sprintf("R MAE: %s\n", paste(sprintf("%s = %.6f", names(mae), mae), collapse = "; ")))
mj <- read.csv("out/fig4_mae_js.csv"); mj <- setNames(mj$mae, paste(mj$schedule, mj$model))
check("MAE identical across R and JavaScript", max(abs(mj[names(mae)] - mae)) < 1e-6, sprintf("max diff %.1e", max(abs(mj[names(mae)] - mae))))

# what the figure drew
d <- read.csv("out/fig4_plotted_values.csv", stringsAsFactors = FALSE)
for (sch in c("season", "monthly")) {
  o <- d[d$schedule == sch & d$series == "obs", ]
  check(sprintf("%s panel: 36 observed points equal the workbook (Jan 2022 - Dec 2024)", sch), identical(o$month, mk[test]) && max(abs(o$y - obs)) == 0)
  for (M in c("M0", "M1", "M2")) { g <- d[d$schedule == sch & d$series == M, ]
    check(sprintf("%s panel: drawn %s line equals R predictions", sch, M), identical(g$month, mk[test]) && rel(g$y, R[[paste(sch, M)]]) < 1e-8) }
}
check("no plotted value above the 40,000 axis limit (nothing clipped)", max(d$y) < 40000, sprintf("max %.0f", max(d$y)))
tx <- read.csv("out/fig4_text.csv", stringsAsFactors = FALSE); tv <- setNames(tx$text, tx$item)
fm <- function(v) formatC(round(v), format = "d", big.mark = ",")
for (p in list(c("panel_a", "season", "Refitted once per season (primary)"), c("panel_b", "monthly", "Refitted every month (post hoc)")))
  check(sprintf("%s text states the R-computed MAEs", p[1]), tv[[p[1]]] == sprintf("%s | Mean absolute error: M0 %s, M1 %s, M2 %s", p[3],
        fm(mae[paste(p[2], "M0")]), fm(mae[paste(p[2], "M1")]), fm(mae[paste(p[2], "M2")])), tv[[p[1]]])
check("key lists the four series", tv[["key"]] == "Observed | M0, last month carried forward | M1, seasonal terms + last month | M2, M1 + current-month searches")

# logged primary pipeline (statsmodels formulas)
if (file.exists("out/final/final_results.xlsx") && file.exists("out/15_nowcast_monthly_refit.xlsx")) {
  Ls <- as.data.frame(read_excel("out/final/final_results.xlsx", sheet = "Q3_rolling_predictions"))
  Lm <- as.data.frame(read_excel("out/15_nowcast_monthly_refit.xlsx", sheet = "predictions"))
  dl <- max(sapply(c("M0", "M1", "M2"), function(M) max(rel(Ls[[M]], R[[paste("season", M)]]), rel(Lm[[M]], R[[paste("monthly", M)]]))))
  check("all 216 equal the logged primary pipeline (statsmodels)", dl < 1e-8, sprintf("max relative difference %.1e", dl))
} else cat("[SKIP] logged pipeline outputs not present (repository copy)\n")

# facts for the caption
sea <- function(sch, M, y) mean(abs(R[[paste(sch, M)]][YR[test] == y] - obs[YR[test] == y]))
for (sch in c("season", "monthly")) cat(sprintf("seasons where M2 beats M1 (%s): %d of 3 [%s]\n", sch, sum(sapply(2022:2024, function(y) sea(sch, "M2", y) < sea(sch, "M1", y))),
    paste(sprintf("%d: M1 %.0f vs M2 %.0f", 2022:2024, sapply(2022:2024, function(y) sea(sch, "M1", y)), sapply(2022:2024, function(y) sea(sch, "M2", y))), collapse = "; ")))
i9 <- which(mk[test] == "2023-09"); cat(sprintf("Sep 2023: observed %d; M1 %.0f (season), %.0f (monthly); M2 %.0f (season), %.0f (monthly)\n",
    obs[i9], R[["season M1"]][i9], R[["monthly M1"]][i9], R[["season M2"]][i9], R[["monthly M2"]][i9]))
cat(sprintf("RESULT: %d failed check(s)\n", fails)); quit(status = if (fails == 0) 0 else 1)
