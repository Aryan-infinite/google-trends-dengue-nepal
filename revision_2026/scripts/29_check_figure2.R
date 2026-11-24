# Independent check of Figure 2 (R; no Python-derived input). Rebuilds both panels from the original files,
# recomputes r and Cook's distance with base R, and compares with the points read back from the figure
# (out/fig2_plotted_values.csv) and with the logged primary result (out/final/run_log.txt, if present).
# Run from the project folder: Rscript scripts/29_check_figure2.R
suppressMessages(library(readxl))
fails <- 0
check <- function(label, ok) {
  cat(sprintf("[%s] %s\n", if (isTRUE(ok)) "PASS" else "FAIL", label)); if (!isTRUE(ok)) fails <<- fails + 1
}
tol <- 1e-9
months72 <- format(seq(as.Date("2019-01-01"), by = "month", length.out = 72), "%Y-%m")

xlsx <- Sys.getenv("PREPRINT_XLSX", "data_raw/paper_finalized_data_RP.xlsx")
x <- as.data.frame(read_excel(xlsx, col_names = FALSE, .name_repair = "minimal"))
hr <- which(apply(x, 1, function(r) any(r %in% "Confirmed Cases"))); hdr <- as.character(unlist(x[hr, ]))
cases72 <- as.numeric(x[(hr + 1):nrow(x), match("Confirmed Cases", hdr)]); cases72 <- cases72[!is.na(cases72)]
check("workbook gives 72 case values", length(cases72) == 72)
rd <- function(f) { d <- read.csv(f, skip = 2, check.names = FALSE, stringsAsFactors = FALSE)
  v <- d[[2]]; v[v == "<1"] <- "0.5"; stopifnot(identical(d[[1]], months72)); as.numeric(v) }
search72 <- rowMeans(sapply(c("data_raw/gt_dengue_term.csv", "data_raw/gt_dengue_term_day2_2026-11-22.csv",
                              "data_raw/gt_dengue_term_day3_2026-11-23.csv"), rd))

w <- months72 >= "2020-01"                                   # primary window 2020-01..2024-12
m <- months72[w]; cs <- cases72[w]; ss <- search72[w]; mon <- as.integer(substr(m, 6, 7))
check("primary window has 60 months", length(m) == 60)
Y <- log(cs + 1); X <- log(ss + 0.5)
Ya <- Y - ave(Y, mon); Xa <- X - ave(X, mon)               # subtract each calendar month's 2020-2024 mean
r_raw <- cor(X, Y); r_adj <- cor(Xa, Ya)

p <- read.csv("out/fig2_plotted_values.csv", stringsAsFactors = FALSE)
pr <- p[p$panel == "raw", ]; pa <- p[p$panel == "adj", ]
pr <- pr[order(pr$month), ]; pa <- pa[order(pa$month), ]
check("figure draws 60 points in each panel", nrow(pr) == 60 && nrow(pa) == 60)
check("panel months are exactly 2020-01..2024-12", identical(pr$month, m) && identical(pa$month, m))
check("panel (a) points equal log(search+0.5), log(cases+1)", max(abs(pr$x - X)) < tol && max(abs(pr$y - Y)) < tol)
check("panel (b) points equal month-adjusted values", max(abs(pa$x - Xa)) < tol && max(abs(pa$y - Ya)) < tol)
grp <- ifelse(substr(m, 1, 4) >= "2022", "2022–2024", "2020–2021")
check("group assignment (2020-21 vs 2022-24) correct", identical(pr$group, grp) && identical(pa$group, grp))

st <- read.csv("out/fig2_stats.csv", stringsAsFactors = FALSE); sv <- setNames(st$value, st$quantity)
check("r printed on panel (a) = R value", abs(as.numeric(sv["r_raw"]) - r_raw) < tol)
check("r printed on panel (b) = R value", abs(as.numeric(sv["r_adj"]) - r_adj) < tol)
check("text drawn on panel (a) is exactly 'r = 0.87'", sv["text_a"] == sprintf("r = %.2f", r_raw) && sv["text_a"] == "r = 0.87")
check("text drawn on panel (b) is 'r = 0.81' plus the two flagged-month labels",
      identical(strsplit(sv[["text_b"]], "|", fixed = TRUE)[[1]], sort(c("Jun 2020", "Sep 2022", sprintf("r = %.2f", r_adj)))))
if (file.exists("out/final/run_log.txt")) {
  lg <- readLines("out/final/run_log.txt")
  grab <- function(k) as.numeric(sub('.*: ', '', sub(',$', '', grep(k, lg, fixed = TRUE, value = TRUE)[1])))
  check("panel (b) r equals logged primary Q1 r (0.8102...)", abs(grab('"Q1_primary.primary_adjusted_log.r"') - r_adj) < tol)
  check("panel (a) r equals logged unadjusted r (0.8727...)", abs(grab('"Q1_primary.unadjusted_log.r"') - r_raw) < tol)
} else cat("[SKIP] run_log.txt not present (repository copy); r checked against R recomputation only\n")

# Cook's distance from the same regression, fitted with base R lm()
d <- data.frame(y = log1p(cs), y1 = c(NA, head(log1p(cs), -1)), x = X,
                s1 = sin(2 * pi * mon / 12), c1 = cos(2 * pi * mon / 12), s2 = sin(4 * pi * mon / 12), c2 = cos(4 * pi * mon / 12))
keep <- !is.na(d$y1); fit <- lm(y ~ s1 + c1 + s2 + c2 + y1 + x, data = d[keep, ])
cd <- cooks.distance(fit); flag <- m[keep][cd > 4 / sum(keep)]
check("Cook's distance flags exactly 2020-06 and 2022-09 (n = 59)",
      identical(flag, c("2020-06", "2022-09")) && sum(keep) == 59 && sv["flagged"] == "2020-06;2022-09")

# Descriptive diagnostic for the caption/Discussion (post hoc; not a test): where does r_adj come from?
g <- substr(m, 1, 4) >= "2022"
cat(sprintf("DIAG within-group r (panel b): 2020-21 = %.4f (n=%d); 2022-24 = %.4f (n=%d)\n",
            cor(Xa[!g], Ya[!g]), sum(!g), cor(Xa[g], Ya[g]), sum(g)))
yr <- substr(m, 1, 4)
cat(sprintf("DIAG r after also removing each year's mean = %.4f\n", cor(Xa - ave(Xa, yr), Ya - ave(Ya, yr))))
cat(sprintf("DIAG group means (x, y): 2020-21 = (%.3f, %.3f); 2022-24 = (%.3f, %.3f)\n",
            mean(Xa[!g]), mean(Ya[!g]), mean(Xa[g]), mean(Ya[g])))
cat(sprintf("RESULT: %d failed check(s)\n", fails))
quit(status = if (fails == 0) 0 else 1)
