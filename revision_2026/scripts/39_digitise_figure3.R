# Image-level check of the DELIVERED Figure 3 TIFF, second implementation (R; ImageMagick decode; base R).
# Different estimators from scripts/38: bar value = median top-to-bottom extent of blue per column across the middle
# half of the bar (a pixel count would be short where the dashed threshold line is drawn across the lag 0 bar), signed by position relative to the zero line; dashed thresholds = rows
# where the INK2 row-profile peaks inside each bar's span. Expected values recomputed in R from the raw files.
suppressMessages(library(readxl))
fails <- 0
check <- function(l, ok, d = "") { cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", l, if (nzchar(d)) paste0("  (", d, ")") else "")); if (!isTRUE(ok)) fails <<- fails + 1 }
runs <- function(i) { if (!length(i)) return(matrix(numeric(0), ncol = 2)); b <- which(diff(i) > 1); cbind(c(i[1], i[b + 1]), c(i[b], i[length(i)])) }
ppm <- tempfile(fileext = ".ppm"); stopifnot(system2("convert", c("manuscript/figures/Fig3.tif", "-colorspace", "sRGB", "-depth", "8", ppm)) == 0)
con <- file(ppm, "rb"); hdr <- character(0)
while (length(hdr) < 4) { ch <- character(0); repeat { s <- rawToChar(readBin(con, "raw", 1)); if (s %in% c(" ", "\n")) break; ch <- c(ch, s) }; if (length(ch)) hdr <- c(hdr, paste(ch, collapse = "")) }
w <- as.integer(hdr[2]); h <- as.integer(hdr[3]); v <- readBin(con, "raw", w * h * 3); close(con); unlink(ppm)
ch <- function(k) matrix(as.integer(v[seq(k, length(v), 3)]), nrow = h, byrow = TRUE); R_ <- ch(1); G_ <- ch(2); B_ <- ch(3); rm(v)
near <- function(rgb, t) abs(R_ - rgb[1]) + abs(G_ - rgb[2]) + abs(B_ - rgb[3]) <= t
ink2 <- near(c(48, 62, 59), 30); grid <- near(c(208, 221, 211), 18); blue <- near(c(52, 125, 255), 60)

# expected values
mk <- format(seq(as.Date("2019-01-01"), by = "month", length.out = 72), "%Y-%m")
x <- as.data.frame(read_excel(Sys.getenv("PREPRINT_XLSX", "data_raw/paper_finalized_data_RP.xlsx"), col_names = FALSE, .name_repair = "minimal"))
hr <- which(apply(x, 1, function(r) any(r %in% "Confirmed Cases"))); hd <- as.character(unlist(x[hr, ]))
cases <- as.numeric(x[(hr + 1):nrow(x), match("Confirmed Cases", hd)]); cases <- cases[!is.na(cases)]
gt <- function(f) { d <- read.csv(f, skip = 2); vv <- d[[2]]; vv[vv == "<1"] <- "0.5"; as.numeric(vv) }
s <- (gt("data_raw/gt_dengue_term.csv") + gt("data_raw/gt_dengue_term_day2_2026-11-22.csv") + gt("data_raw/gt_dengue_term_day3_2026-11-23.csv")) / 3
ww <- 13:72; mon <- as.integer(substr(mk[ww], 6, 7)); a <- log(s[ww] + .5); a <- a - ave(a, mon); b <- log(cases[ww] + 1); b <- b - ave(b, mon); n <- 60; p <- 2
Xp <- sapply(1:p, function(j) a[(p + 1 - j):(n - j)]); cf <- coef(lm(a[(p + 1):n] ~ Xp))
flt <- function(z) { e <- rep(NA, n); for (t in (p + 1):n) e[t] <- z[t] - cf[1] - sum(cf[-1] * z[t - (1:p)]); e }; ea <- flt(a); eb <- flt(b)
E <- t(sapply(-3:3, function(k) { i <- (p + 1):n; i <- i[i + k >= p + 1 & i + k <= n]; m <- length(i); tq <- qt(.975, m - 2); c(cor(ea[i], eb[i + k]), tq / sqrt(m - 2 + tq^2)) }))

cs <- colSums(ink2); sp <- runs(which(cs > 0.5 * max(cs))); sp <- sp[sp[, 2] - sp[, 1] + 1 >= 9, , drop = FALSE]
check("one left spine", nrow(sp) == 1)
x0 <- sp[1, 2] + 1; yr <- runs(which(rowSums(ink2[, sp[1, 1]:sp[1, 2]]) == sp[1, 2] - sp[1, 1] + 1)); yr <- yr[yr[, 2] - yr[, 1] > 500, , drop = FALSE]
y0 <- yr[1, 1]; y1 <- yr[1, 2]; cc <- (x0 + 10):(w - 20)
hg <- (y0 - 1) + rowMeans(runs(which(rowSums(grid[y0:y1, cc]) > 0.5 * length(cc))))
zr <- runs(which(rowSums(ink2[y0:y1, cc]) > 0.8 * length(cc))); zr <- zr[(y0 - 1 + zr[, 2]) < y1 - 30, , drop = FALSE]
check("one zero line", nrow(zr) == 1); zy <- y0 - 1 + mean(zr[1, ])
fy <- lm(c(0.6, 0.4, 0.2, -0.2, 0) ~ c(sort(hg), zy)); sc <- coef(fy)[[2]]
check("y scale fits within 2 px", length(hg) == 4 && max(abs(resid(fy) / sc)) < 2, sprintf("%.2f px", max(abs(resid(fy) / sc))))
bm <- blue; bm[-(y0:y1), ] <- FALSE; bm[, 1:(x0 + 5)] <- FALSE
br <- runs(which(colSums(bm) > 0)); br <- br[br[, 2] - br[, 1] > 40, , drop = FALSE]
check("seven bars", nrow(br) == 7)
got <- sapply(1:nrow(br), function(k) { cols <- round(br[k, 1] + 0.25 * (br[k, 2] - br[k, 1])):round(br[k, 1] + 0.75 * (br[k, 2] - br[k, 1]))
  hgt <- median(sapply(cols, function(cc_) { r_ <- which(bm[, cc_]); max(r_) - min(r_) + 1 }))   # extent: unaffected by the dashed line drawn across a bar
  mid <- mean(which(bm[, round(mean(br[k, ]))])); sign <- if (mid < zy) 1 else -1
  sign * hgt * abs(sc) })
err <- abs(got - E[, 1]) / abs(sc)
check("bar heights (top-to-bottom extent) within 2 px of recomputed r", max(err) <= 2, sprintf("max %.2f px; read: %s", max(err), paste(sprintf("%.3f", got), collapse = ", ")))
unit <- mean(diff(rowMeans(br)))
derr <- c()
for (k in 1:7) { c0 <- mean(br[k, ]); cols <- round(c0 - 0.38 * unit):round(c0 + 0.38 * unit)
  pr <- rowSums(ink2[y0:y1, cols]) / length(cols); pr[(round(zy) - y0 + 1) + (-12:12)] <- 0; pr[(y1 - y0 + 1) - (0:25)] <- 0
  rr <- (y0 - 1) + rowMeans(runs(which(pr > 0.3)))
  exp_rows <- (c(E[k, 2], -E[k, 2]) - coef(fy)[[1]]) / sc
  derr <- c(derr, if (length(rr) == 2) abs(sort(rr) - sort(exp_rows)) else Inf) }
check("dashed thresholds within 2 px of +/- r_crit", max(derr) <= 2, sprintf("max %.2f px", max(derr)))
# no blue ink outside the seven expected bars (expected rectangles from the recomputed r, +3 px)
blue2 <- blue; blue2[-(y0:y1), ] <- FALSE; blue2[, 1:(x0 + 5)] <- FALSE
cxs <- rowMeans(br); ok <- matrix(FALSE, h, w)
for (k in 1:7) { yy <- sort(c((0 - coef(fy)[[1]]) / sc, (E[k, 1] - coef(fy)[[1]]) / sc))
  ok[max(1, floor(yy[1]) - 3):min(h, ceiling(yy[2]) + 3), max(1, floor(cxs[k] - 0.3 * unit) - 3):min(w, ceiling(cxs[k] + 0.3 * unit) + 3)] <- TRUE }
check("no blue ink outside the expected bars", !any(blue2 & !ok), sprintf("%d stray px", sum(blue2 & !ok)))
cat(sprintf("RESULT: %d failed check(s)\n", fails)); quit(status = if (fails == 0) 0 else 1)
