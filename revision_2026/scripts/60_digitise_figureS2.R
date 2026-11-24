# Image-level check of the DELIVERED Figure S2 TIFF, second implementation (R; ImageMagick decode; base R).
# Different estimators from scripts/59: in each region's row band, the CI ends are the outermost black columns and the dot
# is the middle of the columns where the black ink is tallest; markers in panel (b) by bounding-box centre of the ink in
# each row band. Expected values recomputed in R from the raw district workbook (own zone grouping as in scripts/58).
suppressMessages(library(readxl))
fails <- 0
check <- function(l, ok, d = "") { cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", l, if (nzchar(d)) paste0("  (", d, ")") else "")); if (!isTRUE(ok)) fails <<- fails + 1 }
runs <- function(i) { if (!length(i)) return(matrix(numeric(0), ncol = 2)); b <- which(diff(i) > 1); cbind(c(i[1], i[b + 1]), c(i[b], i[length(i)])) }
src <- readLines("scripts/58_check_figureS2.R"); eval(parse(text = src[grep("^ZONE <- list", src):grep("^print\\(format\\(Rv", src) - 1]))
ppm <- tempfile(fileext = ".ppm"); stopifnot(system2("convert", c("manuscript/figures/FigS2.tif", "-colorspace", "sRGB", "-depth", "8", ppm)) == 0)
con <- file(ppm, "rb"); hd <- character(0)
while (length(hd) < 4) { ch <- character(0); repeat { z <- rawToChar(readBin(con, "raw", 1)); if (z %in% c(" ", "\n")) break; ch <- c(ch, z) }; if (length(ch)) hd <- c(hd, paste(ch, collapse = "")) }
w <- as.integer(hd[2]); h <- as.integer(hd[3]); v <- readBin(con, "raw", w * h * 3); close(con); unlink(ppm)
k3 <- function(i) matrix(as.integer(v[seq(i, length(v), 3)]), nrow = h, byrow = TRUE); R_ <- k3(1); G_ <- k3(2); B_ <- k3(3); rm(v)
near <- function(c, t) abs(R_ - c[1]) + abs(G_ - c[2]) + abs(B_ - c[3]) <= t
ink2 <- near(c(48, 62, 59), 30); grid <- near(c(208, 221, 211), 18); blk <- R_ + G_ + B_ <= 60; blu <- near(c(52, 125, 255), 60)
cs <- colSums(ink2); bd <- runs(which(cs > 0.5 * max(cs))); bd <- bd[bd[, 2] - bd[, 1] + 1 >= 9, , drop = FALSE]
check("two spines and one reference line", nrow(bd) == 3)
calib <- function(s1, right, y0, yb, extra) { g <- (s1 + 4) + rowMeans(runs(which(colSums(grid[y0:(yb - 10), (s1 + 5):right]) > 0.5 * (yb - 10 - y0)))); sort(c(extra, g)) }
panel <- function(k) { s0 <- bd[k, 1]; s1 <- bd[k, 2]; yr <- runs(which(rowSums(ink2[, s0:s1]) == s1 - s0 + 1)); yr <- yr[yr[, 2] - yr[, 1] > 300, , drop = FALSE]
  y0 <- yr[1, 1]; bt <- runs(which(rowSums(ink2[y0:(yr[1, 2] + 30), (s1 + 20):(s1 + 1200)]) > 1000)); yb <- round(y0 - 1 + mean(bt[nrow(bt), ])); list(s0 = s0, s1 = s1, y0 = y0, yb = yb) }
# panel (a)
P <- panel(1); right <- bd[2, 1] - 800
xa <- calib(P$s1, right, P$y0, P$yb, mean(bd[1, ])); fa <- lm(c(0, .25, .5, .75, 1) ~ xa)
check("(a) x scale residual < 2 px", length(xa) == 5 && max(abs(resid(fa) / coef(fa)[2])) < 2)
m <- blk; m[-(P$y0:(P$yb - 8)), ] <- FALSE; m[, -((P$s1 + 8):right)] <- FALSE
rb <- runs(which(rowSums(m) > 0)); check("(a) three row bands of black ink", nrow(rb) == 3)
toX <- function(v) (v - coef(fa)[1]) / coef(fa)[2]; err <- c()
for (j in 1:nrow(rb)) { r <- Rv[Rv$region == c("Central", "Eastern", "Western")[j], ]; sub <- m[rb[j, 1]:rb[j, 2], , drop = FALSE]
  cols <- which(colSums(sub) > 0); hgt <- colSums(sub); dc <- which(hgt >= max(hgt) - 2)
  err <- c(err, abs(min(cols) - toX(r$ci_lo)), abs(max(cols) - toX(r$ci_hi)), abs((min(dc) + max(dc)) / 2 - toX(r$r))) }
check("(a) CI ends and dots within 2 px", max(err) <= 2, sprintf("max %.2f px", max(err)))
# panel (b)
P <- panel(2); xb <- calib(P$s1, w - 5, P$y0, P$yb, mean(bd[3, ])); fb <- lm(log10(c(.1, .3, 1, 3, 10)) ~ xb)
check("(b) log x scale residual < 2 px", length(xb) == 5 && max(abs(resid(fb) / coef(fb)[2])) < 2)
toXb <- function(v) (log10(v) - coef(fb)[1]) / coef(fb)[2]; err <- c()
for (cl in c("cases", "search")) { m <- if (cl == "cases") blk else blu; m[-(P$y0:(P$yb - 8)), ] <- FALSE; m[, -((P$s1 + 8):(w - 5))] <- FALSE
  rb <- runs(which(rowSums(m) > 0)); check(sprintf("(b) three %s markers (one per row band)", cl), nrow(rb) == 3)
  for (j in 1:nrow(rb)) { cols <- which(colSums(m[rb[j, 1]:rb[j, 2], , drop = FALSE]) > 0); r <- Rv[Rv$region == c("Central", "Eastern", "Western")[j], ]
    err <- c(err, abs((min(cols) + max(cols)) / 2 - toXb(if (cl == "cases") r$ratio_cases else r$ratio_search))) } }
check("(b) six markers within 2 px of the peak ratios (bounding-box centre)", max(err) <= 2, sprintf("max %.2f px", max(err)))
cat(sprintf("RESULT: %d failed check(s)\n", fails)); quit(status = if (fails == 0) 0 else 1)
