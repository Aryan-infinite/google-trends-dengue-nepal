# Image-level check of the DELIVERED Figure S1 TIFF, second implementation (R; ImageMagick decode; base R).
# Different estimators from scripts/53: bar value = median top-to-bottom extent of blue over the middle half of the bar,
# signed by which side of the zero line the bar lies; y scale from a least-squares fit of the gridline rows; bar centres
# matched to lags by their order. Expected values recomputed in R (stats::acf and lm).
suppressMessages(library(readxl))
fails <- 0
check <- function(l, ok, d = "") { cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", l, if (nzchar(d)) paste0("  (", d, ")") else "")); if (!isTRUE(ok)) fails <<- fails + 1 }
runs <- function(i) { if (!length(i)) return(matrix(numeric(0), ncol = 2)); b <- which(diff(i) > 1); cbind(c(i[1], i[b + 1]), c(i[b], i[length(i)])) }
mk <- format(seq(as.Date("2019-01-01"), by = "month", length.out = 72), "%Y-%m")
x <- as.data.frame(read_excel(Sys.getenv("PREPRINT_XLSX", "data_raw/paper_finalized_data_RP.xlsx"), col_names = FALSE, .name_repair = "minimal"))
hr <- which(apply(x, 1, function(r) any(r %in% "Confirmed Cases"))); hd <- as.character(unlist(x[hr, ]))
cs <- as.numeric(x[(hr + 1):nrow(x), match("Confirmed Cases", hd)]); cs <- cs[!is.na(cs)][13:72]
gt <- function(f) { d <- read.csv(f, skip = 2); v <- d[[2]]; v[v == "<1"] <- "0.5"; as.numeric(v) }
s <- ((gt("data_raw/gt_dengue_term.csv") + gt("data_raw/gt_dengue_term_day2_2026-11-22.csv") + gt("data_raw/gt_dengue_term_day3_2026-11-23.csv")) / 3)[13:72]
mo <- rep(1:12, 5); n <- 60; ly <- log(cs + 1); lx <- log(s + .5)
S <- list(ly, ly - ave(ly, mo), lx, lx - ave(lx, mo))
EXP <- lapply(S, function(v) list(as.numeric(acf(v, lag.max = 20, plot = FALSE)$acf)[-1],
  sapply(1:20, function(L) coef(lm(v[(L + 1):n] ~ sapply(1:L, function(j) v[(L + 1 - j):(n - j)])))[[L + 1]])))
ppm <- tempfile(fileext = ".ppm"); stopifnot(system2("convert", c("manuscript/figures/FigS1.tif", "-colorspace", "sRGB", "-depth", "8", ppm)) == 0)
con <- file(ppm, "rb"); hd <- character(0)
while (length(hd) < 4) { ch <- character(0); repeat { z <- rawToChar(readBin(con, "raw", 1)); if (z %in% c(" ", "\n")) break; ch <- c(ch, z) }; if (length(ch)) hd <- c(hd, paste(ch, collapse = "")) }
w <- as.integer(hd[2]); h <- as.integer(hd[3]); v <- readBin(con, "raw", w * h * 3); close(con); unlink(ppm)
k3 <- function(i) matrix(as.integer(v[seq(i, length(v), 3)]), nrow = h, byrow = TRUE); R_ <- k3(1); G_ <- k3(2); B_ <- k3(3); rm(v); invisible(gc())
near <- function(c, t) abs(R_ - c[1]) + abs(G_ - c[2]) + abs(B_ - c[3]) <= t
ink2 <- near(c(48, 62, 59), 30); grid <- near(c(208, 221, 211), 18); blue <- near(c(52, 125, 255), 60)
cs_ <- colSums(ink2); sp <- runs(which(cs_ > 0.4 * max(cs_))); sp <- sp[sp[, 2] - sp[, 1] + 1 >= 9, , drop = FALSE]
check("two columns of panels", nrow(sp) == 2)
worst <- 0; nread <- 0; bandw <- 0
for (cc in 1:2) { s1 <- sp[cc, 2]; pr <- runs(which(rowSums(ink2[, sp[cc, 1]:s1]) == s1 - sp[cc, 1] + 1)); pr <- pr[pr[, 2] - pr[, 1] > 500, , drop = FALSE]
  check(sprintf("column %d: four panels", cc), nrow(pr) == 4)
  for (ri in 1:4) { lab <- letters[2 * (ri - 1) + cc]; y0 <- pr[ri, 1]; y1 <- pr[ri, 2]
    bt <- runs(which(rowSums(ink2[y0:(y1 + 40), (s1 + 20):(s1 + 1500)]) > 1200)); yb <- y0 - 1 + mean(bt[nrow(bt), ])
    cl <- which(colSums(ink2[round(yb) + (-3:3), (s1 + 1):w]) > 0); br <- which(diff(cl) > 60); xr <- s1 + (if (length(br)) cl[br[1]] else cl[length(cl)])
    cols <- (s1 + 20):(xr - 10)
    hg <- (y0 - 13) + rowMeans(runs(which(rowSums(grid[(y0 - 12):round(yb), cols]) > 0.5 * length(cols))))
    tri <- (y0 - 1) + rowMeans(runs(which(rowSums(ink2[y0:(round(yb) - 30), cols]) > 0.25 * length(cols))))
    check(sprintf("(%s) gridlines and three dark-grey bands found", lab), length(tri) == 3 && length(hg) %in% 2:3)
    yv <- c(if (length(hg) == 3) c(1, 0.5, -0.5) else c(0.5, -0.5), 0, -1); fy <- lm(yv ~ c(sort(hg), tri[2], yb)); sc <- coef(fy)[[2]]
    check(sprintf("(%s) y scale residual < 2 px", lab), max(abs(resid(fy) / sc)) < 2, sprintf("%.2f px", max(abs(resid(fy) / sc))))
    zy <- tri[2]; ex <- EXP[[ri]][[cc]]
    bm <- blue; bm[-(y0:round(yb)), ] <- FALSE; bm[, -((s1 + 5):(xr - 1))] <- FALSE
    bars <- runs(which(colSums(bm) > 0)); bars <- bars[bars[, 2] - bars[, 1] > 20, , drop = FALSE]
    vis <- which(abs(ex) / abs(sc) >= 1)
    check(sprintf("(%s) one visible bar per lag at least 1 px tall", lab), nrow(bars) == length(vis), sprintf("%d found, %d expected", nrow(bars), length(vis)))
    if (nrow(bars) == length(vis)) for (j in seq_along(vis)) { k <- vis[j]; if (abs(ex[k]) / abs(sc) < 3) next
      cl2 <- round(bars[j, 1] + 0.25 * (bars[j, 2] - bars[j, 1])):round(bars[j, 1] + 0.75 * (bars[j, 2] - bars[j, 1]))
      hpx <- median(sapply(cl2, function(c2) { r <- which(bm[, c2]); max(r) - min(r) + 1 })); mid <- mean(which(bm[, round(mean(bars[j, ]))]))
      val <- (if (mid < zy) 1 else -1) * hpx * abs(sc); e <- abs(val - ex[k]) / abs(sc); worst <- max(worst, e); nread <- nread + 1 }
    bw <- max(abs(tri[1] - (1.96 / sqrt(n) - coef(fy)[[1]]) / sc), abs(tri[3] - (-1.96 / sqrt(n) - coef(fy)[[1]]) / sc)); bandw <- max(bandw, bw)
    unit <- mean(diff(rowMeans(bars))); ok <- matrix(FALSE, h, w)
    for (j in seq_along(vis)) { yy <- sort(c(zy, (ex[vis[j]] - coef(fy)[[1]]) / sc)); cx <- mean(bars[j, ])
      ok[max(1, floor(yy[1]) - 3):(ceiling(yy[2]) + 3), floor(cx - 0.3 * unit - 3):ceiling(cx + 0.3 * unit + 3)] <- TRUE }   # parentheses: ':' binds tighter than '+'
    check(sprintf("(%s) no blue ink outside the expected bars", lab), !any(bm & !ok), sprintf("%d px", sum(bm & !ok))) } }
check("bar extents within 2 px of recomputed values", worst <= 2, sprintf("%d bars read; max error %.2f px", nread, worst))
check("dashed bands within 2 px of +/- 1.96/sqrt(60) in all panels", bandw <= 2, sprintf("max %.2f px", bandw))
cat(sprintf("RESULT: %d failed check(s)\n", fails)); quit(status = if (fails == 0) 0 else 1)
