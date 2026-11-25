# Image-level check of the DELIVERED Supplementary Figure S2 TIFF, second implementation (R; ImageMagick decode; base R).
# Different estimators from scripts/65:
#   - y scale from the bottom spine (0) and the gridline at 1; x scale from the left spine (0) and the 12 tick marks;
#   - line values by COLUMN CENTRES: in every column of every shallow segment (|slope| < 1 px/px), away from vertices and
#     from the other line, the centre of the ink run must lie within 2 px of the expected line;
#   - stray ink: every blue / black pixel within half the line width + 2.5 px of its expected line (exact segment distance,
#     ends extended by the projecting caps);
#   - events: one column-run per marker; shape from fill ratio and the centre pixel; position by bounding box
#     (triangles by their flat top edge).
# Expected values recomputed in R from the raw files; event dates and categories from the R table in scripts/63.
suppressMessages(library(readxl))
fails <- 0
check <- function(l, ok, d = "") { cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", l, if (nzchar(d)) paste0("  (", d, ")") else "")); if (!isTRUE(ok)) fails <<- fails + 1 }
runs <- function(i) { if (!length(i)) return(matrix(numeric(0), ncol = 2)); b <- which(diff(i) > 1); cbind(c(i[1], i[b + 1]), c(i[b], i[length(i)])) }
src <- readLines("scripts/63_check_figureS2.R")
eval(parse(text = src[grep("^mk <- format", src):grep("^Rv <- data.frame", src)]))                      # cases, s, cc, crel, srel
eval(parse(text = src[grep("^CATR <- c", src):(grep("^CATR <- c", src) + 1)]))
L6 <- read.csv("out/06_2024_reporting_lag.csv", stringsAsFactors = FALSE)
ev <- rbind(data.frame(date = names(CATR), category = unname(CATR)), data.frame(date = L6$as_of[L6$as_of != "2024-12-31"], category = "count"))
dd <- as.Date(ev$date); ev$year <- as.integer(format(dd, "%Y")); mth <- as.integer(format(dd, "%m"))
ev$x <- (mth - 1) + (as.integer(format(dd, "%d")) - 0.5) / as.integer(format(as.Date(sprintf("%d-%02d-01", ev$year + (mth == 12), mth %% 12 + 1)) - 1, "%d"))

tif <- Sys.getenv("FIGS3_TIF", "manuscript/figures/FigS2.tif")
ppm <- tempfile(fileext = ".ppm"); stopifnot(system2("convert", c(tif, "-colorspace", "sRGB", "-depth", "8", ppm)) == 0)
con <- file(ppm, "rb"); hd <- character(0)
while (length(hd) < 4) { ch <- character(0); repeat { z <- rawToChar(readBin(con, "raw", 1)); if (z %in% c(" ", "\n")) break; ch <- c(ch, z) }; if (length(ch)) hd <- c(hd, paste(ch, collapse = "")) }
w <- as.integer(hd[2]); h <- as.integer(hd[3]); v <- readBin(con, "raw", w * h * 3); close(con); unlink(ppm)
k3 <- function(i) matrix(as.integer(v[seq(i, length(v), 3)]), nrow = h, byrow = TRUE); R_ <- k3(1); G_ <- k3(2); B_ <- k3(3); rm(v); invisible(gc())
near <- function(c, t) abs(R_ - c[1]) + abs(G_ - c[2]) + abs(B_ - c[3]) <= t
ink2 <- near(c(48, 62, 59), 30); grid <- near(c(208, 221, 211), 18); blu <- near(c(52, 125, 255), 60); blk <- R_ + G_ + B_ <= 60
rm(R_, G_, B_); invisible(gc())
HB <- 1200 / 72 * 1.0 / 2; HK <- 1200 / 72 * 1.2 / 2

cs <- colSums(ink2); sp <- runs(which(cs > 0.4 * max(cs))); sp <- sp[sp[, 2] - sp[, 1] + 1 >= 8, , drop = FALSE]
check("two left spines", nrow(sp) == 2)
P <- NULL
for (j in 1:2) { pr <- runs(which(rowSums(ink2[, sp[j, 1]:sp[j, 2]]) == sp[j, 2] - sp[j, 1] + 1)); pr <- pr[pr[, 2] - pr[, 1] > 800, , drop = FALSE]
  P <- rbind(P, cbind(sp[j, 1], sp[j, 2], pr)) }
P <- P[order(P[, 3], P[, 1]), ]; check("six panels", nrow(P) == 6)
segd <- function(px, py, X, Y) { d <- rep(Inf, length(px))                                          # exact distance to a polyline
  for (i in seq_len(length(X) - 1)) { vx <- X[i + 1] - X[i]; vy <- Y[i + 1] - Y[i]; t <- pmin(pmax(((px - X[i]) * vx + (py - Y[i]) * vy) / (vx^2 + vy^2), 0), 1)
    d <- pmin(d, sqrt((px - X[i] - t * vx)^2 + (py - Y[i] - t * vy)^2)) }; d }
capx <- function(X, Y, H) { n <- length(X); u0 <- c(X[1] - X[2], Y[1] - Y[2]); u0 <- u0 / sqrt(sum(u0^2)); u1 <- c(X[n] - X[n - 1], Y[n] - Y[n - 1]); u1 <- u1 / sqrt(sum(u1^2))
  list(X = c(X[1] + H * u0[1], X[2:(n - 1)], X[n] + H * u1[1]), Y = c(Y[1] + H * u0[2], Y[2:(n - 1)], Y[n] + H * u1[2])) }
colerr <- c(); ncol_ <- 0; stray <- c(); everr <- c()
for (k in 1:6) {
  yr <- 2018 + k; s0 <- P[k, 1]; s1 <- P[k, 2]; top <- P[k, 3]; bot <- P[k, 4]
  right <- if (s0 > w / 3) w - 40 else sp[2, 1] - 80
  bs <- runs(which(rowMeans(ink2[top:(bot + 40), (s1 + 30):right]) > 0.5)); yb <- top - 1 + mean(bs[nrow(bs), ])
  g1 <- runs(which(rowSums(grid[top:(top + 500), (s1 + 5):right]) > 0.1 * (right - s1))); y1 <- top - 1 + mean(g1[1, ])
  tk <- runs(which(colSums(ink2[(round(yb) + 12):(round(yb) + 40), (s1 + 2):right]) > 0)); tx <- s1 + 1 + rowMeans(tk)
  fx <- lm(tx ~ I(0:11 + 0.5)); x0 <- mean(c(s0, s1))
  check(sprintf("%d: 12 ticks; x scale agrees with the left spine (< 1.5 px)", yr), length(tx) == 12 && abs(coef(fx)[1] - x0) < 1.5, sprintf("%.2f px", abs(coef(fx)[1] - x0)))
  X <- function(m) coef(fx)[1] + coef(fx)[2] * m; Y <- function(u) yb + (y1 - yb) * u
  ytop <- floor(Y(1.25)) + 2; ylo <- floor(yb) - 7
  E <- list(search = srel[YR == yr], cases = crel[YR == yr]); M <- list(search = blu, cases = blk); H <- c(search = HB, cases = HK)
  ex <- lapply(E, function(u) { ok <- !is.na(u); list(X = X(which(ok) - 0.5), Y = Y(u[ok])) })
  for (nm in names(E)) {
    oth <- ex[[setdiff(names(E), nm)]]; e <- ex[[nm]]
    for (i in seq_len(length(e$X) - 1)) { sl <- (e$Y[i + 1] - e$Y[i]) / (e$X[i + 1] - e$X[i]); if (abs(sl) >= 1) next
      for (cx in seq(ceiling(e$X[i] + 25), floor(e$X[i + 1] - 25), by = 4)) {
        ey <- e$Y[i] + sl * (cx - e$X[i]); if (ey + H[nm] + 3 > ylo || ey - H[nm] - 3 < ytop) next
        if (min(segd(cx, ey, oth$X, oth$Y)) < 40) next
        cr <- runs(which(M[[nm]][ytop:ylo, cx])) + ytop - 1; if (!nrow(cr)) { colerr <- c(colerr, 99); next }
        cc_ <- rowMeans(cr); colerr <- c(colerr, min(abs(cc_ - ey))); ncol_ <- ncol_ + 1 } }
    cp <- capx(e$X, e$Y, H[nm]); idx <- which(M[[nm]][ytop:ylo, (s1 + 7):right], arr.ind = TRUE)
    if (nrow(idx)) stray <- c(stray, max(segd(idx[, 2] + s1 + 6, idx[, 1] + ytop - 1, cp$X, cp$Y)) - H[nm]) }
  # events
  want <- ev[ev$year == yr, ]; want <- want[order(want$x), ]
  sub <- ink2[ytop:floor(Y(1.05)), (s1 + 7):right]
  cr <- runs(which(colSums(sub) > 0))
  if (yr == 2021) { check("2021: no markers; strip ink only inside x 2.5-9.5 months (the 'No national alert found' text)",
                          nrow(want) == 0 && min(cr) + s1 + 6 > X(2.5) && max(cr) + s1 + 6 < X(9.5)); next }
  cr <- cr[cr[, 2] - cr[, 1] >= 30, , drop = FALSE]; got <- NULL
  for (j in seq_len(nrow(cr))) { b <- sub[, cr[j, 1]:cr[j, 2], drop = FALSE]; rr <- range(which(rowSums(b) > 0))
    bb <- b[rr[1]:rr[2], , drop = FALSE]; ctr <- bb[round(nrow(bb) / 2) + (-4:4), round(ncol(bb) / 2) + (-4:4)]
    kind <- if (mean(ctr) < 0.2) "count" else if (mean(bb) > 0.85) "report" else "alert"
    got <- rbind(got, data.frame(xc = s1 + 6 + mean(cr[j, ]), yc = ytop - 1 + mean(rr), ytopedge = ytop - 1 + rr[1] - 0.5, kind = kind)) }
  ok <- !is.null(got) && nrow(got) == nrow(want) && all(got$kind == want$category)
  check(sprintf("%d: %d event markers, shapes = categories", yr, nrow(want)), ok, paste(got$kind, collapse = " "))
  if (ok) everr <- c(everr, abs(got$xc - X(want$x)), ifelse(got$kind == "alert", abs(got$ytopedge - (Y(1.15) - (0.5 * 4.2 + 0.5) * 1200 / 72)), abs(got$yc - Y(1.15)))) else everr <- c(everr, 99)
}
check(sprintf("column centres: %d columns on shallow segments within 2 px", ncol_), ncol_ >= 500 && max(colerr) <= 2, sprintf("max %.2f px", max(colerr)))
check("stray ink: every blue / black data pixel within half line width + 2.5 px of its expected line", max(stray) <= 2.5, sprintf("max excess %.2f px", max(stray)))
check(sprintf("events: %d of 15 x positions and 15 y positions within 2 px", sum(seq_along(everr) <= length(everr) / 2)), length(everr) == 30 && max(everr) <= 2, sprintf("max %.2f px", max(everr)))
cat(sprintf("RESULT: %d failed check(s)\n", fails)); quit(status = if (fails == 0) 0 else 1)
