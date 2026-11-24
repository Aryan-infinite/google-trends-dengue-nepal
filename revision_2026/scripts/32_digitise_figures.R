# Image-level check of the DELIVERED TIFFs, second implementation (R; base R only).
# Decoding: ImageMagick converts each CMYK TIFF to an 8-bit PPM with its own decoder; R reads the bytes.
# Methods differ from scripts/31_digitise_figures.py on purpose:
#   Figure 1: each month's value = where the fitted centrelines of the two neighbouring line segments
#             reach that month (straight middle parts of segments only; markers are not used).
#   Figure 2: isolated markers read from their bounding box (circles) or flat base (triangles);
#             leader lines tested by looking for leader ink in a ring around every one of the 60 points.
# Run from the project folder: Rscript scripts/32_digitise_figures.R
fails <- 0
check <- function(label, ok, detail = "") {
  cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", label, if (nzchar(detail)) paste0("  (", detail, ")") else ""))
  if (!isTRUE(ok)) fails <<- fails + 1
}
read_ppm <- function(tif) {
  ppm <- tempfile(fileext = ".ppm")
  stopifnot(system2("convert", c(shQuote(tif), "-colorspace", "sRGB", "-depth", "8", shQuote(ppm))) == 0)
  con <- file(ppm, "rb"); on.exit({close(con); unlink(ppm)})
  tok <- character(0)
  while (length(tok) < 4) {                                   # P6, width, height, maxval
    ch <- character(0)
    repeat { b <- readBin(con, "raw", 1); s <- rawToChar(b); if (s %in% c(" ", "\n", "\r", "\t")) break; ch <- c(ch, s) }
    if (length(ch)) tok <- c(tok, paste(ch, collapse = ""))
  }
  stopifnot(tok[1] == "P6", tok[4] == "255")
  w <- as.integer(tok[2]); h <- as.integer(tok[3])
  v <- readBin(con, "raw", w * h * 3)
  ch <- function(k) matrix(as.integer(v[seq(k, length(v), 3)]), nrow = h, ncol = w, byrow = TRUE)
  list(r = ch(1), g = ch(2), b = ch(3), w = w, h = h)
}
near <- function(im, rgb, tol) abs(im$r - rgb[1]) + abs(im$g - rgb[2]) + abs(im$b - rgb[3]) <= tol
runs <- function(idx) {
  if (!length(idx)) return(matrix(numeric(0), ncol = 2))
  br <- which(diff(idx) > 1); cbind(c(idx[1], idx[br + 1]), c(idx[br], idx[length(idx)]))
}
fitres <- function(px, val) { f <- lm(val ~ px); list(a = coef(f)[[1]], b = coef(f)[[2]], res = max(abs(resid(f) / coef(f)[[2]]))) }
match_ticks <- function(pxs, ticks, desc) {
  pxs <- sort(pxs); tv <- sort(ticks); if (desc) tv <- rev(tv); best <- NULL
  for (k in 0:(length(tv) - length(pxs))) { f <- fitres(pxs, tv[(k + 1):(k + length(pxs))]); if (is.null(best) || f$res < best$res) best <- f }
  best
}

# ---------- expected values from the raw files ----------
suppressMessages(library(readxl))
months <- seq(as.Date("2019-01-01"), by = "month", length.out = 72); mk <- format(months, "%Y-%m")
x <- as.data.frame(read_excel(Sys.getenv("PREPRINT_XLSX", "data_raw/paper_finalized_data_RP.xlsx"), col_names = FALSE, .name_repair = "minimal"))
hr <- which(apply(x, 1, function(r) any(r %in% "Confirmed Cases"))); hd <- as.character(unlist(x[hr, ]))
body <- x[(hr + 1):nrow(x), ]; body <- body[!is.na(body[[match("Confirmed Cases", hd)]]), ]
cases <- as.numeric(body[[match("Confirmed Cases", hd)]]); fever <- as.numeric(body[[match("Fever", hd)]])
gt <- function(f) { d <- read.csv(f, skip = 2, stringsAsFactors = FALSE); v <- d[[2]]; v[v == "<1"] <- "0.5"; as.numeric(v) }
search <- (gt("data_raw/gt_dengue_term.csv") + gt("data_raw/gt_dengue_term_day2_2026-11-22.csv") + gt("data_raw/gt_dengue_term_day3_2026-11-23.csv")) / 3
casesA <- cases; casesA[1:5] <- NA
day <- function(d) as.numeric(d - as.Date("1970-01-01"))

# ================= FIGURE 1 =================
im <- read_ppm("manuscript/figures/Fig1.tif")
INK2 <- c(48, 62, 59); GRID <- c(208, 221, 211)
ink2 <- near(im, INK2, 30); grid <- near(im, GRID, 18)
cs <- colSums(ink2); sp <- runs(which(cs > 0.5 * max(cs)))
check("F1: one left spine", nrow(sp) == 1)
xl <- mean(sp[1, ]); spc <- sp[1, 1]:sp[1, 2]
pr <- runs(which(rowSums(ink2[, spc]) == length(spc))); pr <- pr[pr[, 2] - pr[, 1] > 500, , drop = FALSE]
check("F1: three panels", nrow(pr) == 3)
gc <- colSums(grid[pr[2, 1]:pr[2, 2], ]); vx <- rowMeans(runs(which(gc > 0.6 * (pr[2, 2] - pr[2, 1]))))
check("F1: five year gridlines", length(vx) == 5)
fx <- fitres(c(xl, vx), day(as.Date(paste0(2019:2024, "-01-01"))))
check("F1: x scale residual < 2 px", fx$res < 2, sprintf("%.2f px", fx$res))
d2px <- function(d) (day(d) - fx$a) / fx$b
vxm <- d2px(months + 14)
cols <- list(c(0, 0, 0), c(52, 125, 255), c(248, 71, 30))
vals <- list(casesA, search, fever); ticks <- list(c(10, 100, 1000, 10000), c(25, 50, 75, 100), c(10, 20, 30, 40))
base <- c(1, 0, 0); logy <- c(TRUE, FALSE, FALSE)
for (p in 1:3) {
  lab <- letters[p]; y0 <- pr[p, 1]; y1 <- pr[p, 2]
  rs <- rowSums(ink2[y0:(y1 + 40), (round(xl) + 20):(round(xl) + 2000)]); bt <- runs(which(rs > 1500))
  ybot <- y0 - 1 + mean(bt[nrow(bt), ])
  rg <- rowSums(grid[(y0 - 10):y1, (round(xl) + 20):(im$w - 20)])
  hy <- (y0 - 11) + rowMeans(runs(which(rg > 0.5 * (im$w - 40 - xl))))
  tv <- if (logy[p]) log10(c(ticks[[p]], base[p])) else c(ticks[[p]], base[p])
  fy <- fitres(c(sort(hy), ybot), c(sort(tv[1:4], decreasing = TRUE), tv[5]))
  check(sprintf("F1(%s): 4 gridlines, y scale residual < 2 px", lab), length(hy) == 4 && fy$res < 2, sprintf("%.2f px", fy$res))
  m <- near(im, cols[[p]], 60); m[1:(y0 - 30), ] <- FALSE; m[(round(ybot) + 30):im$h, ] <- FALSE; m[, 1:(round(xl) + 5)] <- FALSE
  v <- vals[[p]]; ok <- which(!is.na(v))
  ey <- ((if (logy[p]) log10(v) else v) - fy$a) / fy$b
  # segment centrelines from the middle 40% of each segment (straight parts only; markers not used)
  L <- list()                                   # fitted line through segment k: point (x0, y0) and direction (1, b)
  for (k in seq_along(ok)[-length(ok)]) {
    i <- ok[k]; j <- ok[k + 1]; xa <- vxm[i]; xb <- vxm[j]
    cc <- ceiling(xa + 0.3 * (xb - xa)):floor(xa + 0.7 * (xb - xa)); ctr <- c(); xs <- c()
    for (cx in cc) { r <- runs(which(m[, cx])); if (nrow(r) == 1) { ctr <- c(ctr, mean(r[1, ])); xs <- c(xs, cx) } }
    L[[k]] <- if (length(xs) >= 10) { f <- lm(ctr ~ xs); c(a = coef(f)[[1]], b = coef(f)[[2]]) } else NULL
  }
  perp <- function(l, x0, y0) abs(l[["a"]] + l[["b"]] * x0 - y0) / sqrt(1 + l[["b"]]^2)
  err <- sapply(seq_along(ok), function(k) {
    i <- ok[k]; x0 <- vxm[i]; y0e <- ey[i]
    l1 <- if (k > 1) L[[k - 1]] else NULL; l2 <- if (k < length(ok)) L[[k]] else NULL
    if (is.null(l1) && is.null(l2)) return(Inf)
    if (!is.null(l1) && !is.null(l2)) {
      ang <- abs(atan(l1[["b"]]) - atan(l2[["b"]])) * 180 / pi
      if (ang >= 15) {                          # well-conditioned: intersect the two centrelines
        xi <- (l2[["a"]] - l1[["a"]]) / (l1[["b"]] - l2[["b"]]); yi <- l1[["a"]] + l1[["b"]] * xi
        return(sqrt((xi - x0)^2 + (yi - y0e)^2))
      }
      return(max(perp(l1, x0, y0e), perp(l2, x0, y0e)))
    }
    perp(if (is.null(l1)) l2 else l1, x0, y0e)
  })
  check(sprintf("F1(%s): every month's position from segment centrelines within 2 px (2-D)", lab), max(err) <= 2,
        sprintf("%d points; max error %.2f px; median %.2f px", length(ok), max(err), median(err)))
  if (p == 1) {
    empty <- all(sapply(1:5, function(i) !any(m[y0:round(ybot), round(vxm[i]) + (-2:2)])))
    check("F1(a): no ink in Jan-May 2019 columns", empty)
  }
  # no extra ink: every series pixel within 19 px of the expected polyline
  pix <- which(m, arr.ind = TRUE); P <- cbind(vxm[ok], ey[ok]); dmin <- rep(Inf, nrow(pix))
  for (k in 1:(nrow(P) - 1)) {
    dx <- P[k + 1, 1] - P[k, 1]; dy <- P[k + 1, 2] - P[k, 2]
    t <- pmin(pmax(((pix[, 2] - P[k, 1]) * dx + (pix[, 1] - P[k, 2]) * dy) / (dx^2 + dy^2), 0), 1)
    dmin <- pmin(dmin, sqrt((pix[, 2] - P[k, 1] - t * dx)^2 + (pix[, 1] - P[k, 2] - t * dy)^2))
  }
  samp <- do.call(rbind, lapply(1:(nrow(P) - 1), function(k) { n <- max(2, floor(sqrt(sum((P[k + 1, ] - P[k, ])^2)) / 2))
    t <- seq(0, 1, length.out = n); cbind(P[k, 1] + t * (P[k + 1, 1] - P[k, 1]), P[k, 2] + t * (P[k + 1, 2] - P[k, 2])) }))
  inked <- apply(round(samp), 1, function(q) any(m[q[2] + (-3:3), q[1] + (-3:3)]))
  check(sprintf("F1(%s): expected line fully inked (sampled every 2 px)", lab), all(inked), sprintf("%d samples; %d without ink", length(inked), sum(!inked)))
  check(sprintf("F1(%s): no series pixel off the expected line", lab), max(dmin) <= 19, sprintf("%d px; max %.1f", nrow(pix), max(dmin)))
}
rm(im, ink2, grid, m); invisible(gc())

# ================= FIGURE 2 =================
im <- read_ppm("manuscript/figures/Fig2.tif")
ink2 <- near(im, INK2, 30); grid <- near(im, GRID, 18)
cs <- colSums(ink2); sp <- runs(which(cs > 0.6 * max(cs))); sp <- sp[sp[, 2] - sp[, 1] + 1 >= 9, , drop = FALSE]
check("F2: two left spines", nrow(sp) == 2)
w <- 13:72; cs60 <- cases[w]; ss <- search[w]; mon <- as.integer(format(months[w], "%m")); late <- months[w] >= as.Date("2022-01-01")
Y <- log(cs60 + 1); X <- log(ss + 0.5); Ya <- Y - ave(Y, mon); Xa <- X - ave(X, mon)
d <- data.frame(y = log1p(cs60), y1 = c(NA, head(log1p(cs60), -1)), x = X, s1 = sin(2*pi*mon/12), c1 = cos(2*pi*mon/12), s2 = sin(4*pi*mon/12), c2 = cos(4*pi*mon/12))
kp <- !is.na(d$y1); cdist <- cooks.distance(lm(y ~ s1 + c1 + s2 + c2 + y1 + x, d[kp, ]))
flagged <- which(kp)[cdist > 4 / sum(kp)]
check("F2: Cook flags 2020-06 and 2022-09", identical(format(months[w][flagged], "%Y-%m"), c("2020-06", "2022-09")))
R <- sqrt(14) / 2 * 1200 / 72; E <- 0.5 / 2 * 1200 / 72                 # marker radius; half white-edge width
chroma <- pmax(im$r, im$g, im$b) - pmin(im$r, im$g, im$b)    # stacked translucent markers darken; keep any saturated ink
blue <- chroma > 60 & im$b > im$r + 20; orng <- chroma > 60 & im$r > im$b + 20
spec <- list(list("a", X, Y, log(c(1, 2, 5, 10, 25, 100) + 0.5), log(c(1, 10, 100, 1000, 10000) + 1)),
             list("b", Xa, Ya, -2:2, -3:2))
for (p in 1:2) {
  lab <- spec[[p]][[1]]; xv <- spec[[p]][[2]]; yv <- spec[[p]][[3]]
  bx1 <- sp[p, 2]; spc <- sp[p, 1]:sp[p, 2]
  yr <- runs(which(rowSums(ink2[, spc]) == length(spc))); yr <- yr[yr[, 2] - yr[, 1] > 500, , drop = FALSE]
  y0 <- yr[1, 1]; y1 <- yr[1, 2]; right <- if (p == 1) sp[2, 1] - 400 else im$w - 20
  cols_ <- (bx1 + 5):right
  vxs <- (bx1 + 4) + rowMeans(runs(which(colSums(grid[y0:y1, cols_]) > 0.5 * (y1 - y0))))
  hys <- (y0 - 1) + rowMeans(runs(which(rowSums(grid[y0:y1, cols_]) > 0.5 * length(cols_))))
  if (p == 2) {
    zr <- runs(which(rowSums(ink2[y0:y1, (bx1 + 40):right]) > 0.8 * (right - bx1 - 40))); zr <- zr[(y0 - 1 + zr[, 2]) < y1 - 30, , drop = FALSE]
    zc <- runs(which(colSums(ink2[(y0 + 20):(y1 - 40), cols_]) > 0.8 * (y1 - y0 - 60))); zc <- zc[zc[, 2] - zc[, 1] + 1 <= 9, , drop = FALSE]
    check("F2(b): one zero line each way", nrow(zr) == 1 && nrow(zc) == 1)
    zy <- y0 - 1 + mean(zr[1, ]); zx <- bx1 + 4 + mean(zc[1, ])
    vxs <- sort(c(vxs, zx)); hys <- sort(c(hys, zy))
  }
  fxx <- match_ticks(vxs, spec[[p]][[4]], FALSE); fyy <- match_ticks(hys, spec[[p]][[5]], TRUE)
  check(sprintf("F2(%s): gridlines matched to ticks, residual < 2 px", lab), fxx$res < 2 && fyy$res < 2, sprintf("%.2f / %.2f px", fxx$res, fyy$res))
  px <- (xv - fxx$a) / fxx$b; py <- (yv - fyy$a) / fyy$b
  inside <- matrix(FALSE, im$h, im$w); inside[y0:y1, cols_] <- TRUE
  bm <- blue & inside; om <- orng & inside
  crowd <- sapply(1:60, function(i) sum(sqrt((px - px[i])^2 + (py - py[i])^2) < 2.2 * R) - 1)
  frac <- function(mk, i) { g <- expand.grid(dy = -16:16, dx = -16:16); g <- g[sqrt(g$dx^2 + g$dy^2) <= 0.5 * R, ]
    mean(mk[cbind(round(py[i]) + g$dy, round(px[i]) + g$dx)]) }
  ownf <- sapply(1:60, function(i) frac(if (late[i]) om else bm, i)); anyf <- sapply(1:60, function(i) frac(om | bm, i))
  present <- ifelse(crowd == 0, ownf > 0.5, anyf > 0.2)
  check(sprintf("F2(%s): marker ink at all 60 positions", lab), all(present), sprintf("%d isolated, %d overlapping", sum(crowd == 0), sum(crowd > 0)))
  for (cl in c("blue", "orange")) {
    mk <- if (cl == "blue") bm else om; g <- if (cl == "blue") !late else late
    pix <- which(mk, arr.ind = TRUE); dd <- rep(Inf, nrow(pix))
    for (i in which(g)) dd <- pmin(dd, sqrt((pix[, 2] - px[i])^2 + (pix[, 1] - py[i])^2))
    check(sprintf("F2(%s): no %s pixel outside expected markers", lab, cl), max(dd) <= 1.3 * R, sprintf("%d px; max %.1f", nrow(pix), max(dd)))
  }
  # residual ink (geometric): non-white pixels other than gridlines, zero lines, dark label text/leaders and their fringes
  # must lie within marker reach of an expected point (independent re-implementation of the rule in scripts/31)
  nw <- (abs(im$r - 255) + abs(im$g - 255) + abs(im$b - 255)) > 45
  dn <- (abs(im$r - 48) + abs(im$g - 62) + abs(im$b - 59) <= 30) | (im$r + im$g + im$b <= 60)   # axis grey and black only (not any grey)
  sk <- dn; for (it in 1:3) { s2 <- sk; s2[-1, ] <- s2[-1, ] | sk[-nrow(sk), ]; s2[-nrow(sk), ] <- s2[-nrow(sk), ] | sk[-1, ]
    s2[, -1] <- s2[, -1] | sk[, -ncol(sk)]; s2[, -ncol(sk)] <- s2[, -ncol(sk)] | sk[, -1]; sk <- s2 }
  for (gy in hys) sk[round(gy) + (-7:7), ] <- TRUE; for (gx in vxs) sk[, round(gx) + (-7:7)] <- TRUE
  cand <- nw & inside & !sk; cand[(y1 - 15):y1, ] <- FALSE
  nb <- matrix(0L, nrow(cand), ncol(cand))                              # drop isolated single pixels (< 2 neighbours)
  for (dy in -1:1) for (dx in -1:1) if (dy | dx) { ys <- max(1, 1 + dy):min(nrow(cand), nrow(cand) + dy); xs <- max(1, 1 + dx):min(ncol(cand), ncol(cand) + dx)
    nb[ys - dy, xs - dx] <- nb[ys - dy, xs - dx] + cand[ys, xs] }
  cand <- cand & nb >= 2
  pix <- which(cand, arr.ind = TRUE); dd <- rep(Inf, nrow(pix)); for (i in 1:60) dd <- pmin(dd, sqrt((pix[, 2] - px[i])^2 + (pix[, 1] - py[i])^2))
  check(sprintf("F2(%s): every other non-white pixel within marker reach of an expected point", lab), nrow(pix) == 0 || max(dd) <= 1.3 * R + 6,
        sprintf("%d px; max distance %.1f (limit %.0f)", nrow(pix), if (nrow(pix)) max(dd) else 0, 1.3 * R + 6))
  # isolated markers: bounding box (circles) or flat base (triangles)
  er <- c()
  far <- sapply(1:60, function(i) sum(sqrt((px - px[i])^2 + (py - py[i])^2) < 3 * R) - 1)
  for (i in which(far == 0)) {
    mk <- if (late[i]) om else bm; rr <- round(py[i]) + (-45:45); cc <- round(px[i]) + (-45:45)
    sub <- which(mk[rr, cc], arr.ind = TRUE); yy <- rr[sub[, 1]]; xx <- cc[sub[, 2]]
    ex <- (min(xx) + max(xx)) / 2; eyy <- if (late[i]) max(yy) - (R - E) else (min(yy) + max(yy)) / 2
    er <- c(er, sqrt((ex - px[i])^2 + (eyy - py[i])^2))
  }
  check(sprintf("F2(%s): isolated markers read back within 2 px (bounding box / base)", lab), max(er) <= 2,
        sprintf("%d points; max error %.2f px", length(er), max(er)))
  if (p == 2) {
    lm_ <- ink2 & inside; lm_[round(zy) + (-12:12), ] <- FALSE; lm_[, round(zx) + (-12:12)] <- FALSE
    ring <- sapply(1:60, function(i) { g <- expand.grid(dy = -50:50, dx = -50:50); dd <- sqrt(g$dx^2 + g$dy^2)
      g <- g[dd >= 30 & dd <= 50, ]; sum(lm_[cbind(round(py[i]) + g$dy, round(px[i]) + g$dx)]) })
    check("F2(b): leader ink next to exactly the two flagged months and no other point",
          identical(sort(which(ring > 0)), sort(flagged)),
          sprintf("leader ink near: %s", paste(format(months[w][which(ring > 0)], "%Y-%m"), collapse = ", ")))
  }
}
cat(sprintf("RESULT: %d failed check(s)\n", fails))
quit(status = if (fails == 0) 0 else 1)
