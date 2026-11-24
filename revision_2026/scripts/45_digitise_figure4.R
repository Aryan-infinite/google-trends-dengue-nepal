# Image-level check of the DELIVERED Figure 4 TIFF, second implementation (R; ImageMagick decode; base R).
# Methods differ from scripts/44: anti-aliasing fringes removed by a neighbour-count filter (>= 6 of 8 neighbours
# the same colour); stroke positions checked by fitting the centreline of every segment that runs clear of the
# other lines (middle 40% of the segment) and measuring its perpendicular distance from the expected segment at
# both ends; tube and coverage tests as distance computations. Predictions recomputed here with lm.fit.
suppressMessages(library(readxl))
fails <- 0
check <- function(l, ok, d = "") { cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", l, if (nzchar(d)) paste0("  (", d, ")") else "")); if (!isTRUE(ok)) fails <<- fails + 1 }
runs <- function(i) { if (!length(i)) return(matrix(numeric(0), ncol = 2)); b <- which(diff(i) > 1); cbind(c(i[1], i[b + 1]), c(i[b], i[length(i)])) }
ppm <- tempfile(fileext = ".ppm"); stopifnot(system2("convert", c("manuscript/figures/Fig4.tif", "-colorspace", "sRGB", "-depth", "8", ppm)) == 0)
con <- file(ppm, "rb"); hdr <- character(0)
while (length(hdr) < 4) { ch <- character(0); repeat { s_ <- rawToChar(readBin(con, "raw", 1)); if (s_ %in% c(" ", "\n")) break; ch <- c(ch, s_) }; if (length(ch)) hdr <- c(hdr, paste(ch, collapse = "")) }
w <- as.integer(hdr[2]); h <- as.integer(hdr[3]); v <- readBin(con, "raw", w * h * 3); close(con); unlink(ppm)
chn <- function(k) matrix(as.integer(v[seq(k, length(v), 3)]), nrow = h, byrow = TRUE); R_ <- chn(1); G_ <- chn(2); B_ <- chn(3); rm(v); invisible(gc())
near <- function(rgb, t) abs(R_ - rgb[1]) + abs(G_ - rgb[2]) + abs(B_ - rgb[3]) <= t
dense <- function(m) { n <- matrix(0L, h, w)                       # count same-colour neighbours
  for (dy in -1:1) for (dx in -1:1) if (dy | dx) { ys <- max(1, 1 + dy):min(h, h + dy); xs <- max(1, 1 + dx):min(w, w + dx)
    n[ys - dy, xs - dx] <- n[ys - dy, xs - dx] + m[ys, xs] }
  m & n >= 6 }

# expected values (lm.fit)
mk <- format(seq(as.Date("2019-01-01"), by = "month", length.out = 72), "%Y-%m")
x <- as.data.frame(read_excel(Sys.getenv("PREPRINT_XLSX", "data_raw/paper_finalized_data_RP.xlsx"), col_names = FALSE, .name_repair = "minimal"))
hr <- which(apply(x, 1, function(r) any(r %in% "Confirmed Cases"))); hd <- as.character(unlist(x[hr, ]))
cases <- as.numeric(x[(hr + 1):nrow(x), match("Confirmed Cases", hd)]); cases <- cases[!is.na(cases)]
gt <- function(f) { d <- read.csv(f, skip = 2); vv <- d[[2]]; vv[vv == "<1"] <- "0.5"; as.numeric(vv) }
s <- (gt("data_raw/gt_dengue_term.csv") + gt("data_raw/gt_dengue_term_day2_2026-11-22.csv") + gt("data_raw/gt_dengue_term_day3_2026-11-23.csv")) / 3
LY <- log(cases + 1); LX <- log(s + .5); MO <- rep(1:12, 6)
D <- function(r, M) { X <- cbind(outer(MO[r], 1:12, "==") * 1, LY[r - 1]); if (M == "M2") cbind(X, LX[r]) else X }
pr <- function(tr, te, M) pmax(expm1(D(te, M) %*% lm.fit(D(tr, M), LY[tr])$coefficients), 0)
test <- 37:72; E <- list()
for (sch in c("season", "monthly")) { E[[paste(sch, "M0")]] <- cases[test - 1]; E[[paste(sch, "obs")]] <- cases[test]
  for (M in c("M1", "M2")) E[[paste(sch, M)]] <- if (sch == "season") c(pr(7:36, 37:48, M), pr(7:48, 49:60, M), pr(7:60, 61:72, M)) else sapply(test, function(t) pr(7:(t - 1), t, M)) }

COL <- list(M0 = c(134, 155, 148), M1 = c(52, 125, 255), M2 = c(248, 71, 30), obs = c(0, 0, 0)); ORD <- c("M0", "M1", "M2", "obs")
HALF <- c(M0 = 1200 / 144, M1 = 1200 / 144, M2 = 1200 / 144, obs = 1.2 * 1200 / 144); DISC <- 2.2 * 1200 / 144
raw <- lapply(COL, function(cc) near(cc, 60))
clean <- raw; clean$M0 <- dense(raw$M0); clean$obs <- dense(raw$obs)   # neighbour filter only where fringes mimic the colour (grey M0, black)
ink2 <- near(c(48, 62, 59), 30); grid <- near(c(208, 221, 211), 18)
cs <- colSums(ink2); sp <- runs(which(cs > 0.5 * max(cs))); sp <- sp[sp[, 2] - sp[, 1] + 1 >= 9, , drop = FALSE]
check("one left spine", nrow(sp) == 1)
pn <- runs(which(rowSums(ink2[, sp[1, 1]:sp[1, 2]]) == sp[1, 2] - sp[1, 1] + 1)); pn <- pn[pn[, 2] - pn[, 1] > 500, , drop = FALSE]
check("two panels", nrow(pn) == 2)
dpoly <- function(px, py, Q) { d <- rep(Inf, length(px))
  for (k in 1:(nrow(Q) - 1)) { dx <- Q[k + 1, 1] - Q[k, 1]; dy <- Q[k + 1, 2] - Q[k, 2]
    t <- pmin(pmax(((px - Q[k, 1]) * dx + (py - Q[k, 2]) * dy) / (dx^2 + dy^2), 0), 1); d <- pmin(d, sqrt((px - Q[k, 1] - t * dx)^2 + (py - Q[k, 2] - t * dy)^2)) }
  d }
for (pi_ in 1:2) {
  lab <- letters[pi_]; sch <- c("season", "monthly")[pi_]; y0 <- pn[pi_, 1]; y1 <- pn[pi_, 2]; xs1 <- sp[1, 2]
  bt <- runs(which(rowSums(ink2[y0:(y1 + 40), (xs1 + 20):(xs1 + 2000)]) > 1500)); ybot <- y0 - 1 + mean(bt[nrow(bt), ])
  tk <- colSums(ink2[(round(ybot) + 8):(round(ybot) + 50), (xs1 + 3):(w - 2)]) == 43
  tx <- c(mean(sp[1, ]), xs1 + 2 + rowMeans(runs(which(tk))))
  days <- as.numeric(seq(as.Date("2022-01-01"), as.Date("2024-10-01"), by = "quarter") - as.Date("1970-01-01"))
  check(sprintf("(%s) 12 x calibration points", lab), length(tx) == 12)
  fx <- lm(days ~ tx); hy <- (y0 - 13) + rowMeans(runs(which(rowSums(grid[(y0 - 12):y1, (xs1 + 20):(w - 20)]) > 0.5 * (w - 40 - xs1))))
  fy <- lm(c(40000, 30000, 20000, 10000, 0) ~ c(sort(hy), ybot))
  check(sprintf("(%s) scales fit within 2 px", lab), max(abs(resid(fx) / coef(fx)[2])) < 2 && max(abs(resid(fy) / coef(fy)[2])) < 2,
        sprintf("x %.2f, y %.2f px", max(abs(resid(fx) / coef(fx)[2])), max(abs(resid(fy) / coef(fy)[2]))))
  xm <- (as.numeric(seq(as.Date("2022-01-15"), by = "month", length.out = 36) - as.Date("1970-01-01")) - coef(fx)[1]) / coef(fx)[2]
  xm <- (as.numeric(as.Date(paste0(mk[test], "-01")) + 14 - as.Date("1970-01-01")) - coef(fx)[1]) / coef(fx)[2]
  P <- lapply(setNames(ORD, ORD), function(k) cbind(xm, (E[[paste(sch, k)]] - coef(fy)[1]) / coef(fy)[2]))
  inside <- matrix(FALSE, h, w); inside[(y0 - 30):(round(ybot) + 25), (xs1 + 6):(w - 5)] <- TRUE
  for (k in ORD) {
    pix <- which(clean[[k]] & inside, arr.ind = TRUE); dm <- dpoly(pix[, 2], pix[, 1], P[[k]])
    lim <- (if (k == "obs") DISC else HALF[[k]]) + 2
    if (k == "obs") { tb <- pix[, 1] < y0 + 0.2 * (ybot - y0) & pix[, 2] < xs1 + 0.55 * (w - xs1)
      check(sprintf("(%s) obs: black ink on the observed line or in the text box only", lab), max(dm[!tb]) <= lim, sprintf("max %.1f px outside text box (limit %.1f)", max(dm[!tb]), lim))
    } else check(sprintf("(%s) %s: stroke within 2 px of its expected line", lab, k), max(dm) <= lim, sprintf("%d px; max %.1f (limit %.1f)", nrow(pix), max(dm), lim))
    later <- Reduce(`|`, raw[ORD[seq_along(ORD) > match(k, ORD)]], matrix(FALSE, h, w)); cov <- raw[[k]] | later
    Q <- P[[k]]; samp <- do.call(rbind, lapply(1:35, function(j) { n <- max(2, floor(sqrt(sum((Q[j + 1, ] - Q[j, ])^2)) / 2)); t <- seq(0, 1, length.out = n)
      cbind(Q[j, 1] + t * (Q[j + 1, 1] - Q[j, 1]), Q[j, 2] + t * (Q[j + 1, 2] - Q[j, 2])) }))
    inked <- apply(round(samp), 1, function(q) any(cov[q[2] + (-3:3), q[1] + (-3:3)]))
    check(sprintf("(%s) %s: expected line fully inked", lab, k), all(inked), sprintf("%d samples; %d without ink", length(inked), sum(!inked)))
  }
  # centrelines of segments that run clear of the other series
  nseg <- 0; worst <- 0
  for (k in ORD) for (j in 1:35) {
    a <- P[[k]][j, ]; b <- P[[k]][j + 1, ]; cc <- ceiling(a[1] + 0.3 * (b[1] - a[1])):floor(a[1] + 0.7 * (b[1] - a[1]))
    mid <- cbind(cc, a[2] + (cc - a[1]) / (b[1] - a[1]) * (b[2] - a[2]))
    if (min(sapply(setdiff(ORD, k), function(k2) min(dpoly(mid[, 1], mid[, 2], P[[k2]])))) <= 35) next
    ctr <- c(); xs <- c()
    for (cx in cc) { r <- runs(which(raw[[k]][, cx] & inside[, cx])); if (nrow(r) == 1) { ctr <- c(ctr, mean(r[1, ])); xs <- c(xs, cx) } }
    if (length(xs) < 10) next
    f <- coef(lm(ctr ~ xs)); nseg <- nseg + 1
    perp <- function(px, py) abs(f[1] + f[2] * px - py) / sqrt(1 + f[2]^2)
    worst <- max(worst, perp(a[1], a[2]), perp(b[1], b[2]))
  }
  check(sprintf("(%s) clear segments: fitted centreline within 2 px of the expected segment at both ends", lab), nseg > 0 && worst <= 2,
        sprintf("%d segments; max perpendicular distance %.2f px", nseg, worst))
}
cat(sprintf("RESULT: %d failed check(s)\n", fails)); quit(status = if (fails == 0) 0 else 1)
