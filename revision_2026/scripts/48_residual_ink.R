# Residual-ink check, second implementation (R; ImageMagick decode; base R). Same rules as scripts/47:
# inside each plot area, (1) no ink of an unknown colour (non-white pixels that are not a figure colour and lie more
# than 3 px from one); (2) gridline grey only on gridlines; (3) dark grey / black only where the figure draws them,
# including Fig 2(b): exactly 16 dark-grey shapes besides the zero lines (2 leaders + 14 characters), counted by flood fill.
# Fig 2 unknown-colour and gridline tests are geometric and live in scripts/31 and 32.
fails <- 0
check <- function(l, ok, d = "") { cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", l, if (nzchar(d)) paste0("  (", d, ")") else "")); if (!isTRUE(ok)) fails <<- fails + 1 }
runs <- function(i) { if (!length(i)) return(matrix(numeric(0), ncol = 2)); b <- which(diff(i) > 1); cbind(c(i[1], i[b + 1]), c(i[b], i[length(i)])) }
load <- function(f) { ppm <- tempfile(fileext = ".ppm"); stopifnot(system2("convert", c(f, "-colorspace", "sRGB", "-depth", "8", ppm)) == 0)
  con <- file(ppm, "rb"); hd <- character(0)
  while (length(hd) < 4) { ch <- character(0); repeat { s <- rawToChar(readBin(con, "raw", 1)); if (s %in% c(" ", "\n")) break; ch <- c(ch, s) }; if (length(ch)) hd <- c(hd, paste(ch, collapse = "")) }
  w <- as.integer(hd[2]); h <- as.integer(hd[3]); v <- readBin(con, "raw", w * h * 3); close(con); unlink(ppm)
  k <- function(i) matrix(as.integer(v[seq(i, length(v), 3)]), nrow = h, byrow = TRUE); list(r = k(1), g = k(2), b = k(3), w = w, h = h) }
crop <- function(im, x0, y0, x1, y1) lapply(im[c("r", "g", "b")], function(m) m[y0:y1, x0:x1])
nearc <- function(s, c, t) abs(s$r - c[1]) + abs(s$g - c[2]) + abs(s$b - c[3]) <= t
dil <- function(m, n) { for (it in seq_len(n)) { d <- m; d[-1, ] <- d[-1, ] | m[-nrow(m), ]; d[-nrow(m), ] <- d[-nrow(m), ] | m[-1, ]
  d[, -1] <- d[, -1] | m[, -ncol(m)]; d[, -ncol(m)] <- d[, -ncol(m)] | m[, -1]; m <- d }; m }
nbcount <- function(m) { n <- matrix(0L, nrow(m), ncol(m)); for (dy in -1:1) for (dx in -1:1) if (dy | dx) {
  ys <- max(1, 1 + dy):min(nrow(m), nrow(m) + dy); xs <- max(1, 1 + dx):min(ncol(m), ncol(m) + dx); n[ys - dy, xs - dx] <- n[ys - dy, xs - dx] + m[ys, xs] }; n }
floodcount <- function(m) {                     # number of 8-connected components (explicit queue)
  lab <- matrix(0L, nrow(m), ncol(m)); n <- 0L; idx <- which(m)
  for (s in idx) { if (lab[s]) next; n <- n + 1L; q <- s; lab[s] <- n
    while (length(q)) { cur <- q[1]; q <- q[-1]; r <- (cur - 1) %% nrow(m) + 1; c <- (cur - 1) %/% nrow(m) + 1
      for (dy in -1:1) for (dx in -1:1) { rr <- r + dy; cc <- c + dx
        if (rr >= 1 && rr <= nrow(m) && cc >= 1 && cc <= ncol(m)) { k <- (cc - 1) * nrow(m) + rr; if (m[k] && !lab[k]) { lab[k] <- n; q <- c(q, k) } } } } }
  n }
WHITE <- c(255, 255, 255); BLACK <- c(0, 0, 0); INK2 <- c(48, 62, 59); GRID <- c(208, 221, 211)
PAL <- list(Fig1 = list(BLACK, INK2, GRID, c(232, 237, 230), c(52, 125, 255), c(248, 71, 30)), Fig2 = list(BLACK, INK2, GRID),
            Fig3 = list(BLACK, INK2, GRID, c(52, 125, 255)), Fig4 = list(BLACK, INK2, GRID, c(134, 155, 148), c(52, 125, 255), c(248, 71, 30)),
            FigS1 = list(BLACK, INK2, GRID, c(52, 125, 255)), FigS3 = list(BLACK, INK2, GRID, c(52, 125, 255)), FigS2 = list(BLACK, INK2, GRID, c(52, 125, 255)))
boxes <- function(im) { ink2 <- nearc(im, INK2, 30); cs <- colSums(ink2); b <- runs(which(cs > 0.5 * max(cs))); b <- b[b[, 2] - b[, 1] + 1 >= 9, , drop = FALSE]; out <- list()
  for (k in seq_len(nrow(b))) { e <- b[k, 2]; yr <- runs(which(rowSums(ink2[, b[k, 1]:e]) == e - b[k, 1] + 1)); yr <- yr[yr[, 2] - yr[, 1] > 500, , drop = FALSE]
    for (j in seq_len(nrow(yr))) { y0 <- yr[j, 1]; y1 <- yr[j, 2]; rr <- which(rowSums(ink2[y0:y1, (e + 20):(e + 1500)]) > 1200); yb <- y0 - 1 + runs(rr)[nrow(runs(rr)), 1]
      cl <- which(colSums(ink2[yb:(yb + 5), (e + 1):im$w]) > 0); br <- which(diff(cl) > 60); xr <- e + (if (length(br)) cl[br[1]] else cl[length(cl)])
      out[[length(out) + 1]] <- c(e + 1, y0, xr, yb) } }
  out }
for (nm in c("Fig1", "Fig2", "Fig3", "Fig4", "FigS1", "FigS3", "FigS2")) {
  if (nm == "FigS3" && !file.exists("manuscript/figures/FigS3.tif")) { cat("[SKIP] FigS3: not built (needs the EDCD district workbook, which is not public)\n"); next }
  im <- load(sprintf("manuscript/figures/%s.tif", nm)); bx <- boxes(im)
  for (k in seq_along(bx)) { lab <- letters[k]; B <- bx[[k]]
    if (nm == "FigS2") lab <- as.character(2019 + 2 * (match(B[2], sort(unique(sapply(bx, `[`, 2)))) - 1) + (B[1] > im$w / 2)); s <- crop(im, B[1], B[2], B[3] - 1, B[4] - 1); H <- nrow(s$r); Wd <- ncol(s$r)
    chroma <- pmax(s$r, s$g, s$b) - pmin(s$r, s$g, s$b); lum <- (s$r + s$g + s$b) / 3
    darkn <- dil((chroma <= 20 & lum < 170) | nearc(s, c(134, 155, 148), 40), 3)   # fringes of dark neutral ink and of the grey M0 line (decoded with a tint)
    blk <- nearc(s, BLACK, 30); ink2 <- nearc(s, INK2, 30) & !dil(nearc(s, BLACK, 60), 3)
    if (nm != "Fig2") {
      nw <- !nearc(s, WHITE, 45); pure <- Reduce(`|`, lapply(PAL[[nm]], function(cc) nearc(s, cc, 30)))
      bad <- nw & !pure & !dil(pure & nw, 3); bad <- bad & nbcount(bad) >= 2
      check(sprintf("%s(%s) no ink of an unknown colour", nm, lab), !any(bad), sprintf("%d px", sum(bad)))
      g <- nearc(s, GRID, 18) & !darkn; ok <- matrix(FALSE, H, Wd)   # S3: a gridline need show on 15% of the row (2021: 0.5 line mostly under the flat search line)
      gr <- runs(which(rowSums(nearc(s, GRID, 18)) > (if (nm == "FigS2") 0.15 else 0.5) * Wd)); gc <- runs(which(colSums(nearc(s, GRID, 18)) > 0.5 * H))
      for (i in seq_len(nrow(gr))) ok[max(1, gr[i, 1] - 7):min(H, gr[i, 2] + 7), ] <- TRUE
      for (i in seq_len(nrow(gc))) ok[, max(1, gc[i, 1] - 7):min(Wd, gc[i, 2] + 7)] <- TRUE
      g <- g & !ok; g <- g & nbcount(g) >= 2
      check(sprintf("%s(%s) gridline grey only on gridlines", nm, lab), !any(g), sprintf("%d px off gridlines", sum(g)))
    }
    o <- ink2; o[(H - 15):H, ] <- FALSE
    if (nm == "Fig1") {
      if (lab == "a") { bc <- which(colSums(nearc(s, c(232, 237, 230), 30)) > 0.3 * H); o[, 1:max(bc)] <- FALSE
        check("Fig1(a) dark grey only inside the Not recorded band", !any(o), sprintf("%d px outside", sum(o))) }
      else check(sprintf("Fig1(%s) no dark grey inside the panel", lab), !any(o), sprintf("%d px", sum(o))) }
    if (nm == "Fig2") {
      tl <- matrix(FALSE, H, Wd); tl[1:floor(0.15 * H), 1:floor(0.3 * Wd)] <- TRUE
      check(sprintf("Fig2(%s) black only in the top-left r label area", lab), !any(blk & !tl), sprintf("%d px elsewhere", sum(blk & !tl)))
      if (lab == "a") check("Fig2(a) no dark grey inside the panel", !any(o), sprintf("%d px", sum(o))) else {
        zr <- runs(which(rowSums(o) > 0.8 * Wd)); zc <- runs(which(colSums(o[21:(H - 40), ]) > 0.8 * (H - 60)))
        for (i in seq_len(nrow(zr))) o[(zr[i, 1] - 12):(zr[i, 2] + 12), ] <- FALSE
        for (i in seq_len(nrow(zc))) o[, (zc[i, 1] - 12):(zc[i, 2] + 12)] <- FALSE
        o <- o & nbcount(o) >= 3                                             # drop 1-px crumbs at cut edges
        rr <- range(which(rowSums(o) > 0)); cc <- range(which(colSums(o) > 0))
        n <- floodcount(o[rr[1]:rr[2], cc[1]:cc[2]])
        check("Fig2(b) dark-grey shapes besides the zero lines = 2 leaders + 14 characters", n == 16, sprintf("%d shapes", n)) } }
    if (nm == "Fig3") {
      zr <- runs(which(rowSums(o) > 0.5 * Wd)); for (i in seq_len(nrow(zr))) o[(zr[i, 1] - 12):(zr[i, 2] + 12), ] <- FALSE
      dr <- runs(which(rowSums(o) > 0.15 * Wd)); for (i in seq_len(nrow(dr))) o[(dr[i, 1] - 6):(dr[i, 2] + 6), ] <- FALSE
      check("Fig3 dark grey only on the zero line and the dashed threshold rows", !any(o) && nrow(dr) == 2, sprintf("%d dashed bands; %d px elsewhere", nrow(dr), sum(o))) }
    if (nm == "FigS1") { tri <- runs(which(rowSums(o) > 0.25 * Wd)); for (i in seq_len(nrow(tri))) o[max(1, tri[i, 1] - 8):min(H, tri[i, 2] + 8), ] <- FALSE
      check(sprintf("FigS1(%s) dark grey only on the zero line and the two dashed bands", lab), nrow(tri) == 3 && !any(o), sprintf("%d bands; %d px elsewhere", nrow(tri), sum(o))) }
    if (nm == "FigS3") { vl <- runs(which(colSums(o) > 0.8 * H)); for (i in seq_len(nrow(vl))) o[, max(1, vl[i, 1] - 8):min(Wd, vl[i, 2] + 8)] <- FALSE
      check(sprintf("FigS3(%s) dark grey only on the reference line (panel b)", lab), nrow(vl) == (if (lab == "b") 1 else 0) && !any(o), sprintf("%d lines; %d px elsewhere", nrow(vl), sum(o))) }
    if (nm == "FigS2") { st <- floor(0.17 * H); o[1:st, ] <- FALSE; bl <- nearc(s, c(52, 125, 255), 60)
      check(sprintf("FigS2(%s) dark grey only in the event strip; no line ink in the strip", lab), !any(o) && !any(blk[1:st, ]) && !any(bl[1:st, ]),
            sprintf("%d grey px below; %d black, %d blue px in strip", sum(o), sum(blk[1:st, ]), sum(bl[1:st, ]))) }
    if (nm == "Fig4") check(sprintf("Fig4(%s) no dark grey inside the panel", lab), !any(o), sprintf("%d px", sum(o)))
  }
  rm(im); invisible(gc())
}
cat(sprintf("RESULT: %d failed check(s)\n", fails)); quit(status = if (fails == 0) 0 else 1)
