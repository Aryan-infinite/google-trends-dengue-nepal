# Table 5, independent implementation (R; qr.solve). Nowcast metrics for every schedule (per-season refit = primary;
# monthly refit = post hoc; leave-one-season-out = secondary) and every search series (dengue 3-day mean, preprint
# Topic, day-1 download, "Fever" control). Conventions as in verification/VERIFICATION_SPEC.md section 6, extended to all
# models by season and to the leave-one-season-out season count. Bootstrap limits here use R's own RNG (a Monte Carlo
# check of the pipeline's limits, not the reported values). Output: out/table5_R.csv (key, value).
rd <- function(f) { L <- readLines(paste0("data_raw/", f), encoding = "UTF-8"); L <- L[grepl("^20[0-9]{2}-[0-9]{2},", L)]
  v <- trimws(sub("^[^,]*,", "", L)); v[v == "<1"] <- "0.5"; as.numeric(v) }
f3 <- c("gt_dengue_term.csv", "gt_dengue_term_day2_2026-11-22.csv", "gt_dengue_term_day3_2026-11-23.csv")
pt_ <- read.csv("out/paper_72month_table.csv"); cases <- pt_$cases
SER <- list(dengue = Reduce(`+`, lapply(f3, rd)) / 3, topic = pt_$DF, day1 = rd(f3[1]), fever = pt_$Fever)
YR <- rep(2019:2024, each = 12); MO <- rep(1:12, 6); YM <- YR * 100 + MO; W <- function(lo, hi) which(YM >= lo & YM <= hi)
LY <- log(cases + 1); K <- list()
put <- function(k, v) K[[k]] <<- v
for (sn in names(SER)) { LX <- log(SER[[sn]] + 0.5)
  des <- function(r, M) { X <- cbind(outer(MO[r], 1:12, "==") * 1, LY[r - 1]); if (M == "M2") X <- cbind(X, LX[r]); if (M == "M3") X <- cbind(X, LX[r - 1]); X }
  pred <- function(tr, te, M) pmax(expm1(des(te, M) %*% qr.solve(des(tr, M), LY[tr])), 0)
  ii <- W(201907, 202412); tm <- W(202201, 202412); lo <- W(202002, 202412)
  sch <- list(ps = list(lapply(2022:2024, function(y) tm[YR[tm] == y]), function(g) ii[ii < min(g)]),
              mo = list(as.list(tm), function(g) ii[ii < g]),
              loso = list(lapply(2020:2024, function(y) lo[YR[lo] == y]), function(g) lo[YR[lo] != YR[g[1]]]))
  for (s in names(sch)) { P <- list(M0 = c(), M1 = c(), M2 = c(), M3 = c()); obs <- c(); yrs <- c()
    for (g in sch[[s]][[1]]) { tr <- sch[[s]][[2]](g); P$M0 <- c(P$M0, cases[g - 1]); obs <- c(obs, cases[g]); yrs <- c(yrs, YR[g])
      for (M in c("M1", "M2", "M3")) P[[M]] <- c(P[[M]], pred(tr, g, M)) }
    for (M in names(P)) { e <- P[[M]] - obs; key <- function(z) sprintf("%s_%s_%s_%s", sn, s, z, M)
      put(key("mae"), mean(abs(e))); put(key("rmse"), sqrt(mean(e^2))); put(key("maelog"), mean(abs(log1p(P[[M]]) - log1p(obs))))
      for (y in unique(yrs)) put(sprintf("%s_%s_mae%d_%s", sn, s, y, M), mean(abs(e[yrs == y])))
      o <- obs > 1000; f <- P[[M]] > 1000; put(key("sens1000"), sum(o & f) / sum(o)); put(key("spec1000"), sum(!o & !f) / sum(!o)) }
    d <- abs(P$M2 - obs) - abs(P$M1 - obs); put(sprintf("%s_%s_diff", sn, s), mean(d))
    put(sprintf("%s_%s_better", sn, s), sum(sapply(unique(yrs), function(y) mean(abs(P$M2 - obs)[yrs == y]) < mean(abs(P$M1 - obs)[yrs == y]))))
    put(sprintf("%s_%s_nseasons", sn, s), length(unique(yrs))); put(sprintf("%s_%s_ntest", sn, s), length(obs))
    set.seed(4000); bm <- replicate(4000, mean(d[sample.int(length(d), replace = TRUE)]))     # R's own RNG: Monte Carlo check only
    put(sprintf("%s_%s_bootlo_Rrng", sn, s), quantile(bm, .025, names = FALSE)); put(sprintf("%s_%s_boothi_Rrng", sn, s), quantile(bm, .975, names = FALSE)) } }
write.csv(data.frame(key = names(K), value = sprintf("%.17g", unlist(K))), "out/table5_R.csv", row.names = FALSE)
cat(length(K), "keys\n"); for (k in c("dengue_ps_mae_M2", "dengue_ps_diff", "dengue_mo_diff", "dengue_loso_diff", "dengue_loso_better", "dengue_ps_mae2022_M3")) cat(k, K[[k]], "\n")
