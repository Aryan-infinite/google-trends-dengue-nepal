# Independent check of Supplementary Figure S2 (R) and three-program comparison.
# R reads the RAW district workbook (readxl, all yearly sheets), groups the 77 districts into development regions by zone
# membership, and recomputes per region: adjusted correlation, n_eff, Fisher 95% CI, zero months, 2023:2022 peak ratios.
# Compared with Python (out/figS2_values.csv), JavaScript (out/figS2_values_js.csv), the drawn values
# (out/figS2_plotted_values.csv) and, where present, the logged Q5 results and the four-way verification table.
suppressMessages(library(readxl))
fails <- 0
check <- function(l, ok, d = "") { cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", l, if (nzchar(d)) paste0("  (", d, ")") else "")); if (!isTRUE(ok)) fails <<- fails + 1 }
ZONE <- list(   # districts by former zone -> development region (both Nawalparasi and both Rukum districts lie within one region)
  Eastern = c("TAPLEJUNG", "PANCHTHAR", "ILAM", "JHAPA",  "SANKHUWASABHA", "TERHATHUM", "DHANKUTA", "BHOJPUR", "MORANG", "SUNSARI",
              "SOLUKHUMBU", "OKHALDHUNGA", "KHOTANG", "UDAYAPUR", "SAPTARI", "SIRAHA"),
  Central = c("DOLAKHA", "RAMECHHAP", "SINDHULI", "DHANUSHA", "DHANUSA", "MAHOTTARI", "SARLAHI", "SINDHUPALCHOK", "KAVREPALANCHOK", "LALITPUR",
              "BHAKTAPUR", "KATHMANDU", "NUWAKOT", "RASUWA", "DHADING", "MAKAWANPUR", "MAKWANPUR", "RAUTAHAT", "BARA", "PARSA", "CHITAWAN", "CHITWAN"),
  Western = c("GORKHA", "LAMJUNG", "TANAHU", "TANAHUN", "SYANGJA", "KASKI", "MANANG", "MUSTANG", "MYAGDI", "PARBAT", "BAGLUNG", "GULMI",
              "PALPA", "NAWALPARASI", "RUPANDEHI", "KAPILBASTU", "KAPILVASTU", "ARGHAKHANCHI"),
  MidWestern = c("PYUTHAN", "ROLPA", "RUKUM", "SALYAN", "DANG", "DAILEKH", "JAJARKOT", "SURKHET", "BANKE", "BARDIYA", "DOLPA", "JUMLA", "KALIKOT", "MUGU", "HUMLA"),
  FarWestern = c("BAJURA", "BAJHANG", "ACHHAM", "DOTI", "KAILALI", "DARCHULA", "BAITADI", "DADELDHURA", "KANCHANPUR"))
reg_of <- function(nm) { w <- strsplit(toupper(gsub("[^A-Za-z ]", " ", nm)), " +")[[1]]; names(ZONE)[sapply(ZONE, function(z) any(z %in% w))] }
f <- "data_raw/edcd_district_monthly_2019_2023.xlsx"; RC <- list(Central = c(), Eastern = c(), Western = c()); bad <- 0
for (y in 2020:2023) { d <- as.data.frame(read_excel(f, sheet = sprintf("%d Monthly", y), col_names = FALSE, .name_repair = "minimal"))
  h <- which(apply(d, 1, function(r) any(r %in% "District with Code"))); d <- d[(h + 1):nrow(d), ]
  d <- d[grepl("^[0-9]{3} ", d[[2]]), ]; rg <- lapply(d[[2]], reg_of); bad <- bad + sum(lengths(rg) != 1) + (nrow(d) != 77)
  m <- sapply(4:15, function(j) as.numeric(d[[j]])); m[is.na(m)] <- 0
  for (r in names(RC)) RC[[r]] <- c(RC[[r]], colSums(m[sapply(rg, function(z) identical(z, r)), , drop = FALSE])) }
check("77 districts per year, each in exactly one region (2020-2023)", bad == 0)
gt <- function(r) { L <- readLines(sprintf("data_raw/gt_dengue_region_%s.csv", tolower(r)), warn = FALSE); L <- L[grepl("^202[0-3]-[0-9]{2},", L)]
  v <- sub("^[^,]*,", "", L); v[v == "<1"] <- "0.5"; as.numeric(v) }
mo <- rep(1:12, 4); lagr <- function(v, k) cor(head(v, -k), tail(v, -k))
Rv <- do.call(rbind, lapply(names(RC), function(r) { s <- gt(r); cc <- RC[[r]]; a <- log(s + .5) - ave(log(s + .5), mo); b <- log(cc + 1) - ave(log(cc + 1), mo)
  rr <- cor(a, b); ne <- max(48 / (1 + 2 * sum(sapply(1:12, function(k) lagr(a, k) * lagr(b, k)))), 4); ci <- tanh(atanh(rr) + c(-1, 1) * 1.96 / sqrt(ne - 3))
  data.frame(region = r, r = rr, n_eff = ne, ci_lo = ci[1], ci_hi = ci[2], zeros = sum(s == 0), ratio_cases = max(cc[37:48]) / max(cc[25:36]), ratio_search = max(s[37:48]) / max(s[25:36])) }))
print(format(Rv, digits = 6), row.names = FALSE)
py <- read.csv("out/figS2_values.csv"); js <- read.csv("out/figS2_values_js.csv"); cols <- c("r", "n_eff", "ci_lo", "ci_hi", "zeros", "ratio_cases", "ratio_search")
dif <- function(A) max(abs(as.matrix(A[match(Rv$region, A$region), cols]) - as.matrix(Rv[, cols])))
check("Python = R (all values, 3 regions)", dif(py) < 1e-10, sprintf("max diff %.1e", dif(py)))
check("JavaScript = R (all values, 3 regions)", dif(js) < 1e-10, sprintf("max diff %.1e", dif(js)))
check("plotted CI equals the Fisher interval (bootstrap interval inside it, so the union is the Fisher interval)",
      all(py$boot_lo >= Rv$ci_lo - 1e-12 & py$boot_hi <= Rv$ci_hi + 1e-12))
D <- read.csv("out/figS2_plotted_values.csv", stringsAsFactors = FALSE); D <- D[match(Rv$region, D$region), ]
check("drawn r, CI ends and peak ratios equal R", max(abs(c(D$r_drawn - Rv$r, D$lo_drawn - Rv$ci_lo, D$hi_drawn - Rv$ci_hi, D$cases_drawn - Rv$ratio_cases, D$search_drawn - Rv$ratio_search))) < 1e-10)
D0 <- read.csv("out/figS2_plotted_values.csv", stringsAsFactors = FALSE)        # as written: one row per tick, in drawing order
check("rows ordered Central, Eastern, Western top to bottom, same labels in both panels",
      identical(D0$region[order(-D0$y)], c("Central", "Eastern", "Western")) && identical(D0$ytick_a, D0$region) && identical(D0$ytick_b, D0$region),
      sprintf("y = %s; labels a = %s; b = %s", paste(D0$y, collapse = ","), paste(D0$ytick_a, collapse = ","), paste(D0$ytick_b, collapse = ",")))
raw_lines <- readLines("out/figS2_plotted_values.csv", warn = FALSE)              # compare bytes: R may run in a C locale
check("axis titles exactly as intended (including the division sign, compared as UTF-8 bytes)",
      any(grepl("Correlation, seasonally adjusted (95% CI)", raw_lines, fixed = TRUE, useBytes = TRUE)) &&
      any(grepl("2023 peak \xc3\xb7 2022 peak (log scale)", raw_lines, fixed = TRUE, useBytes = TRUE)))
check("all drawn ratios inside the 0.08-15 axis", all(c(D$cases_drawn, D$search_drawn) > 0.08 & c(D$cases_drawn, D$search_drawn) < 15))
if (file.exists("out/final/final_results.xlsx")) { L <- as.data.frame(read_excel("out/final/final_results.xlsx", sheet = "Q5_regional")); L <- L[match(Rv$region, L$region), ]
  dl <- max(abs(c(L$adj_log_r - Rv$r, L$adj_log_n_eff - Rv$n_eff, L$adj_log_ci_lo - Rv$ci_lo, L$adj_log_ci_hi - Rv$ci_hi, L$cases_2023_vs_2022_peak - Rv$ratio_cases, L$search_2023_vs_2022_peak - Rv$ratio_search)))
  check("equal to the logged Q5 results", dl < 1e-10, sprintf("max diff %.1e", dl)) } else cat("[SKIP] logged results not present (repository copy)\n")
if (file.exists("verification/comparison_table.csv")) { V <- read.csv("verification/comparison_table.csv", stringsAsFactors = FALSE)
  g <- function(k) as.numeric(V$impl2_R[V$key == k]); dv <- max(abs(c(sapply(Rv$region, function(r) g(sprintf("q5_%s_r", r))) - Rv$r, sapply(Rv$region, function(r) g(sprintf("q5_%s_flo", r))) - Rv$ci_lo)))
  check("equal to the earlier four-way verification (R implementation)", dv < 1e-9, sprintf("max diff %.1e", dv)) }
cat(sprintf("RESULT: %d failed check(s)\n", fails)); quit(status = if (fails == 0) 0 else 1)
