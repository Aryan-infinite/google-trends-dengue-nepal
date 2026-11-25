# Independent check of manuscript/Tables_4-5.docx (R). Reads the Word file itself (unzip + XML text of every cell and
# note), formats R's own recomputation (scripts/68 for Table 4, scripts/70 for Table 5) with R's own formatting rules,
# and requires every data cell to match. Also checks against the pipeline (scripts/09, 15 outputs), the shared-resample
# verification table, the decision-rule arithmetic, and the facts stated in the notes.
suppressMessages(library(readxl)); invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))   # UTF-8 strings (minus signs, dashes) compared as text
fails <- 0
check <- function(l, ok, d = "") { cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", l, if (nzchar(d)) paste0("  (", d, ")") else "")); if (!isTRUE(ok)) fails <<- fails + 1 }
hex <- function(v) vapply(v, function(z) paste(charToRaw(enc2utf8(z)), collapse = ""), "", USE.NAMES = FALSE)   # locale-proof comparison
eq <- function(a, b) length(a) == length(b) && all(hex(a) == hex(b))
MINUS <- "−"
num <- function(v, d) gsub("-", MINUS, formatC(v, format = "f", digits = d, big.mark = ","), fixed = TRUE)
civ <- function(lo, hi, d) paste0(num(lo, d), " to ", num(hi, d))
pv <- function(p) if (p < 0.001) "<0.001" else if (p < 0.01) sprintf("%.3f", p) else sprintf("%.2f", p)

# ---- read the document ----
f <- Sys.getenv("TABLES45_DOCX", "manuscript/Tables_4-5.docx"); td <- tempfile(); unzip(f, "word/document.xml", exdir = td)
x <- paste(readLines(file.path(td, "word/document.xml"), warn = FALSE, encoding = "UTF-8"), collapse = "")
txt <- function(s) { t <- regmatches(s, gregexpr("<w:t(?: [^>]*)?>[^<]*</w:t>", s, perl = TRUE))[[1]]
  u <- paste(sub("^<w:t(?: [^>]*)?>", "", sub("</w:t>$", "", t)), collapse = ""); u <- gsub("&lt;", "<", gsub("&gt;", ">", gsub("&quot;", "\"", gsub("&amp;", "&", u, fixed = TRUE), fixed = TRUE), fixed = TRUE), fixed = TRUE)
  Encoding(u) <- "UTF-8"; trimws(u) }
tabs <- regmatches(x, gregexpr("<w:tbl>.*?</w:tbl>", x, perl = TRUE))[[1]]
check("two tables in the document", length(tabs) == 2)
cells <- function(tb) lapply(regmatches(tb, gregexpr("<w:tr>.*?</w:tr>|<w:tr .*?</w:tr>", tb, perl = TRUE))[[1]],
                             function(r) vapply(regmatches(r, gregexpr("<w:tc>.*?</w:tc>", r, perl = TRUE))[[1]], txt, ""))
body <- gsub("<w:tbl>.*?</w:tbl>", "", x, perl = TRUE)
paras <- vapply(regmatches(body, gregexpr("<w:p[ >].*?</w:p>", body, perl = TRUE))[[1]], txt, "", USE.NAMES = FALSE); paras <- paras[nzchar(paras)]

# ---- Table 4 ----
T4 <- read.csv("out/table4_runs_R.csv", stringsAsFactors = FALSE); N4 <- read.csv("out/table4_runs_np.csv", stringsAsFactors = FALSE)
check("Table 4 inputs: R (68) = NumPy (69), every column, 1e-9", all(sapply(names(T4)[-1], function(k) all(abs(T4[[k]] - N4[[k]]) <= 1e-9 * pmax(1, abs(T4[[k]]))))))
G <- as.data.frame(read_excel("out/final/final_results.xlsx", "Robustness_grid"))
gl <- c("primary (2020\u201324, 3-day mean)", "A: 2019 Jun\u20132024", "B: full 72 months", "C: excl. 2020\u201321 (COVID), 2022\u201324", "D: excl. 2022",
        "Topic instead of term", "Day-1 download only", "\"<1\" coded 0", "Control: \"Fever\"")
G <- G[match(hex(gl), hex(G$run)), ]
check("Table 4: pipeline grid = R for r, lags, Holm lag +1, Wald P (1e-8)", !anyNA(G$n) &&
      max(abs(c(G$adj_r - T4$r, G$lag_minus1 - T4$lag_m1, G$lag0 - T4$lag_0, G$lag_plus1 - T4$lag_p1, G$lag_plus1_holm - T4$holm_p1, G$wald_x1_p - T4$fwd_p, G$ci_lo - T4$ci_lo))) < 1e-8)
CT <- read.csv("verification/comparison_table.csv", stringsAsFactors = FALSE); ct <- function(k) as.numeric(CT$impl2_R[CT$key == k])
hi <- ifelse(abs(G$ci_hi - T4$ci_hi) < 1e-8, T4$ci_hi, G$ci_hi)
check("Table 4: CI upper limit = Fisher in 8 analyses; Fever = pipeline bootstrap limit, within Monte Carlo 0.03 of the shared-resample value",
      sum(abs(G$ci_hi - T4$ci_hi) < 1e-8) == 8 && T4$run[abs(G$ci_hi - T4$ci_hi) >= 1e-8] == "fever" && abs(G$ci_hi[9] - ct("q1_fever_bhi")) < 0.03,
      sprintf("pipeline %.4f, shared %.4f", G$ci_hi[9], ct("q1_fever_bhi")))
lab4 <- c("Primary analysis, January 2020–December 2024", "(i) Including June–December 2019", "(ii) All 72 months, as in the preprint",
          "(iii) 2022–2024 only", "(iv) Excluding 2022", "(v) Preprint Topic series", "(vi) First retrieval only", "(vii) Values \"<1\" coded 0", "Control series \"Fever\"")
met <- T4$holm_p1 < 0.05 & T4$lag_p1 > T4$lag_m1 & T4$fwd_p < 0.05
want4 <- lapply(1:9, function(i) with(T4[i, ], c(lab4[i], as.character(n), paste0(num(r, 2), " (", civ(ci_lo, hi[i], 2), ")"), num(n_eff, 1),
  paste0(num(lag_m1, 2), " (", pv(holm_m1), ")"), paste0(num(lag_0, 2), " (", pv(holm_0), ")"), paste0(num(lag_p1, 2), " (", pv(holm_p1), ")"),
  paste0(num(fwd_coef, 2), " (", pv(fwd_p), ")"), if (met[i]) "Yes" else "No")))
C4 <- cells(tabs[1]); D4 <- C4[sapply(C4, length) == 9][-1]                    # data rows (9 cells), header row 2 removed
check("Table 4: 9 data rows, every cell = R's own value and formatting", length(D4) == 9 && all(mapply(eq, D4, want4)),
      if (length(D4) == 9) paste(which(!mapply(eq, D4, want4)), collapse = " ") else sprintf("%d rows", length(D4)))
check("Table 4: no analysis meets the lead rule (R arithmetic)", !any(met))
check("Table 4 header labels", eq(C4[[2]], c("Analysis", "Months", "r (95% CI)ᵃ", "neff", "Lag −1", "Lag 0", "Lag +1",
      "Lagged search → cases, coefficient (P)ᶜ", "Lead rule metᵈ")) && grepl("Prewhitened cross-correlation", C4[[1]][2]))
# notes
nD <- paras[grepl("^Analysis \\(iv\\)", paras)]
check("Table 4 note: run (iv) coefficient and P as computed", length(nD) == 1 && grepl(paste0("(", num(T4$fwd_coef[5], 2), "; P = ", pv(T4$fwd_p[5]), ")"), nD, fixed = TRUE) && T4$fwd_coef[5] < 0)
lt <- sapply(c("gt_dengue_term.csv", "gt_dengue_term_day2_2026-11-22.csv", "gt_dengue_term_day3_2026-11-23.csv"), function(fn) {
  L <- readLines(paste0("data_raw/", fn), encoding = "UTF-8"); L <- L[grepl("^20[0-9]{2}-[0-9]{2},", L)]; paste(sub(",.*", "", L[trimws(sub("^[^,]*,", "", L)) == "<1"]), collapse = " ") })
check("Table 4 note: the only \"<1\" value in the primary series is April 2019 (all three downloads)", all(lt %in% c("2019-04", "")) && any(lt == "2019-04"), paste(lt, collapse = " | "))
nA <- paras[grepl("^ᵃ Pearson", paras)]
check("Table 4 note a: Fever upper limit text = table value", length(nA) == 1 && grepl(paste0("upper limit (", num(hi[9], 2), ")"), nA, fixed = TRUE))
check("Table 4: (vii) row identical to primary row", eq(D4[[8]][-1], D4[[1]][-1]))

# ---- Table 5 ----
K <- read.csv("out/table5_R.csv", stringsAsFactors = FALSE); K <- setNames(as.numeric(K$value), K$key)
NP <- read.csv("out/table5_np.csv", stringsAsFactors = FALSE); NP <- setNames(as.numeric(NP$value), NP$key)
cmn <- intersect(names(K), names(NP))
check("Table 5 inputs: R (70) = NumPy (71) on all shared keys (1e-8)", length(cmn) == 584 &&   # 464 + 120 keys for the 500-case threshold (2026-11-24)
      all(abs(K[cmn] - NP[cmn]) <= 1e-8 * pmax(1, abs(K[cmn])) | (is.na(K[cmn]) & is.na(NP[cmn]))))
ps <- as.data.frame(read_excel("out/final/final_results.xlsx", "Q3_robustness")); mo <- as.data.frame(read_excel("out/15_nowcast_monthly_refit.xlsx", "robustness"))
S <- as.data.frame(read_excel("out/final/final_results.xlsx", "Summary")); sv <- function(k) as.numeric(S$value[S$item == k])
ser <- c("dengue", "topic", "day1", "fever"); lab <- c("\"dengue\" (primary series)", "Preprint Topic series", "First retrieval only", "Control series \"Fever\"")
CIs <- rbind(cbind(ps$lo, ps$hi), cbind(mo$lo, mo$hi), c(sv("Q3_loso_pairdiff_M2_minus_M1.lo"), sv("Q3_loso_pairdiff_M2_minus_M1.hi")))
sch <- c(rep("ps", 4), rep("mo", 4), "loso"); sr <- c(ser, ser, "dengue")
mc <- sapply(1:9, function(i) { w <- CIs[i, 2] - CIs[i, 1]; max(abs(CIs[i, 1] - K[[sprintf("%s_%s_bootlo_Rrng", sr[i], sch[i])]]), abs(CIs[i, 2] - K[[sprintf("%s_%s_boothi_Rrng", sr[i], sch[i])]])) / w })
check("Table 5: pipeline 95% CIs within Monte Carlo error of R's own-RNG bootstrap (5% of interval width)", max(mc) <= 0.05, sprintf("max %.1f%% of width", 100 * max(mc)))
check("Table 5: pipeline differences and season counts = R", max(abs(c(ps$pairdiff, mo$pairdiff, sv("Q3_loso_pairdiff_M2_minus_M1.mean")) - sapply(1:9, function(i) K[[sprintf("%s_%s_diff", sr[i], sch[i])]]))) < 1e-6 &&
      all(c(ps$seasons_M2_better, mo$seasons_M2_better, sv("Q3_loso_pairdiff_M2_minus_M1.seasons_M2_better")) == sapply(1:9, function(i) K[[sprintf("%s_%s_better", sr[i], sch[i])]])))
want5 <- lapply(1:9, function(i) { s <- sr[i]; h <- sch[i]; ns <- K[[sprintf("%s_%s_nseasons", s, h)]]; b <- K[[sprintf("%s_%s_better", s, h)]]
  rule <- if (h == "loso") "No rule" else if (b >= 2 && CIs[i, 2] < 0) "Met" else "Not met"
  c(lab[match(s, ser)], as.character(K[[sprintf("%s_%s_ntest", s, h)]]), sapply(c("M0", "M1", "M2", "M3"), function(M) num(K[[sprintf("%s_%s_mae_%s", s, h, M)]], 0)),
    paste0(num(K[[sprintf("%s_%s_diff", s, h)]], 0), " (", civ(CIs[i, 1], CIs[i, 2], 0), ")"), sprintf("%d of %d", as.integer(b), as.integer(ns)), rule) })
C5 <- cells(tabs[2]); D5 <- C5[sapply(C5, length) == 9][-1]
check("Table 5: 9 data rows, every cell = R's own value and formatting", length(D5) == 9 && all(mapply(eq, D5, want5)),
      if (length(D5) == 9) paste(which(!mapply(eq, D5, want5)), collapse = " ") else sprintf("%d rows", length(D5)))
grp <- unlist(C5[sapply(C5, length) == 1])
check("Table 5: group rows in order (per season, every month post hoc, leave one season out)",
      eq(grp, c("Refitted once per season (primary schedule)", "Refitted every month (post hoc) †", "Leave one season out, 2020–2024 (secondary)")))
check("Table 5: primary row decides Question 3 and does not meet the rule", want5[[1]][9] == "Not met" && D5[[1]][9] == "Not met")
nB <- paras[grepl("^ᵇ Mean over test months", paras)]
bt <- sapply(c("ps", "mo"), function(s) { k <- paste0("q3_", s, "_diff_b3", c("lo", "hi")); v <- sapply(k, ct)
  ok <- all(CT$ok[CT$key %in% k] == "True") && all(abs(v - as.numeric(CT$impl3_numpy[match(k, CT$key)])) < 1e-6); if (ok) civ(v[1], v[2], 0) else "MISMATCH" })
check("Table 5 note b: 3-month block intervals = shared-resample values (R = NumPy)", length(nB) == 1 && all(bt != "MISMATCH") &&
      grepl(paste0("were ", bt[["ps"]], " (once per season) and ", bt[["mo"]], " (every month)"), nB, fixed = TRUE), paste(bt, collapse = "; "))
check("Table 5: leave-one-season-out has 59 test months in 5 seasons", K[["dengue_loso_ntest"]] == 59 && K[["dengue_loso_nseasons"]] == 5)
check("Table 5 note a: M0 and M1 errors identical across series within each schedule", all(sapply(c("ps", "mo"), function(h) all(sapply(c("M0", "M1"), function(M)
      length(unique(round(sapply(ser, function(s) K[[sprintf("%s_%s_mae_%s", s, h, M)]]), 8))) == 1)))))
cat(sprintf("RESULT: %d failed check(s)\n", fails)); quit(status = if (fails == 0) 0 else 1)
