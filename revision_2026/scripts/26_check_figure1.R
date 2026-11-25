# Independent check of Figure 1 (R; does not use any Python-derived input file).
# Rebuilds the three plotted series from the original files and compares them with the
# values read back from the drawn lines (out/fig1_plotted_values.csv, written by scripts/25_figure1.py).
# Run from the project root: Rscript scripts/26_check_figure1.R
suppressMessages(library(readxl))
fails <- 0
check <- function(label, ok) {
  cat(sprintf("[%s] %s\n", if (isTRUE(ok)) "PASS" else "FAIL", label))
  if (!isTRUE(ok)) fails <<- fails + 1
}

# 1. Preprint workbook (PREPRINT_XLSX; default data_raw/paper_finalized_data_RP.xlsx), located by header names
xlsx <- Sys.getenv("PREPRINT_XLSX", "data_raw/paper_finalized_data_RP.xlsx")
x <- as.data.frame(read_excel(xlsx, col_names = FALSE, .name_repair = "minimal"))
hr <- which(apply(x, 1, function(r) any(r %in% "Confirmed Cases")))
check("header row found once in workbook", length(hr) == 1)
hdr <- as.character(unlist(x[hr, ]))
x <- x[(hr + 1):nrow(x), c(match("Fever", hdr), match("Confirmed Cases", hdr))]
names(x) <- c("Fever", "cases"); x <- x[rowSums(is.na(x)) < 2, ]
check("workbook has 72 monthly rows", nrow(x) == 72)
cases <- as.numeric(x$cases); fever <- as.numeric(x$Fever)
months <- format(seq(as.Date("2019-01-01"), by = "month", length.out = 72), "%Y-%m")
cases[months < "2019-06"] <- NA

# 2. Three raw Google Trends retrievals: skip the 2 header lines; "<1" -> 0.5; mean of the three
rd <- function(f) {
  d <- read.csv(f, skip = 2, check.names = FALSE, stringsAsFactors = FALSE)
  v <- d[[2]]; v[v == "<1"] <- "0.5"; stopifnot(identical(d[[1]], months)); as.numeric(v)
}
files <- c("data_raw/gt_dengue_term.csv", "data_raw/gt_dengue_term_day2_2026-11-22.csv",
           "data_raw/gt_dengue_term_day3_2026-11-23.csv")
search <- rowMeans(sapply(files, rd))

# 3. Values read back from the figure
p <- read.csv("out/fig1_plotted_values.csv", stringsAsFactors = FALSE)
check("figure has 72 points per panel, Jan 2019 to Dec 2024", identical(p$month, months))
check("panel (a): missing months are exactly Jan-May 2019",
      identical(p$month[is.na(p$cases_plotted)], months[1:5]))
tol <- 1e-9
same <- function(a, b) length(a) == length(b) && identical(is.na(a), is.na(b)) &&
  sum(!is.na(a)) > 0 && max(abs(a - b), na.rm = TRUE) < tol
check("panel (a): 67 plotted case counts equal the workbook",
      same(p$cases_plotted, cases) && sum(!is.na(p$cases_plotted)) == 67)
check("panel (b): 72 search values equal the mean of the 3 raw retrievals",
      same(p$dengue_search_plotted, search))
check("panel (c): 72 Fever values equal the workbook", same(p$fever_search_plotted, fever))

# 4. Facts the caption will state, recomputed here
cat(sprintf("caption facts: case peak %d (%s); search peak %.0f (%s); Fever peak %.0f (%s); cases 2019 Jun-Dec %d\n",
    as.integer(max(cases, na.rm = TRUE)), months[which.max(cases)], max(search), months[which.max(search)],
    max(fever), months[which.max(fever)], as.integer(sum(cases[months >= "2019-06" & months <= "2019-12"]))))
yr <- substr(months, 1, 4)
cat("annual reported cases:", paste(names(tapply(cases, yr, sum, na.rm = TRUE)),
    tapply(cases, yr, sum, na.rm = TRUE), sep = "=", collapse = "  "), "\n")
cat(sprintf("RESULT: %d failed check(s)\n", fails))
quit(status = if (fails == 0) 0 else 1)
