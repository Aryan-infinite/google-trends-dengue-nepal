# Independent check of Supplementary Figure S2 (R) and three-program comparison.
# R recomputes every plotted value and every event position from the raw files and compares:
#   (1) R vs Python (out/figS2_values.csv) vs JavaScript (out/figS2_values_js.csv, out/figS2_events_js.csv);
#   (2) what was drawn (out/figS2_plotted_values.csv, out/figS2_events_drawn.csv, out/figS2_text.csv);
#   (3) the check E timing metrics logged earlier (out/23_announcement_timing.json, out/24_announcement_timing.json);
#   (4) the delivered TIFF's mode, resolution and size (ImageMagick identify).
# Event categories here come from a table keyed by DATE, written from the statements in the events file
# (Python keys on event_type through a lookup table; JavaScript uses keyword rules).
suppressMessages(library(readxl))
fails <- 0
check <- function(l, ok, d = "") { cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", l, if (nzchar(d)) paste0("  (", d, ")") else "")); if (!isTRUE(ok)) fails <<- fails + 1 }
same <- function(a, b, tol = 1e-12) length(a) == length(b) && length(a) > 0 && identical(is.na(a), is.na(b)) && all(abs(a[!is.na(a)] - b[!is.na(b)]) <= tol)
mk <- format(seq(as.Date("2019-01-01"), by = "month", length.out = 72), "%Y-%m")
x <- as.data.frame(read_excel(Sys.getenv("PREPRINT_XLSX", "data_raw/paper_finalized_data_RP.xlsx"), col_names = FALSE, .name_repair = "minimal"))
hr <- which(apply(x, 1, function(r) any(r %in% "Confirmed Cases"))); hd <- as.character(unlist(x[hr, ]))
cases <- as.numeric(x[(hr + 1):nrow(x), match("Confirmed Cases", hd)]); cases <- cases[!is.na(cases)]
gt <- function(f) { d <- read.csv(f, skip = 2, stringsAsFactors = FALSE); stopifnot(identical(d[[1]], mk)); v <- d[[2]]; v[v == "<1"] <- "0.5"; as.numeric(v) }
s <- (gt("data_raw/gt_dengue_term.csv") + gt("data_raw/gt_dengue_term_day2_2026-11-22.csv") + gt("data_raw/gt_dengue_term_day3_2026-11-23.csv")) / 3
check("72 months of cases and search", length(cases) == 72 && length(s) == 72)
YR <- rep(2019:2024, each = 12); MO <- rep(1:12, 6)
cc <- cases; cc[YR == 2019 & MO <= 5] <- NA
crel <- ave(cc, YR, FUN = function(v) v / max(v, na.rm = TRUE)); srel <- ave(s, YR, FUN = function(v) v / max(v))
Rv <- data.frame(year = YR, month = MO, cases_rel = crel, search_rel = srel)

P <- read.csv("out/figS2_values.csv"); J <- read.csv("out/figS2_values_js.csv")
ok_key <- function(d) identical(as.integer(d$year), YR) && identical(as.integer(d$month), MO)
check("(1) R = Python, 144 relative values (1e-12; same missing months)", ok_key(P) && same(P$cases_rel, Rv$cases_rel) && same(P$search_rel, Rv$search_rel))
check("(1) R = JavaScript, 144 relative values", ok_key(J) && same(J$cases_rel, Rv$cases_rel) && same(J$search_rel, Rv$search_rel))
check("(1) Python raw monthly values = raw files", same(P$cases[!is.na(cc)], cc[!is.na(cc)], 0) && same(P$search, s))

D <- read.csv("out/figS2_plotted_values.csv")
dv <- function(nm) { d <- D[D$series == nm, ]; d[order(d$year, d$x), ] }
dc <- dv("cases"); ds <- dv("search")
check("(2) drawn x = month centres 0.5..11.5 in every panel", nrow(dc) == 72 && nrow(ds) == 72 && same(dc$x, MO - 0.5) && same(ds$x, MO - 0.5) && identical(as.integer(dc$year), YR))
check("(2) drawn case line = R values (Jan-May 2019 not drawn)", same(dc$y, Rv$cases_rel))
check("(2) drawn search line = R values", same(ds$y, Rv$search_rel))

CATR <- c("2019-06-05" = "alert", "2019-09-21" = "count", "2020-02-26" = "alert", "2022-05-06" = "count", "2022-08-05" = "count",
          "2022-08-31" = "report", "2022-10-10" = "report", "2023-05-10" = "alert", "2023-10-11" = "count", "2024-05-14" = "alert")
E <- read.csv("governance/announcement_events.csv", stringsAsFactors = FALSE, encoding = "UTF-8")
L6 <- read.csv("out/06_2024_reporting_lag.csv", stringsAsFactors = FALSE)
check("events file: every dated single event has an R category; the rest are the two known non-events",
      setequal(E$date[grepl("^20[0-9]{2}-[0-9]{2}-[0-9]{2}$", E$date)], c(names(CATR), "2021-11-03")) && sum(!grepl("^20[0-9]{2}-[0-9]{2}-[0-9]{2}$", E$date)) == 1)
ev <- rbind(data.frame(date = names(CATR), category = unname(CATR)), data.frame(date = L6$as_of[L6$as_of != "2024-12-31"], category = "count"))
dd <- as.Date(ev$date); y <- as.integer(format(dd, "%Y")); m <- as.integer(format(dd, "%m")); d <- as.integer(format(dd, "%d"))
nd <- as.integer(format(as.Date(sprintf("%d-%02d-01", y + (m == 12), m %% 12 + 1)) - 1, "%d"))
ev$year <- y; ev$x <- (m - 1) + (d - 0.5) / nd; ev <- ev[order(ev$date), ]
ED <- read.csv("out/figS2_events_drawn.csv", stringsAsFactors = FALSE); ED <- ED[order(ED$date), ]
EJ <- read.csv("out/figS2_events_js.csv", stringsAsFactors = FALSE); EJ <- EJ[order(EJ$date), ]
check("(2) 15 drawn events = R events (date, year, category, x)", nrow(ED) == 15 && identical(ED$date, ev$date) && identical(ED$category, ev$category) &&
      identical(as.integer(ED$year), ev$year) && same(ED$x, ev$x))
check("(2) marker shape follows category (v alert, o count, s report)", identical(ED$marker, unname(c(alert = "v", count = "o", report = "s")[ED$category])))
check("(1) JavaScript events = R events", nrow(EJ) == 15 && identical(EJ$date, ev$date) && identical(EJ$category, ev$category) && same(EJ$x, ev$x))
check("(2) no event drawn in 2021; all event markers in the strip (y = 1.15)", !any(ED$year == 2021) && all(ED$y == 1.15))

half <- function(v) which(v >= max(v, na.rm = TRUE) / 2)[1]
mine <- unlist(lapply(2019:2024, function(yy) { a <- cc[YR == yy]; b <- s[YR == yy]
  setNames(c(half(b), which.max(b), half(a), which.max(a)), sprintf("E_%d_%s", yy, c("search_half_month", "search_peak_month", "case_half_month", "case_peak_month"))) }))
rd <- function(f) { L <- readLines(f, warn = FALSE); L <- L[grepl('"E_', L)]; setNames(as.integer(sub('.*: *([0-9]+).*', "\\1", L)), sub('^ *"([^"]+)".*', "\\1", L)) }
j23 <- rd("out/23_announcement_timing.json"); j24 <- rd("out/24_announcement_timing.json")
check("(3) half-maximum and peak months from the drawn data = check E (scripts/23 and /24), 24 values",
      length(j23) == 24 && identical(mine[names(j23)], j23) && identical(j24[names(j23)], j23))

T <- read.csv("out/figS2_text.csv", stringsAsFactors = FALSE, encoding = "UTF-8")
want <- c(sprintf("(%s) %d", letters[1:6], 2019:2024), "Cases not recorded Jan–May", "Search index ≤ 3 all year", "Search index ≤ 2 all year",
          "No national alert found", "Cases", "Search interest", "National alert or EDCD statement", "Count announced", "Daily reports began / WHO report")
hex <- function(v) vapply(v, function(z) paste(charToRaw(z), collapse = ""), "", USE.NAMES = FALSE)   # byte comparison: locale-proof
check("(2) all drawn text exactly as intended (UTF-8 bytes)", nrow(T) == length(want) && setequal(hex(T$text), hex(want)))
check("(2) notes sit on the right panels", identical(T$year[T$kind == "note"], c(2019L, 2020L, 2021L)) && identical(T$year[T$kind == "strip"], 2021L))

id <- system2("identify", c("-format", shQuote("%[colorspace] %x %y %w %h %[compression]"), "manuscript/figures/FigS2.tif"), stdout = TRUE)
p <- strsplit(id, " ")[[1]]
check("(4) TIFF: CMYK, 1200 dpi, 170 x 150 mm, LZW", p[1] == "CMYK" && as.numeric(p[2]) == 1200 && as.numeric(p[3]) == 1200 &&
      abs(as.numeric(p[4]) / 1200 * 25.4 - 170) < 0.05 && abs(as.numeric(p[5]) / 1200 * 25.4 - 150) < 0.05 && p[6] == "LZW", id)
cat(sprintf("RESULT: %d failed check(s)\n", fails)); quit(status = fails > 0)
