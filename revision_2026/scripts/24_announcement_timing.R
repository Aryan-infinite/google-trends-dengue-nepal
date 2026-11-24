# Check E timing metrics, reimplemented in R from the addendum (not from scripts/23).
suppressMessages(library(readxl))   # run from the project root (as run_all.sh does)
gt <- function(f) { L <- readLines(paste0("data_raw/", f), encoding = "UTF-8"); L <- L[grepl("^20[0-9]{2}-[0-9]{2},", L)]; v <- sub("^[^,]*,", "", L); v[trimws(v) == "<1"] <- "0.5"; as.numeric(v) }
x <- (gt("gt_dengue_term.csv") + gt("gt_dengue_term_day2_2026-11-22.csv") + gt("gt_dengue_term_day3_2026-11-23.csv")) / 3
xl <- Sys.getenv("PREPRINT_XLSX", "data_raw/paper_finalized_data_RP.xlsx")          # preprint table; columns found by header name
pp <- as.data.frame(read_excel(xl, col_names = FALSE, .name_repair = "minimal"))
hr <- which(apply(pp, 1, function(r) any(r %in% "Confirmed Cases"))); hd <- as.character(unlist(pp[hr, ])); pp <- pp[(hr + 1):nrow(pp), ]
mo_ <- suppressWarnings(as.numeric(pp[[match("Months", hd)]])); cs_ <- suppressWarnings(as.numeric(pp[[match("Confirmed Cases", hd)]]))
cs_all <- cs_[!is.na(mo_)][order(mo_[!is.na(mo_)])]; stopifnot(length(cs_all) == 72, !anyNA(cs_all))
out <- c()
for (i in 0:5) { y <- 2019 + i; xs <- x[(i * 12 + 1):(i * 12 + 12)]; cs <- cs_all[(i * 12 + 1):(i * 12 + 12)]; if (y == 2019) cs[1:5] <- NA
  v <- c(search_half_month = which(xs >= max(xs) / 2)[1], search_peak_month = which.max(xs), case_half_month = which(cs >= max(cs, na.rm = TRUE) / 2)[1], case_peak_month = which.max(cs))
  for (k in names(v)) out[sprintf("E_%d_%s", y, k)] <- v[[k]] }
writeLines(paste0("{\n", paste(sprintf('"%s": %d', names(out), as.integer(out)), collapse = ",\n"), "\n}"), "out/24_announcement_timing.json")
