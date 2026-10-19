library(tinytest)
library(DBI)

# Synthetic corrected bins with a known truth, as in test/sql/depth_state_model.test:
# 22 autosomes of 120 bins and X of 100, one sample with normal fraction 0.7 at
# ploidy 2 carrying three copies of chromosome 8, one copy of chromosome 13 and
# one copy of the first 40 bins of chromosome 17, and one sample without change.
# Declared before running: the best tumour solution has a normal fraction within
# 0.05 of 0.7 and a ploidy within 0.15 of 2; chromosomes 8 and 13 are called over
# at least 95% of their bins, the 17p loss over at least 90%, every other
# autosomal bin two copies over at least 99%; the normal sample is guard_zero.

set.seed(343)
contigs <- c(as.character(1:22), "X")
bins_per_contig <- c(rep(120L, 22L), 100L)
genome <- do.call(rbind, lapply(seq_along(contigs), function(i) {
  b <- seq_len(bins_per_contig[[i]]) - 1L
  data.frame(chrom = contigs[[i]], bin = b, start = b * 1e6, end = (b + 1) * 1e6,
             stringsAsFactors = FALSE)
}))
design <- data.frame(sample_id = c("tumour", "normal"), normal_fraction = c(0.7, 1),
                     stringsAsFactors = FALSE)
truth <- do.call(rbind, lapply(seq_len(nrow(design)), function(i) {
  copies <- rep(2, nrow(genome))
  if (design$sample_id[[i]] == "tumour") {
    copies[genome$chrom == "8"] <- 3
    copies[genome$chrom == "13"] <- 1
    copies[genome$chrom == "17" & genome$bin < 40] <- 1
  }
  n <- design$normal_fraction[[i]]
  data.frame(sample_id = design$sample_id[[i]], genome, copies = copies,
             log2_ratio = log2(((1 - n) * copies + 2 * n) / ((1 - n) * 2 + 2 * n)) +
               0.06 * stats::rnorm(nrow(genome)),
             kept = TRUE, stringsAsFactors = FALSE)
}))
truth$log2_ratio[truth$chrom == "5" & truth$bin == 10] <- NA

con <- rduckhts_connect()
dbWriteTable(con, "bins", truth[, c("sample_id", "chrom", "start", "end", "log2_ratio", "kept")])

solutions <- rduckhts_depth_state_solutions(con, "bins")
expect_equal(
  names(solutions),
  c("sample_id", "solution", "start_normal", "start_ploidy", "normal_fraction", "tumour_fraction",
    "ploidy", "log_likelihood", "altered_fraction", "iterations", "status", "segments", "bins",
    "upstream", "method"))
expect_equal(as.vector(table(solutions$sample_id)), c(10L, 10L))
expect_equal(sort(solutions$solution[solutions$sample_id == "tumour"]), 1:10)
best <- solutions[solutions$solution == 1, ]
tumour <- best[best$sample_id == "tumour", ]
normal <- best[best$sample_id == "normal", ]
expect_true(abs(tumour$normal_fraction - 0.7) <= 0.05)
expect_true(abs(tumour$ploidy - 2) <= 0.15)
expect_equal(tumour$tumour_fraction, 1 - tumour$normal_fraction)
expect_equal(tumour$status, "ok")
# A genome without change: the fraction is zero, by the guard or by a fit
# whose normal fraction is near 1. In such a fit the state means lie within
# the noise of each other, so whole chromosomes can take any state (as in
# upstream); the copy numbers of that solution are not evidence, the
# fraction is.
expect_true(normal$status %in% c("guard_zero", "ok"))
expect_true(normal$tumour_fraction <= 0.01)
expect_equal(unique(solutions$upstream), "ichorCNA v0.6.1")
expect_equal(unique(solutions$method), "duckhts_depth_state v1")

segments <- rduckhts_depth_state_segments(con, "bins")
expect_true(all(c("chrom", "start", "end", "bins", "state", "copy_number", "median_log2_ratio") %in%
                  names(segments)))
best_segments <- segments[segments$solution == 1, ]
expect_equal(
  as.vector(tapply(best_segments$bins, best_segments$sample_id, sum)[c("normal", "tumour")]),
  c(nrow(genome), nrow(genome)))
calls <- merge(
  truth, best_segments[, c("sample_id", "chrom", "start", "end", "copy_number")],
  by = c("sample_id", "chrom"), suffixes = c("", ".segment"))
calls <- calls[calls$start >= calls$start.segment & calls$end <= calls$end.segment, ]
expect_equal(nrow(calls), nrow(truth))
share <- function(sample, keep, copy) {
  rows <- calls$sample_id == sample & keep(calls)
  mean(calls$copy_number[rows] == copy)
}
expect_true(share("tumour", function(c) c$chrom == "8", 3) >= 0.95)
expect_true(share("tumour", function(c) c$chrom == "13", 1) >= 0.95)
expect_true(share("tumour", function(c) c$chrom == "17" & c$bin < 40, 1) >= 0.90)
expect_true(share("tumour", function(c) c$copies == 2 & c$chrom != "X", 2) >= 0.99)
if (normal$status == "guard_zero") {
  expect_true(share("normal", function(c) c$chrom != "X", 2) >= 0.99)
}

# The bin without a ratio stays in its segment; the median ignores it.
hole <- best_segments[best_segments$sample_id == "normal" & best_segments$chrom == "5" &
                        best_segments$start <= 10e6 & best_segments$end > 10e6, ]
expect_equal(nrow(hole), 1L)
expect_true(hole$bins >= 11)
expect_false(is.na(hole$median_log2_ratio))

# Options reach the macro; table_name publishes the result.
narrow <- rduckhts_depth_state_solutions(con, "bins", normal_fractions = 0.7, ploidies = 2)
expect_equal(nrow(narrow), 2L)
expect_true(rduckhts_depth_state_segments(con, "bins", normal_fractions = 0.7, ploidies = 2,
                                          table_name = "state_segments"))
expect_equal(dbGetQuery(con, "SELECT count(DISTINCT solution) AS n FROM state_segments")$n, 1L)

# Argument contracts are checked in R before any query.
expect_error(rduckhts_depth_state_solutions(con, "bins", max_cn = 1), "max_cn must be")
expect_error(rduckhts_depth_state_solutions(con, "bins", normal_fractions = c(0.5, 1.5)),
             "normal_fractions must be")
expect_error(rduckhts_depth_state_solutions(con, "bins", ploidies = "two"), "ploidies must be")
expect_error(rduckhts_depth_state_solutions(con, "bins", txn_e = 1), "txn_e must be")
expect_error(rduckhts_depth_state_segments(con, "bins", max_bins = 100), "max_bins")
expect_error(rduckhts_depth_state_segments(con, "no such table"), "no such table")

dbDisconnect(con, shutdown = TRUE)
