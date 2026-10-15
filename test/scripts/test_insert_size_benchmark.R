#!/usr/bin/env Rscript
# The sparse comparator must retain every actual row, independent of key order.
root <- normalizePath(commandArgs(trailingOnly = TRUE)[[1L]], mustWork = TRUE)
source(file.path(root, "benchmarks/insert_size_driver.R"))
raw <- data.frame(insert_size = c("1", "2"), fr_count = c("2", "0"),
                  rf_count = c("0", "3"), tandem_count = c("0", "0"))
counts <- data.frame(pair_orientation = c("RF", "FR"), insert_size = c("2", "1"),
                     read_pairs = c("3", "2"))
insert_size_check_raw(counts, raw)
changed <- counts
changed$read_pairs[[1L]] <- "4"
zero <- data.frame(pair_orientation = "TANDEM", insert_size = "2", read_pairs = "0")
extra <- data.frame(pair_orientation = "FR", insert_size = "3", read_pairs = "1")
for (damaged in list(rbind(counts, counts[1L, ]), counts[-1L, ], changed,
                     rbind(counts, extra), rbind(counts, zero))) {
  failed <- tryCatch({
    insert_size_check_raw(damaged, raw)
    FALSE
  }, error = function(error) TRUE)
  stopifnot(failed)
}
message("PASS: keyed counts; duplicates, missing/extra rows, changed counts and zero rows rejected")
