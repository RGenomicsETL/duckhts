#!/usr/bin/env Rscript
# The receipt must see every output field, not just the uploaded variant ID.
args <- commandArgs(TRUE)
script <- if (length(args)) args[[1L]] else "benchmarks/benchmark_duckvep_fastvep_receipt.R"
local({
  work <- tempfile("duckhts-fastvep-receipt-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  fields <- c("Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type",
    "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
    "Codons", "Existing_variation", "IMPACT", "DISTANCE", "STRAND", "FLAGS")
  data <- as.data.frame(setNames(lapply(seq_along(fields), function(i) paste0("v", i, "-", 1:3)), fields))
  input <- file.path(work, "output.tab")
  output <- file.path(work, "receipt.csv")
  fingerprint <- function(data) {
    write.table(data, input, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    status <- system2("Rscript", c(shQuote(script), "--input", shQuote(input),
      "--tool", "duckvep", "--output-contract", "matched_tab_17_fields",
      "--output", shQuote(output)), stdout = file.path(work, "receipt.log"), stderr = "")
    stopifnot(status == 0L)
    receipt <- read.csv(output, colClasses = "character")
    receipt[c("row_count", "xor_hash", "low32_sum", "high32_sum")]
  }
  original <- fingerprint(data)
  stopifnot(identical(original, fingerprint(data[3:1, ])))
  for (field in fields) {
    changed <- data
    changed[[field]][2L] <- "changed"
    if (identical(original, fingerprint(changed))) stop("receipt ignores field: ", field)
  }
  stopifnot(!identical(original, fingerprint(rbind(data, data[1L, ]))))
  missing <- data
  missing$FLAGS[2L] <- NA_character_
  stopifnot(!identical(original, fingerprint(missing)))
  cat("FastVEP output receipt: all 17 fields, missing values, duplicate multiplicity and row-order invariance: OK\n")
})
