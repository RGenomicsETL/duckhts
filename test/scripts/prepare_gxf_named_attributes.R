#!/usr/bin/env Rscript
# Reconstruct BGZF/tabix fixtures from the committed, authored text inputs.
# Run from the repository root with bgzip and tabix (HTSlib 1.19) on PATH.
inputs <- c("gff_named_attributes.gff3", "gtf_named_attributes.gtf")
for (name in inputs) {
  path <- file.path("test", "data", name)
  output <- paste0(path, ".gz")
  status <- system2("bgzip", c("-c", shQuote(path)), stdout = output)
  stopifnot(status == 0L)
  status <- system2("tabix", c("-f", "-p", "gff", shQuote(output)))
  stopifnot(status == 0L)
  files <- c(path, output, paste0(output, ".tbi"))
  stopifnot(all(file.copy(files, "r/Rduckhts/inst/extdata", overwrite = TRUE)))
}
