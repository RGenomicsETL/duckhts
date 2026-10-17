source("benchmarks/count_mixture_lowpass_genome/common.R")
args <- commandArgs(TRUE)
output <- args[[1L]]
checked <- parallel::mclapply(genome_settings$chromosomes, function(chrom) {
  id <- paste0("roh_ancestry_chr", chrom, "_children_bcf")
  path <- duckhts_bench_artifact_path(id)
  duckhts_bench_validate_identity(id, path)
  temporary <- tempfile()
  on.exit(unlink(temporary))
  status <- system2(RBCFTools::bcftools_path(), c(
    "query", "-i", shQuote('GT~"/"'),
    "-f", shQuote("%CHROM\t%POS\n"), shQuote(path)
  ), stdout = temporary)
  if (status != 0L) stop("Phase check failed: ", id)
  counterexamples <- readLines(temporary, n = 10L)
  if (length(counterexamples)) stop("Unphased children GTs: ", paste(counterexamples, collapse = "; "))
  data.frame(chrom = chrom, artifact = id, unphased_records = 0L)
}, mc.cores = 2L)
if (any(vapply(checked, inherits, logical(1L), "try-error"))) stop("Children phase audit failed")
write.csv(do.call(rbind, checked), output, row.names = FALSE)
