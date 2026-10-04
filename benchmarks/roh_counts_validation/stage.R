# Stage the registered genotype BCF of the read-count validation into the cache.
# Run from the repository root: Rscript benchmarks/roh_counts_validation/stage.R
source("benchmarks/roh_counts_validation/stage_genotypes.R")
staged <- stage_roh_validation_genotypes_from_registry(
  "roh_counts_validation_chr20_genotypes",
  bcftools = Sys.getenv("BCFTOOLS", "/usr/local/bin/bcftools"))
cat("staged", staged$output, "records", staged$records, "samples", staged$samples,
    "record SHA-256", staged$records_sha256, "\n")
