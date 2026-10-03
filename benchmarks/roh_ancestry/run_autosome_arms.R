# Run every genome-wide ROH arm, each in a fresh R process under /usr/bin/time.
# Usage, from the repository root: Rscript benchmarks/roh_ancestry/run_autosome_arms.R [REPETITIONS]
# The accuracy arms are deterministic, so one repetition (the default) answers the
# criterion; summarise_genome.R requires any further repetitions to match the first.
# A completed repetition (output, timing receipt and .ok marker) is reused only
# when its marker records the same extension and proportions SHA-256 as now.
args <- commandArgs(trailingOnly = TRUE)
repetitions <- if (length(args)) as.integer(args[[1L]]) else 1L
if (length(repetitions) != 1L || is.na(repetitions) || repetitions < 1L || repetitions > 3L) {
  stop("REPETITIONS must be 1, 2 or 3", call. = FALSE)
}
cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-ancestry")
script <- normalizePath("benchmarks/roh_ancestry/run_autosome_arm.R", mustWork = TRUE)
arms <- c("pooled", "single_AFR", "single_AMR", "single_EUR", "single_EAS", "ancestry")
sha256 <- function(path) strsplit(system2("sha256sum", shQuote(path), stdout = TRUE), " ")[[1L]][[1L]]
identity <- c(
  extension_sha256 = sha256(normalizePath("build/release/duckhts.duckdb_extension", mustWork = TRUE)),
  proportions_sha256 = sha256(file.path(cache, "autosome.q.csv")))
identity_lines <- paste0(names(identity), "=", identity)

for (arm in arms) {
  for (repetition in seq_len(repetitions)) {
    stem <- file.path(cache, sprintf("autosome.arm_%s_%d", arm, repetition))
    output <- paste0(stem, ".csv")
    timing <- paste0(stem, ".time.txt")
    success <- paste0(stem, ".ok")
    if (file.exists(success) && file.size(output) > 0 && file.size(timing) > 0 &&
        identical(readLines(success), identity_lines)) {
      message("reuse ", arm, " repetition ", repetition)
      next
    }
    unlink(c(output, timing, success))
    message("run ", arm, " repetition ", repetition)
    status <- system2("/usr/bin/time",
      c("-f", shQuote("%e\t%M"), "-o", shQuote(timing), "Rscript", shQuote(script),
        arm, repetition),
      stdout = paste0(stem, ".log"), stderr = paste0(stem, ".log"))
    if (!identical(status, 0L)) {
      stop("arm ", arm, " repetition ", repetition, " failed with status ", status,
           call. = FALSE)
    }
    writeLines(identity_lines, success)
  }
}
