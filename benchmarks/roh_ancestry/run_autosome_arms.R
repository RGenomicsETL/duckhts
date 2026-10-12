# Run every genome-wide ROH arm, each in a fresh R process under /usr/bin/time.
# Usage, from the repository root: Rscript benchmarks/roh_ancestry/run_autosome_arms.R [REPETITIONS]
# The accuracy arms are deterministic, so one repetition (the default) answers the
# criterion; summarise_genome.R requires any further repetitions to match the first.
# A completed repetition (output, timing receipt and .ok marker) is reused only
# when its marker records the same SHA-256 for every input the arms read (the
# extension, the arm script, the shared helpers, the pedigree, the children
# list, the proportions, and each autosome's children BCF, AF sites and long
# reference), the same chromosome set and duckdb version, and the output's
# current SHA-256. The arm processes inherit ROH_GENOME_CHROMOSOMES.
source("benchmarks/roh_ancestry/genome_common.R")
args <- commandArgs(trailingOnly = TRUE)
repetitions <- if (length(args)) as.integer(args[[1L]]) else 1L
if (length(repetitions) != 1L || is.na(repetitions) || repetitions < 1L || repetitions > 3L) {
  stop("REPETITIONS must be 1, 2 or 3", call. = FALSE)
}
chromosomes <- genome_chromosomes()
script <- normalizePath("benchmarks/roh_ancestry/run_autosome_arm.R", mustWork = TRUE)
arms <- c("pooled", "single_AFR", "single_AMR", "single_EUR", "single_EAS", "ancestry")
per_chromosome <- unlist(lapply(chromosomes, function(chromosome) {
  stats::setNames(c(genome_children_bcf(chromosome), genome_af_sites(chromosome),
                    genome_reference_long(chromosome)),
                  paste0(c("bcf_", "af_sites_", "reference_"), chromosome))
}))
identity <- c(genome_set_line(chromosomes), genome_identity(
  c(extension = normalizePath("build/release/duckhts.duckdb_extension", mustWork = TRUE),
    arm_script = script,
    helpers = normalizePath("benchmarks/roh_ancestry/genome_common.R", mustWork = TRUE),
    pedigree = duckhtsbench::duckhts_bench_artifact_path("roh_ancestry_pedigree"),
    children = genome_children_list,
    proportions = file.path(genome_cache, "autosome.q.csv"),
    per_chromosome),
  c(duckdb_r = as.character(utils::packageVersion("duckdb")))))
genome_check_output(file.path(genome_cache, "autosome.q.csv"), chromosomes)

for (arm in arms) {
  for (repetition in seq_len(repetitions)) {
    stem <- file.path(genome_cache, sprintf("autosome.arm_%s_%d", arm, repetition))
    output <- paste0(stem, ".csv")
    timing <- paste0(stem, ".time.txt")
    success <- paste0(stem, ".ok")
    if (genome_current(output, identity, success) && isTRUE(file.size(timing) > 0)) {
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
    genome_record(output, identity, success)
  }
}
