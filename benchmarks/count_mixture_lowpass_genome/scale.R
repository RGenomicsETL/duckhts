source("benchmarks/count_mixture_lowpass_genome/common.R")
source("benchmarks/count_mixture_lowpass_genome/collect.R")
args <- commandArgs(TRUE)
extension <- normalizePath(args[[1L]])
manifest_path <- normalizePath(args[[2L]])
root <- args[[3L]]
destination <- args[[4L]]
dir.create(root, recursive = TRUE, showWarnings = FALSE)
manifest <- readRDS(manifest_path)
stopifnot(identical(manifest$extension_sha256, digest::digest(file = extension, algo = "sha256")))
worker <- normalizePath("benchmarks/count_mixture_lowpass_genome/worker.R")
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
rows <- list()
for (depth in genome_declaration$primary_depths) {
  for (stride in genome_declaration$scaling_strides) {
    for (repetition in 1:3) {
      output <- file.path(root, sprintf("depth%g-stride%d-rep%d.rds", depth, stride, repetition))
      log <- sub("\\.rds$", ".log", output)
      status <- system2("/usr/bin/time", c(
        "-v", "Rscript", shQuote(worker), shQuote(extension),
        shQuote(output), depth, .1, "same", .0001, genome_declaration$seed, "FALSE", 1L,
        stride, shQuote(manifest_path)
      ), stdout = log, stderr = log)
      if (status != 0L) stop("Scaling worker failed: ", log)
      saved <- readRDS(output)
      measured <- genome_log_metrics(log)
      sites <- max(saved$fit$sites)
      ceiling <- genome_declaration$memory_overhead_bytes + genome_declaration$live_memory_multiplier *
        (genome_declaration$sql_live_bytes_per_site * sites + saved$live_bytes +
          genome_declaration$max_state_bytes)
      rows[[length(rows) + 1L]] <- data.frame(
        depth = depth, stride = stride,
        repetition = repetition, sites = sites, output_rows = nrow(saved$fit),
        wall_seconds = measured[["wall_seconds"]], simulation_seconds = saved$simulation_seconds,
        fit_seconds = saved$fit_seconds, rss_bytes = measured[["rss_bytes"]],
        peak_buffer_bytes = saved$peak_buffer_bytes, spill_bytes = saved$spill_bytes,
        live_R_bytes = saved$live_bytes, memory_ceiling_bytes = ceiling,
        status = paste(saved$fit$status, collapse = ";"), threads = saved$threads,
        manifest_sha256 = saved$manifest_sha256
      )
    }
  }
}
write.csv(do.call(rbind, rows), destination, row.names = FALSE)
