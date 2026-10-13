source("benchmarks/count_mixture_lowpass/common.R")
args <- commandArgs(TRUE)
extension <- normalizePath(args[[1L]])
scratch <- normalizePath(args[[2L]], mustWork = FALSE)
dir.create(scratch, recursive = TRUE, showWarnings = FALSE)

jobs <- rbind(data.frame(kind = "synthetic", id = lowpass_strata()$id),
              data.frame(kind = "thinning", id = lowpass_thinning_strata()$id))
run_one <- function(i) {
  kind <- jobs$kind[[i]]
  id <- jobs$id[[i]]
  output <- file.path(scratch, sprintf("%s-%d.rds", kind, id))
  log <- paste0(output, ".log")
  status <- system2("Rscript", shQuote(c("benchmarks/count_mixture_lowpass/worker.R", kind, id, output, extension)),
                    stdout = log, stderr = log)
  if (status != 0L) stop("Worker failed: ", log)
  output
}
# Recovery has no elapsed-time criterion. Three workers share the assigned
# six cores; scaling below is strictly sequential in fresh processes.
outputs <- parallel::mclapply(seq_len(nrow(jobs)), run_one, mc.cores = 3L, mc.preschedule = FALSE)
if (any(vapply(outputs, inherits, logical(1L), "try-error"))) stop("Campaign workers failed")

for (axis in c("sites", "groups")) {
  for (factor in c(1L, 2L, 4L)) {
    for (repeat_id in seq_len(lowpass_declaration$scaling_repeats)) {
      size <- lowpass_declaration$scaling_sites[[1L]] * if (axis == "sites") factor else 1L
      groups <- if (axis == "groups") factor else 1L
      output <- file.path(scratch, sprintf("scale-%s-%d-%d.rds", axis, factor, repeat_id))
      rss <- paste0(output, ".rss")
      log <- paste0(output, ".log")
      status <- system2("/usr/bin/time", c("-f", shQuote("%M"), "-o", shQuote(rss), "Rscript",
        shQuote(c("benchmarks/count_mixture_lowpass/memory_worker.R", "unordered", size,
                 lowpass_declaration$threads, output, extension, groups))), stdout = log, stderr = log)
      if (status != 0L) stop("Scaling worker failed: ", log)
    }
  }
}
metadata <- list(source_revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  extension_sha256 = system2("sha256sum", extension, stdout = TRUE),
  dirty = system2("git", c("status", "--porcelain"), stdout = TRUE),
  R = R.version.string, duckdb = as.character(packageVersion("duckdb")),
  cpu = Sys.info(), cores = "0-5", fresh_processes = TRUE,
  declaration_commit = "6e23123c", cache = "warm registered inputs; each DuckDB process fresh")
saveRDS(metadata, file.path(scratch, "metadata.rds"))
