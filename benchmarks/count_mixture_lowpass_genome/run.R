source("benchmarks/count_mixture_lowpass_genome/manifest.R")
source("benchmarks/count_mixture_lowpass_genome/declaration.R")
args <- commandArgs(TRUE)
extension <- normalizePath(args[[1L]])
root <- args[[2L]]
dir.create(root, recursive = TRUE, showWarnings = FALSE)
manifest_path <- file.path(root, "manifest.rds")
if (!file.exists(manifest_path)) genome_manifest(extension, manifest_path)
manifest <- readRDS(manifest_path)
stopifnot(
  identical(manifest$input_identity, genome_input_identity()),
  identical(manifest$extension_sha256, digest::digest(file = extension, algo = "sha256"))
)
for (path in names(manifest$driver_sha256)) {
  stopifnot(identical(manifest$driver_sha256[[path]], digest::digest(file = path, algo = "sha256")))
}
key <- digest::digest(file = manifest_path, algo = "sha256")
run_root <- file.path(root, key)
dir.create(run_root, showWarnings = FALSE)
worker <- normalizePath("benchmarks/count_mixture_lowpass_genome/batch_worker.R")
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")

run_worker <- function(ids, label) {
  available <- as.numeric(trimws(tail(system2("df", c("-B1", "--output=avail", "/"), stdout = TRUE), 1L)))
  if (available < 8 * 1024^3 + 64 * 1024^2) stop("Disk reserve below 8 GiB")
  log <- file.path(root, paste0(label, ".log"))
  status <- system2("/usr/bin/time", c(
    "-v", "Rscript", shQuote(worker), shQuote(extension), shQuote(manifest_path),
    shQuote(run_root), paste(ids, collapse = ",")
  ), stdout = log, stderr = log)
  if (status != 0L) stop("Worker failed: inspect ", log)
  invisible(log)
}

# Retain the first complete stratum as a throughput pilot before concurrency.
started <- proc.time()[["elapsed"]]
run_worker(4L, "pilot-stratum04")
pilot_seconds <- proc.time()[["elapsed"]] - started
expected <- 18000 * pilot_seconds / 200 / genome_declaration$workers
line <- sprintf("200-replicate depth .1 throughput pilot: %.1f s, %.3f s/replicate. Depth-weighted probe predicts about 1-1.5 hours on 14 workers; pilot-only lower-depth projection %.2f hours.",
  pilot_seconds, pilot_seconds / 200, expected / 3600)
cat(line, "\n")
progress <- Sys.getenv("SOL_D_PROGRESS")
if (nzchar(progress)) cat(line, "\n", file = progress, append = TRUE)

ids <- setdiff(genome_strata()$stratum, 4L)
groups <- split(ids, rep(seq_len(genome_declaration$workers), length.out = length(ids)))
started <- proc.time()[["elapsed"]]
result <- parallel::mclapply(seq_along(groups), function(i) {
  run_worker(groups[[i]], sprintf("worker%02d", i))
}, mc.cores = genome_declaration$workers, mc.preschedule = TRUE)
if (any(vapply(result, inherits, logical(1L), "try-error"))) stop("Campaign worker failed")
line <- sprintf("Completed campaign at %s; concurrent phase %.1f seconds; 90 strata, 18,000 replicates.",
  run_root, proc.time()[["elapsed"]] - started)
cat(line, "\n")
if (nzchar(progress)) cat(line, "\n", file = progress, append = TRUE)
