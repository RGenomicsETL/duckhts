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
progress <- Sys.getenv("SOL_D_PROGRESS")
strata <- genome_strata()
worker <- normalizePath("benchmarks/count_mixture_lowpass_genome/worker.R")
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")

run_one <- function(setting, replicate) {
  seed <- genome_declaration$seed + setting$stratum * 1000L + replicate
  output <- file.path(run_root, sprintf("stratum%02d-seed%d.rds", setting$stratum, seed))
  if (file.exists(output)) {
    saved <- readRDS(output)
    stopifnot(
      identical(saved$manifest_sha256, key), saved$seed == seed,
      identical(saved$setting, as.list(setting[c("depth", "share", "relationship", "error")]))
    )
    return(invisible(output))
  }
  available <- as.numeric(strsplit(trimws(tail(system2("df", c("-B1", "--output=avail", "/"),
    stdout = TRUE
  ), 1L)), "[[:space:]]+")[[1L]][[1L]])
  if (available < 8 * 1024^3 + 16 * 1024^2) stop("Disk reserve below 8 GiB")
  log <- sub("\\.rds$", ".log", output)
  command <- c(
    "-v", "Rscript", shQuote(worker), shQuote(extension), shQuote(output),
    setting$depth, setting$share, setting$relationship, setting$error, seed,
    replicate == 1L, genome_declaration$threads, 1L, shQuote(manifest_path)
  )
  status <- system2("/usr/bin/time", command, stdout = log, stderr = log)
  if (status != 0L || !file.exists(output)) stop("Replicate failed: ", output, "; inspect ", log)
  invisible(output)
}

for (i in seq_len(nrow(strata))) {
  setting <- strata[i, ]
  result <- parallel::mclapply(seq_len(setting$replicates), function(replicate) run_one(setting, replicate),
    mc.cores = genome_declaration$workers, mc.preschedule = FALSE
  )
  if (any(vapply(result, inherits, logical(1L), "try-error"))) stop("Stratum failed: ", i)
  line <- sprintf(
    "Completed genome stratum %d/%d: share=%g depth=%g relationship=%s; 200 seeded outputs retained.",
    i, nrow(strata), setting$share, setting$depth, setting$relationship
  )
  cat(line, "\n")
  if (nzchar(progress)) cat(line, "\n", file = progress, append = TRUE)
}
cat("Completed campaign at", run_root, "\n")
