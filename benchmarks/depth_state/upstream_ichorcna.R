# The upstream reference for benchmarks/benchmark_depth_state_model.Rmd: the
# ichorCNA v0.6.1 model run from its own, unmodified R sources on corrected
# log2 ratios, with the HMMcopy C passes it calls. The sources are the pinned
# copy under .sync/ichorCNA-v0.6.1-dl, checked against the tarball SHA-256
# below. This is the function-level entry (getDefaultParameters, HMMsegment,
# correctIntegerCN) that run_ichorCNA calls after its own read-count
# correction, with run_ichorCNA's default arguments and its guard; the CLI
# cannot take corrected bins directly.

upstream_tarball_sha256 <- "ebcaec4b0a4dbf8e1248c83a82e2ea317df0febeb05d49c7629982ca23c17cda"

load_upstream_ichorcna <- function(repo_root = ".") {
  stopifnot(requireNamespace("HMMcopy", quietly = TRUE),
            requireNamespace("GenomicRanges", quietly = TRUE),
            requireNamespace("foreach", quietly = TRUE),
            requireNamespace("digest", quietly = TRUE))
  mirror <- file.path(repo_root, ".sync", "ichorCNA-v0.6.1-dl")
  tarball <- file.path(mirror, "ichorCNA-v0.6.1.tar.gz")
  if (!file.exists(tarball)) stop("pinned ichorCNA tarball is missing: ", tarball)
  if (digest::digest(file = tarball, algo = "sha256") != upstream_tarball_sha256) {
    stop("pinned ichorCNA tarball does not match its SHA-256")
  }
  source_dir <- file.path(mirror, "ichorCNA-0.6.1", "R")
  env <- new.env(parent = globalenv())
  suppressPackageStartupMessages({
    library(HMMcopy)
    library(GenomicRanges)
    library(foreach)
  })
  for (file in c("segmentation.R", "EM.R", "utils.R")) {
    sys.source(file.path(source_dir, file), envir = env)
  }
  env
}

# Fits one sample's bins (data.frame: chrom, start, end, log2_ratio; ordered by
# chrom then start; log2_ratio may be NA) from every start of the grid with
# run_ichorCNA's defaults, and returns the solutions with upstream's guard
# applied, the Viterbi state of every bin for every start, and the corrected
# segments of every start.
upstream_state_model <- function(env, bins, normal_starts = c(0.5, 0.6, 0.7, 0.8, 0.9),
                                 ploidy_starts = c(2, 3), max_cn = 7) {
  ncbi <- sub("^chr", "", bins$chrom)
  chr_train <- as.character(1:22)
  granges <- GenomicRanges::GRanges(
    seqnames = ncbi, ranges = IRanges::IRanges(start = bins$start + 1, end = bins$end),
    copy = bins$log2_ratio)
  tumour_copy <- list(sample = granges)
  valid <- rep(TRUE, nrow(bins))
  chr_ind <- ncbi %in% chr_train
  log_ratio <- data.frame(copy = bins$log2_ratio)
  solutions <- list()
  states <- list()
  segments <- list()
  counter <- 0L
  for (p in ploidy_starts) {
    for (n in normal_starts) {
      if (n == 0.95 && p != 2) next
      param <- env$getDefaultParameters(
        log_ratio[valid & chr_ind, , drop = FALSE], n_0 = n, maxCN = max_cn,
        includeHOMD = FALSE, ct.sc = NULL, normal2IgnoreSC = 1.0, ploidy_0 = floor(p),
        e = 0.9999999, e.subclone = 0.1, e.sameState = 1, strength = 1e7, likModel = "t")
      results <- suppressMessages(env$HMMsegment(
        tumour_copy, valid, dataType = "copy", param = param, chrTrain = chr_train,
        maxiter = 20, estimateNormal = TRUE, estimatePloidy = TRUE, estimatePrecision = TRUE,
        estimateVar = TRUE, estimateSubclone = TRUE, estimateTransition = TRUE,
        estimateInitDist = TRUE, likChangeConvergence = 1e-4, verbose = FALSE))
      iter <- results$results$iter
      n_est <- results$results$n[1, iter]
      phi_est <- results$results$phi[1, iter]
      sp_est <- results$results$sp[1, iter]
      corrected <- env$correctIntegerCN(
        cn = results$cna[[1]], segs = results$results$segs[[1]], purity = 1 - n_est,
        ploidy = phi_est, cellPrev = 1 - sp_est, maxCNtoCorrect.autosomes = max_cn,
        maxCNtoCorrect.X = max_cn, minPurityToCorrect = 0.05, gender = "female",
        chrs = c(1:22, "X"), correctHOMD = FALSE)
      segs <- corrected$segs
      cna <- corrected$cn
      # The guard of run_ichorCNA (R/runIchorCNA.R, the block after HMMsegment).
      segs_train <- segs[segs$chr %in% chr_train, ]
      alt <- which(segs_train$event != "NEUT")
      max_bin_length <- -Inf
      if (length(alt) > 0) {
        longest <- alt[which.max(segs_train$end[alt] - segs_train$start[alt] + 1)]
        region <- GenomicRanges::GRanges(
          seqnames = segs_train$chr[longest],
          ranges = IRanges::IRanges(start = segs_train$start[longest], end = segs_train$end[longest]))
        hits <- GenomicRanges::findOverlaps(query = region, subject = tumour_copy[[1]][valid, ])
        max_bin_length <- length(S4Vectors::subjectHits(hits))
      }
      neutral <- cna[cna$chr %in% chr_train, "event"] == "NEUT"
      alt_fraction <- sum(!neutral, na.rm = TRUE) / length(neutral)
      guard_zero <- (max_bin_length <= 50) && (alt_fraction <= 0.05)
      if (guard_zero) n_est <- 1.0
      counter <- counter + 1L
      solutions[[counter]] <- data.frame(
        start_normal = n, start_ploidy = p, normal_fraction = n_est, ploidy = phi_est,
        log_likelihood = results$results$loglik[iter], iterations = iter - 1L,
        altered_fraction = alt_fraction, guard_zero = guard_zero, stringsAsFactors = FALSE)
      states[[counter]] <- data.frame(
        start_normal = n, start_ploidy = p, chrom = bins$chrom, start = bins$start,
        state = param$jointCNstates[results$viterbiResults$states, 1], stringsAsFactors = FALSE)
      segments[[counter]] <- data.frame(
        start_normal = n, start_ploidy = p, chrom = paste0("chr", segs$chr), start = segs$start - 1,
        end = segs$end, state = segs$copy.number, copy_number = segs$Corrected_Copy_Number,
        median_log2_ratio = segs$median, stringsAsFactors = FALSE)
    }
  }
  solutions <- do.call(rbind, solutions)
  solutions <- solutions[order(-solutions$log_likelihood, seq_len(nrow(solutions))), ]
  solutions$solution <- seq_len(nrow(solutions))
  list(solutions = solutions, states = do.call(rbind, states), segments = do.call(rbind, segments))
}
