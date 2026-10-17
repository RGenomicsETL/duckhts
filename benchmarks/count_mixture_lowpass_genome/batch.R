source("benchmarks/count_mixture_lowpass_genome/common.R")
source("benchmarks/count_mixture_lowpass_genome/declaration.R")

genome_load_panels <- function() {
  panels <- lapply(genome_settings$chromosomes, function(chrom) {
    panel <- readRDS(duckhts_bench_artifact_path(genome_id(chrom)))
    stopifnot(identical(panel$settings, genome_settings))
    panel
  })
  stopifnot(all(vapply(panels, function(panel) {
    identical(panel$selected, panels[[1L]]$selected)
  }, logical(1L))))
  panels
}

genome_simulate_covered <- function(panel, setting, identity) {
  sites <- panel$sites
  n <- rpois(nrow(sites), setting$depth)
  keep <- which(n > 0)
  h <- as.integer(panel$genotypes[[identity$sample]][keep])
  donor <- as.integer(panel$genotypes[[identity$donor]][keep])
  a <- (h %/% 2L + h %% 2L) / 2
  if (setting$relationship == "ibd1") {
    transmitted <- if (identity$sample_haplotype == 1L) h %/% 2L else h %% 2L
    independent <- if (identity$donor_haplotype == 1L) donor %/% 2L else donor %% 2L
    b <- (transmitted + independent) / 2
  } else {
    b <- (donor %/% 2L + donor %% 2L) / 2
  }
  p <- setting$error + (1 - 2 * setting$error) * ((1 - setting$share) * a + setting$share * b)
  alt <- rbinom(length(keep), n[keep], p)
  artefact <- runif(length(n)) < genome_settings$artefact_weight
  artefact_uniform <- numeric(length(n))
  artefact_uniform[artefact] <- runif(sum(artefact))
  covered_artefact <- artefact[keep]
  alt[covered_artefact] <- floor(artefact_uniform[keep][covered_artefact] *
    (n[keep][covered_artefact] + 1L))
  list(counts = data.frame(
    chrom = sites$chrom[keep], pos = sites$pos[keep], ref_count = n[keep] - alt,
    alt_count = alt, other_count = 0L, sample_af = sites$f[keep],
    candidate_af = if (setting$relationship == "other") sites$g[keep] else sites$f[keep]
  ), qc = c(sites = length(n), zero_depth_sites = sum(n == 0),
            excluded_sites = sum(n > genome_settings$max_depth), unavailable_sites = 0))
}

genome_counts <- function(panels, setting, replicate, covered_only = TRUE) {
  seed <- genome_declaration$seed + setting$stratum * 1000L + replicate
  set.seed(seed)
  identity <- genome_draw_identity(panels[[1L]]$selected, setting$relationship)
  qc <- c(sites = 0, zero_depth_sites = 0, excluded_sites = 0, unavailable_sites = 0)
  counts <- lapply(panels, function(panel) {
    if (covered_only) {
      simulated <- genome_simulate_covered(panel, setting, identity)
      qc <<- qc + simulated$qc
      simulated$counts
    } else {
      x <- genome_simulate(panel, setting, identity)
      depth <- x$ref_count + x$alt_count + x$other_count
      qc <<- qc + c(nrow(x), sum(depth == 0), sum(depth > genome_settings$max_depth), 0)
      x
    }
  })
  counts <- do.call(rbind, counts)
  rownames(counts) <- NULL
  counts$sample_id <- sprintf("%d-%d", setting$stratum, replicate)
  counts$candidate <- if (setting$relationship == "other") "AFR" else "matched"
  list(counts = counts, identity = identity, qc = qc, seed = seed)
}

genome_bulk_counts <- function(con, counts) {
  duckdb_register(con, "counts_frame", counts)
  on.exit(duckdb_unregister(con, "counts_frame"))
  dbExecute(con, "CREATE OR REPLACE TEMP TABLE counts AS SELECT * FROM counts_frame")
}

genome_batch_fit <- function(con, counts, error, threads, materialize = TRUE) {
  dbExecute(con, paste("SET threads =", threads))
  started <- proc.time()[["elapsed"]]
  if (materialize) {
    genome_bulk_counts(con, counts)
    on.exit(dbExecute(con, "DROP TABLE counts"))
  } else {
    duckdb_register(con, "counts", counts)
    on.exit(duckdb_unregister(con, "counts"))
  }
  loading_seconds <- proc.time()[["elapsed"]] - started
  profile <- tempfile(fileext = ".json")
  on.exit(unlink(profile), add = TRUE)
  dbExecute(con, "PRAGMA enable_profiling='json'")
  dbExecute(con, paste("SET profiling_output =", dbQuoteString(con, profile)))
  started <- proc.time()[["elapsed"]]
  fit <- genome_fit(con, error)
  seconds <- proc.time()[["elapsed"]] - started
  dbExecute(con, "PRAGMA disable_profiling")
  metrics <- jsonlite::fromJSON(profile, simplifyVector = FALSE)
  stopifnot(metrics$system_peak_temp_dir_size == 0)
  list(fit = fit, seconds = seconds, loading_seconds = loading_seconds,
       peak_buffer_bytes = metrics$system_peak_buffer_memory,
       spill_bytes = metrics$system_peak_temp_dir_size)
}

genome_exact_fields <- c(
  "share", "share_lower", "share_upper", "share_se", "singles_share", "pair_residual",
  "status", "relation_diagnostic", "blocks", "pair_blocks"
)

genome_order_fit <- function(fit) {
  fit <- fit[order(fit$sample, fit$candidate, fit$relationship_model), ]
  rownames(fit) <- NULL
  fit
}
