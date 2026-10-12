# The read-mixture evaluation: reads of a second 30x sample mixed into a
# receiver's chr20 slice at declared shares, counted at the panel sites as the
# staged counts are, fitted with duckhts_count_error_fit and decoded with
# duckhts_roh_counts_fitted with and without the contamination term.
#
# Run from the repository root after sourcing common.R.

roh_mixture_samtools <- function(arguments) {
  value <- system2(bam_mismatch_tool("samtools"), shQuote(arguments), stdout = TRUE,
                   env = bam_mismatch_tool_environment())
  if (!is.null(attr(value, "status"))) {
    stop("samtools failed: ", paste(arguments, collapse = " "), call. = FALSE)
  }
  value
}

roh_mixture_free_bytes <- function(directory) {
  value <- system2("df", c("-B1", "--output=avail", shQuote(directory)), stdout = TRUE)
  if (!is.null(attr(value, "status")) || length(value) != 2L) {
    stop("could not read the free space of ", directory, call. = FALSE)
  }
  as.numeric(value[[2L]])
}

roh_mixture_label <- function(alpha) sprintf("alpha_%.2f", alpha)

roh_mixture_run <- function(extension, cli, output_dir, size = c("full", "reduced"),
                            declaration = roh_mixture_declaration) {
  size <- match.arg(size)
  alphas <- if (size == "full") declaration$alphas else declaration$reduced_alphas
  extension <- normalizePath(extension, mustWork = TRUE)
  work <- tempfile("roh-mixtures-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  con <- roh_truth_connect(extension)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  q <- roh_truth_quote
  receiver <- declaration$receiver
  donor <- declaration$donor

  # Inputs and their registered identities.
  slices <- vapply(declaration$slice_artifacts, duckhtsbench::duckhts_bench_artifact_path, "")
  slice_identities <- do.call(rbind, lapply(names(slices), function(sample) {
    data.frame(sample = sample, artifact = declaration$slice_artifacts[[sample]],
               bam_mismatch_check_staged(declaration$slice_artifacts[[sample]], extension))
  }))
  sites_path <- roh_truth_check_file(declaration$sites_artifact)
  reference <- roh_truth_check_file(declaration$reference_artifact)
  settings <- roh_truth_registered(declaration$count_artifacts[[receiver]])
  donor_settings <- roh_truth_registered(declaration$count_artifacts[[donor]])
  setting_names <- c("assembly", "min_mapq", "min_baseq", "exclude_flags", "overlap_policy")
  if (!identical(settings[setting_names], donor_settings[setting_names])) {
    stop("the two registered count artifacts use different settings", call. = FALSE)
  }
  panel <- file.path(work, "panel.parquet")
  DBI::dbExecute(con, roh_counts_panel_sql(sites_path, settings[["assembly"]], panel))
  shape <- DBI::dbGetQuery(con, sprintf(
    "SELECT count(*) AS sites, min(pos) AS first_pos, max(pos) AS last_pos, any_value(chrom) AS chrom,
            count(DISTINCT chrom) AS chroms FROM read_parquet(%s)", q(sites_path)))
  if (shape$chroms != 1) stop("the panel must hold one chromosome", call. = FALSE)
  sites <- as.numeric(shape$sites)
  span_bases <- as.numeric(shape$last_pos) - as.numeric(shape$first_pos) + 1
  fai <- utils::read.delim(paste0(reference, ".fai"), header = FALSE, stringsAsFactors = FALSE)
  contig_bytes <- as.numeric(fai$V2[fai$V1 == shape$chrom])
  if (length(contig_bytes) != 1L) stop("the reference index lacks ", shape$chrom, call. = FALSE)

  recorder <- roh_truth_recorder(cli, extension, work)
  overhead_decode <- recorder$overhead(declaration$threads)
  overhead_count <- recorder$overhead(declaration$count_threads)
  ceiling <- roh_truth_ceiling_bytes

  count_alignments <- function(name, source, sample_id) {
    output <- file.path(work, sprintf("%s.counts.parquet", name))
    recorder$record(sprintf("count_%s", name), "read_counting", character(0),
      roh_counts_query_sql(source, sample_id, reference, panel, sites_path,
                           as.integer(settings[["min_mapq"]]), as.integer(settings[["min_baseq"]]),
                           as.integer(settings[["exclude_flags"]]), settings[["overlap_policy"]],
                           declaration$worker_count),
      declaration$count_threads, output, "parquet",
      live_bytes = ceiling$read_counting(declaration$worker_count, contig_bytes, sites),
      overhead_mib = overhead_count)
    output
  }

  # Each slice alone reproduces its registered counts.
  alone <- list()
  identities <- list()
  for (sample in c(receiver, donor)) {
    alone[[sample]] <- count_alignments(sample, slices[[sample]], sample)
    registered <- roh_truth_registered(declaration$count_artifacts[[sample]])
    observed <- roh_counts_identity(con, alone[[sample]])
    identities[[sample]] <- data.frame(sample = sample, rows = observed$rows,
                                       counts_sha256 = observed$counts_sha256,
                                       registered_sha256 = registered[["counts_sha256"]],
                                       identical = identical(observed$counts_sha256,
                                                             registered[["counts_sha256"]]))
    if (!identities[[sample]]$identical || observed$rows != as.numeric(registered[["rows"]])) {
      stop("the counts of the ", sample, " slice differ from ", declaration$count_artifacts[[sample]],
           call. = FALSE)
    }
  }
  identities <- do.call(rbind, identities)

  # The mixtures, written and counted one at a time.
  slice_bytes <- sum(file.size(slices))
  mixture_counts <- list()
  reads <- list()
  for (alpha in alphas) {
    free <- roh_mixture_free_bytes(work)
    if (free - 2 * slice_bytes < declaration$min_free_bytes) {
      stop("writing a mixture would leave less than ", declaration$min_free_bytes / 2^30,
           " GiB free", call. = FALSE)
    }
    kept_receiver <- file.path(work, "receiver.cram")
    kept_donor <- file.path(work, "donor.cram")
    mixture <- file.path(work, "mixture.cram")
    subsample <- function(source, share, output) {
      roh_mixture_samtools(c("view", "--no-PG", "-C", "-T", reference, "-@", "2",
                             "--subsample", format(share, digits = 15),
                             "--subsample-seed", format(declaration$subsample_seed),
                             "-o", output, source))
      as.numeric(roh_mixture_samtools(c("view", "-c", "-@", "2", "-T", reference, output)))
    }
    receiver_alignments <- subsample(slices[[receiver]], 1 - alpha, kept_receiver)
    donor_alignments <- subsample(slices[[donor]], alpha, kept_donor)
    roh_mixture_samtools(c("merge", "--no-PG", "-f", "-@", "2", "--reference", reference,
                           "-O", "cram", "--write-index", "-o", mixture, kept_receiver, kept_donor))
    label <- roh_mixture_label(alpha)
    mixture_counts[[label]] <- count_alignments(label, mixture, label)
    reads[[label]] <- data.frame(
      alpha = alpha, receiver_alignments = receiver_alignments, donor_alignments = donor_alignments,
      read_share = donor_alignments / (receiver_alignments + donor_alignments),
      mixture_bytes = file.size(mixture))
    unlink(c(kept_receiver, kept_donor, mixture, paste0(mixture, ".crai")), force = TRUE)
  }
  reads <- do.call(rbind, reads)

  # Fit and decode alpha 0 and every mixture together.
  counts_path <- file.path(work, "mixtures.parquet")
  parts <- c(sprintf("SELECT %s AS sample_id, chrom, pos, ref_count, alt_count, af FROM read_parquet(%s)",
                     q(roh_mixture_label(0)), q(alone[[receiver]])),
             vapply(names(mixture_counts), function(label) {
               sprintf("SELECT sample_id, chrom, pos, ref_count, alt_count, af FROM read_parquet(%s)",
                       q(mixture_counts[[label]]))
             }, ""))
  DBI::dbExecute(con, sprintf("COPY (%s) TO %s (FORMAT parquet)",
                              paste(parts, collapse = " UNION ALL "), q(counts_path)))
  samples <- c(roh_mixture_label(0), names(mixture_counts))
  denominators <- DBI::dbGetQuery(con, sprintf(paste(
    "WITH r AS (SELECT sample_id, chrom, (pos - 1) // 10000000 AS block_index,",
    "floor(af * 64) AS bin, ref_count + alt_count AS depth, alt_count AS alt FROM read_parquet(%s)",
    "WHERE ref_count IS NOT NULL AND alt_count IS NOT NULL AND af > 0 AND af < 1",
    "AND ref_count + alt_count BETWEEN 1 AND 1000)",
    "SELECT (SELECT count(*) FROM (SELECT DISTINCT sample_id, chrom, block_index, bin, depth, alt FROM r)) AS cells,",
    "(SELECT count(*) FROM (SELECT DISTINCT sample_id, depth, alt FROM r)) AS pairs"), q(counts_path)))
  counts_setup <- sprintf("CREATE TEMP VIEW counts AS SELECT * FROM read_parquet(%s);", q(counts_path))
  fit_path <- file.path(work, "fit.parquet")
  recorder$record("mixture_fit", "fit", counts_setup,
    "SELECT * FROM duckhts_count_error_fit('counts') ORDER BY sample", declaration$threads, fit_path,
    "parquet", live_bytes = ceiling$fit(denominators$cells, denominators$pairs, length(samples)),
    overhead_mib = overhead_decode)
  fit <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s) ORDER BY sample", q(fit_path)))
  numbers <- c("seq_error", "contamination", "allele_balance", "artefact_weight")
  if (!setequal(fit$sample, samples) || anyNA(fit[numbers])) {
    stop("the mixture fit lacks a sample or a fitted number", call. = FALSE)
  }
  fit$alpha <- c(0, alphas)[match(fit$sample, samples)]
  fit <- fit[order(fit$alpha), , drop = FALSE]

  runs_columns <- "sample, chrom, start, \"end\", length, n_markers, quality"
  decode <- function(name, parameters) {
    fit_file <- file.path(work, sprintf("fit_%s.parquet", name))
    DBI::dbWriteTable(con, "decode_fit", parameters, overwrite = TRUE)
    DBI::dbExecute(con, sprintf("COPY decode_fit TO %s (FORMAT parquet)", q(fit_file)))
    output <- file.path(work, sprintf("%s.csv", name))
    recorder$record(name, "count_decode", c(counts_setup,
      sprintf("CREATE TEMP VIEW decode_fit AS SELECT * FROM read_parquet(%s);", q(fit_file))),
      sprintf("SELECT %s FROM duckhts_roh_counts_fitted('counts', 'decode_fit') ORDER BY sample, chrom, start",
              runs_columns), declaration$threads, output,
      live_bytes = ceiling$count_decode(length(samples) * sites, sites, declaration$threads),
      overhead_mib = overhead_decode)
    data.frame(decode = name, roh_truth_read_runs(output))
  }
  with_term <- data.frame(fit[c("sample", numbers)], artefact_weight_sd = 0)
  without_term <- with_term
  without_term$contamination <- 0
  runs <- rbind(decode("with_term", with_term), decode("without_term", without_term))

  observations <- recorder$table()
  roh_truth_gate(observations)

  reference_runs <- runs[runs$decode == "with_term" & runs$sample == roh_mixture_label(0), , drop = FALSE]
  froh <- do.call(rbind, lapply(c("with_term", "without_term"), function(name) {
    do.call(rbind, lapply(samples, function(sample) {
      test <- runs[runs$decode == name & runs$sample == sample, , drop = FALSE]
      comparison <- roh_run_comparison(test, reference_runs, span_bases)
      data.frame(decode = name, sample = sample, alpha = c(0, alphas)[match(sample, samples)],
                 runs = nrow(test), froh = comparison$test_froh,
                 erosion = comparison$froh_difference, jaccard = comparison$jaccard)
    }))
  }))
  contamination <- fit$contamination[order(fit$alpha)]
  m1_pass <- all(diff(contamination) > 0)
  m2 <- do.call(rbind, lapply(alphas, function(alpha) {
    with_erosion <- froh$erosion[froh$decode == "with_term" & froh$alpha == alpha]
    without_erosion <- froh$erosion[froh$decode == "without_term" & froh$alpha == alpha]
    data.frame(alpha = alpha, erosion_with_term = with_erosion, erosion_without_term = without_erosion,
               pass = abs(with_erosion) <= abs(without_erosion) + declaration$m2_tolerance)
  }))
  criteria <- data.frame(criterion = c("M1", "M2"), pass = c(m1_pass, all(m2$pass)))

  metadata <- data.frame(
    key = c("receiver", "donor", "alphas", "size", "sites", "span_bases", "samtools",
            "rbcftools", "fit_cells", "fit_pairs", "overhead_decode_mib", "overhead_count_mib"),
    value = c(receiver, donor, paste(alphas, collapse = ","), size, sites,
              format(span_bases, scientific = FALSE), RBCFTools::samtools_version(),
              as.character(utils::packageVersion("RBCFTools")), format(denominators$cells),
              format(denominators$pairs), sprintf("%.1f", overhead_decode), sprintf("%.1f", overhead_count)),
    stringsAsFactors = FALSE)
  write <- function(data, name) utils::write.csv(data, file.path(output_dir, name), row.names = FALSE)
  write(metadata, "mixture_metadata.csv")
  write(slice_identities, "mixture_slice_identities.csv")
  write(identities, "mixture_count_identities.csv")
  write(reads, "mixture_reads.csv")
  write(fit, "mixture_fit.csv")
  write(runs, "mixture_runs.csv")
  write(froh, "mixture_froh.csv")
  write(m2, "mixture_erosion.csv")
  write(criteria, "mixture_criteria.csv")
  write(observations, "mixture_observations.csv")
  list(metadata = metadata, slice_identities = slice_identities, identities = identities, reads = reads,
       fit = fit, froh = froh, m2 = m2, criteria = criteria, observations = observations)
}
