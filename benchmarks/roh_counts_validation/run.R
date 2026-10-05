# Decode the read-count validation and write every observation as CSV.
#
# (a) Agreement: duckhts_roh_af_table on the staged genotypes (the reference)
#     against duckhts_roh_counts on the staged read counts, same sites and af.
# (b) Titration: the counts of a receiver and a contaminant are thinned and
#     added per allele, then decoded without and with the contamination term.
#     The reference is the decode of the receiver's unmixed counts.
#
# Inputs come from the duckhtsbench registry. Run from the repository root.

source("benchmarks/roh_counts_validation/declaration.R")
source("benchmarks/roh_counts_validation/metrics.R")
source("benchmarks/roh_counts_validation/inputs.R")

roh_validation_connect <- function(extension, threads) {
  driver <- duckdb::duckdb(dbdir = ":memory:", shared_home = FALSE,
    allow_extensions = TRUE, config = list(threads = as.character(threads),
      allow_unsigned_extensions = "true", autoinstall_known_extensions = "false",
      autoload_known_extensions = "false"))
  con <- DBI::dbConnect(driver)
  DBI::dbExecute(con, sprintf("LOAD %s", DBI::dbQuoteString(con, normalizePath(extension))))
  con
}

# Thin and add the counts of one allele: Binomial(receiver, 1 - alpha) reads
# stay and Binomial(contaminant, alpha) reads arrive.
roh_validation_mix_counts <- function(receiver, contaminant, alpha) {
  kept <- stats::rbinom(length(receiver), receiver, 1 - alpha)
  arrived <- stats::rbinom(length(contaminant), contaminant, alpha)
  list(mixed = kept + arrived, kept = sum(kept), arrived = sum(arrived))
}

roh_validation_runs_of <- function(runs, sample) {
  runs[runs$sample == sample, , drop = FALSE]
}

roh_validation_run <- function(extension, output_dir, threads = 2L,
                               declaration = roh_validation_declaration) {
  samples <- declaration$samples
  con <- roh_validation_connect(extension, threads)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  q <- function(value) as.character(DBI::dbQuoteString(con, value))
  count_paths <- vapply(declaration$count_artifacts[samples],
                        duckhtsbench::duckhts_bench_artifact_path, character(1L))
  genotype_path <- duckhtsbench::duckhts_bench_artifact_path(declaration$genotype_artifact)
  for (path in c(count_paths, genotype_path)) {
    if (!file.exists(path)) stop("missing staged input: ", path, call. = FALSE)
  }
  roh_validation_check_inputs(con, declaration)

  DBI::dbExecute(con, sprintf(paste0(
    "CREATE TABLE counts AS SELECT sample_id, chrom, pos, ref, alt, ref_count, alt_count, af ",
    "FROM read_parquet([%s])"), paste(q(count_paths), collapse = ", ")))
  # Every sample must have the same sites and frequencies, so one af relation
  # serves both paths.
  shape <- DBI::dbGetQuery(con, paste0(
    "SELECT count(*) AS rows, count(DISTINCT sample_id) AS samples, ",
    "count(DISTINCT (chrom, pos, ref, alt, af)) AS sites, count(DISTINCT chrom) AS chroms, ",
    "min(pos) AS first_pos, max(pos) AS last_pos, ",
    "count(*) FILTER (WHERE ref_count IS NULL OR alt_count IS NULL) AS null_counts FROM counts"))
  if (shape$samples != length(samples) || shape$chroms != 1 || shape$null_counts != 0 ||
      shape$rows != shape$samples * shape$sites) {
    stop("the count artifacts do not share one chromosome's sites and frequencies",
         call. = FALSE)
  }
  span_bases <- as.numeric(shape$last_pos) - as.numeric(shape$first_pos) + 1
  DBI::dbExecute(con, sprintf(
    "CREATE VIEW af_sites AS SELECT chrom, pos, ref, alt, af FROM counts WHERE sample_id = %s",
    q(samples[[1L]])))

  # (a) Agreement.
  vcf_runs <- DBI::dbGetQuery(con, sprintf(paste0(
    "SELECT * FROM duckhts_roh_af_table(%s, 'af_sites', gt_error := %s) ",
    "ORDER BY sample, chrom, start"), q(genotype_path), format(declaration$gt_error)))
  count_runs <- DBI::dbGetQuery(con, paste0(
    "SELECT * FROM duckhts_roh_counts('counts', contamination := 0.0) ",
    "ORDER BY sample, chrom, start"))
  if (length(setdiff(c(vcf_runs$sample, count_runs$sample), samples))) {
    stop("a decode returned a sample outside the declared ones", call. = FALSE)
  }
  agreement <- do.call(rbind, lapply(samples, function(sample) {
    data.frame(sample = sample, roh_run_verdict(
      roh_run_comparison_by_class(roh_validation_runs_of(count_runs, sample),
                                  roh_validation_runs_of(vcf_runs, sample),
                                  span_bases, declaration$long_run_bases), declaration))
  }))
  agreement_runs <- rbind(data.frame(path = "vcf_gt", vcf_runs),
                          data.frame(path = "read_counts", count_runs))

  # (b) Titration.
  counts <- DBI::dbGetQuery(con,
    "SELECT sample_id, chrom, pos, ref_count, alt_count, af FROM counts ORDER BY sample_id, pos")
  by_sample <- split(counts, counts$sample_id)
  positions <- by_sample[[samples[[1L]]]]$pos
  for (sample in samples) {
    if (!identical(by_sample[[sample]]$pos, positions)) {
      stop("the count artifacts are not aligned by position", call. = FALSE)
    }
  }
  grid <- expand.grid(alpha = declaration$alphas, contaminant = samples, receiver = samples,
                      stringsAsFactors = FALSE)
  grid <- grid[grid$receiver != grid$contaminant, c("receiver", "contaminant", "alpha")]
  grid$cell <- seq_len(nrow(grid))
  grid$seed <- declaration$seed + grid$cell
  grid$mixture <- sprintf("%s+%s@%.2f", grid$receiver, grid$contaminant, grid$alpha)
  grid$receiver_reads <- NA_real_
  grid$contaminant_reads <- NA_real_
  mixed <- vector("list", nrow(grid))
  for (index in seq_len(nrow(grid))) {
    receiver <- by_sample[[grid$receiver[[index]]]]
    contaminant <- by_sample[[grid$contaminant[[index]]]]
    set.seed(grid$seed[[index]])
    ref <- roh_validation_mix_counts(receiver$ref_count, contaminant$ref_count,
                                     grid$alpha[[index]])
    alt <- roh_validation_mix_counts(receiver$alt_count, contaminant$alt_count,
                                     grid$alpha[[index]])
    mixed[[index]] <- data.frame(sample_id = grid$mixture[[index]], chrom = receiver$chrom,
                                 pos = receiver$pos, ref_count = ref$mixed,
                                 alt_count = alt$mixed, af = receiver$af,
                                 stringsAsFactors = FALSE)
    grid$receiver_reads[index] <- ref$kept + alt$kept
    grid$contaminant_reads[index] <- ref$arrived + alt$arrived
  }
  grid$realised_alpha <- grid$contaminant_reads / (grid$receiver_reads + grid$contaminant_reads)
  grid$mean_depth <- (grid$receiver_reads + grid$contaminant_reads) / length(positions)

  titration_runs <- list()
  titration <- list()
  for (alpha in declaration$alphas) {
    cells <- grid[grid$alpha == alpha, , drop = FALSE]
    DBI::dbWriteTable(con, "mixed_counts", do.call(rbind, mixed[cells$cell]), overwrite = TRUE)
    for (contamination in c(0, alpha)) {
      runs <- DBI::dbGetQuery(con, sprintf(paste0(
        "SELECT * FROM duckhts_roh_counts('mixed_counts', contamination := %s) ",
        "ORDER BY sample, chrom, start"), sprintf("%.2f", contamination)))
      term <- if (contamination == 0) "none" else "alpha"
      for (index in seq_len(nrow(cells))) {
        cell_runs <- roh_validation_runs_of(runs, cells$mixture[[index]])
        titration[[length(titration) + 1L]] <- data.frame(
          cells[index, c("cell", "receiver", "contaminant", "alpha", "seed",
                         "realised_alpha", "mean_depth")],
          contamination_term = term, contamination = contamination,
          roh_run_verdict(roh_run_comparison_by_class(
            cell_runs, roh_validation_runs_of(count_runs, cells$receiver[[index]]),
            span_bases, declaration$long_run_bases), declaration), row.names = NULL)
        if (nrow(cell_runs)) {
          titration_runs[[length(titration_runs) + 1L]] <- data.frame(
            cells[index, c("cell", "receiver", "contaminant", "alpha")],
            contamination_term = term, contamination = contamination,
            cell_runs[c("chrom", "start", "end", "length", "n_markers", "quality")],
            row.names = NULL)
        }
      }
    }
  }
  titration <- do.call(rbind, titration)
  titration_runs <- do.call(rbind, titration_runs)

  revision <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
  dirty <- length(system2("git", c("status", "--porcelain", "--untracked-files=no"),
                          stdout = TRUE)) > 0L
  extension_version <- DBI::dbGetQuery(con, paste0(
    "SELECT extension_version FROM duckdb_extensions() WHERE extension_name = 'duckhts'"))
  metadata <- data.frame(
    key = c("revision", "tracked_changes", "input_identities", "extension_version",
            "duckdb_version", "r_version",
            "rng_kind", "seed", "threads", "sites", "first_pos", "last_pos", "span_bases",
            "gt_error", "long_run_bases", "max_froh_difference", "min_jaccard_long",
            "min_jaccard_all"),
    value = c(revision, if (dirty) "yes" else "no", "match the registry",
              extension_version$extension_version[[1L]],
              DBI::dbGetQuery(con, "SELECT version() AS v")$v, R.version.string,
              paste(RNGkind(), collapse = "/"), declaration$seed, threads,
              format(shape$sites), format(shape$first_pos), format(shape$last_pos),
              format(span_bases, scientific = FALSE), declaration$gt_error,
              format(declaration$long_run_bases, scientific = FALSE),
              declaration$max_froh_difference, declaration$min_jaccard_long,
              declaration$min_jaccard_all),
    stringsAsFactors = FALSE)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write <- function(data, name) {
    utils::write.csv(data, file.path(output_dir, name), row.names = FALSE)
  }
  write(metadata, "run_metadata.csv")
  write(agreement_runs, "agreement_runs.csv")
  write(agreement, "agreement_metrics.csv")
  write(titration_runs, "titration_runs.csv")
  write(titration, "titration_metrics.csv")
  list(metadata = metadata, agreement = agreement, titration = titration)
}
