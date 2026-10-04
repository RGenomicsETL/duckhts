# Stage synthetic children with planted autozygous segments, and their truth.
#
# The source is a phased, biallelic, single-contig BCF. For every sample the
# plan draws one segment of each declared length. Inside a planted segment
# haplotype 1 is copied onto haplotype 2, so every genotype there is homozygous
# for the allele of haplotype 1. Everything else is copied: sites, INFO, sample
# names and the genotypes outside the segments. INFO is not recomputed, so
# INFO/AC and INFO/AN still describe the source genotypes.
#
# The truth is the plan. It depends on the seed, the sample order and the
# record positions of the source, and on no caller's output.
#
# Coordinates are 1-based and inclusive: a segment changes the records with
# start <= POS <= end.

# One segment of each length in `lengths_bp` per sample, drawn in the order of
# `samples` and of `lengths_bp`. A start is uniform over the positions that keep
# the segment between the first and the last record. A draw is rejected when the
# segment overlaps a window of `window_bp` bases holding fewer than
# `min_window_records` records, or lies closer than `min_gap_bp` bases to a
# segment already planted in the same sample.
roh_synthetic_plan <- function(samples, positions, seed, lengths_bp, min_gap_bp,
                               min_window_records, window_bp, max_draws = 10000L) {
  settings <- c(seed, lengths_bp, min_gap_bp, min_window_records, window_bp)
  if (anyNA(settings) || any(settings < 0) || any(settings != round(settings)) ||
      !length(lengths_bp) || any(lengths_bp < 1) || window_bp < 1) {
    stop("planting settings must be non-negative integers", call. = FALSE)
  }
  if (!length(samples) || anyDuplicated(samples) || !length(positions) ||
      is.unsorted(positions)) {
    stop("planting needs unique samples and sorted record positions", call. = FALSE)
  }
  first <- positions[[1L]]
  last <- positions[[length(positions)]]
  if (any(last - first + 1 < lengths_bp)) {
    stop("a planted segment is longer than the span of the records", call. = FALSE)
  }
  window_records <- tabulate((positions - 1) %/% window_bp + 1,
                             nbins = (last - 1) %/% window_bp + 1)
  blocked <- window_records < min_window_records

  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  plan <- vector("list", length(samples))
  for (sample_index in seq_along(samples)) {
    starts <- numeric()
    ends <- numeric()
    for (length_bp in lengths_bp) {
      accepted <- FALSE
      for (draw in seq_len(max_draws)) {
        start <- first + sample.int(last - length_bp - first + 2, 1L) - 1
        end <- start + length_bp - 1
        windows <- seq.int((start - 1) %/% window_bp + 1, (end - 1) %/% window_bp + 1)
        apart <- all(start > ends + min_gap_bp | end < starts - min_gap_bp)
        if (!any(blocked[windows]) && apart) {
          accepted <- TRUE
          break
        }
      }
      if (!accepted) stop("could not place a planted segment", call. = FALSE)
      starts <- c(starts, start)
      ends <- c(ends, end)
    }
    plan[[sample_index]] <- data.frame(
      sample_id = samples[[sample_index]], segment = seq_along(lengths_bp),
      start = starts, end = ends, length_bp = lengths_bp, stringsAsFactors = FALSE)
  }
  do.call(rbind, plan)
}

roh_synthetic_bcftools <- function(bcftools, arguments) {
  value <- system2(bcftools, shQuote(arguments), stdout = TRUE, stderr = FALSE)
  if (!is.null(attr(value, "status"))) {
    stop("bcftools command failed: ", paste(arguments, collapse = " "), call. = FALSE)
  }
  value
}

roh_synthetic_close <- function(connection, what) {
  status <- close(connection)
  if (!is.null(status) && status != 0L) stop(what, " failed", call. = FALSE)
}

# `expected`, when given, is list(records = , samples = , records_sha256 = ,
# truth_rows = , truth_sha256 = ); any mismatch is an error and publishes
# nothing. records_sha256 is the SHA-256 of the record lines (`bcftools view -H`)
# and truth_sha256 that of the truth CSV. The BCF's own bytes are not compared,
# because bcftools writes its command line and date into the header.
stage_roh_synthetic <- function(source_bcf, output_bcf, output_truth, seed, lengths_bp,
                                min_gap_bp, min_window_records, window_bp,
                                bcftools = "/usr/local/bin/bcftools", expected = NULL,
                                chunk_records = 20000L) {
  if (!file.exists(source_bcf)) stop("missing staging input: ", source_bcf, call. = FALSE)
  if (!file.exists(bcftools)) stop("bcftools is required: ", bcftools, call. = FALSE)
  samples <- roh_synthetic_bcftools(bcftools, c("query", "-l", source_bcf))
  sites <- roh_synthetic_bcftools(bcftools, c("query", "-f", "%CHROM\\t%POS\\n", source_bcf))
  chrom <- unique(sub("\t.*$", "", sites))
  if (length(chrom) != 1L) stop("the source BCF must hold one contig", call. = FALSE)
  positions <- as.numeric(sub("^.*\t", "", sites))
  plan <- roh_synthetic_plan(samples, positions, seed, lengths_bp, min_gap_bp,
                             min_window_records, window_bp)
  plan$records <- 0
  plan$source_heterozygous <- 0
  plan_column <- match(plan$sample_id, samples)

  dir.create(dirname(output_bcf), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(output_truth), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(output_bcf, ".partial-", Sys.getpid(), ".bcf")
  temporary_index <- paste0(temporary, ".csi")
  temporary_truth <- paste0(output_truth, ".partial-", Sys.getpid())
  on.exit(unlink(c(temporary, temporary_index, temporary_truth), force = TRUE), add = TRUE)

  header <- roh_synthetic_bcftools(bcftools, c("view", "-h", source_bcf))
  planting <- sprintf(paste0(
    "##roh_synthetic_planting=seed=%d;lengths_bp=%s;min_gap_bp=%d;",
    "min_window_records=%d;window_bp=%d"), seed, paste(sprintf("%.0f", lengths_bp), collapse = ","),
    min_gap_bp, min_window_records, window_bp)
  header <- c(header[-length(header)], planting, header[[length(header)]])

  site_lines <- pipe(paste(shQuote(bcftools), "view -H -G", shQuote(source_bcf)), "r")
  genotype_lines <- pipe(paste(shQuote(bcftools), "query -f '[\\t%GT]\\n'",
                               shQuote(source_bcf)), "r")
  writer <- pipe(paste(shQuote(bcftools), "view --no-update -Ob -o", shQuote(temporary), "-"), "w")
  open_connections <- TRUE
  on.exit(if (open_connections) {
    for (connection in list(site_lines, genotype_lines, writer)) try(close(connection), silent = TRUE)
  }, add = TRUE)
  writeLines(header, writer)

  # Every genotype is tab, allele, bar, allele: four characters per sample.
  genotype_width <- 4L
  genotype_pattern <- sprintf("^(\t[01][|][01]){%d}$", length(samples))
  done <- 0
  repeat {
    fixed <- readLines(site_lines, n = chunk_records)
    genotypes <- readLines(genotype_lines, n = chunk_records)
    if (length(fixed) != length(genotypes)) {
      stop("site and genotype streams differ in length", call. = FALSE)
    }
    if (!length(fixed)) break
    rows <- done + seq_along(fixed)
    if (rows[[length(rows)]] > length(positions) ||
        !all(startsWith(fixed, sprintf("%s\t%.0f\t", chrom, positions[rows])))) {
      stop("site stream differs from the planned record positions", call. = FALSE)
    }
    if (!all(grepl(genotype_pattern, genotypes, perl = TRUE))) {
      stop("every genotype must be phased, diploid, biallelic and called", call. = FALSE)
    }
    chunk_positions <- positions[rows]
    touched <- which(plan$start <= chunk_positions[[length(rows)]] &
                       plan$end >= chunk_positions[[1L]])
    for (index in touched) {
      inside <- which(chunk_positions >= plan$start[[index]] &
                        chunk_positions <= plan$end[[index]])
      if (!length(inside)) next
      offset <- genotype_width * (plan_column[[index]] - 1L)
      first_allele <- substr(genotypes[inside], offset + 2L, offset + 2L)
      second_allele <- substr(genotypes[inside], offset + 4L, offset + 4L)
      plan$records[[index]] <- plan$records[[index]] + length(inside)
      plan$source_heterozygous[[index]] <- plan$source_heterozygous[[index]] +
        sum(first_allele != second_allele)
      substr(genotypes[inside], offset + 4L, offset + 4L) <- first_allele
    }
    writeLines(paste0(fixed, "\tGT", genotypes), writer)
    done <- done + length(fixed)
  }
  open_connections <- FALSE
  roh_synthetic_close(site_lines, "reading the source sites")
  roh_synthetic_close(genotype_lines, "reading the source genotypes")
  roh_synthetic_close(writer, "writing the synthetic BCF")
  if (done != length(positions)) stop("the source records were not all copied", call. = FALSE)

  status <- system2(bcftools, shQuote(c("index", "-f", temporary)))
  if (status != 0L || !file.exists(temporary_index)) {
    stop("bcftools failed to index the synthetic BCF", call. = FALSE)
  }
  if (!identical(roh_synthetic_bcftools(bcftools, c("query", "-l", temporary)), samples)) {
    stop("synthetic BCF samples differ from the source", call. = FALSE)
  }
  records <- as.numeric(roh_synthetic_bcftools(bcftools, c("index", "-n", temporary)))
  records_sha256 <- system2("/bin/bash", c("-o", "pipefail", "-c", shQuote(paste(
    shQuote(bcftools), "view -H", shQuote(temporary), "| sha256sum"))), stdout = TRUE)
  if (!is.null(attr(records_sha256, "status")) || length(records_sha256) != 1L) {
    stop("could not digest the synthetic BCF records", call. = FALSE)
  }
  records_sha256 <- sub(" .*$", "", records_sha256)

  truth <- data.frame(sample_id = plan$sample_id, segment = plan$segment, chrom = chrom,
                      plan[c("start", "end", "length_bp", "records", "source_heterozygous")],
                      stringsAsFactors = FALSE)
  truth <- truth[order(match(truth$sample_id, samples), truth$start), , drop = FALSE]
  utils::write.table(format(truth, scientific = FALSE, trim = TRUE), temporary_truth,
                     sep = ",", quote = FALSE, row.names = FALSE)
  identity <- list(records = records, samples = length(samples),
                   records_sha256 = records_sha256, truth_rows = nrow(truth),
                   truth_sha256 = unname(digest::digest(file = temporary_truth, algo = "sha256")))
  if (!is.null(expected)) {
    counts <- c("records", "samples", "truth_rows")
    digests <- c("records_sha256", "truth_sha256")
    same_counts <- unlist(identity[counts]) == as.numeric(unlist(expected[counts]))
    same_digests <- unlist(identity[digests]) == unlist(expected[digests])
    if (length(same_counts) != 3L || length(same_digests) != 2L ||
        !all(same_counts, same_digests)) {
      stop("synthetic children differ from their registered identity: ",
           sprintf("records=%.0f;samples=%.0f;records_sha256=%s;truth_rows=%.0f;truth_sha256=%s",
                   identity$records, identity$samples, identity$records_sha256,
                   identity$truth_rows, identity$truth_sha256),
           call. = FALSE)
    }
  }
  if (!file.rename(temporary, output_bcf) ||
      !file.rename(temporary_index, paste0(output_bcf, ".csi")) ||
      !file.rename(temporary_truth, output_truth)) {
    unlink(c(output_bcf, paste0(output_bcf, ".csi"), output_truth), force = TRUE)
    stop("could not publish the synthetic children", call. = FALSE)
  }
  c(identity, list(output_bcf = output_bcf, output_truth = output_truth))
}

# The synthetic BCF artifact's locator names its source BCF. Its supplier
# identity holds the planting settings and the expected records, samples and
# records_sha256. The truth artifact names the same source and holds the
# expected rows and sha256 of the truth CSV.
stage_roh_synthetic_from_registry <- function(bcf_id, truth_id,
                                              bcftools = "/usr/local/bin/bcftools") {
  registry <- duckhtsbench::duckhts_bench_registry()
  row_of <- function(id) {
    row <- registry[registry$id == id, , drop = FALSE]
    if (nrow(row) != 1L) stop("unknown or non-unique benchmark artifact: ", id, call. = FALSE)
    row
  }
  bcf_row <- row_of(bcf_id)
  truth_row <- row_of(truth_id)
  source_id <- sub("^artifact:", "", bcf_row$locator[[1L]])
  if (!source_id %in% registry$id || !identical(truth_row$locator[[1L]], bcf_row$locator[[1L]])) {
    stop("synthetic BCF and truth must name one registered source BCF", call. = FALSE)
  }
  fields <- duckhtsbench:::duckhts_bench_identity_fields(bcf_row$supplier_identity[[1L]])
  truth_fields <- duckhtsbench:::duckhts_bench_identity_fields(truth_row$supplier_identity[[1L]])
  required <- c("seed", "lengths_bp", "min_gap_bp", "min_window_records", "window_bp",
                "records", "samples", "records_sha256")
  if (!all(required %in% names(fields)) || !all(c("rows", "sha256") %in% names(truth_fields))) {
    stop("synthetic artifact identities lack planting settings or expected digests",
         call. = FALSE)
  }
  stage_roh_synthetic(
    source_bcf = duckhtsbench::duckhts_bench_artifact_path(source_id),
    output_bcf = duckhtsbench::duckhts_bench_artifact_path(bcf_id),
    output_truth = duckhtsbench::duckhts_bench_artifact_path(truth_id),
    seed = as.numeric(fields[["seed"]]),
    lengths_bp = as.numeric(strsplit(fields[["lengths_bp"]], ",", fixed = TRUE)[[1L]]),
    min_gap_bp = as.numeric(fields[["min_gap_bp"]]),
    min_window_records = as.numeric(fields[["min_window_records"]]),
    window_bp = as.numeric(fields[["window_bp"]]), bcftools = bcftools,
    expected = list(records = fields[["records"]], samples = fields[["samples"]],
                    records_sha256 = fields[["records_sha256"]],
                    truth_rows = truth_fields[["rows"]], truth_sha256 = truth_fields[["sha256"]]))
}
