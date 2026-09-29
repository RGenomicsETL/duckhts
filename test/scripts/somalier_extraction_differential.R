#!/usr/bin/env Rscript

# Differential panel-site extraction against Somalier v0.3.4 and an
# independent SAM/CIGAR oracle. The generated fixture is deliberately small,
# but retains every compared site and exercises hileup v0.1.0 filtering before
# encounter-ordered overlapping-mate suppression.
# Stage the third argument with
# duckhtsbench:::duckhts_bench_stage_somalier_v034().

SOMALIER_COMMIT <- "ff58fdade8f4f8293d904f10e0a4a13f1fac808d"
SOMALIER_LINUX_SHA256 <-
  "18717c205a9c4b65d479f1d2cf069a30b047a5378edf544d403d4081f06a3a78"
HILEUP_COMMIT <- "3133320d9660620a4a7f4da5647c85bf2b2738a6"

fail <- function(...) stop(..., call. = FALSE)

run <- function(command, args, label, stdout = TRUE, stderr = TRUE) {
  output <- suppressWarnings(system2(command, args, stdout = stdout,
    stderr = stderr))
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    fail(label, " exited ", status, if (length(output))
      paste0(":\n", paste(output, collapse = "\n")) else "")
  }
  output
}

sha256 <- function(path) {
  command <- Sys.which("sha256sum")
  if (!nzchar(command)) fail("sha256sum is required")
  fields <- strsplit(run(command, shQuote(path), "sha256sum"),
    "[[:space:]]+")[[1L]]
  fields[[1L]]
}

replace_base <- function(sequence, index1, base) {
  substr(sequence, index1, index1) <- base
  sequence
}

sam_record <- function(qname, flag, pos, mapq, cigar, sequence, quality,
                       rnext = "*", pnext = 0L, tlen = 0L) {
  paste(qname, flag, "chr1", pos, mapq, cigar, rnext, pnext, tlen,
    sequence, quality, sep = "\t")
}

write_fixture <- function(directory, samtools) {
  reference <- file.path(directory, "reference.fa")
  writeLines(c(">chr1", strwrap(paste(rep("A", 400L), collapse = ""), 60L)),
    reference)
  run(samtools, c("faidx", shQuote(reference)), "samtools faidx")

  bases20 <- paste(rep("A", 20L), collapse = "")
  qualities20 <- paste(rep("I", 20L), collapse = "")
  records <- c(
    sam_record("basic-a", 0L, 45L, 60L, "20M", bases20, qualities20),
    sam_record("basic-c", 0L, 45L, 60L, "20M",
      replace_base(bases20, 6L, "C"), qualities20),
    sam_record("basic-g", 0L, 45L, 60L, "20M",
      replace_base(bases20, 6L, "G"), qualities20),
    sam_record("mapq-zero", 0L, 45L, 0L, "20M",
      replace_base(bases20, 6L, "T"), qualities20),
    sam_record("duplicate", 1024L, 45L, 60L, "20M",
      replace_base(bases20, 6L, "C"), qualities20),
    sam_record("secondary", 256L, 45L, 60L, "20M",
      replace_base(bases20, 6L, "C"), qualities20),
    sam_record("qc-fail", 512L, 45L, 60L, "20M",
      replace_base(bases20, 6L, "C"), qualities20),
    sam_record("supplementary", 2048L, 45L, 60L, "20M",
      replace_base(bases20, 6L, "C"), qualities20),
    sam_record("deletion", 0L, 75L, 60L, "5M1D15M",
      paste(rep("C", 20L), collapse = ""), qualities20),
    sam_record("refskip", 0L, 75L, 60L, "5M1N15M",
      paste(rep("C", 20L), collapse = ""), qualities20),
    sam_record("insertion", 0L, 75L, 60L, "5M1I15M",
      replace_base(paste(rep("A", 21L), collapse = ""), 7L, "C"),
      paste(rep("I", 21L), collapse = "")),
    sam_record("softclip", 0L, 75L, 60L, "5S20M",
      replace_base(paste(rep("A", 25L), collapse = ""), 11L, "G"),
      paste(rep("I", 25L), collapse = "")),
    sam_record("missing-qual", 0L, 95L, 60L, "20M",
      replace_base(bases20, 6L, "C"), "*"),
    sam_record("ordinary-overlap", 99L, 130L, 60L, "40M",
      paste(rep("A", 40L), collapse = ""), paste(rep("I", 40L), collapse = ""),
      "=", 140L, 50L),
    sam_record("ordinary-overlap", 147L, 140L, 60L, "40M",
      replace_base(paste(rep("A", 40L), collapse = ""), 11L, "C"),
      paste(rep("I", 40L), collapse = ""), "=", 130L, -50L),
    sam_record("filtered-overlap", 99L, 160L, 0L, "40M",
      paste(rep("A", 40L), collapse = ""), paste(rep("I", 40L), collapse = ""),
      "=", 170L, 50L),
    sam_record("filtered-overlap", 147L, 170L, 60L, "40M",
      replace_base(paste(rep("A", 40L), collapse = ""), 11L, "C"),
      paste(rep("I", 40L), collapse = ""), "=", 160L, -50L),
    sam_record("ambiguous", 0L, 195L, 60L, "20M",
      replace_base(bases20, 6L, "N"), qualities20),
    sam_record("reverse-other", 16L, 195L, 60L, "20M",
      replace_base(bases20, 6L, "T"), qualities20)
  )
  sam <- file.path(directory, "reads.sam")
  header_lines <- c("@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:300",
    "@RG\tID:rg1\tSM:sample-1")
  header <- file.path(directory, "header.sam")
  writeLines(header_lines, header)
  writeLines(c(header_lines, records), sam)
  bam <- file.path(directory, "reads.bam")
  run(samtools, c("view", "-b", "-o", shQuote(bam), shQuote(sam)),
    "samtools BAM conversion")
  run(samtools, c("index", shQuote(bam)), "samtools BAM index")
  cram <- file.path(directory, "reads.cram")
  run(samtools, c("view", "-C", "-T", shQuote(reference), "-o",
    shQuote(cram), shQuote(sam)), "samtools CRAM conversion")
  run(samtools, c("reheader", "-P", "-i", shQuote(header), shQuote(cram)),
    "samtools CRAM reheader")
  run(samtools, c("index", shQuote(cram)), "samtools CRAM index")

  panel <- data.frame(
    site_index = 0:6,
    region = "chr1",
    position = c(50L, 80L, 100L, 150L, 180L, 200L, 350L),
    allele_a = "A",
    allele_b = "C",
    relation = c(rep("equal", 6L), "alignment_position_unavailable"),
    stringsAsFactors = FALSE)
  sites <- file.path(directory, "sites.vcf")
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=chr1,length=300>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    paste(panel$region, panel$position, paste0("site", panel$site_index),
      panel$allele_a, panel$allele_b, ".", "PASS", ".", sep = "\t")), sites)
  list(reference = reference, bam = bam, cram = cram, panel = panel,
    sites = sites)
}

parse_cigar <- function(cigar) {
  starts <- gregexpr("[0-9]+[MIDNSHP=X]", cigar, perl = TRUE)[[1L]]
  if (identical(starts, -1L)) fail("invalid CIGAR: ", cigar)
  tokens <- regmatches(cigar, list(starts))[[1L]]
  if (paste(tokens, collapse = "") != cigar) fail("invalid CIGAR: ", cigar)
  data.frame(length = as.integer(sub("[MIDNSHP=X]$", "", tokens)),
    operation = sub("^[0-9]+", "", tokens), stringsAsFactors = FALSE)
}

reference_span <- function(cigar) {
  operations <- parse_cigar(cigar)
  sum(operations$length[operations$operation %in% c("M", "D", "N", "=", "X")])
}

query_index <- function(pos1, cigar, site1) {
  operations <- parse_cigar(cigar)
  reference <- pos1
  query <- 1L
  for (i in seq_len(nrow(operations))) {
    length <- operations$length[[i]]
    operation <- operations$operation[[i]]
    if (operation %in% c("M", "=", "X")) {
      if (site1 >= reference && site1 < reference + length) {
        return(query + site1 - reference)
      }
      reference <- reference + length
      query <- query + length
    } else if (operation %in% c("D", "N")) {
      if (site1 >= reference && site1 < reference + length) return(NA_integer_)
      reference <- reference + length
    } else if (operation %in% c("I", "S")) {
      query <- query + length
    } else if (!operation %in% c("H", "P")) {
      fail("unsupported CIGAR operation: ", operation)
    }
  }
  NA_integer_
}

read_alignments <- function(path, reference, samtools) {
  lines <- run(samtools, c("view", "-T", shQuote(reference), shQuote(path)),
    "samtools view", stderr = FALSE)
  fields <- strsplit(lines, "\t", fixed = TRUE)
  if (any(lengths(fields) < 11L)) fail("samtools emitted a truncated SAM record")
  data.frame(
    qname = vapply(fields, `[[`, character(1), 1L),
    flag = as.integer(vapply(fields, `[[`, character(1), 2L)),
    region = vapply(fields, `[[`, character(1), 3L),
    position = as.integer(vapply(fields, `[[`, character(1), 4L)),
    mapq = as.integer(vapply(fields, `[[`, character(1), 5L)),
    cigar = vapply(fields, `[[`, character(1), 6L),
    mate_region = vapply(fields, `[[`, character(1), 7L),
    mate_position = as.integer(vapply(fields, `[[`, character(1), 8L)),
    sequence = vapply(fields, `[[`, character(1), 10L),
    quality = vapply(fields, `[[`, character(1), 11L),
    stringsAsFactors = FALSE)
}

oracle_counts <- function(alignments, panel) {
  output <- panel[c("site_index", "region", "position", "allele_a", "allele_b")]
  output$a <- output$b <- output$other <- integer(nrow(output))
  for (site_row in seq_len(nrow(panel))) {
    site <- panel[site_row, ]
    spans <- vapply(alignments$cigar, reference_span, integer(1))
    candidates <- alignments$region == site$region &
      alignments$position <= site$position &
      alignments$position + spans > site$position
    reads <- alignments[candidates, , drop = FALSE]
    tracked <- character()
    for (read_row in seq_len(nrow(reads))) {
      read <- reads[read_row, ]
      eligible <- read$mapq >= 1L && bitwAnd(read$flag, 1796L) == 0L
      if (!eligible) next
      index <- query_index(read$position, read$cigar, site$position)
      if (is.na(index)) next
      if (read$qname %in% tracked) {
        tracked <- tracked[tracked != read$qname]
        next
      }
      base <- substr(read$sequence, index, index)
      if (base == site$allele_a) {
        output$a[[site_row]] <- output$a[[site_row]] + 1L
      } else if (base == site$allele_b) {
        output$b[[site_row]] <- output$b[[site_row]] + 1L
      } else {
        output$other[[site_row]] <- output$other[[site_row]] + 1L
      }
      primary <- bitwAnd(read$flag, bitwOr(256L, 2048L)) == 0L
      same_contig <- read$mate_region == "=" || read$mate_region == read$region
      start0 <- read$position - 1L
      mate0 <- read$mate_position - 1L
      end0 <- start0 + reference_span(read$cigar)
      if (primary && same_contig && read$mate_position > 0L &&
          end0 > mate0 && mate0 <= site$position - 1L && start0 <= mate0) {
        tracked <- union(tracked, read$qname)
      }
    }
  }
  output
}

read_somalier <- function(path, panel) {
  connection <- file(path, "rb")
  on.exit(close(connection))
  version <- readBin(connection, integer(), 1L, size = 1L, signed = FALSE)
  name_length <- readBin(connection, integer(), 1L, size = 1L, signed = FALSE)
  sample <- rawToChar(readBin(connection, raw(), name_length))
  site_counts <- readBin(connection, integer(), 3L, size = 2L,
    signed = FALSE, endian = "little")
  if (version != 2L || sample != "sample-1" ||
      !identical(site_counts, c(nrow(panel), 0L, 0L))) {
    fail("unexpected Somalier digest header")
  }
  values <- readBin(connection, integer(), 3L * nrow(panel), size = 4L,
    endian = "little")
  if (length(values) != 3L * nrow(panel) ||
      length(readBin(connection, raw(), 1L)) != 0L) {
    fail("unexpected Somalier digest length")
  }
  counts <- matrix(values, ncol = 3L, byrow = TRUE)
  data.frame(panel[c("site_index", "region", "position", "allele_a", "allele_b")],
    a = counts[, 1L], b = counts[, 2L], other = counts[, 3L],
    stringsAsFactors = FALSE)
}

upstream_counts <- function(somalier, source, fixture, directory) {
  out <- file.path(directory, paste0("somalier-", basename(source)))
  dir.create(out)
  run(somalier, c("extract", "--sites", shQuote(fixture$sites), "--fasta",
    shQuote(fixture$reference), "--out-dir", shQuote(out), shQuote(source)),
    paste("Somalier extraction for", basename(source)))
  digests <- list.files(out, pattern = "[.]somalier$", full.names = TRUE)
  if (length(digests) != 1L) fail("Somalier did not emit exactly one digest")
  read_somalier(digests[[1L]], fixture$panel)
}

sql_quote <- function(value) paste0("'", gsub("'", "''", value, fixed = TRUE), "'")

duckhts_sql <- function(extension, source, index, fixture, output,
                        max_depth = 100000L, worker_count = 1L, copy = TRUE,
                        reference_index = paste0(fixture$reference, ".fai")) {
  panel_values <- paste(sprintf("('test', %d, '%s', %d, '%s', '%s')",
    fixture$panel$site_index, fixture$panel$region, fixture$panel$position,
    fixture$panel$allele_a, fixture$panel$allele_b), collapse = ",")
  panel_file <- if (nzchar(output)) {
    paste0(output, ".panel.parquet")
  } else {
    tempfile(fileext = ".parquet")
  }
  query <- paste0(
    "SELECT site_index, region, position, allele_a, allele_b, a, b, other, status ",
    "FROM duckhts_somalier_bam_counts(", sql_quote(source),
    ", NULL, 'sample-1', ", sql_quote(fixture$reference),
    ", panel_parquet := ", sql_quote(panel_file),
    ", index_path := ", sql_quote(index),
    if (is.null(reference_index)) "" else paste0(
      ", reference_index_path := ", sql_quote(reference_index)
    ),
    ", worker_count := ", worker_count,
    ", max_depth := ", max_depth, ") ORDER BY site_index")
  statement <- if (copy) paste0("COPY (", query, ") TO ", sql_quote(output),
    " (FORMAT CSV, HEADER, DELIMITER '\t');") else query
  paste0(
    "LOAD ", sql_quote(extension), ";",
    "SET threads = ", max(1L, worker_count), ";",
    "CREATE TABLE panel(assembly VARCHAR, site_index UBIGINT, region VARCHAR, ",
    "position UBIGINT, allele_a VARCHAR, allele_b VARCHAR);",
    "INSERT INTO panel VALUES ", panel_values, ";",
    "COPY panel TO ", sql_quote(panel_file), " (FORMAT parquet);", statement, ";")
}

duckhts_counts <- function(duckdb, extension, source, index, fixture, directory,
                           reference_index = paste0(fixture$reference, ".fai"),
                           output_tag = basename(source), worker_count = 1L) {
  output <- file.path(directory, paste0("duckhts-", output_tag, ".tsv"))
  sql <- duckhts_sql(extension, source, index, fixture, output,
    reference_index = reference_index, worker_count = worker_count)
  run(duckdb, c("-unsigned", "-c", shQuote(sql)),
    paste("DuckHTS extraction for", basename(source)))
  observed <- utils::read.delim(output, stringsAsFactors = FALSE,
    check.names = FALSE, na.strings = "")
  if (nrow(observed) != nrow(fixture$panel)) {
    fail("DuckHTS did not return one row for every panel site")
  }
  observed
}

assert_explicit_cram_fai <- function(duckdb, extension, fixture, directory,
                                     expected) {
  default_index <- paste0(fixture$reference, ".fai")
  custom_directory <- file.path(directory, "custom-reference-index")
  dir.create(custom_directory)
  custom_index <- file.path(custom_directory, "reference.custom.fai")
  if (!file.copy(default_index, custom_index)) {
    fail("could not copy the reference FAI to its explicit non-colocated path")
  }
  restored <- FALSE
  on.exit({
    if (!restored && !file.exists(default_index)) {
      file.copy(custom_index, default_index)
    }
  }, add = TRUE)
  if (unlink(default_index) != 0L || file.exists(default_index)) {
    fail("could not remove the default reference FAI before CRAM extraction")
  }

  missing_index_sql <- duckhts_sql(
    extension, fixture$cram, paste0(fixture$cram, ".crai"), fixture,
    output = "", reference_index = NULL, copy = FALSE
  )
  missing_index_output <- suppressWarnings(system2(
    duckdb, c("-unsigned", "-c", shQuote(missing_index_sql)),
    stdout = TRUE, stderr = TRUE
  ))
  missing_index_status <- attr(missing_index_output, "status")
  if (is.null(missing_index_status) || missing_index_status == 0L ||
      !any(grepl("failed to open reference and existing faidx",
                 missing_index_output, fixed = TRUE))) {
    fail("CRAM extraction without an existing inferred FAI did not fail explicitly")
  }
  if (file.exists(default_index)) {
    fail("CRAM extraction created a default reference FAI sidecar")
  }

  observed <- duckhts_counts(duckdb, extension, fixture$cram,
    paste0(fixture$cram, ".crai"), fixture, directory,
    reference_index = custom_index, output_tag = "explicit-fai.cram")
  assert_identical_counts(observed, expected,
    "CRAM extraction with an explicit non-colocated FAI changed counts")
  if (!identical(observed$status, expected$status)) {
    fail("CRAM extraction with an explicit non-colocated FAI changed statuses")
  }
  if (file.exists(default_index)) {
    fail("CRAM extraction recreated the default reference FAI sidecar")
  }
  if (!file.copy(custom_index, default_index)) {
    fail("could not restore the default reference FAI after CRAM extraction")
  }
  restored <- TRUE
}

assert_depth_limit <- function(duckdb, extension, fixture) {
  sql <- duckhts_sql(extension, fixture$bam, paste0(fixture$bam, ".bai"),
    fixture, output = "", max_depth = 3L, copy = FALSE)
  output <- suppressWarnings(system2(duckdb,
    c("-unsigned", "-c", shQuote(sql)), stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  if (is.null(status) || status == 0L ||
      !any(grepl("max_depth", output, fixed = TRUE))) {
    fail("max_depth below same-start passing depth did not fail explicitly")
  }
}

assert_identical_counts <- function(left, right, label) {
  columns <- c("site_index", "region", "position", "allele_a", "allele_b",
    "a", "b", "other")
  if (!identical(unname(as.matrix(left[columns])),
      unname(as.matrix(right[columns])))) {
    print(merge(left[columns], right[columns], by = columns[1:5], all = TRUE,
      suffixes = c(".left", ".right")))
    fail(label)
  }
}

main <- function(args) {
  if (length(args) != 3L) {
    fail("usage: Rscript somalier_extraction_differential.R ",
      "<duckdb-cli> <duckhts-extension> <somalier-v0.3.4>")
  }
  duckdb <- normalizePath(args[[1L]], mustWork = TRUE)
  extension <- normalizePath(args[[2L]], mustWork = TRUE)
  somalier <- normalizePath(args[[3L]], mustWork = TRUE)
  samtools <- Sys.which("samtools")
  if (!nzchar(samtools)) fail("samtools is required")
  if (sha256(somalier) != SOMALIER_LINUX_SHA256) {
    fail("Somalier executable is not the pinned v0.3.4 Linux release asset")
  }
  help <- run(somalier, "--help", "Somalier --help")
  if (!any(grepl("somalier version: 0.3.4", help, fixed = TRUE))) {
    fail("Somalier did not report version 0.3.4")
  }

  directory <- tempfile("somalier-extraction-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  fixture <- write_fixture(directory, samtools)
  oracle_bam <- oracle_counts(read_alignments(fixture$bam,
    fixture$reference, samtools), fixture$panel)
  oracle_cram <- oracle_counts(read_alignments(fixture$cram,
    fixture$reference, samtools), fixture$panel)
  assert_identical_counts(oracle_bam, oracle_cram,
    "independent BAM and CRAM oracles differ")

  upstream_bam <- upstream_counts(somalier, fixture$bam, fixture, directory)
  upstream_cram <- upstream_counts(somalier, fixture$cram, fixture, directory)
  assert_identical_counts(upstream_bam, upstream_cram,
    "pinned Somalier BAM and CRAM extraction differ")

  duckhts_bam <- duckhts_counts(duckdb, extension, fixture$bam,
    paste0(fixture$bam, ".bai"), fixture, directory)
  duckhts_cram <- duckhts_counts(duckdb, extension, fixture$cram,
    paste0(fixture$cram, ".crai"), fixture, directory)
  duckhts_bam_workers <- duckhts_counts(duckdb, extension, fixture$bam,
    paste0(fixture$bam, ".bai"), fixture, directory,
    output_tag = "bam-workers", worker_count = 4L)
  duckhts_cram_workers <- duckhts_counts(duckdb, extension, fixture$cram,
    paste0(fixture$cram, ".crai"), fixture, directory,
    output_tag = "cram-workers", worker_count = 4L)
  assert_identical_counts(duckhts_bam, duckhts_cram,
    "DuckHTS BAM and CRAM extraction differ")
  assert_identical_counts(duckhts_bam, duckhts_bam_workers,
    "DuckHTS one/four-worker BAM extraction differs")
  assert_identical_counts(duckhts_cram, duckhts_cram_workers,
    "DuckHTS one/four-worker CRAM extraction differs")
  if (!identical(duckhts_bam$status, duckhts_cram$status)) {
    fail("DuckHTS BAM and CRAM availability statuses differ")
  }

  available <- fixture$panel$relation != "alignment_position_unavailable"
  equal <- fixture$panel$relation == "equal"
  if (!identical(equal, available)) {
    fail("fixture comparison domains are incomplete")
  }
  assert_identical_counts(upstream_bam, oracle_bam,
    "Somalier release differs from the independent oracle")
  assert_identical_counts(duckhts_bam[available, , drop = FALSE],
    oracle_bam[available, , drop = FALSE],
    "DuckHTS differs from Somalier/oracle inside the supported subset")
  unavailable <- fixture$panel$relation == "alignment_position_unavailable"
  triplet <- function(data, rows) {
    unname(as.integer(data[rows, c("a", "b", "other")]))
  }
  if (sum(unavailable) != 1L ||
      duckhts_bam$status[unavailable] != "unavailable_alignment_position" ||
      duckhts_cram$status[unavailable] != "unavailable_alignment_position" ||
      any(!is.na(duckhts_bam[unavailable, c("a", "b", "other")])) ||
      any(!is.na(duckhts_cram[unavailable, c("a", "b", "other")]))) {
    fail("reference-valid coordinate beyond @SQ length was not unavailable: BAM status=",
      duckhts_bam$status[unavailable], ", CRAM status=",
      duckhts_cram$status[unavailable], ", BAM counts=",
      paste(triplet(duckhts_bam, unavailable), collapse = "/"),
      ", CRAM counts=", paste(triplet(duckhts_cram, unavailable), collapse = "/"))
  }
  assert_explicit_cram_fai(duckdb, extension, fixture, directory,
    duckhts_cram)
  assert_depth_limit(duckdb, extension, fixture)

  cat("Somalier extraction differential: OK (", nrow(fixture$panel),
    " sites x BAM/CRAM; ", sum(equal),
    " exact count matches; 1/4 workers; 1 unavailable alignment-position row; ",
    "explicit CRAM FAI)\n",
    sep = "")
  cat("Pinned sources: Somalier ", SOMALIER_COMMIT, "; hileup ",
    HILEUP_COMMIT, "\n", sep = "")
}

if (sys.nframe() == 0L) main(commandArgs(trailingOnly = TRUE))
