#!/usr/bin/env Rscript

# Differential test of duckhts_somalier_sex against Somalier v0.3.4
# `extract` + `relate` on a deterministic synthetic cohort VCF with FORMAT/AD:
# autosomal, X and Y sites, an XX, an XY, an ambiguous X dosage, an XX sample
# with Y reads, a sample with no usable X sites and a sample with zero
# autosomal depth. Integer counts must agree exactly; depth ratios are
# compared with the full-precision values Somalier embeds in its HTML report.
# Stage the third argument with
# duckhtsbench:::duckhts_bench_stage_somalier_v034().

SOMALIER_LINUX_SHA256 <-
  "18717c205a9c4b65d479f1d2cf069a30b047a5378edf544d403d4081f06a3a78"
RATIO_TOLERANCE <- 1e-12

fail <- function(...) stop(..., call. = FALSE)

run <- function(command, args, label) {
  output <- suppressWarnings(system2(command, args, stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    fail(label, " exited ", status, ":\n", paste(output, collapse = "\n"))
  }
  output
}

sql_quote <- function(value) paste0("'", gsub("'", "''", value, fixed = TRUE), "'")

# Sites in coordinate order: 60 autosomal (chr1), 30 X, 12 Y. Every third
# autosomal site has REF > ALT so A/B orientation is exercised.
site_table <- function() {
  autosomal <- data.frame(chrom = "chr1", pos = 1000L + 100L * seq_len(60L),
    kind = "A", stringsAsFactors = FALSE)
  x <- data.frame(chrom = "chrX", pos = 5000000L + 1000L * seq_len(30L),
    kind = "X", stringsAsFactors = FALSE)
  y <- data.frame(chrom = "chrY", pos = 7000000L + 1000L * seq_len(12L),
    kind = "Y", stringsAsFactors = FALSE)
  sites <- rbind(autosomal, x, y)
  flipped <- seq_len(nrow(sites)) %% 3L == 0L
  sites$ref <- ifelse(flipped, "G", "A")
  sites$alt <- ifelse(flipped, "A", "G")
  sites$af <- round(0.05 + 0.9 * ((seq_len(nrow(sites)) * 37L) %% 100L) / 100, 2)
  sites
}

# One sample profile is a function of site kind and index within that kind.
# Genotypes: 0 hom-ref, 1 het, 2 hom-alt; `depth` is the A+B read total.
autosomal_counts <- function(n, depth, rng_seed) {
  set.seed(rng_seed)
  genotype <- sample(0:2, n, replace = TRUE, prob = c(0.45, 0.4, 0.15))
  total <- pmax(0L, depth + sample(-5:5, n, replace = TRUE))
  total[c(7L, 23L)] <- 3L
  alt <- ifelse(genotype == 0L, 0L,
    ifelse(genotype == 2L, total, rbinom(n, total, 0.5)))
  cbind(total - alt, alt)
}

x_counts <- function(profile, depth) {
  n <- 30L
  het <- switch(profile, female = 12L, ambiguous = 3L, male = 0L, none = 0L)
  homalt <- switch(profile, female = 5L, ambiguous = 12L, male = 6L, none = 0L)
  genotype <- c(rep(1L, het), rep(2L, homalt), rep(0L, n - het - homalt))
  genotype <- genotype[order((seq_len(n) * 7L) %% n)]
  total <- rep(depth, n)
  if (profile == "none") total[] <- 4L
  # A hom-alt-like site below the depth gate and a middling-balance site.
  total[3L] <- 5L
  alt <- ifelse(genotype == 0L, 0L, ifelse(genotype == 2L, total, total %/% 2L))
  if (profile != "none") alt[9L] <- as.integer(round(total[9L] * 0.15))
  cbind(total - alt, alt)
}

y_counts <- function(depth) {
  n <- 12L
  total <- rep(depth, n)
  alt <- ifelse(seq_len(n) %% 4L == 0L, total, 0L)
  cbind(total - alt, alt)
}

sample_profiles <- function(kind) {
  if (kind == "main") {
    list(
      list(id = "xx", auto = 30L, x = "female", xd = 30L, yd = 0L),
      list(id = "xy", auto = 30L, x = "male", xd = 15L, yd = 15L),
      list(id = "xy_b", auto = 40L, x = "male", xd = 20L, yd = 20L),
      list(id = "xy_c", auto = 25L, x = "male", xd = 12L, yd = 13L),
      list(id = "no_x", auto = 30L, x = "none", xd = 30L, yd = 0L),
      list(id = "zero_auto", auto = 0L, x = "female", xd = 30L, yd = 0L),
      list(id = "ambiguous", auto = 30L, x = "ambiguous", xd = 28L, yd = 0L),
      list(id = "xx_with_y", auto = 30L, x = "female", xd = 30L, yd = 15L)
    )
  } else {
    c(lapply(sprintf("xx%02d", 1:11), function(id) {
      list(id = id, auto = 30L, x = "female", xd = 30L, yd = 0L)
    }), list(list(id = "xx_with_y", auto = 30L, x = "female", xd = 30L, yd = 15L)))
  }
}

cell <- function(counts) {
  alt <- counts[[2L]]
  ref <- counts[[1L]]
  gt <- if (ref + alt == 0L) "./." else if (alt == 0L) "0/0" else if (ref == 0L) "1/1" else "0/1"
  paste0(gt, ":", ref, ",", alt)
}

write_vcf <- function(path, sites, profiles) {
  columns <- lapply(seq_along(profiles), function(i) {
    profile <- profiles[[i]]
    counts <- list(
      A = autosomal_counts(sum(sites$kind == "A"), profile$auto, 100L + i),
      X = x_counts(profile$x, profile$xd), Y = y_counts(profile$yd))
    index <- c(A = 0L, X = 0L, Y = 0L)
    vapply(seq_len(nrow(sites)), function(row) {
      kind <- sites$kind[[row]]
      index[[kind]] <<- index[[kind]] + 1L
      pair <- counts[[kind]][index[[kind]], ]
      # Somalier and DuckHTS both order counts by REF, ALT in the VCF.
      cell(pair)
    }, character(1))
  })
  filter <- ifelse(sites$kind == "X" & seq_len(nrow(sites)) %% 11L == 0L, "lowq",
    ifelse(sites$kind == "A" & seq_len(nrow(sites)) %% 29L == 0L, "lowq", "PASS"))
  body <- do.call(paste, c(list(sites$chrom, sites$pos, ".", sites$ref, sites$alt,
    ".", filter, paste0("AF=", sites$af), "GT:AD"), columns, sep = "\t"))
  writeLines(c("##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=100000>", "##contig=<ID=chrX,length=156040895>",
    "##contig=<ID=chrY,length=57227415>",
    "##FILTER=<ID=lowq,Description=\"low quality\">",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate allele frequency\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
      "FORMAT", vapply(profiles, `[[`, character(1), "id")), collapse = "\t"),
    body), path)
  run("bgzip", c("-f", shQuote(path)), "bgzip")
  run("tabix", c("-f", "-p", "vcf", shQuote(paste0(path, ".gz"))), "tabix")
  paste0(path, ".gz")
}

# X/Y sites keep their positions as records in the sites file used by both tools.
write_sites <- function(path, sites) {
  writeLines(c("##fileformat=VCFv4.2",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate allele frequency\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    paste(sites$chrom, sites$pos, ".", sites$ref, sites$alt, ".", "PASS",
      paste0("AF=", sites$af), sep = "\t")), path)
  path
}

# Somalier requires --fasta even for VCF extraction; it is not read for counts.
write_reference <- function(path) {
  lengths <- c(chr1 = 100000L, chrX = 5100000L, chrY = 7100000L)
  lines <- unlist(lapply(names(lengths), function(name) {
    c(paste0(">", name), rep(strrep("A", 60L), lengths[[name]] %/% 60L),
      strrep("A", lengths[[name]] %% 60L))
  }))
  writeLines(lines[nzchar(lines)], path)
  run("samtools", c("faidx", shQuote(path)), "samtools faidx")
  path
}

upstream_stats <- function(somalier, sites_path, reference, vcf, directory, tag) {
  out <- file.path(directory, paste0("extract-", tag))
  dir.create(out)
  run(somalier, c("extract", "--sites", shQuote(sites_path), "--fasta",
    shQuote(reference), "--out-dir", shQuote(out), shQuote(vcf)),
    "somalier extract")
  digests <- list.files(out, pattern = "[.]somalier$", full.names = TRUE)
  prefix <- file.path(directory, paste0("relate-", tag))
  relate_output <- run(somalier, c("relate", "--infer", "--output-prefix", shQuote(prefix),
    shQuote(digests)), "somalier relate")
  if (Sys.getenv("SEX_DEBUG") != "") print(relate_output)
  samples <- utils::read.delim(paste0(prefix, ".samples.tsv"),
    check.names = FALSE, stringsAsFactors = FALSE)
  names(samples)[1L] <- "family_id"
  html <- readLines(paste0(prefix, ".html"), warn = FALSE)
  json <- grep("\"x_depth_mean\"", html, value = TRUE)
  if (!length(json)) fail("Somalier HTML report has no sample statistics")
  list(samples = samples, html = paste(html, collapse = "\n"))
}

json_number <- function(html, sample, field) {
  pattern <- sprintf("\"sample\":\"%s\".*?\"%s\":([^,}]+)", sample, field)
  hit <- regmatches(html, regexec(pattern, html, perl = TRUE))[[1L]]
  if (length(hit) != 2L) fail("no ", field, " for ", sample, " in the Somalier report")
  value <- gsub("\"", "", hit[[2L]], fixed = TRUE)
  if (value %in% c("nan", "NaN")) NA_real_ else if (value %in% c("inf", "Infinity")) Inf
  else as.numeric(value)
}

duckhts_sex <- function(duckdb, extension, sites_path, vcf, y_gate, output) {
  sql <- paste0(
    "LOAD ", sql_quote(extension), ";",
    "CREATE TABLE panel AS SELECT * FROM duckhts_somalier_import_sites(",
    sql_quote(sites_path), ", 'synthetic');",
    "CREATE TABLE counts AS SELECT * FROM duckhts_somalier_vcf_counts(",
    sql_quote(vcf), ", 'panel');",
    "COPY (SELECT sample_id, autosomal_genotyped_sites, autosomal_depth_mean, ",
    "x_usable_sites, x_depth_ratio, x_hom_ref, x_het, x_hom_alt, y_usable_sites, ",
    "y_depth_ratio, cohort_has_y, y_signal, inferred_sex, status ",
    "FROM duckhts_somalier_sex('counts', 'panel', y_gate := ", sql_quote(y_gate),
    ") ORDER BY sample_id) TO ", sql_quote(output),
    " (FORMAT CSV, HEADER, DELIMITER '\t');")
  run(duckdb, c("-unsigned", "-c", shQuote(sql)), "DuckHTS sex")
  utils::read.delim(output, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = "")
}

close_enough <- function(observed, expected) {
  if (is.na(expected)) return(is.na(observed))
  if (is.infinite(expected)) return(identical(observed, expected))
  !is.na(observed) && abs(observed - expected) <= RATIO_TOLERANCE * max(1, abs(expected))
}

compare_cohort <- function(duckdb, extension, somalier, sites_path, reference,
                           sites, profiles, directory, tag, y_gate) {
  vcf <- write_vcf(file.path(directory, paste0(tag, ".vcf")), sites, profiles)
  upstream <- upstream_stats(somalier, sites_path, reference, vcf, directory, tag)
  observed <- duckhts_sex(duckdb, extension, sites_path, vcf, y_gate,
    file.path(directory, paste0(tag, "-duckhts.tsv")))
  rows <- list()
  for (id in observed$sample_id) {
    u <- upstream$samples[upstream$samples$sample_id == id, , drop = FALSE]
    o <- observed[observed$sample_id == id, , drop = FALSE]
    integer_pairs <- list(
      c("x_usable_sites", "X_n"), c("x_hom_ref", "X_hom_ref"),
      c("x_het", "X_het"), c("x_hom_alt", "X_hom_alt"),
      c("y_usable_sites", "Y_n"))
    for (pair in integer_pairs) {
      if (!identical(as.integer(o[[pair[[1L]]]]), as.integer(u[[pair[[2L]]]]))) {
        fail(tag, "/", id, ": ", pair[[1L]], " ", o[[pair[[1L]]]], " != upstream ",
          u[[pair[[2L]]]])
      }
    }
    # Somalier reports 2 * 0 / gt_mean = 0 when no site is usable and NaN or
    # infinity when there is no autosomal depth; DuckHTS reports NULL for both.
    for (pair in list(c("x_depth_ratio", "x_depth_mean", "X_n"),
                      c("y_depth_ratio", "y_depth_mean", "Y_n"))) {
      expected <- json_number(upstream$html, id, pair[[2L]])
      if (u[[pair[[3L]]]] == 0L || u$gt_depth_mean == 0) expected <- NA_real_
      if (!close_enough(o[[pair[[1L]]]], expected)) {
        fail(tag, "/", id, ": ", pair[[1L]], " ", o[[pair[[1L]]]],
          " != upstream ", expected)
      }
    }
    # samples.tsv prints depth means at fixed precision; check that
    # rounding, which bounds the tolerance justified for the ratios.
    if (!is.na(o$autosomal_depth_mean) &&
        abs(o$autosomal_depth_mean - u$gt_depth_mean) > 0.05 + 1e-9) {
      fail(tag, "/", id, ": autosomal depth mean differs beyond .1f rounding")
    }
    call <- c("1" = "XY", "2" = "XX", "-2" = "ambiguous")[as.character(u$sex)]
    rows[[id]] <- data.frame(cohort = tag, sample = id,
      x_usable = o$x_usable_sites, x_het = o$x_het, x_hom_alt = o$x_hom_alt,
      y_usable = o$y_usable_sites,
      x_ratio = o$x_depth_ratio, upstream_x_ratio =
        json_number(upstream$html, id, "x_depth_mean"),
      call = o$inferred_sex, upstream_call = if (is.na(call)) "unset" else unname(call),
      status = o$status, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

# The BAM/CRAM path counts every panel site in one pass. A small alignment
# fixture with chr1, chrX and chrY sites must give Somalier's autosomal, X and
# Y counts exactly.
write_bam_fixture <- function(directory) {
  contigs <- c("chr1", "chrX", "chrY")
  reference <- file.path(directory, "bam-reference.fa")
  writeLines(unlist(lapply(contigs, function(name) {
    c(paste0(">", name), rep(strrep("A", 60L), 5L))
  })), reference)
  run("samtools", c("faidx", shQuote(reference)), "samtools faidx")
  sites <- data.frame(chrom = c("chr1", "chr1", "chrX", "chrX", "chrY"),
    pos = c(50L, 150L, 60L, 160L, 70L), ref = "A", alt = "C", af = 0.3,
    stringsAsFactors = FALSE)
  sites_path <- write_sites(file.path(directory, "bam-sites.vcf"), sites)
  read <- function(name, chrom, site, base, flag = 0L) {
    bases <- replace_base(strrep("A", 40L), site - (site - 10L) + 1L, base)
    paste(name, flag, chrom, site - 10L, 60L, "40M", "*", 0L, 0L, bases,
      strrep("I", 40L), sep = "\t")
  }
  bases <- list(chr1 = list("50" = c("A", "A", "C", "G"), "150" = c("C", "C", "C")),
    chrX = list("60" = c("A", "C", "A", "A", "T"), "160" = c("C", "C", "C", "C")),
    chrY = list("70" = c("A", "A", "C")))
  records <- character()
  for (chrom in names(bases)) {
    for (site in names(bases[[chrom]])) {
      for (i in seq_along(bases[[chrom]][[site]])) {
        records <- c(records, read(sprintf("%s-%s-%d", chrom, site, i), chrom,
          as.integer(site), bases[[chrom]][[site]][[i]]))
      }
    }
  }
  header <- c("@HD\tVN:1.6\tSO:coordinate", sprintf("@SQ\tSN:%s\tLN:300", contigs),
    "@RG\tID:rg1\tSM:bam-sample")
  sam <- file.path(directory, "reads.sam")
  writeLines(c(header, records[order(match(sub("-.*", "", records), contigs),
    as.integer(vapply(strsplit(records, "\t"), `[[`, "", 4L)))]), sam)
  bam <- file.path(directory, "reads.bam")
  run("samtools", c("sort", "-o", shQuote(bam), shQuote(sam)), "samtools sort")
  run("samtools", c("index", shQuote(bam)), "samtools index")
  list(reference = reference, sites = sites_path, bam = bam)
}

replace_base <- function(sequence, index1, base) {
  substr(sequence, index1, index1) <- base
  sequence
}

read_digest <- function(path) {
  connection <- file(path, "rb")
  on.exit(close(connection))
  readBin(connection, integer(), 1L, size = 1L, signed = FALSE)
  name_length <- readBin(connection, integer(), 1L, size = 1L, signed = FALSE)
  readBin(connection, raw(), name_length)
  sizes <- readBin(connection, integer(), 3L, size = 2L, signed = FALSE,
    endian = "little")
  values <- readBin(connection, integer(), 3L * sum(sizes), size = 4L,
    endian = "little")
  matrix(values, ncol = 3L, byrow = TRUE)
}

compare_bam <- function(duckdb, extension, somalier, directory) {
  fixture <- write_bam_fixture(directory)
  out <- file.path(directory, "bam-extract")
  dir.create(out)
  run(somalier, c("extract", "--sites", shQuote(fixture$sites), "--fasta",
    shQuote(fixture$reference), "--out-dir", shQuote(out), shQuote(fixture$bam)),
    "somalier extract BAM")
  expected <- read_digest(list.files(out, pattern = "[.]somalier$", full.names = TRUE))
  panel <- file.path(directory, "bam-panel.parquet")
  output <- file.path(directory, "bam-duckhts.tsv")
  sql <- paste0(
    "LOAD ", sql_quote(extension), ";",
    "COPY (SELECT * FROM duckhts_somalier_import_sites(", sql_quote(fixture$sites),
    ", 'synthetic')) TO ", sql_quote(panel), " (FORMAT parquet);",
    "COPY (SELECT site_index, region, a, b, other FROM duckhts_somalier_bam_counts(",
    sql_quote(fixture$bam), ", NULL, 'bam-sample', ", sql_quote(fixture$reference),
    ", panel_parquet := ", sql_quote(panel), ") ORDER BY site_index) TO ",
    sql_quote(output), " (FORMAT CSV, HEADER, DELIMITER '\t');")
  run(duckdb, c("-unsigned", "-c", shQuote(sql)), "DuckHTS BAM counts")
  observed <- utils::read.delim(output, stringsAsFactors = FALSE)
  if (!identical(unname(as.matrix(observed[c("a", "b", "other")])), expected) ||
      !identical(observed$region, c("chr1", "chr1", "chrX", "chrX", "chrY"))) {
    print(observed)
    print(expected)
    fail("BAM counts differ from Somalier on autosomal, X and Y sites")
  }
  nrow(observed)
}

main <- function(args) {
  if (length(args) != 3L) {
    fail("usage: Rscript somalier_sex_differential.R <duckdb-cli> ",
      "<duckhts-extension> <somalier-v0.3.4>")
  }
  duckdb <- normalizePath(args[[1L]], mustWork = TRUE)
  extension <- normalizePath(args[[2L]], mustWork = TRUE)
  somalier <- normalizePath(args[[3L]], mustWork = TRUE)
  for (tool in c("bgzip", "tabix", "samtools")) {
    if (!nzchar(Sys.which(tool))) fail(tool, " is required")
  }
  checksum <- strsplit(run(Sys.which("sha256sum"), shQuote(somalier), "sha256sum"),
    " ")[[1L]][[1L]]
  if (checksum != SOMALIER_LINUX_SHA256) {
    fail("Somalier executable is not the pinned v0.3.4 Linux release asset")
  }
  directory <- tempfile("somalier-sex-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  sites <- site_table()
  sites_path <- write_sites(file.path(directory, "sites.vcf"), sites)

  reference <- write_reference(file.path(directory, "reference.fa"))
  main_rows <- compare_cohort(duckdb, extension, somalier, sites_path, reference, sites,
    sample_profiles("main"), directory, "main", "cohort")
  # Eleven XX samples and one XX sample with Y reads: Somalier's cohort gate
  # (n > 5 or n / N > 0.1 samples with Y depth) is closed, so it keeps XX,
  # whereas the per-sample gate reports the apparent Y.
  gate_rows <- compare_cohort(duckdb, extension, somalier, sites_path, reference, sites,
    sample_profiles("gate"), directory, "gate", "cohort")
  sample_gate <- duckhts_sex(duckdb, extension, sites_path,
    file.path(directory, "gate.vcf.gz"), "sample",
    file.path(directory, "gate-sample.tsv"))

  bam_sites <- compare_bam(duckdb, extension, somalier, directory)

  calls <- rbind(main_rows, gate_rows)
  # Somalier leaves the sex unset (-9) for samples that fail its autosomal
  # quality gate and for X dosage between its two thresholds; DuckHTS reports
  # `ambiguous` for the latter and evaluates the X dosage regardless of the
  # autosomal gate, reporting no_autosomal_depth.
  decided <- calls$upstream_call != "unset"
  no_autosomal <- calls$status == "no_autosomal_depth"
  if (any(calls$call[decided] != calls$upstream_call[decided]) ||
      any(calls$call[!decided & !no_autosomal] != "ambiguous")) {
    print(calls)
    fail("DuckHTS and Somalier sex calls differ")
  }
  if (sum(decided) < 10L) fail("too few Somalier sex calls to compare")
  expected <- c(xx = "XX", xy = "XY", xy_b = "XY", xy_c = "XY", no_x = "ambiguous",
    zero_auto = "XX", ambiguous = "ambiguous", xx_with_y = "ambiguous")
  main_calls <- stats::setNames(main_rows$call, main_rows$sample)
  if (!identical(unname(main_calls[names(expected)]), unname(expected))) {
    print(main_rows)
    fail("fixture calls are not the designed XX, XY and edge-case calls")
  }
  if (gate_rows$call[gate_rows$sample == "xx_with_y"] != "XX" ||
      sample_gate$inferred_sex[sample_gate$sample_id == "xx_with_y"] != "ambiguous") {
    fail("cohort and per-sample Y gates did not diverge as designed")
  }
  print(calls, row.names = FALSE)
  cat("Somalier sex differential passed:", nrow(calls), "sample comparisons,",
    bam_sites, "BAM sites\n")
}

if (!interactive() && Sys.getenv("SEX_LIB") == "") main(commandArgs(trailingOnly = TRUE))
