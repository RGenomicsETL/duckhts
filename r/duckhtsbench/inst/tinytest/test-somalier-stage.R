library(tinytest)

test_somalier_stage <- function() {
  previous <- Sys.getenv(c("DUCKHTSBENCH_REGISTRY", "DUCKHTS_CACHE_DIR"), unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(previous[name]))
  })
  registry <- Sys.getenv("DUCKHTSBENCH_REGISTRY", unset = "")
  if (!nzchar(registry)) {
    registry <- system.file("benchmark_registry.tsv", package = "duckhtsbench")
  }
  cache <- tempfile("somalier-stage-")
  dir.create(cache)
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry, DUCKHTS_CACHE_DIR = cache)

  upstream <- duckhts_bench_stage_plan("somalier-upstream")
  expect_equal(upstream$id, c(
    "somalier_v034_source_archive", "somalier_v034_linux_x86_64"
  ))
  expect_true(all(upstream$transform == "direct_download"))
  expect_match(upstream$supplier_identity[[1L]],
    "sha256=acd2dc11be6051c80d15628703a9419965cea3fe7a0563b47e14743e7ac6339e",
    fixed = TRUE)
  expect_match(upstream$supplier_identity[[1L]],
    "commit=ff58fdade8f4f8293d904f10e0a4a13f1fac808d", fixed = TRUE)
  expect_match(upstream$supplier_identity[[2L]],
    "sha256=18717c205a9c4b65d479f1d2cf069a30b047a5378edf544d403d4081f06a3a78",
    fixed = TRUE)

  system <- Sys.info()
  if (identical(unname(system[["sysname"]]), "Linux") &&
      unname(system[["machine"]]) %in% c("x86_64", "amd64")) {
    fake <- file.path(cache, "somalier-fixture")
    writeLines(c("#!/bin/sh", "printf 'somalier version: 0.3.4\\n'"), fake)
    fake_sha256 <- digest::digest(file = fake, algo = "sha256")
    fake_registry <- file.path(cache, "somalier-registry.tsv")
    utils::write.table(data.frame(
      id = "somalier_v034_linux_x86_64", workload = "somalier-upstream",
      role = "pinned_linux_executable", release = "Somalier_v0.3.4_test",
      locator = paste0("file://", normalizePath(fake, winslash = "/")),
      access = "public", cache_relpath = "upstream/somalier/v0.3.4/somalier",
      transform = "direct_download",
      consumer = "test/scripts/somalier_extraction_differential.R",
      stage_order = 1L,
      supplier_identity = paste0(
        "sha256=", fake_sha256, ";bytes=", file.info(fake)$size
      ), stringsAsFactors = FALSE
    ), fake_registry, sep = "\t", row.names = FALSE, quote = FALSE)
    Sys.setenv(DUCKHTSBENCH_REGISTRY = fake_registry,
               DUCKHTS_CACHE_DIR = file.path(cache, "somalier-binary"))
    executable <- duckhtsbench:::duckhts_bench_stage_somalier_v034()
    expect_true(file.exists(executable))
    expect_equal(as.character(file.info(executable)$mode), "755")
    expect_true(file.exists(paste0(executable, ".provenance.tsv")))
    Sys.setenv(DUCKHTSBENCH_REGISTRY = registry, DUCKHTS_CACHE_DIR = cache)
  }

  plan <- duckhts_bench_stage_plan("somalier-synthetic")
  expect_equal(plan$id, c(
    "somalier_synthetic_panel", "somalier_synthetic_frequency",
    "somalier_synthetic_evidence", "somalier_synthetic_pairs"
  ))
  expect_true(all(plan$access == "local_derived"))
  expect_true(all(grepl("^algorithm:|^artifact:", plan$locator)))

  paths <- duckhts_bench_stage_somalier(4L)
  expect_true(all(file.exists(paths)))
  expect_true(all(file.exists(paste0(paths, ".provenance.tsv"))))
  expect_match(paths[["evidence"]], "samples-4/evidence\\.parquet$")
  expect_match(paths[["pairs"]], "samples-4/selected-3/pairs-2/pairs\\.parquet$")
  hashes <- vapply(paths, digest::digest, character(1L), file = TRUE, algo = "sha256")
  mtimes <- file.info(c(paths, paste0(paths, ".provenance.tsv")))$mtime
  expect_identical(duckhts_bench_stage_somalier(4L), paths)
  expect_identical(vapply(paths, digest::digest, character(1L),
    file = TRUE,
    algo = "sha256"
  ), hashes)
  expect_identical(file.info(c(paths, paste0(paths, ".provenance.tsv")))$mtime, mtimes)
  second_cache <- file.path(cache, "independent-regeneration")
  Sys.setenv(DUCKHTS_CACHE_DIR = second_cache)
  second_paths <- duckhts_bench_stage_somalier(4L)
  expect_identical(unname(vapply(second_paths, digest::digest, character(1L),
    file = TRUE,
    algo = "sha256"
  )), unname(hashes))
  collision_cache <- file.path(cache, "samples-250")
  Sys.setenv(DUCKHTS_CACHE_DIR = collision_cache)
  collision_paths <- duckhts_bench_stage_somalier(4L)
  expect_match(collision_paths[["evidence"]],
    "samples-250/benchmarks/somalier/synthetic-v2/samples-4/evidence\\.parquet$")
  expect_match(collision_paths[["pairs"]],
    "samples-250/benchmarks/somalier/synthetic-v2/samples-4/selected-3/pairs-2/pairs\\.parquet$")

  # Registry paths use canonical forward slashes on every platform. Shadowing
  # .Platform in a copied closure reproduces the Windows separator lookup.
  windows_platform <- .Platform
  windows_platform$file.sep <- "\\"
  windows_stage <- duckhts_bench_stage_somalier
  environment(windows_stage) <- list2env(
    list(.Platform = windows_platform),
    parent = environment(windows_stage)
  )
  Sys.setenv(DUCKHTS_CACHE_DIR = file.path(cache, "windows-separator"))
  windows_paths <- windows_stage(4L)
  expect_true(all(file.exists(windows_paths)))
  Sys.setenv(DUCKHTS_CACHE_DIR = cache)

  alignment_plan <- duckhts_bench_stage_plan("somalier-site-extraction")
  expect_equal(alignment_plan$id, c(
    "somalier_site_panel", "somalier_site_reference",
    "somalier_site_reference_fai", "somalier_site_bam",
    "somalier_site_bam_bai", "somalier_site_cram",
    "somalier_site_cram_crai"
  ))
  expect_true(all(alignment_plan$access %in% c("local_generated", "local_derived")))
  expect_true(all(grepl("^algorithm:|^artifact:", alignment_plan$locator)))
  expect_true(all(grepl("benchmark_somalier_site_extraction.Rmd",
                        alignment_plan$consumer, fixed = TRUE)))
  stage_alignment <- duckhtsbench:::duckhts_bench_stage_somalier_site_extraction
  alignment_paths <- stage_alignment()
  expect_true(all(file.exists(alignment_paths)))
  expect_true(all(file.exists(paste0(alignment_paths, ".provenance.tsv"))))
  expect_true(all(file.info(alignment_paths)$size > 0))
  expect_false(any(grepl("\\.sam$", list.files(dirname(alignment_paths[["bam"]])))))
  alignment_mtimes <- file.info(c(
    alignment_paths, paste0(alignment_paths, ".provenance.tsv")
  ))$mtime
  expect_identical(stage_alignment(), alignment_paths)
  expect_identical(file.info(c(
    alignment_paths, paste0(alignment_paths, ".provenance.tsv")
  ))$mtime, alignment_mtimes)

  bam_bytes <- readBin(alignment_paths[["bam"]], "raw",
                       n = file.info(alignment_paths[["bam"]])$size)
  corrupted_bam <- bam_bytes
  last_byte <- length(corrupted_bam)
  corrupted_bam[[last_byte]] <- as.raw(bitwXor(
    as.integer(corrupted_bam[[last_byte]]), 1L
  ))
  writeBin(corrupted_bam, alignment_paths[["bam"]])
  expect_error(
    stage_alignment(), "existing Somalier site-extraction provenance is invalid"
  )
  writeBin(bam_bytes, alignment_paths[["bam"]])
  expect_identical(stage_alignment(), alignment_paths)

  samtools <- Sys.which("samtools")
  Sys.setenv(DUCKHTS_CACHE_DIR = file.path(cache, "alignment-regeneration"))
  regenerated_alignment_paths <- stage_alignment()
  expect_identical(
    readLines(regenerated_alignment_paths[["reference"]], warn = FALSE),
    readLines(alignment_paths[["reference"]], warn = FALSE)
  )
  for (format in c("bam", "cram")) {
    original_arguments <- c("view")
    regenerated_arguments <- c("view")
    if (format == "cram") {
      original_arguments <- c(
        original_arguments, "-T", shQuote(alignment_paths[["reference"]])
      )
      regenerated_arguments <- c(
        regenerated_arguments, "-T",
        shQuote(regenerated_alignment_paths[["reference"]])
      )
    }
    original_records <- system2(
      samtools, c(original_arguments, shQuote(alignment_paths[[format]])),
      stdout = TRUE
    )
    regenerated_records <- system2(
      samtools,
      c(regenerated_arguments, shQuote(regenerated_alignment_paths[[format]])),
      stdout = TRUE
    )
    expect_identical(regenerated_records, original_records)
  }
  Sys.setenv(DUCKHTS_CACHE_DIR = cache)

  for (format in c("bam", "cram")) {
    arguments <- c("view", "-c")
    if (format == "cram") {
      arguments <- c(arguments, "-T", shQuote(alignment_paths[["reference"]]))
    }
    count <- suppressWarnings(system2(
      samtools, c(arguments, shQuote(alignment_paths[[format]])),
      stdout = TRUE, stderr = FALSE
    ))
    expect_equal(as.numeric(count), 17000)
  }
  reference_index <- utils::read.delim(
    alignment_paths[["reference_fai"]], header = FALSE,
    colClasses = "character", check.names = FALSE
  )
  expect_equal(reference_index[[1L]], paste0("chr", seq_len(22L)))
  alignment_receipt <- utils::read.delim(
    paste0(alignment_paths[["bam"]], ".provenance.tsv"),
    colClasses = "character", check.names = FALSE
  )
  alignment_fields <- stats::setNames(alignment_receipt$value, alignment_receipt$field)
  expect_equal(alignment_fields[["panel_sites"]], "17000")
  expect_equal(alignment_fields[["source_reads"]], "17000")
  expect_equal(alignment_fields[["read_length"]], "101")
  expect_equal(alignment_fields[["artifact_sha256"]],
               digest::digest(file = alignment_paths[["bam"]], algo = "sha256"))

  con <- DBI::dbConnect(duckdb::duckdb())
  table <- function(path) sprintf("read_parquet(%s)", as.character(DBI::dbQuoteString(con, path)))
  observed <- DBI::dbGetQuery(con, sprintf(
    "SELECT
  (SELECT count(*) FROM %s)::DOUBLE AS sites,
  (SELECT count(*) FROM %s)::DOUBLE AS frequencies,
  (SELECT count(*) FROM %s)::DOUBLE AS evidence,
  (SELECT count(DISTINCT sample_id) FROM %s)::DOUBLE AS samples,
  (SELECT count(*) FROM %s)::DOUBLE AS pairs,
  (SELECT count(*) FROM %s WHERE a IS NULL AND b IS NULL AND other IS NULL)::DOUBLE
    AS unavailable,
  (SELECT count(*) FROM %s WHERE a IS NOT NULL AND a + b = 30 AND other = 0)::DOUBLE
    AS measured",
    table(paths[["panel"]]), table(paths[["frequency"]]), table(paths[["evidence"]]),
    table(paths[["evidence"]]), table(paths[["pairs"]]), table(paths[["evidence"]]),
    table(paths[["evidence"]])
  ))
  expect_equal(
    unname(as.numeric(observed[1, c(
      "sites", "frequencies", "evidence",
      "samples", "pairs"
    )])),
    c(17000, 17000, 68000, 4, 2)
  )
  expect_equal(observed$unavailable + observed$measured, observed$evidence)
  canonical <- DBI::dbGetQuery(con, sprintf(
    "SELECT count(*) AS sites,
  count(*) FILTER (WHERE allele_a < allele_b) AS canonical FROM %s",
    table(paths[["panel"]])
  ))
  expect_equal(canonical$canonical, canonical$sites)
  alignment_panel <- DBI::dbGetQuery(con, sprintf(
    "SELECT count(*) AS sites, count(DISTINCT (region, position)) AS positions,
  min(site_index) AS first_site, max(site_index) AS last_site FROM %s",
    table(alignment_paths[["panel"]])
  ))
  expect_equal(
    unname(as.numeric(alignment_panel[1, c(
      "sites", "positions", "first_site", "last_site"
    )])),
    c(17000, 17000, 0, 16999)
  )
  panel_differences <- DBI::dbGetQuery(con, sprintf(
    "SELECT count(*) AS n FROM (
  (SELECT * FROM %s EXCEPT ALL SELECT * FROM %s)
  UNION ALL
  (SELECT * FROM %s EXCEPT ALL SELECT * FROM %s))",
    table(alignment_paths[["panel"]]),
    table(regenerated_alignment_paths[["panel"]]),
    table(regenerated_alignment_paths[["panel"]]),
    table(alignment_paths[["panel"]])
  ))$n[[1L]]
  expect_equal(panel_differences, 0)
  topology <- duckhts_bench_stage_somalier(4L,
    selected_samples = 4L,
    selected_pairs = 6L
  )
  expect_match(topology[["pairs"]], "samples-4/selected-4/pairs-6/pairs\\.parquet$")
  topology_rows <- DBI::dbGetQuery(con, sprintf(
    "SELECT
  (SELECT count(*) FROM %s) AS pairs,
  (SELECT count(DISTINCT struct_pack(receiver_id := receiver_id,
    anchor_id := anchor_id)) FROM %s) AS distinct_pairs,
  (SELECT count(DISTINCT sample_id) FROM (
    SELECT receiver_id AS sample_id FROM %s
    UNION ALL SELECT anchor_id AS sample_id FROM %s)) AS selected_samples",
    table(topology[["pairs"]]), table(topology[["pairs"]]),
    table(topology[["pairs"]]), table(topology[["pairs"]])
  ))
  expect_equal(topology_rows$pairs, 6)
  expect_equal(topology_rows$distinct_pairs, 6)
  expect_equal(topology_rows$selected_samples, 4)
  DBI::dbDisconnect(con, shutdown = TRUE)

  receipt <- paste0(paths[["evidence"]], ".provenance.tsv")
  fields <- utils::read.delim(receipt, colClasses = "character", check.names = FALSE)
  expect_equal(fields$value[fields$field == "samples"], "4")
  expect_equal(fields$value[fields$field == "sites"], "17000")
  expect_equal(fields$value[fields$field == "rows"], "68000")
  expect_equal(fields$value[fields$field == "instance_key"], "samples=4")
  expect_true(all(c(
    "generator_sha256", "duckdb_version", "r_version",
    "writer_options"
  ) %in% fields$field))
  fields$value[fields$field == "sha256"] <- paste0(fields$value[fields$field == "sha256"], "0")
  utils::write.table(fields, receipt, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(duckhts_bench_stage_somalier(4L), "receipt does not match bytes")

  expect_error(duckhts_bench_stage_somalier(1L), "from 2 through 1000")
  expect_error(duckhts_bench_stage_somalier(1001L), "from 2 through 1000")
  expect_error(duckhts_bench_stage_somalier(3.5), "from 2 through 1000")
  expect_error(duckhts_bench_stage_somalier(4L, 5L), "selected_samples")
  expect_error(duckhts_bench_stage_somalier(4L, 4L, 2L), "selected_pairs")
  expect_error(duckhts_bench_stage_somalier(4L, 4L, 13L), "selected_pairs")
}

test_somalier_stage()
