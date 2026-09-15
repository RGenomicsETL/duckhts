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
  expect_equal(upstream$id, "somalier_v034_source_archive")
  expect_equal(upstream$transform, "direct_download")
  expect_match(upstream$supplier_identity,
    "sha256=acd2dc11be6051c80d15628703a9419965cea3fe7a0563b47e14743e7ac6339e",
    fixed = TRUE)
  expect_match(upstream$supplier_identity,
    "commit=ff58fdade8f4f8293d904f10e0a4a13f1fac808d", fixed = TRUE)

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
