library(tinytest)

if (requireNamespace("DBI", quietly = TRUE) &&
    requireNamespace("duckdb", quietly = TRUE)) local({
  root <- tempfile("ancestry-parquet-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  plan <- duckhts_bench_stage_plan("ancestry-reference")
  plan <- plan[plan$id %in% c("ancestry_ref_freqs", "ancestry_projection",
                               "ancestry_reference_parquet"), , drop = FALSE]
  plan$cache_relpath <- basename(plan$cache_relpath)
  keys <- data.frame(chr = c(22L, 1L), pos = c(200L, 100L),
                     rsid = c("rs2", "rs1"), a0 = c("A", "C"), a1 = c("G", "T"))
  sources <- list(ancestry_ref_freqs = cbind(keys, as.data.frame(
    setNames(rep(list(c(0.2, 0.8)), 21L), paste0("group", seq_len(21L))))),
    ancestry_projection = cbind(keys, as.data.frame(
      setNames(rep(list(c(0.5, -0.5)), 16L), paste0("PC", seq_len(16L))))))
  for (id in names(sources)) {
    path <- file.path(root, plan$cache_relpath[plan$id == id])
    con <- gzfile(path, "wt")
    utils::write.csv(sources[[id]], con, row.names = FALSE)
    close(con)
    plan$supplier_identity[plan$id == id] <- paste0(
      "bytes=", file.info(path)$size, ";sha256=",
      digest::digest(file = path, algo = "sha256"))
  }
  registry <- file.path(root, "registry.tsv")
  utils::write.table(plan, registry, sep = "\t", row.names = FALSE, quote = FALSE)
  previous <- Sys.getenv(c("DUCKHTSBENCH_REGISTRY", "DUCKHTS_CACHE_DIR"), unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(previous[name]))
  }, add = TRUE)
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry, DUCKHTS_CACHE_DIR = root)
  path <- duckhts_bench_stage_ancestry_parquet()
  receipt <- utils::read.delim(paste0(path, ".sources.tsv"), colClasses = "character")
  expect_equal(receipt$sha256[receipt$id == "validation"], "frequency_unit_unique_v2")
  expect_equal(receipt$sha256[receipt$id == "derivation"],
               digest::digest(plan$transform[plan$id == "ancestry_reference_parquet"],
                              algo = "sha256", serialize = FALSE))
  expect_identical(duckhts_bench_stage_ancestry_parquet(), path)
  changed <- plan
  changed$transform[changed$id == "ancestry_reference_parquet"] <- paste0(
    changed$transform[changed$id == "ancestry_reference_parquet"], ";revision=other")
  utils::write.table(changed, registry, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(duckhts_bench_stage_ancestry_parquet(), pattern = "identity")
  utils::write.table(plan, registry, sep = "\t", row.names = FALSE, quote = FALSE)
  driver_args <- list(dbdir = ":memory:")
  if ("shared_home" %in% names(formals(duckdb::duckdb))) {
    driver_args$shared_home <- FALSE
  }
  con <- DBI::dbConnect(do.call(duckdb::duckdb, driver_args))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  rows <- DBI::dbGetQuery(con, paste0("SELECT * FROM read_parquet(",
                                     as.character(DBI::dbQuoteString(con, path)), ")"))
  expect_equal(rows$chromosome, c(1L, 22L))
  expect_equal(rows$PC1, c(-0.5, 0.5))
  expect_equal(rows$group1, c(0.8, 0.2))
  expect_equal(ncol(rows), 41L)
  writeLines("corrupt", path)
  expect_error(duckhts_bench_stage_ancestry_parquet(), pattern = "identity")
  unlink(c(path, paste0(path, ".sources.tsv")))
  frequencies <- sources[["ancestry_ref_freqs"]]
  frequencies$group1[1L] <- Inf
  source_path <- file.path(root, plan$cache_relpath[plan$id == "ancestry_ref_freqs"])
  source_connection <- gzfile(source_path, "wt")
  utils::write.csv(frequencies, source_connection, row.names = FALSE)
  close(source_connection)
  plan$supplier_identity[plan$id == "ancestry_ref_freqs"] <- paste0(
    "bytes=", file.info(source_path)$size, ";sha256=",
    digest::digest(file = source_path, algo = "sha256"))
  utils::write.table(plan, registry, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(duckhts_bench_stage_ancestry_parquet(), pattern = "complete keyed row per locus")
  expect_false(file.exists(path))
  expect_false(file.exists(paste0(path, ".sources.tsv")))
  for (value in c(-0.01, 1.01)) {
    frequencies <- sources[["ancestry_ref_freqs"]]
    frequencies$group1[1L] <- value
    source_connection <- gzfile(source_path, "wt")
    utils::write.csv(frequencies, source_connection, row.names = FALSE)
    close(source_connection)
    plan$supplier_identity[plan$id == "ancestry_ref_freqs"] <- paste0(
      "bytes=", file.info(source_path)$size, ";sha256=",
      digest::digest(file = source_path, algo = "sha256"))
    utils::write.table(plan, registry, sep = "\t", row.names = FALSE, quote = FALSE)
    expect_error(duckhts_bench_stage_ancestry_parquet(),
                 pattern = "complete keyed row per locus")
  }
  frequencies <- sources[["ancestry_ref_freqs"]]
  frequencies[1L, names(keys)] <- frequencies[2L, names(keys)]
  source_path <- file.path(root, plan$cache_relpath[plan$id == "ancestry_ref_freqs"])
  source_connection <- gzfile(source_path, "wt")
  utils::write.csv(frequencies, source_connection, row.names = FALSE)
  close(source_connection)
  plan$supplier_identity[plan$id == "ancestry_ref_freqs"] <- paste0(
    "bytes=", file.info(source_path)$size, ";sha256=",
    digest::digest(file = source_path, algo = "sha256"))
  utils::write.table(plan, registry, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(duckhts_bench_stage_ancestry_parquet(), pattern = "complete keyed row per locus")
  for (id in names(sources)) {
    duplicate <- sources[[id]][2L, , drop = FALSE]
    duplicate$a1 <- "A"
    staged_source <- rbind(sources[[id]], duplicate)
    source_path <- file.path(root, plan$cache_relpath[plan$id == id])
    source_connection <- gzfile(source_path, "wt")
    utils::write.csv(staged_source, source_connection, row.names = FALSE)
    close(source_connection)
    plan$supplier_identity[plan$id == id] <- paste0(
      "bytes=", file.info(source_path)$size, ";sha256=",
      digest::digest(file = source_path, algo = "sha256"))
  }
  utils::write.table(plan, registry, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(duckhts_bench_stage_ancestry_parquet(), pattern = "one complete keyed row per locus")
})

local({
  plan <- duckhts_bench_stage_plan("ancestry-projection")
  expect_equal(plan$id, "ancestry_correction")
  expect_equal(plan$transform, "copy_committed_fixture")
  root <- tempfile("ancestry-stage-")
  dir.create(file.path(root, "test/data"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE))
  source <- file.path(root, "test/data/correction.tsv")
  writeLines(c("pc\tcoefficient", "1\t1"), source)
  plan$locator <- "repo:test/data/correction.tsv"
  plan$supplier_identity <- paste0("bytes=", file.info(source)$size,
                                   ";md5=", unname(tools::md5sum(source)))
  registry <- file.path(root, "registry.tsv")
  utils::write.table(plan, registry, sep = "\t", row.names = FALSE, quote = FALSE)
  previous <- Sys.getenv(c("DUCKHTSBENCH_REGISTRY", "DUCKHTS_CACHE_DIR"), unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(previous[name]))
  }, add = TRUE)
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry, DUCKHTS_CACHE_DIR = file.path(root, "cache"))
  staged <- duckhts_bench_stage_repository_fixtures(root, "ancestry-projection")
  expect_equal(unname(tools::md5sum(staged)), unname(tools::md5sum(source)))
  expect_equal(duckhts_bench_stage_repository_fixtures(root, "ancestry-projection"), staged)
  writeLines("corrupted", staged)
  expect_error(duckhts_bench_stage_repository_fixtures(root, "ancestry-projection"),
               pattern = "identity does not match")
})

local({
  plan <- rbind(duckhts_bench_stage_plan("ancestry-reference"),
                duckhts_bench_stage_plan("ancestry-1000g-chr22"),
                duckhts_bench_stage_plan("ancestry-30x-cram"))
  expect_equal(plan$id, c("ancestry_ref_freqs", "ancestry_projection",
                          "ancestry_epilepsy", "ancestry_reference_parquet",
                          "ancestry_1000g_chr22", "ancestry_1000g_panel",
                          "ancestry_30x_na18507", "ancestry_30x_hg00403"))
  direct <- plan$transform == "direct_download"
  expect_true(all(grepl("^https://", plan$locator[direct])))
  expect_equal(plan$transform[!direct], "duckhts_bench_stage_ancestry_parquet")
  root <- tempfile("ancestry-reference-stage-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  previous <- Sys.getenv(c("DUCKHTSBENCH_REGISTRY", "DUCKHTS_CACHE_DIR"),
                         unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) {
      Sys.unsetenv(name)
    } else {
      do.call(Sys.setenv, as.list(previous[name]))
    }
  }, add = TRUE)
  plan$cache_relpath <- file.path("benchmarks/ancestry-reference", basename(plan$cache_relpath))
  plan$supplier_identity <- ""
  for (i in which(direct)) {
    path <- file.path(root, plan$cache_relpath[[i]])
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(paste("reference", i), path)
    plan$supplier_identity[[i]] <- paste0(
      "bytes=", file.info(path)$size, ";sha256=",
      digest::digest(file = path, algo = "sha256")
    )
  }
  registry <- file.path(root, "registry.tsv")
  utils::write.table(plan, registry, sep = "\t", row.names = FALSE, quote = FALSE)
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry, DUCKHTS_CACHE_DIR = root)
  for (id in plan$id[direct]) {
    path <- duckhts_bench_artifact_path(id)
    expect_identical(duckhts_bench_fetch(id), path)
    expect_true(file.exists(paste0(path, ".provenance.tsv")))
    expect_equal(duckhts_bench_validate_identity(id, path), TRUE)
  }
  writeLines("changed", duckhts_bench_artifact_path(plan$id[[1L]]))
  expect_error(duckhts_bench_validate_identity(plan$id[[1L]]),
               pattern = "identity does not match")
})
