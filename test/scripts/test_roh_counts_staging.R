# Network-free test of benchmarks/roh_counts_stage.R: a temporary registry and
# cache hold a small BAM fixture, its reference and an AF-site Parquet. Staging
# must orient the counts and clamp af, and must refuse to publish counts whose
# identity differs from the registered one.
test_roh_counts_staging <- function() {
  extension <- normalizePath(Sys.getenv("DUCKHTS_EXTENSION",
                                        "build/release/duckhts.duckdb_extension"), mustWork = TRUE)
  source("benchmarks/roh_counts_stage.R")
  cache <- tempfile("roh-counts-stage-")
  dir.create(file.path(cache, "fixture"), recursive = TRUE)
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)
  old <- Sys.getenv(c("DUCKHTS_CACHE_DIR", "DUCKHTSBENCH_REGISTRY"), unset = NA)
  on.exit({
    for (name in names(old)) {
      if (is.na(old[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(old[name]))
    }
  }, add = TRUE)
  for (file in c("range.bam", "range.bam.bai", "ce.fa", "ce.fa.fai")) {
    stopifnot(file.copy(file.path("test/data", file), file.path(cache, "fixture", file)))
  }

  # Reference bases: CHROMOSOME_I:1 is covered by no read; 914, 917 and 919
  # are covered by one read with bases A, T and G. At 919 the AF allele is the
  # reference base, so the read counts as alt.
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  DBI::dbExecute(con, sprintf(paste(
    "COPY (SELECT * FROM (VALUES",
    "('CHROMOSOME_I', 1::BIGINT, 'A', 'G', [0.2]::FLOAT[]),",
    "('CHROMOSOME_I', 914::BIGINT, 'A', 'C', [0.5]::FLOAT[]),",
    "('CHROMOSOME_I', 917::BIGINT, 'T', 'C', [0.0005]::FLOAT[]),",
    "('CHROMOSOME_I', 919::BIGINT, 'C', 'G', [0.9995]::FLOAT[]))",
    "AS t(chrom, pos, ref, alt, INFO_AF)) TO %s (FORMAT PARQUET)"),
    DBI::dbQuoteString(con, file.path(cache, "fixture", "sites.parquet"))))
  DBI::dbDisconnect(con, shutdown = TRUE)

  header <- c("id", "workload", "role", "release", "locator", "access", "cache_relpath",
              "transform", "consumer", "stage_order", "supplier_identity")
  row <- function(...) paste(c(...), collapse = "\t")
  settings <- "sample=fixture;assembly=WBcel235;min_mapq=1;min_baseq=0;exclude_flags=1796;overlap_policy=hileup_v0.1.0"
  write_registry <- function(identity) {
    path <- file.path(cache, "registry.tsv")
    writeLines(c(paste(header, collapse = "\t"),
      row("fixture_bam", "test", "input_bam", "test", "local", "local_generated",
          "fixture/range.bam", "copy_committed_fixture", "test", "1", ""),
      row("fixture_sites", "test", "af_sites", "test", "local", "local_generated",
          "fixture/sites.parquet", "copy_committed_fixture", "test", "1", ""),
      row("fixture_reference", "test", "reference", "test", "local", "local_generated",
          "fixture/ce.fa", "copy_committed_fixture", "test", "1", ""),
      row("fixture_counts", "test", "site_read_counts", "test",
          "artifact:fixture_bam;artifact:fixture_sites;artifact:fixture_reference",
          "local_derived", "counts/fixture.counts.parquet", "stage_roh_counts_from_registry",
          "test", "2", paste0(settings, ";", identity))), path)
    path
  }
  Sys.setenv(DUCKHTS_CACHE_DIR = cache)

  # The first staging has no expected identity and supplies the one to register.
  first <- stage_roh_counts(
    source = file.path(cache, "fixture", "range.bam"), sample_id = "fixture",
    af_sites = file.path(cache, "fixture", "sites.parquet"),
    reference = file.path(cache, "fixture", "ce.fa"),
    output = file.path(cache, "first.parquet"), extension = extension,
    assembly = "WBcel235", min_mapq = 1L, min_baseq = 0L, exclude_flags = 1796L,
    overlap_policy = "hileup_v0.1.0", worker_count = 1L)
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  counts <- DBI::dbGetQuery(con, sprintf(
    "SELECT pos, ref_count, alt_count, other_count, af, status FROM read_parquet(%s) ORDER BY pos",
    DBI::dbQuoteString(con, first$output)))
  DBI::dbDisconnect(con, shutdown = TRUE)
  stopifnot(first$rows == 4, grepl("^[0-9a-f]{64}$", first$counts_sha256),
            identical(counts$pos, c(1, 914, 917, 919)),
            identical(counts$ref_count, c(0L, 1L, 1L, 0L)),
            identical(counts$alt_count, c(0L, 0L, 0L, 1L)),
            identical(counts$other_count, c(0L, 0L, 0L, 0L)),
            isTRUE(all.equal(counts$af, c(0.2, 0.5, 1e-3, 0.999), tolerance = 1e-6)),
            all(counts$status == "measured"))

  # A wrong registered identity is an error that publishes nothing.
  Sys.setenv(DUCKHTSBENCH_REGISTRY = write_registry(
    paste0("rows=4;counts_sha256=", strrep("0", 64L))))
  destination <- file.path(cache, "counts", "fixture.counts.parquet")
  failed <- tryCatch({
    stage_roh_counts_from_registry("fixture_counts", extension, worker_count = 1L)
    FALSE
  }, error = function(error) grepl("differ from their registered identity", conditionMessage(error)))
  stopifnot(failed, !file.exists(destination))

  # The registered identity is accepted, and staging reproduces it.
  Sys.setenv(DUCKHTSBENCH_REGISTRY = write_registry(
    paste0("rows=4;counts_sha256=", first$counts_sha256)))
  accepted <- stage_roh_counts_from_registry("fixture_counts", extension, worker_count = 1L)
  stopifnot(file.exists(destination), identical(accepted$counts_sha256, first$counts_sha256),
            accepted$rows == 4)
}

test_roh_counts_staging()
