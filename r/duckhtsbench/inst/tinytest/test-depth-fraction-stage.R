library(tinytest)

registry <- duckhts_bench_registry()
transform <- "duckhts_bench_stage_depth_fraction"
rows <- registry[registry$workload == "depth-fraction", ]
expect_equal(sort(rows$id), sort(c("depth_fraction_hg00403_buckets", "depth_fraction_na18507_buckets",
                                   "depth_fraction_hg00188_buckets")))
expect_true(all(rows$transform == transform))
expect_true(all(grepl("^artifact:[a-z0-9_]+;artifact:liftover_grch38_fasta$", rows$locator)))
expect_true(all(sub("^artifact:([^;]+);.*$", "\\1", rows$locator) %in% registry$id))
identity_fields <- getFromNamespace("duckhts_bench_identity_fields", "duckhtsbench")
for (identity in rows$supplier_identity) {
  fields <- identity_fields(identity)
  expect_true(all(c("sample", "regions", "bin_width", "buckets", "kept_buckets", "min_mapq",
                    "exclude_flags", "rows", "reads") %in% names(fields)))
}
for (id in c(sub("^artifact:([^;]+);.*$", "\\1", rows$locator), "liftover_grch38_fasta",
             "ancestry_1000g_panel")) {
  expect_true(grepl("benchmark_depth_fraction.Rmd", registry$consumer[registry$id == id], fixed = TRUE))
}

# Network-free staging of the bundled CRAM fixture of Rduckhts through a
# private registry: the counts of all buckets equal the reads that pass the
# filters, read again independently; a subset of buckets is a subset of
# them; a matching receipt is reused; a registered count that differs is an
# error; and a row without its parameters is rejected before any work.
if (requireNamespace("DBI", quietly = TRUE) && requireNamespace("duckdb", quietly = TRUE) &&
    requireNamespace("Rduckhts", quietly = TRUE) &&
    nzchar(system.file("extdata", "bam_mismatch.cram", package = "Rduckhts"))) local({
  root <- tempfile("depth-fraction-stage-")
  dir.create(file.path(root, "cache", "fixture"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE))
  extdata <- system.file("extdata", package = "Rduckhts")
  for (name in c("bam_mismatch.cram", "bam_mismatch.cram.crai", "bam_mismatch.fa", "bam_mismatch.fa.fai")) {
    file.copy(file.path(extdata, name), file.path(root, "cache", "fixture", name))
  }
  parameters <- "sample=FIXTURE;regions=ref1,ref2;bin_width=100;buckets=3;kept_buckets=%d;min_mapq=0;exclude_flags=3844"
  plan_row <- function(id, kept, extra = "") {
    c(id, "depth-fraction", "read_name_bucket_bin_counts", "fixture",
      "artifact:fixture_cram;artifact:fixture_reference", "local_derived",
      file.path("fixture", paste0(id, ".parquet")), transform, "tinytest", "1",
      paste0(sprintf(parameters, kept), extra))
  }
  plan <- rbind(
    c("fixture_cram", "fixture", "input_cram", "fixture", "local", "local", "fixture/bam_mismatch.cram",
      "direct_download", "tinytest", "1", ""),
    c("fixture_reference", "fixture", "reference", "fixture", "local", "local", "fixture/bam_mismatch.fa",
      "direct_download", "tinytest", "1", ""),
    plan_row("all_buckets", 3L), plan_row("one_bucket", 1L))
  colnames(plan) <- names(registry)
  for (input in c("fixture_cram", "fixture_reference")) {
    source <- file.path(root, "cache", plan[plan[, "id"] == input, "cache_relpath"])
    plan[plan[, "id"] == input, "supplier_identity"] <- paste0(
      "sha256=", unname(digest::digest(file = source, algo = "sha256")))
  }
  registry_path <- file.path(root, "registry.tsv")
  write_registry <- function(plan) {
    utils::write.table(plan, registry_path, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  write_registry(plan)
  previous <- Sys.getenv(c("DUCKHTSBENCH_REGISTRY", "DUCKHTS_CACHE_DIR"), unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) {
      Sys.unsetenv(name)
    } else {
      do.call(Sys.setenv, as.list(previous[name]))
    }
  }, add = TRUE)
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry_path, DUCKHTS_CACHE_DIR = file.path(root, "cache"))

  staged <- duckhts_bench_stage_depth_fraction()
  expect_equal(staged$id, c("all_buckets", "one_bucket"))
  expect_equal(staged$reused, c(FALSE, FALSE))
  con <- Rduckhts::rduckhts_connect()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  passing <- DBI::dbGetQuery(con, sprintf(paste(
    "SELECT RNAME AS chrom, (POS - 1) // 100 * 100 AS start, count(*)::DOUBLE AS reads",
    "FROM read_bam(%s, reference := %s) WHERE (FLAG & 3844) = 0 GROUP BY ALL ORDER BY ALL"),
    DBI::dbQuoteString(con, file.path(root, "cache", "fixture", "bam_mismatch.cram")),
    DBI::dbQuoteString(con, file.path(root, "cache", "fixture", "bam_mismatch.fa"))))
  all_buckets <- DBI::dbGetQuery(con, sprintf(paste(
    "SELECT chrom, start, sum(reads)::DOUBLE AS reads, min(bucket) AS low, max(bucket) AS high,",
    "min(sample_id) AS sample FROM read_parquet(%s) GROUP BY chrom, start ORDER BY ALL"),
    DBI::dbQuoteString(con, staged$path[[1L]])))
  expect_equal(all_buckets[c("chrom", "start", "reads")], passing)
  expect_true(all(all_buckets$low >= 0 & all_buckets$high < 3))
  expect_equal(unique(all_buckets$sample), "FIXTURE")
  expect_equal(staged$reads[[1L]], sum(passing$reads))
  one_bucket <- DBI::dbGetQuery(con, sprintf(
    "SELECT count(*) FILTER (WHERE bucket != 0) AS other, coalesce(sum(reads), 0) AS reads FROM read_parquet(%s)",
    DBI::dbQuoteString(con, staged$path[[2L]])))
  expect_equal(one_bucket$other, 0)
  expect_true(one_bucket$reads <= sum(passing$reads))
  receipt <- utils::read.delim(paste0(staged$path[[1L]], ".sources.tsv"), colClasses = "character")
  expect_equal(receipt$field, c("artifact_id", "supplier_identity", "source", "derivation_sha256",
                                "duckdb_version", "extension_sha256", "alignments_identity",
                                "reference_identity", "output_sha256", "rows", "reads"))
  expect_equal(receipt$value[receipt$field == "output_sha256"],
               unname(digest::digest(file = staged$path[[1L]], algo = "sha256")))

  # Recording the counts in the registry keeps the receipt valid.
  plan[plan[, "id"] == "all_buckets", "supplier_identity"] <- paste0(
    sprintf(parameters, 3L), ";rows=", staged$rows[[1L]], ";reads=", staged$reads[[1L]])
  write_registry(plan)
  expect_equal(duckhts_bench_stage_depth_fraction(ids = "all_buckets")$reused, TRUE)
  # A changed registered identity of either input rebuilds the output.
  for (input in c("fixture_cram", "fixture_reference")) {
    plan[plan[, "id"] == input, "supplier_identity"] <- paste0(
      plan[plan[, "id"] == input, "supplier_identity"], ";revision=2")
    write_registry(plan)
    expect_equal(duckhts_bench_stage_depth_fraction(ids = "all_buckets")$reused, FALSE)
    expect_equal(duckhts_bench_stage_depth_fraction(ids = "all_buckets")$reused, TRUE)
  }
  # A forced rebuild must verify the actual cached inputs before recording
  # their registered identities. Corruption leaves the published output
  # and receipt intact, and no partial artifact is published.
  for (input in c("fixture_cram", "fixture_reference")) {
    path <- file.path(root, "cache", plan[plan[, "id"] == input, "cache_relpath"])
    original <- readBin(path, "raw", n = file.size(path))
    output_digest <- unname(digest::digest(file = staged$path[[1L]], algo = "sha256"))
    receipt_path <- paste0(staged$path[[1L]], ".sources.tsv")
    receipt_digest <- unname(digest::digest(file = receipt_path, algo = "sha256"))
    plan[plan[, "id"] == input, "supplier_identity"] <- paste0(
      plan[plan[, "id"] == input, "supplier_identity"], ";corruption_check=1")
    write_registry(plan)
    damaged <- original
    damaged[[1L]] <- as.raw(bitwXor(as.integer(damaged[[1L]]), 1L))
    writeBin(damaged, path)
    expect_error(duckhts_bench_stage_depth_fraction(ids = "all_buckets"), "supplier SHA-256 identity")
    expect_equal(unname(digest::digest(file = staged$path[[1L]], algo = "sha256")), output_digest)
    expect_equal(unname(digest::digest(file = receipt_path, algo = "sha256")), receipt_digest)
    expect_equal(length(list.files(dirname(staged$path[[1L]]), pattern = "[.]partial-")), 0L)
    writeBin(original, path)
    expect_equal(duckhts_bench_stage_depth_fraction(ids = "all_buckets")$reused, FALSE)
  }
  plan[plan[, "id"] == "all_buckets", "supplier_identity"] <- paste0(
    sprintf(parameters, 3L), ";rows=", staged$rows[[1L]], ";reads=", staged$reads[[1L]] + 1)
  write_registry(plan)
  expect_error(duckhts_bench_stage_depth_fraction(ids = "all_buckets"), "differs from its registered identity")
  plan[plan[, "id"] == "one_bucket", "supplier_identity"] <- "sample=FIXTURE;regions=ref1"
  write_registry(plan)
  expect_error(duckhts_bench_stage_depth_fraction(ids = "one_bucket"), "lacks")
  expect_error(duckhts_bench_stage_depth_fraction(ids = "missing"), "ids must name")
})
