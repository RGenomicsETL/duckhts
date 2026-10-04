# Network-free test of benchmarks/roh_counts_validation/: a temporary registry
# and cache hold a small VCF and an AF-site Parquet. Staging must keep the
# requested samples in order, GT only, and exactly one record per AF site, and
# must refuse to publish a BCF whose identity differs from the registered one.
# The run comparison used by the report is checked on hand-computed intervals.
test_roh_counts_validation_staging <- function() {
  bcftools <- Sys.which("bcftools")
  if (!nzchar(bcftools)) {
    stop("bcftools is required for the network-free staging test", call. = FALSE)
  }
  source("benchmarks/roh_counts_validation/stage_genotypes.R")
  cache <- tempfile("roh-counts-validation-stage-")
  dir.create(file.path(cache, "fixture"), recursive = TRUE)
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)
  old <- Sys.getenv(c("DUCKHTS_CACHE_DIR", "DUCKHTSBENCH_REGISTRY"), unset = NA)
  on.exit({
    for (name in names(old)) {
      if (is.na(old[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(old[name]))
    }
  }, add = TRUE)

  # Position 20 has two records; only G>T is an AF site. Position 30 is not a
  # site. The INFO tag and the DP field must not reach the output.
  samples <- c("S1", "S2", "S3", "S4", "S5")
  record <- function(pos, ref, alt, genotypes) {
    paste(c("chr20", pos, ".", ref, alt, ".", "PASS", "AF=0.25", "GT:DP",
            paste0(genotypes, ":9")), collapse = "\t")
  }
  vcf_path <- file.path(cache, "fixture", "source.vcf")
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=chr20,length=100>",
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
            samples), collapse = "\t"),
    record(10, "A", "C", c("0|0", "0|1", "1|1", "0|0", "1|0")),
    record(20, "G", "A", c("0|1", "0|0", "0|0", "0|0", "0|0")),
    record(20, "G", "T", c("1|1", "0|0", "0|1", "0|0", "0|0")),
    record(30, "C", "T", c("0|1", "0|1", "0|1", "0|1", "0|1")),
    record(40, "T", "C", c("0|0", "1|1", "0|0", "1|0", "1|1"))), vcf_path)

  write_sites <- function(name, values) {
    path <- file.path(cache, "fixture", name)
    con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
    DBI::dbExecute(con, sprintf(
      "COPY (SELECT * FROM (VALUES %s) AS t(chrom, pos, ref, alt)) TO %s (FORMAT PARQUET)",
      values, DBI::dbQuoteString(con, path)))
    path
  }
  sites <- write_sites("sites.parquet", paste(
    "('chr20', 10::BIGINT, 'A', 'C'), ('chr20', 20::BIGINT, 'G', 'T'),",
    "('chr20', 40::BIGINT, 'T', 'C')"))
  uncovered <- write_sites("uncovered.parquet", paste(
    "('chr20', 10::BIGINT, 'A', 'C'), ('chr20', 50::BIGINT, 'G', 'T')"))

  # The first staging has no expected identity and supplies the one to register.
  first <- stage_roh_validation_genotypes(
    source_vcf = vcf_path, af_sites = sites, samples = c("S5", "S1", "S3"),
    output = file.path(cache, "first.bcf"), bcftools = bcftools)
  query <- function(arguments) system2(bcftools, shQuote(arguments), stdout = TRUE)
  stopifnot(first$records == 3L, first$samples == 3L,
            grepl("^[0-9a-f]{64}$", first$records_sha256),
            file.exists(paste0(first$output, ".csi")),
            identical(query(c("query", "-l", first$output)), c("S5", "S1", "S3")),
            identical(query(c("view", "-H", first$output)),
                      c("chr20\t10\tchr20_10_A_C\tA\tC\t.\tPASS\t.\tGT\t1|0\t0|0\t1|1",
                        "chr20\t20\tchr20_20_G_T\tG\tT\t.\tPASS\t.\tGT\t0|0\t1|1\t0|1",
                        "chr20\t40\tchr20_40_T_C\tT\tC\t.\tPASS\t.\tGT\t1|1\t0|0\t0|0")),
            !length(list.files(cache, pattern = "targets|ids|partial")))

  # A site with no source record is an error that publishes nothing.
  missing_output <- file.path(cache, "uncovered.bcf")
  failed <- tryCatch({
    stage_roh_validation_genotypes(vcf_path, uncovered, c("S1", "S2"), missing_output,
                                   bcftools = bcftools)
    FALSE
  }, error = function(error) grepl("each site needs exactly one record", conditionMessage(error)))
  stopifnot(failed, !file.exists(missing_output))

  header <- c("id", "workload", "role", "release", "locator", "access", "cache_relpath",
              "transform", "consumer", "stage_order", "supplier_identity")
  row <- function(...) paste(c(...), collapse = "\t")
  source_sha256 <- unname(digest::digest(file = vcf_path, algo = "sha256"))
  write_registry <- function(source_identity, identity) {
    path <- file.path(cache, "registry.tsv")
    writeLines(c(paste(header, collapse = "\t"),
      row("fixture_source", "test", "source_vcf", "test", "local", "local_generated",
          "fixture/source.vcf", "copy_committed_fixture", "test", "1", source_identity),
      row("fixture_sites", "test", "af_sites", "test", "local", "local_generated",
          "fixture/sites.parquet", "copy_committed_fixture", "test", "1", ""),
      row("fixture_genotypes", "test", "site_genotypes", "test",
          "artifact:fixture_source;artifact:fixture_sites", "local_derived",
          "genotypes/fixture.bcf", "stage_roh_validation_genotypes_from_registry",
          "test", "2", paste0("sample_list=S5,S1,S3;", identity))), path)
    path
  }
  Sys.setenv(DUCKHTS_CACHE_DIR = cache)
  destination <- file.path(cache, "genotypes", "fixture.bcf")
  refused <- function(pattern) {
    failed <- tryCatch({
      stage_roh_validation_genotypes_from_registry("fixture_genotypes", bcftools)
      FALSE
    }, error = function(error) grepl(pattern, conditionMessage(error)))
    failed && !file.exists(destination) && !file.exists(paste0(destination, ".csi"))
  }
  right_identity <- paste0("records=3;samples=3;records_sha256=", first$records_sha256)

  # A wrong registered output identity is an error that publishes nothing.
  Sys.setenv(DUCKHTSBENCH_REGISTRY = write_registry(
    paste0("sha256=", source_sha256),
    paste0("records=3;samples=3;records_sha256=", strrep("0", 64L))))
  stopifnot(refused("differs from its registered identity"))

  # So is a cached source that differs from its registered SHA-256.
  Sys.setenv(DUCKHTSBENCH_REGISTRY = write_registry(
    paste0("sha256=", strrep("0", 64L)), right_identity))
  stopifnot(refused("checksum differs from the expected source identity"))

  # The registered identities are accepted, and staging reproduces the records.
  Sys.setenv(DUCKHTSBENCH_REGISTRY = write_registry(
    paste0("sha256=", source_sha256), right_identity))
  accepted <- stage_roh_validation_genotypes_from_registry("fixture_genotypes", bcftools)
  stopifnot(file.exists(destination), file.exists(paste0(destination, ".csi")),
            identical(accepted$records_sha256, first$records_sha256))
}

# Hand-computed intervals, one-based and inclusive. Test runs cover 100 + 51
# bases and reference runs 101 + 1,000,000; they share 51 + 1 bases.
test_roh_run_comparison <- function() {
  source("benchmarks/roh_counts_validation/metrics.R")
  test <- data.frame(start = c(1, 950), end = c(100, 1000))
  reference <- data.frame(start = c(50, 1000), end = c(150, 1000999))
  none <- test[0L, ]
  compared <- roh_run_comparison(test, reference, span_bases = 2e6)
  stopifnot(roh_run_bases(test) == 151, roh_run_bases(reference) == 1000101,
            roh_run_overlap_bases(test, reference) == 52,
            roh_run_overlap_bases(reference, test) == 52,
            roh_run_overlap_bases(test, none) == 0,
            compared$shared_bases == 52, compared$either_bases == 1000200,
            isTRUE(all.equal(compared$jaccard, 52 / 1000200)),
            isTRUE(all.equal(compared$froh_difference, (151 - 1000101) / 2e6)),
            is.na(roh_run_comparison(none, none, 2e6)$jaccard),
            roh_run_comparison(test, none, 2e6)$jaccard == 0)
  overlapping <- data.frame(start = c(1, 50), end = c(100, 150))
  failed <- tryCatch({
    roh_run_overlap_bases(overlapping, reference)
    FALSE
  }, error = function(error) grepl("must not overlap", conditionMessage(error)))
  stopifnot(failed)

  # Only the reference has a run of at least 1,000,000 bases, so the long
  # class has a Jaccard of 0 and fails; an empty class agrees.
  declaration <- list(max_froh_difference = 0.01, min_jaccard_long = 0.9,
                      min_jaccard_all = 0.7)
  classes <- roh_run_verdict(
    roh_run_comparison_by_class(test, reference, 2e6, long_run_bases = 1e6), declaration)
  stopifnot(identical(classes$run_class, c("all", "long")),
            identical(classes$test_runs, c(2L, 0L)),
            identical(classes$reference_runs, c(2L, 1L)),
            classes$jaccard[[2L]] == 0, !classes$pass[[1L]], !classes$pass[[2L]])
  empty <- roh_run_verdict(
    roh_run_comparison_by_class(test, test, 2e6, long_run_bases = 1e6), declaration)
  stopifnot(empty$jaccard[[1L]] == 1, is.na(empty$jaccard[[2L]]), all(empty$pass))
}

test_roh_counts_validation_staging()
test_roh_run_comparison()
