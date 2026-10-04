# Network-free test of benchmarks/roh_synthetic_stage.R: a temporary registry and
# cache hold a small phased VCF. Staging must plant every segment where the plan
# says, copy haplotype 1 onto haplotype 2 inside it, keep everything else, and
# refuse to publish children whose identity differs from the registered one.
test_roh_synthetic_staging <- function() {
  bcftools <- Sys.which("bcftools")
  if (!nzchar(bcftools)) {
    stop("bcftools is required for the network-free staging test", call. = FALSE)
  }
  source("benchmarks/roh_synthetic_stage.R")
  cache <- tempfile("roh-synthetic-stage-")
  dir.create(file.path(cache, "fixture"), recursive = TRUE)
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)
  old <- Sys.getenv(c("DUCKHTS_CACHE_DIR", "DUCKHTSBENCH_REGISTRY"), unset = NA)
  on.exit({
    for (name in names(old)) {
      if (is.na(old[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(old[name]))
    }
  }, add = TRUE)

  # Records every 5 bases from 98005 to 103000, so a 100-base window holds 20
  # records. The windows of 100001-100400 keep one record each and are closed to
  # planting. Position 100000 is a record: R prints that number as 1e+05 unless
  # it is formatted.
  positions <- seq.int(98005L, 103000L, by = 5L)
  positions <- positions[positions <= 100000L | positions > 100400L | positions %% 100L == 50L]
  samples <- sprintf("child_%d", 1:6)
  set.seed(1L)
  alleles <- function() sample(c("0", "1"), length(positions) * length(samples), replace = TRUE)
  source_genotypes <- matrix(paste0(alleles(), "|", alleles()), nrow = length(positions),
                             dimnames = list(NULL, samples))
  write_vcf <- function(path, genotypes) {
    records <- paste("chrT", positions, paste0("site", positions), "A", "C", ".", "PASS",
                     sprintf("AF=%.3f", (positions - 98000L) / 10000), "GT",
                     apply(genotypes, 1L, paste, collapse = "\t"), sep = "\t")
    writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=chrT,length=110000>",
      '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
      paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
              samples), collapse = "\t"), records), path)
  }
  source_path <- file.path(cache, "fixture", "source.vcf")
  write_vcf(source_path, source_genotypes)

  lengths_bp <- c(600, 300)
  min_gap_bp <- 100
  min_window_records <- 5
  window_bp <- 100
  seed <- 7

  # The plan is a function of the seed; its segments keep the declared rules.
  plan <- roh_synthetic_plan(samples, positions, seed, lengths_bp, min_gap_bp,
                             min_window_records, window_bp)
  again <- roh_synthetic_plan(samples, positions, seed, lengths_bp, min_gap_bp,
                              min_window_records, window_bp)
  other <- roh_synthetic_plan(samples, positions, seed + 1, lengths_bp, min_gap_bp,
                              min_window_records, window_bp)
  stopifnot(identical(plan, again), !identical(plan$start, other$start),
            nrow(plan) == length(samples) * length(lengths_bp),
            identical(plan$length_bp, rep(lengths_bp, length(samples))),
            all(plan$end - plan$start + 1 == plan$length_bp),
            all(plan$start >= min(positions)), all(plan$end <= max(positions)),
            all(plan$end < 100001 | plan$start > 100400))
  for (id in samples) {
    one <- plan[plan$sample_id == id, , drop = FALSE]
    one <- one[order(one$start), , drop = FALSE]
    stopifnot(all(one$start[-1L] - one$end[-nrow(one)] > min_gap_bp))
  }

  # The first staging has no expected identity and supplies the one to register.
  first <- stage_roh_synthetic(
    source_bcf = source_path, output_bcf = file.path(cache, "first.bcf"),
    output_truth = file.path(cache, "first.truth.csv"), seed = seed,
    lengths_bp = lengths_bp, min_gap_bp = min_gap_bp,
    min_window_records = min_window_records, window_bp = window_bp,
    bcftools = bcftools, chunk_records = 64L)
  query <- function(path, arguments) {
    system2(bcftools, shQuote(c(arguments, path)), stdout = TRUE)
  }
  truth <- utils::read.csv(first$output_truth, stringsAsFactors = FALSE)
  ordered <- plan[order(match(plan$sample_id, samples), plan$start), , drop = FALSE]
  stopifnot(identical(truth$sample_id, ordered$sample_id),
            all(truth$segment == ordered$segment), all(truth$start == ordered$start),
            all(truth$end == ordered$end), all(truth$length_bp == ordered$length_bp),
            all(truth$chrom == "chrT"))

  # The expected genotypes are built here from the truth, apart from the staging code.
  expected_genotypes <- source_genotypes
  for (index in seq_len(nrow(truth))) {
    inside <- positions >= truth$start[[index]] & positions <= truth$end[[index]]
    first_allele <- substr(source_genotypes[inside, truth$sample_id[[index]]], 1L, 1L)
    expected_genotypes[inside, truth$sample_id[[index]]] <- paste0(first_allele, "|", first_allele)
    heterozygous <- source_genotypes[inside, truth$sample_id[[index]]] %in% c("0|1", "1|0")
    stopifnot(truth$records[[index]] == sum(inside),
              truth$source_heterozygous[[index]] == sum(heterozygous))
  }
  staged_genotypes <- do.call(rbind, strsplit(
    query(first$output_bcf, c("query", "-f", "[%GT\\t]\\n")), "\t", fixed = TRUE))
  stopifnot(first$records == length(positions), first$samples == length(samples),
            first$truth_rows == nrow(plan),
            grepl("^[0-9a-f]{64}$", first$records_sha256),
            identical(query(first$output_bcf, c("query", "-l")), samples),
            identical(unname(staged_genotypes), unname(expected_genotypes)),
            sum(staged_genotypes != source_genotypes) == sum(truth$source_heterozygous),
            identical(query(first$output_bcf, c("view", "-H", "-G")),
                      query(source_path, c("view", "-H", "-G"))),
            file.exists(paste0(first$output_bcf, ".csi")))

  header <- c("id", "workload", "role", "release", "locator", "access", "cache_relpath",
              "transform", "consumer", "stage_order", "supplier_identity")
  row <- function(...) paste(c(...), collapse = "\t")
  settings <- sprintf("seed=%d;lengths_bp=%s;min_gap_bp=%d;min_window_records=%d;window_bp=%d",
                      seed, paste(lengths_bp, collapse = ","), min_gap_bp,
                      min_window_records, window_bp)
  write_registry <- function(records_sha256, truth_sha256) {
    path <- file.path(cache, "registry.tsv")
    writeLines(c(paste(header, collapse = "\t"),
      row("fixture_source", "test", "phased_source", "test", "local", "local_generated",
          "fixture/source.vcf", "written_by_test", "test", "1", ""),
      row("fixture_synthetic_bcf", "test", "synthetic_children_bcf", "test",
          "artifact:fixture_source", "local_derived", "synthetic/children.bcf",
          "stage_roh_synthetic_from_registry", "test", "2",
          sprintf("%s;records=%d;samples=%d;records_sha256=%s", settings,
                  length(positions), length(samples), records_sha256)),
      row("fixture_synthetic_truth", "test", "planted_truth", "test",
          "artifact:fixture_source", "local_derived", "synthetic/truth.csv",
          "stage_roh_synthetic_from_registry", "test", "2",
          sprintf("rows=%d;sha256=%s", nrow(plan), truth_sha256))), path)
    path
  }
  Sys.setenv(DUCKHTS_CACHE_DIR = cache)
  destinations <- file.path(cache, "synthetic", c("children.bcf", "children.bcf.csi", "truth.csv"))
  stage <- function() {
    stage_roh_synthetic_from_registry("fixture_synthetic_bcf", "fixture_synthetic_truth",
                                      bcftools = bcftools)
  }
  rejected <- function() {
    tryCatch({
      stage()
      FALSE
    }, error = function(error) grepl("differ from their registered identity",
                                     conditionMessage(error)))
  }

  # A wrong record digest or a wrong truth digest is an error that publishes nothing.
  Sys.setenv(DUCKHTSBENCH_REGISTRY = write_registry(strrep("0", 64L), first$truth_sha256))
  stopifnot(rejected(), !any(file.exists(destinations)))
  Sys.setenv(DUCKHTSBENCH_REGISTRY = write_registry(first$records_sha256, strrep("0", 64L)))
  stopifnot(rejected(), !any(file.exists(destinations)))

  # The registered identity is accepted, and staging reproduces it.
  Sys.setenv(DUCKHTSBENCH_REGISTRY = write_registry(first$records_sha256, first$truth_sha256))
  accepted <- stage()
  stopifnot(all(file.exists(destinations)),
            identical(accepted$records_sha256, first$records_sha256),
            identical(accepted$truth_sha256, first$truth_sha256))

  # An unphased genotype cannot be planted: there is no haplotype 1 to copy.
  unphased <- source_genotypes
  unphased[3L, 2L] <- "0/1"
  unphased_path <- file.path(cache, "fixture", "unphased.vcf")
  write_vcf(unphased_path, unphased)
  failed <- tryCatch({
    stage_roh_synthetic(
      source_bcf = unphased_path, output_bcf = file.path(cache, "unphased.bcf"),
      output_truth = file.path(cache, "unphased.truth.csv"), seed = seed,
      lengths_bp = lengths_bp, min_gap_bp = min_gap_bp,
      min_window_records = min_window_records, window_bp = window_bp, bcftools = bcftools)
    FALSE
  }, error = function(error) grepl("must be phased", conditionMessage(error)))
  stopifnot(failed, !file.exists(file.path(cache, "unphased.bcf")),
            !file.exists(file.path(cache, "unphased.truth.csv")))
}

test_roh_synthetic_staging()
