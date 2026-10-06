# Network-free test of benchmarks/bam_mismatch_stage.R: a temporary registry
# and cache hold the committed fixtures of duckhts_bam_mismatch_counts. Staging
# must cut the region and keep the records without genotypes, and must refuse
# to publish a file whose identity differs from the registered one.
test_bam_mismatch_staging <- function() {
  # The direct stagers take the path as given, which is relative by default:
  # the slice stager changes its working directory and must still load it.
  extension_as_given <- Sys.getenv("DUCKHTS_EXTENSION", "build/release/duckhts.duckdb_extension")
  extension <- normalizePath(extension_as_given, mustWork = TRUE)
  source("benchmarks/bam_mismatch_stage.R")
  cache <- tempfile("bam-mismatch-stage-")
  dir.create(file.path(cache, "fixture"), recursive = TRUE)
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)
  old <- Sys.getenv(c("DUCKHTS_CACHE_DIR", "DUCKHTSBENCH_REGISTRY"), unset = NA)
  on.exit({
    for (name in names(old)) {
      if (is.na(old[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(old[name]))
    }
  }, add = TRUE)
  fixtures <- c("bam_mismatch.bam", "bam_mismatch.bam.bai", "bam_mismatch.fa", "bam_mismatch.fa.fai",
                "bam_mismatch.mask.vcf.gz", "bam_mismatch.mask.vcf.gz.tbi")
  stopifnot(all(file.copy(file.path("test/data", fixtures), file.path(cache, "fixture", fixtures))))
  fixture <- function(name) file.path(cache, "fixture", name)

  # Region ref1:1-100 overlaps three alignments: the pair at 11 and 41 (20
  # bases each) and clip_ins at 71 (25 stored bases). The fixture mask has
  # three records, at 173, 180 and 186.
  slice <- stage_bam_mismatch_slice(fixture("bam_mismatch.bam"), "ref1:1-100", fixture("bam_mismatch.fa"),
                                    file.path(cache, "first.cram"), extension_as_given)
  mask <- stage_bam_mismatch_mask(fixture("bam_mismatch.mask.vcf.gz"), "ref1",
                                  file.path(cache, "first.bcf"), extension_as_given)
  stopifnot(slice$alignments == 3, slice$position_sum == 11 + 41 + 71, slice$base_sum == 20 + 20 + 25,
            file.exists(paste0(slice$output, ".crai")),
            mask$records == 3, mask$position_sum == 173 + 180 + 186,
            file.exists(paste0(mask$output, ".csi")))

  header <- c("id", "workload", "role", "release", "locator", "access", "cache_relpath",
              "transform", "consumer", "stage_order", "supplier_identity")
  row <- function(...) paste(c(...), collapse = "\t")
  write_registry <- function(slice_identity, mask_identity) {
    path <- file.path(cache, "registry.tsv")
    writeLines(c(paste(header, collapse = "\t"),
      row("fixture_bam", "test", "input_bam", "test", "local", "local_generated",
          "fixture/bam_mismatch.bam", "copy_committed_fixture", "test", "1", ""),
      row("fixture_reference", "test", "reference", "test", "local", "local_generated",
          "fixture/bam_mismatch.fa", "copy_committed_fixture", "test", "1", ""),
      row("fixture_variants", "test", "variants", "test", "local", "local_generated",
          "fixture/bam_mismatch.mask.vcf.gz", "copy_committed_fixture", "test", "1", ""),
      row("fixture_slice", "test", "alignment_slice", "test", "artifact:fixture_bam;artifact:fixture_reference",
          "local_derived", "staged/slice.cram", "stage_bam_mismatch_slice_from_registry", "test", "2",
          paste0("region=ref1:1-100;", slice_identity)),
      row("fixture_mask", "test", "known_variant_mask", "test", "artifact:fixture_variants",
          "local_derived", "staged/mask.bcf", "stage_bam_mismatch_mask_from_registry", "test", "2",
          paste0("contig=ref1;", mask_identity)),
      # Sources without a cached copy: two with a URL locator, one without.
      row("url_bam", "test", "input_bam", "test", paste0("file://", fixture("bam_mismatch.bam")), "public",
          "absent/alignments.bam", "direct_download", "test", "1", ""),
      row("url_variants", "test", "variants", "test", paste0("file://", fixture("bam_mismatch.mask.vcf.gz")),
          "public", "absent/variants.vcf.gz", "direct_download", "test", "1", ""),
      row("absent_variants", "test", "variants", "test", "local", "local_generated",
          "absent/other.vcf.gz", "copy_committed_fixture", "test", "1", ""),
      row("url_slice", "test", "alignment_slice", "test", "artifact:url_bam;artifact:fixture_reference",
          "local_derived", "staged/url_slice.cram", "stage_bam_mismatch_slice_from_registry", "test", "2",
          paste0("region=ref1:1-100;", slice_identity)),
      row("url_mask", "test", "known_variant_mask", "test", "artifact:url_variants",
          "local_derived", "staged/url_mask.bcf", "stage_bam_mismatch_mask_from_registry", "test", "2",
          paste0("contig=ref1;", mask_identity)),
      row("absent_mask", "test", "known_variant_mask", "test", "artifact:absent_variants",
          "local_derived", "staged/absent_mask.bcf", "stage_bam_mismatch_mask_from_registry", "test", "2",
          paste0("contig=ref1;", mask_identity)),
      # A thinned mask derived from the staged mask: the records at even positions.
      row("half_mask", "test", "known_variant_mask", "test", "artifact:fixture_mask",
          "local_derived", "staged/half_mask.bcf", "stage_bam_mismatch_mask_from_registry", "test", "3",
          "contig=ref1;modulo=2;records=2;position_sum=366")), path)
    path
  }
  Sys.setenv(DUCKHTS_CACHE_DIR = cache)
  slice_identity <- "alignments=3;position_sum=123;base_sum=65"
  mask_identity <- "records=3;position_sum=539"

  # A wrong registered identity is an error that publishes nothing.
  Sys.setenv(DUCKHTSBENCH_REGISTRY = write_registry("alignments=4;position_sum=123;base_sum=65",
                                                    "records=3;position_sum=540"))
  for (staging in list(function() stage_bam_mismatch_slice_from_registry("fixture_slice", extension),
                       function() stage_bam_mismatch_mask_from_registry("fixture_mask", extension))) {
    failure <- tryCatch({ staging(); NULL }, error = function(error) conditionMessage(error))
    stopifnot(is.character(failure), grepl("differs from its registered identity", failure, fixed = TRUE))
  }
  stopifnot(!file.exists(file.path(cache, "staged", "slice.cram")),
            !file.exists(file.path(cache, "staged", "mask.bcf")),
            length(list.files(file.path(cache, "staged"))) == 0L)

  # The right identity publishes both files with their indexes, and the staged
  # files then pass the check that a benchmark runs before it measures.
  Sys.setenv(DUCKHTSBENCH_REGISTRY = write_registry(slice_identity, mask_identity))
  staged_slice <- stage_bam_mismatch_slice_from_registry("fixture_slice", extension)
  staged_mask <- stage_bam_mismatch_mask_from_registry("fixture_mask", extension)
  stopifnot(file.exists(staged_slice$output), file.exists(paste0(staged_slice$output, ".crai")),
            file.exists(staged_mask$output), file.exists(paste0(staged_mask$output, ".csi")))
  bam_mismatch_check_staged("fixture_slice", extension)
  bam_mismatch_check_staged("fixture_mask", extension)

  # A source without a cached copy is read where its URL locator says, and the
  # staged files have the same identities. Without a URL it is an error.
  stage_bam_mismatch_slice_from_registry("url_slice", extension)
  stage_bam_mismatch_mask_from_registry("url_mask", extension)
  bam_mismatch_check_staged("url_slice", extension)
  bam_mismatch_check_staged("url_mask", extension)
  failure <- tryCatch({ stage_bam_mismatch_mask_from_registry("absent_mask", extension); NULL },
                      error = function(error) conditionMessage(error))
  stopifnot(is.character(failure), grepl("is not cached and its locator is not a URL", failure, fixed = TRUE),
            !file.exists(file.path(cache, "staged", "absent_mask.bcf")))

  # The thinned mask keeps the records at 180 and 186, not the one at 173.
  half_mask <- stage_bam_mismatch_mask_from_registry("half_mask", extension)
  stopifnot(half_mask$records == 2, half_mask$position_sum == 180 + 186)
  bam_mismatch_check_staged("half_mask", extension)

  # The staged files give the counts of the three alignments: 20 + 20 bases for
  # the pair and 3 for clip_ins with the default flank, with two mismatches.
  con <- bam_mismatch_connect(extension)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  counts <- DBI::dbGetQuery(con, sprintf(
    "SELECT sum(bases)::INTEGER AS bases, sum(bases) FILTER (WHERE reference_base != read_base)::INTEGER AS mismatches
     FROM duckhts_bam_mismatch_counts(%s, %s, mask := %s)",
    DBI::dbQuoteString(con, staged_slice$output), DBI::dbQuoteString(con, fixture("bam_mismatch.fa")),
    DBI::dbQuoteString(con, staged_mask$output)))
  stopifnot(counts$bases == 43L, counts$mismatches == 2L)
  cat("bam mismatch staging: OK\n")
}

test_bam_mismatch_staging()
