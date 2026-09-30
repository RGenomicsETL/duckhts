library(tinytest)

# The registered workload: one pinned NCBI archive and its gunzip derivation.
plan <- duckhts_bench_stage_plan("genbank-reader")
expect_equal(plan$id, c("genbank_ecoli_k12_gbff_gz", "genbank_ecoli_k12_gbff"))
expect_equal(plan$transform, c("direct_download", "gunzip"))
expect_equal(plan$locator[[2L]], "artifact:genbank_ecoli_k12_gbff_gz")

# The record-count scaling workload: four pinned RefSeq release parts, each
# with NCBI's published MD5, and three derived record files joining 1, 2 and 4
# of them in order.
plasmid_plan <- duckhts_bench_stage_plan("genbank-plasmid")
part_ids <- paste0("genbank_plasmid_part", 1:4, "_gbff_gz")
expect_equal(plasmid_plan$id, c(part_ids, "genbank_plasmid_records_1",
                                "genbank_plasmid_records_2", "genbank_plasmid_records_4"))
expect_equal(plasmid_plan$transform, c(rep("direct_download", 4L), rep("gunzip_concatenate", 3L)))
expect_true(all(grepl("^https://ftp\\.ncbi\\.nlm\\.nih\\.gov/refseq/release/plasmid/plasmid\\.[1-4]\\.genomic\\.gbff\\.gz$",
                      plasmid_plan$locator[1:4])))
expect_equal(plasmid_plan$locator[5:7], c(
  "artifact:genbank_plasmid_part1_gbff_gz",
  paste(paste0("artifact:", part_ids[1:2]), collapse = ";"),
  paste(paste0("artifact:", part_ids[1:4]), collapse = ";")
))
expect_true(all(grepl("benchmark_genbank_named_attributes.Rmd", plasmid_plan$consumer, fixed = TRUE)))
part_identity <- lapply(plasmid_plan$supplier_identity[1:4], duckhtsbench:::duckhts_bench_identity_fields)
expect_true(all(vapply(part_identity, function(x) all(c("md5", "sha256", "bytes") %in% names(x)), logical(1L))))
records_identity <- lapply(plasmid_plan$supplier_identity[5:7], duckhtsbench:::duckhts_bench_identity_fields)
expect_true(all(vapply(records_identity, function(x) all(c("sha256", "bytes", "records") %in% names(x)), logical(1L))))
expect_match(plan$locator[[1L]], "^https://ftp\\.ncbi\\.nlm\\.nih\\.gov/genomes/all/GCF/000/005/845/GCF_000005845\\.2_ASM584v2/")
expect_true(all(grepl("benchmark_genbank_reader.Rmd", plan$consumer, fixed = TRUE)))
expect_true(all(grepl("benchmark_genbank_memory.Rmd", plan$consumer, fixed = TRUE)))
source_identity <- duckhtsbench:::duckhts_bench_identity_fields(plan$supplier_identity[[1L]])
expect_true(all(c("md5", "sha256", "bytes") %in% names(source_identity)))
derived_identity <- duckhtsbench:::duckhts_bench_identity_fields(plan$supplier_identity[[2L]])
expect_true(all(c("sha256", "bytes", "bp") %in% names(derived_identity)))
expect_equal(duckhts_bench_artifact_path("genbank_ecoli_k12_gbff"),
             sub("\\.gz$", "", duckhts_bench_artifact_path("genbank_ecoli_k12_gbff_gz")))

memory_plan <- duckhts_bench_stage_plan("genbank-memory-scaling")
expect_equal(memory_plan$id, c(
  "genbank_memory_phix174", "genbank_memory_phix174_x2", "genbank_memory_lambda"
))
expect_true(all(memory_plan$access == "repository_fixture"))
expect_true(all(memory_plan$transform == "copy_committed_fixture"))
expect_true(all(grepl("^repo:test/data/", memory_plan$locator)))
expect_true(all(memory_plan$consumer == "benchmark_genbank_memory.Rmd"))
expect_equal(
  unname(vapply(memory_plan$supplier_identity, function(x) {
    duckhtsbench:::duckhts_bench_identity_fields(x)[["records"]]
  }, character(1L))),
  c("1", "2", "1")
)

# Network-free derivation against a synthetic archive under a private registry.
test_genbank_derivation <- function() {
  previous <- Sys.getenv(c("DUCKHTSBENCH_REGISTRY", "DUCKHTS_CACHE_DIR"), unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(previous[name]))
  })
  directory <- tempfile("genbank-stage-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)

  plain <- file.path(directory, "record.gbff")
  writeLines(c("LOCUS       PROBE 4 bp DNA linear PHG 01-JAN-2000", "ORIGIN", "        1 acgt", "//"), plain)
  archive <- file.path(directory, "record.gbff.gz")
  handle <- gzfile(archive, open = "wb")
  writeBin(readBin(plain, what = "raw", n = file.info(plain)$size), handle)
  close(handle)

  registry <- duckhts_bench_stage_plan("genbank-reader")
  registry$supplier_identity <- c(
    paste0("bytes=", file.info(archive)$size, ";md5=", unname(tools::md5sum(archive))),
    paste0("bytes=", file.info(plain)$size, ";md5=", unname(tools::md5sum(plain)))
  )
  registry_path <- file.path(directory, "registry.tsv")
  utils::write.table(registry, registry_path, sep = "\t", row.names = FALSE, quote = FALSE)
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry_path, DUCKHTS_CACHE_DIR = file.path(directory, "cache"))

  # Nothing cached yet: rendering must not download.
  expect_error(duckhts_bench_stage_genbank(fetch = FALSE), "not staged")

  source_path <- duckhts_bench_artifact_path("genbank_ecoli_k12_gbff_gz")
  dir.create(dirname(source_path), recursive = TRUE)
  stopifnot(file.copy(archive, source_path))
  paths <- duckhts_bench_stage_genbank(fetch = FALSE)
  expect_equal(names(paths), c("gbff_gz", "gbff"))
  expect_equal(unname(tools::md5sum(paths[["gbff"]])), unname(tools::md5sum(plain)))
  expect_true(file.exists(paste0(paths[["gbff"]], ".provenance.tsv")))
  expect_false(any(grepl("partial", list.files(dirname(paths[["gbff"]])))))

  # A poisoned derived file is rebuilt from the verified source.
  writeLines("poisoned", paths[["gbff"]])
  paths <- duckhts_bench_stage_genbank(fetch = FALSE)
  expect_equal(unname(tools::md5sum(paths[["gbff"]])), unname(tools::md5sum(plain)))

  # A source that no longer matches its registered identity is refused.
  writeLines("not the archive", source_path)
  expect_error(duckhts_bench_stage_genbank(fetch = FALSE), "identity does not match")
}

test_genbank_derivation()


# Network-free record-count staging: four synthetic gzipped parts under a
# private registry, joined 1, 2 and 4 at a time, each derivation verified
# against its registered identity and rebuilt when poisoned.
test_genbank_plasmid_derivation <- function() {
  previous <- Sys.getenv(c("DUCKHTSBENCH_REGISTRY", "DUCKHTS_CACHE_DIR"), unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(previous[name]))
  })
  directory <- tempfile("genbank-plasmid-stage-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)

  record <- function(name) c(sprintf("LOCUS       %s 4 bp DNA linear PHG 01-JAN-2000", name),
                             "ORIGIN", "        1 acgt", "//")
  parts <- vapply(1:4, function(k) {
    plain <- file.path(directory, sprintf("part%d.gbff", k))
    writeLines(record(sprintf("PART%d", k)), plain)
    archive <- file.path(directory, sprintf("part%d.gbff.gz", k))
    handle <- gzfile(archive, open = "wb")
    writeBin(readBin(plain, what = "raw", n = file.info(plain)$size), handle)
    close(handle)
    archive
  }, character(1L))
  joined <- function(n) {
    path <- file.path(directory, sprintf("records_%d.gbff", n))
    writeLines(unlist(lapply(seq_len(n), function(k) record(sprintf("PART%d", k)))), path)
    path
  }
  expected <- list(joined(1L), joined(2L), joined(4L))

  registry <- duckhts_bench_stage_plan("genbank-plasmid")
  registry$supplier_identity <- c(
    vapply(parts, function(p) paste0("bytes=", file.info(p)$size, ";md5=", unname(tools::md5sum(p))), character(1L)),
    vapply(expected, function(p) paste0("bytes=", file.info(p)$size, ";md5=", unname(tools::md5sum(p)),
                                        ";records=", length(grep("^LOCUS", readLines(p)))), character(1L))
  )
  registry_path <- file.path(directory, "registry.tsv")
  utils::write.table(registry, registry_path, sep = "\t", row.names = FALSE, quote = FALSE)
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry_path, DUCKHTS_CACHE_DIR = file.path(directory, "cache"))

  # Nothing cached yet: rendering must not download.
  expect_error(duckhts_bench_stage_genbank_plasmid(fetch = FALSE), "not staged")
  for (k in 1:4) {
    path <- duckhts_bench_artifact_path(sprintf("genbank_plasmid_part%d_gbff_gz", k))
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(parts[[k]], path))
  }
  staged <- duckhts_bench_stage_genbank_plasmid(fetch = FALSE)
  expect_equal(names(staged), c("parts", "records_1", "records_2", "records_4"))
  for (k in 1:3) {
    path <- staged[[c("records_1", "records_2", "records_4")[[k]]]]
    expect_equal(readLines(path), readLines(expected[[k]]))
    expect_true(file.exists(paste0(path, ".provenance.tsv")))
  }
  expect_false(any(grepl("partial", list.files(dirname(staged$records_4)))))

  # A poisoned derived file is rebuilt from the verified parts; a mismatched
  # source list is refused before anything is written.
  writeLines("poisoned", staged$records_4)
  staged <- duckhts_bench_stage_genbank_plasmid(fetch = FALSE)
  expect_equal(readLines(staged$records_4), readLines(expected[[3L]]))
  expect_error(duckhtsbench::duckhts_bench_stage_gunzip_concatenate(
    "genbank_plasmid_records_2", staged$parts[1:3], staged$records_2), "names 2 sources")

  # A part that no longer matches its registered identity is refused.
  writeLines("not the archive", duckhts_bench_artifact_path("genbank_plasmid_part3_gbff_gz"))
  expect_error(duckhts_bench_stage_genbank_plasmid(fetch = FALSE), "identity does not match")
}

test_genbank_plasmid_derivation()
