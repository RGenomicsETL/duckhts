library(tinytest)

# The registered workload: a pinned reference archive, its gunzip derivation, a
# pinned ENA read archive, and the BAM minimap2 and samtools derive from them.
plan <- duckhts_bench_stage_plan("ont-ecoli-k12")
expect_equal(plan$id, c("ont_ecoli_k12_reference_fna_gz", "ont_ecoli_k12_reference_fna",
                        "ont_ecoli_k12_reads_fastq_gz", "ont_ecoli_k12_bam"))
expect_equal(plan$transform, c("direct_download", "gunzip", "direct_download",
                               "minimap2_map_ont;samtools_sort;samtools_index"))
expect_match(plan$locator[[1L]], "^https://ftp\\.ncbi\\.nlm\\.nih\\.gov/genomes/all/GCF/000/005/845/GCF_000005845\\.2_ASM584v2/")
expect_equal(plan$locator[[2L]], "artifact:ont_ecoli_k12_reference_fna_gz")
expect_match(plan$locator[[3L]], "^https://ftp\\.sra\\.ebi\\.ac\\.uk/vol1/fastq/ERR146/055/ERR14686255/")
expect_equal(plan$locator[[4L]], "artifact:ont_ecoli_k12_reference_fna;artifact:ont_ecoli_k12_reads_fastq_gz")
expect_true(all(grepl("benchmark_cigar_aligned_blocks.Rmd", plan$consumer, fixed = TRUE)))
expect_true(all(grepl("benchmark_cigar_validation.Rmd", plan$consumer, fixed = TRUE)))
reference_identity <- duckhtsbench:::duckhts_bench_identity_fields(plan$supplier_identity[[1L]])
expect_true(all(c("md5", "bytes") %in% names(reference_identity)))
derived_identity <- duckhtsbench:::duckhts_bench_identity_fields(plan$supplier_identity[[2L]])
expect_true(all(c("sha256", "bytes", "bp") %in% names(derived_identity)))
reads_identity <- duckhtsbench:::duckhts_bench_identity_fields(plan$supplier_identity[[3L]])
expect_true(all(c("md5", "bytes", "reads") %in% names(reads_identity)))
expect_equal(duckhts_bench_artifact_path("ont_ecoli_k12_reference_fna"),
             sub("\\.gz$", "", duckhts_bench_artifact_path("ont_ecoli_k12_reference_fna_gz")))

# Network-free derivation against synthetic sources under a private registry.
test_ont_ecoli_derivation <- function() {
  previous <- Sys.getenv(c("DUCKHTSBENCH_REGISTRY", "DUCKHTS_CACHE_DIR", "PATH"),
                         unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(previous[name]))
  })
  directory <- tempfile("ont-ecoli-stage-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)

  path_samtools <- unname(Sys.which("samtools"))
  if (nzchar(path_samtools)) {
    expect_identical(duckhtsbench:::duckhts_bench_samtools(), path_samtools)
  }
  bundled_samtools <- ""
  if (requireNamespace("RBCFTools", quietly = TRUE) &&
      package_version(getNamespaceVersion("RBCFTools")) >= "1.24-1.1.0") {
    bundled_samtools <- RBCFTools::samtools_path()
  }
  Sys.setenv(PATH = directory)
  expect_identical(duckhtsbench:::duckhts_bench_samtools(), bundled_samtools)
  Sys.setenv(PATH = previous[["PATH"]])

  gzip_copy <- function(source, destination) {
    handle <- gzfile(destination, open = "wb")
    writeBin(readBin(source, what = "raw", n = file.info(source)$size), handle)
    close(handle)
    destination
  }
  set.seed(4)
  genome <- paste(sample(c("A", "C", "G", "T"), 4000, replace = TRUE), collapse = "")
  reference <- file.path(directory, "reference.fna")
  writeLines(c(">probe", substring(genome, seq(1, 4000, 80), pmin(seq(80, 4000, 80), 4000))), reference)
  reference_gz <- gzip_copy(reference, file.path(directory, "reference.fna.gz"))
  starts <- c(1, 401, 1201, 2001, 3001)
  reads <- file.path(directory, "reads.fastq")
  writeLines(unlist(lapply(seq_along(starts), function(i) {
    sequence <- substr(genome, starts[[i]], starts[[i]] + 899)
    c(paste0("@read", i), sequence, "+", strrep("I", nchar(sequence)))
  })), reads)
  reads_gz <- gzip_copy(reads, file.path(directory, "reads.fastq.gz"))

  registry <- duckhts_bench_stage_plan("ont-ecoli-k12")
  registry$supplier_identity <- c(
    paste0("bytes=", file.info(reference_gz)$size, ";md5=", unname(tools::md5sum(reference_gz))),
    paste0("bytes=", file.info(reference)$size, ";md5=", unname(tools::md5sum(reference))),
    paste0("bytes=", file.info(reads_gz)$size, ";md5=", unname(tools::md5sum(reads_gz))),
    "aligner=minimap2;preset=map-ont;sort=coordinate"
  )
  registry_path <- file.path(directory, "registry.tsv")
  utils::write.table(registry, registry_path, sep = "\t", row.names = FALSE, quote = FALSE)
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry_path, DUCKHTS_CACHE_DIR = file.path(directory, "cache"))

  samtools <- duckhtsbench:::duckhts_bench_samtools()
  minimap2 <- unname(Sys.which("minimap2"))
  expect_error(duckhts_bench_stage_ont_ecoli(fetch = FALSE, samtools = ""),
               "samtools is required")
  if (!nzchar(samtools)) return(invisible(NULL))

  # Nothing cached yet: rendering must not download.
  expect_error(duckhts_bench_stage_ont_ecoli(fetch = FALSE, samtools = samtools), "not staged")

  for (id in c("ont_ecoli_k12_reference_fna_gz", "ont_ecoli_k12_reads_fastq_gz")) {
    path <- duckhts_bench_artifact_path(id)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(if (id == "ont_ecoli_k12_reference_fna_gz") reference_gz else reads_gz, path))
  }
  expect_error(duckhts_bench_stage_ont_ecoli(fetch = FALSE, samtools = samtools, minimap2 = ""),
               "minimap2 is required")

  test_bam_receipt_identity <- function() {
    ids <- c("ont_ecoli_k12_reference_fna_gz", "ont_ecoli_k12_reference_fna",
             "ont_ecoli_k12_reads_fastq_gz", "ont_ecoli_k12_bam")
    paths <- stats::setNames(vapply(ids, duckhts_bench_artifact_path, character(1L)),
                             c("reference_gz", "reference", "reads", "bam"))
    bam <- paths[["bam"]]
    index <- paste0(bam, ".bai")
    receipt <- paste0(bam, ".provenance.tsv")
    dir.create(dirname(bam), recursive = TRUE, showWarnings = FALSE)

    fasta <- readLines(paths[["reference"]], warn = FALSE)
    contig_line <- which(startsWith(fasta, ">"))[[1L]]
    contig <- sub("^>", "", fasta[[contig_line]])
    reference_length <- sum(nchar(fasta[!startsWith(fasta, ">")]))
    fastq_connection <- gzfile(paths[["reads"]])
    fastq <- readLines(fastq_connection, warn = FALSE)
    close(fastq_connection)
    read_starts <- seq.int(1L, length(fastq), by = 4L)
    alignments <- vapply(seq_along(read_starts), function(i) {
      line <- read_starts[[i]]
      paste(sub("^@", "", fastq[[line]]), "0", contig,
            as.character((i - 1L) * 100L + 1L), "60",
            paste0(nchar(fastq[[line + 1L]]), "M"), "*", "0", "0",
            fastq[[line + 1L]], fastq[[line + 3L]], sep = "\t")
    }, character(1L))
    sam <- file.path(directory, "receipt-fixture.sam")
    writeLines(c("@HD\tVN:1.6\tSO:coordinate",
                 paste0("@SQ\tSN:", contig, "\tLN:", reference_length),
                 alignments), sam)

    rebuild_minimap2 <- minimap2
    if (!nzchar(rebuild_minimap2) && .Platform$OS.type == "unix") {
      rebuild_minimap2 <- file.path(directory, "minimap2-mock")
      mock <- c(
        "#!/bin/sh",
        "if [ \"$1\" = \"--version\" ]; then",
        "  printf '%s\\n' 'minimap2 receipt test mock'",
        "  exit 0",
        "fi",
        "output=",
        "while [ \"$#\" -gt 0 ]; do",
        "  if [ \"$1\" = \"-o\" ]; then",
        "    shift",
        "    output=$1",
        "  fi",
        "  shift",
        "done",
        paste0("cp ", shQuote(sam, type = "sh"), " \"$output\"")
      )
      writeLines(mock, rebuild_minimap2)
      Sys.chmod(rebuild_minimap2, mode = "0755")
    }
    if (!nzchar(rebuild_minimap2)) return(invisible(NULL))

    duckhts_bench_stage_ont_ecoli(fetch = FALSE, threads = 1L,
                                 samtools = samtools, minimap2 = rebuild_minimap2)
    expect_identical(duckhts_bench_stage_ont_ecoli(fetch = FALSE, samtools = samtools,
                                                   minimap2 = ""),
                     paths[c("reference", "reads", "bam")])
    expect_receipt_identity <- function() {
      fields <- utils::read.delim(receipt, colClasses = "character")
      expect_true(all(c("bam_sha256", "bam_bytes") %in% fields$field))
      expect_identical(fields$value[match(c("bam_sha256", "bam_bytes"), fields$field)],
                       c(digest::digest(file = bam, algo = "sha256"),
                         as.character(file.info(bam)$size)))
    }
    expect_receipt_identity()
    original_bam <- readBin(bam, what = "raw", n = file.info(bam)$size)
    original_eof <- tail(original_bam, 28L)

    header_path <- file.path(directory, "changed-header.sam")
    header <- system2(samtools, c("view", "-H", shQuote(bam)), stdout = TRUE)
    writeLines(c(header, "@CO\tchanged BAM content"), header_path)
    changed_bam <- file.path(directory, "changed-header.bam")
    expect_equal(system2(samtools, c("reheader", shQuote(header_path), shQuote(bam)),
                         stdout = changed_bam), 0L)
    changed_bytes <- readBin(changed_bam, what = "raw", n = file.info(changed_bam)$size)
    expect_identical(tail(changed_bytes, 28L), original_eof)
    expect_false(identical(digest::digest(file = changed_bam, algo = "sha256"),
                           digest::digest(file = bam, algo = "sha256")))
    writeBin(changed_bytes, bam)
    expect_equal(system2(samtools, c("quickcheck", shQuote(bam))), 0L)
    expect_error(duckhts_bench_stage_ont_ecoli(fetch = FALSE, samtools = samtools,
                                               minimap2 = ""), "minimap2 is required")
    duckhts_bench_stage_ont_ecoli(fetch = FALSE, threads = 1L,
                                 samtools = samtools, minimap2 = rebuild_minimap2)
    expect_receipt_identity()

    fields <- utils::read.delim(receipt, colClasses = "character")
    legacy <- fields[!fields$field %in% c("bam_sha256", "bam_bytes"), ]
    utils::write.table(legacy, receipt, sep = "\t", row.names = FALSE, quote = FALSE)
    expect_error(duckhts_bench_stage_ont_ecoli(fetch = FALSE, samtools = samtools,
                                               minimap2 = ""), "minimap2 is required")
    duckhts_bench_stage_ont_ecoli(fetch = FALSE, threads = 1L,
                                 samtools = samtools, minimap2 = rebuild_minimap2)
    expect_receipt_identity()
    unlink(c(bam, index, receipt), force = TRUE)
  }

  test_bam_receipt_identity()
  if (!nzchar(minimap2)) return(invisible(NULL))

  paths <- duckhts_bench_stage_ont_ecoli(fetch = FALSE, threads = 1L,
                                       samtools = samtools, minimap2 = minimap2)
  expect_equal(names(paths), c("reference", "reads", "bam"))
  expect_equal(unname(tools::md5sum(paths[["reference"]])), unname(tools::md5sum(reference)))
  expect_true(file.exists(paste0(paths[["bam"]], ".bai")))
  receipt <- paste0(paths[["bam"]], ".provenance.tsv")
  expect_true(file.exists(receipt))
  fields <- utils::read.delim(receipt, colClasses = "character")
  expect_true(all(c("aligner", "aligner_preset", "sorter", "index_sha256", "index_bytes",
                    "bam_sha256", "bam_bytes") %in% fields$field))
  expect_equal(fields$value[fields$field == "aligner_preset"], "map-ont")
  index <- paste0(paths[["bam"]], ".bai")
  expect_identical(fields$value[fields$field == "index_sha256"],
                   digest::digest(file = index, algo = "sha256"))
  expect_identical(fields$value[fields$field == "index_bytes"],
                   as.character(file.info(index)$size))
  expect_identical(fields$value[fields$field == "bam_sha256"],
                   digest::digest(file = paths[["bam"]], algo = "sha256"))
  expect_identical(fields$value[fields$field == "bam_bytes"],
                   as.character(file.info(paths[["bam"]])$size))
  aligned <- system2(samtools, c("view", "-c", "-F", "4", shQuote(paths[["bam"]])), stdout = TRUE)
  expect_equal(as.integer(aligned), length(starts))
  expect_false(any(grepl("partial", list.files(dirname(paths[["bam"]])))))
  identities <- c("reference_sha256", "reads_sha256")
  hashes <- vapply(paths[c("reference", "reads")], function(path) {
    digest::digest(file = path, algo = "sha256")
  }, character(1L))
  expect_identical(fields$value[match(identities, fields$field)], unname(hashes))
  expect_identical(duckhts_bench_stage_ont_ecoli(fetch = FALSE, samtools = samtools,
                                               minimap2 = ""), paths)

  # Receipts without source identities require an explicit derivation.
  unidentified <- fields[!fields$field %in% identities, ]
  utils::write.table(unidentified, receipt, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(duckhts_bench_stage_ont_ecoli(fetch = FALSE, samtools = samtools, minimap2 = ""),
               "minimap2 is required")
  paths <- duckhts_bench_stage_ont_ecoli(fetch = FALSE, threads = 1L,
                                       samtools = samtools, minimap2 = minimap2)

  # A changed read identity at the same cache path requires derivation.
  fastq <- readLines(reads)
  fastq[[1L]] <- "@changed-read"
  writeLines(fastq, reads)
  gzip_copy(reads, paths[["reads"]])
  registry$supplier_identity[[3L]] <- paste0(
    "bytes=", file.info(paths[["reads"]])$size, ";md5=", unname(tools::md5sum(paths[["reads"]])))
  utils::write.table(registry, registry_path, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(duckhts_bench_stage_ont_ecoli(fetch = FALSE, samtools = samtools, minimap2 = ""),
               "minimap2 is required")
  paths <- duckhts_bench_stage_ont_ecoli(fetch = FALSE, threads = 1L,
                                       samtools = samtools, minimap2 = minimap2)
  alignments <- system2(samtools, c("view", shQuote(paths[["bam"]])), stdout = TRUE)
  expect_true(any(startsWith(alignments, "changed-read\t")))

  # A changed reference identity at the same cache path requires derivation.
  fasta <- readLines(reference)
  fasta[[1L]] <- ">changed-reference"
  writeLines(fasta, reference)
  archive <- duckhts_bench_artifact_path("ont_ecoli_k12_reference_fna_gz")
  gzip_copy(reference, archive)
  registry$supplier_identity[1:2] <- c(
    paste0("bytes=", file.info(archive)$size, ";md5=", unname(tools::md5sum(archive))),
    paste0("bytes=", file.info(reference)$size, ";md5=", unname(tools::md5sum(reference)))
  )
  utils::write.table(registry, registry_path, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(duckhts_bench_stage_ont_ecoli(fetch = FALSE, samtools = samtools, minimap2 = ""),
               "minimap2 is required")
  paths <- duckhts_bench_stage_ont_ecoli(fetch = FALSE, threads = 1L,
                                       samtools = samtools, minimap2 = minimap2)
  header <- system2(samtools, c("view", "-H", shQuote(paths[["bam"]])), stdout = TRUE)
  expect_true(any(grepl("SN:changed-reference\t", header, fixed = TRUE)))
  fields <- utils::read.delim(receipt, colClasses = "character")
  hashes <- vapply(paths[c("reference", "reads")], function(path) {
    digest::digest(file = path, algo = "sha256")
  }, character(1L))
  expect_identical(fields$value[match(identities, fields$field)], unname(hashes))
  expect_identical(duckhts_bench_stage_ont_ecoli(fetch = FALSE, samtools = samtools,
                                               minimap2 = ""), paths)

  # A poisoned BAM fails quickcheck and is rebuilt from the verified sources.
  writeLines("poisoned", paths[["bam"]])
  paths <- duckhts_bench_stage_ont_ecoli(fetch = FALSE, threads = 1L,
                                       samtools = samtools, minimap2 = minimap2)
  aligned <- system2(samtools, c("view", "-c", "-F", "4", shQuote(paths[["bam"]])), stdout = TRUE)
  expect_equal(as.integer(aligned), length(starts))

  # A same-size replacement of the index cannot be reused with the old receipt.
  index_bytes <- readBin(index, what = "raw", n = file.info(index)$size)
  index_bytes[[1L]] <- as.raw(bitwXor(as.integer(index_bytes[[1L]]), 1L))
  writeBin(index_bytes, index)
  expect_error(duckhts_bench_stage_ont_ecoli(fetch = FALSE, samtools = samtools,
                                             minimap2 = ""), "minimap2 is required")
  paths <- duckhts_bench_stage_ont_ecoli(fetch = FALSE, threads = 1L,
                                       samtools = samtools, minimap2 = minimap2)
  expect_identical(readBin(paste0(paths[["bam"]], ".bai"), what = "raw", n = 4L),
                   charToRaw("BAI\001"))
  fields <- utils::read.delim(receipt, colClasses = "character")
  expect_identical(fields$value[fields$field == "index_sha256"],
                   digest::digest(file = index, algo = "sha256"))
  aligned <- system2(samtools, c("view", "-c", "-F", "4", shQuote(paths[["bam"]])), stdout = TRUE)
  expect_equal(as.integer(aligned), length(starts))

  # A source that no longer matches its registered identity is refused.
  writeLines("not the archive", duckhts_bench_artifact_path("ont_ecoli_k12_reads_fastq_gz"))
  expect_error(duckhts_bench_stage_ont_ecoli(fetch = FALSE, samtools = samtools),
               "identity does not match")
}

test_ont_ecoli_derivation()
