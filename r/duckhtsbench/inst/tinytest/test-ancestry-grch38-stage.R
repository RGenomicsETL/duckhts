library(tinytest)

if (requireNamespace("DBI", quietly = TRUE) && requireNamespace("duckdb", quietly = TRUE) &&
    requireNamespace("Rduckhts", quietly = TRUE)) local({
  root <- tempfile("ancestry-grch38-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  ids <- c("ancestry_ref_freqs", "ancestry_projection", "ancestry_reference_parquet",
           "ancestry_reference_grch38_parquet", "liftover_grch37_grch38_chain",
           "liftover_grch37_fasta", "liftover_grch38_fasta")
  registry <- duckhts_bench_registry()
  plan <- registry[match(ids, registry$id), , drop = FALSE]
  expect_equal(plan$transform[plan$id == "ancestry_reference_grch38_parquet"],
               "duckhts_bench_stage_ancestry_grch38_parquet")
  plan$cache_relpath <- file.path("cache", basename(plan$cache_relpath))
  dir.create(file.path(root, "cache"))
  path_of <- function(id) file.path(root, plan$cache_relpath[plan$id == id])
  identify <- function(id, fields) {
    plan$supplier_identity[plan$id == id] <<- paste0("bytes=", file.info(path_of(id))$size,
                                                     fields(path_of(id)))
  }
  sha <- function(path) paste0(";sha256=", digest::digest(file = path, algo = "sha256"))
  md5 <- function(path) paste0(";md5=", unname(tools::md5sum(path)))

  # Source loci (GRCh37 spelling); each one exercises one liftover outcome.
  loci <- data.frame(
    chr = c(1L, 1L, 1L, 3L, 1L, 1L, 1L, 2L, 2L, 4L),
    pos = c(10L, 20L, 35L, 35L, 50L, 90L, 60L, 30L, 40L, 10L),
    a0 = c("A", "A", "C", "C", "C", "A", "AT", "A", "A", "A"),
    a1 = c("G", "G", "T", "T", "T", "G", "A", "C", "C", "G"),
    stringsAsFactors = FALSE)
  loci$rsid <- paste0("rs", seq_len(nrow(loci)))
  frequency <- matrix(seq_len(nrow(loci) * 21L) / (nrow(loci) * 21L + 1), nrow(loci))
  loading <- matrix(seq_len(nrow(loci) * 16L) - 80, nrow(loci)) / 10
  sources <- list(
    ancestry_ref_freqs = cbind(loci, setNames(as.data.frame(frequency), paste0("group", 1:21))),
    ancestry_projection = cbind(loci, setNames(as.data.frame(loading), paste0("PC", 1:16))))
  for (id in names(sources)) {
    connection <- gzfile(path_of(id), "wt")
    utils::write.csv(sources[[id]], connection, row.names = FALSE)
    close(connection)
    identify(id, sha)
  }

  # Source contigs 1-4 and destination contigs chr1, chr2, chrX, 100 bases each.
  fasta <- function(path, bases) {
    lines <- character()
    offsets <- integer()
    total <- 0L
    for (name in names(bases)) {
      header <- paste0(">", name)
      lines <- c(lines, header, paste(bases[[name]], collapse = ""))
      offsets <- c(offsets, total + nchar(header) + 1L)
      total <- total + nchar(header) + 1L + 101L
    }
    writeLines(lines, path)
    writeLines(paste(names(bases), lengths(bases), offsets, 100L, 101L, sep = "\t"),
               paste0(path, ".fai"))
  }
  filler <- function() rep("A", 100L)
  source_bases <- list("1" = filler(), "2" = filler(), "3" = filler(), "4" = filler())
  source_bases[["1"]][61L] <- "T"
  destination_bases <- list(chr1 = filler(), chr2 = filler(), chrX = filler())
  destination_bases$chr1[c(20L, 35L, 50L)] <- c("G", "C", "G")
  destination_bases$chr2[c(71L, 61L)] <- c("T", "G")
  fasta(path_of("liftover_grch37_fasta"), source_bases)
  fasta(path_of("liftover_grch38_fasta"), destination_bases)
  identify("liftover_grch37_fasta", function(path) "")
  identify("liftover_grch38_fasta", md5)
  chain <- c("chain 800 1 100 + 0 80 chr1 100 + 0 80 1", "80", "",
             "chain 1000 2 100 + 0 100 chr2 100 - 0 100 2", "100", "",
             "chain 1000 3 100 + 0 100 chr1 100 + 0 100 3", "100", "",
             "chain 1000 4 100 + 0 100 chrX 100 + 0 100 4", "100", "")
  connection <- gzfile(path_of("liftover_grch37_grch38_chain"), "wt")
  writeLines(chain, connection)
  close(connection)
  identify("liftover_grch37_grch38_chain", md5)
  registry_path <- file.path(root, "registry.tsv")
  utils::write.table(plan, registry_path, sep = "\t", row.names = FALSE, quote = FALSE)
  previous <- Sys.getenv(c("DUCKHTSBENCH_REGISTRY", "DUCKHTS_CACHE_DIR"), unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(previous[name]))
  }, add = TRUE)
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry_path, DUCKHTS_CACHE_DIR = root)

  output <- duckhts_bench_stage_ancestry_grch38_parquet()
  expect_equal(output, path_of("ancestry_reference_grch38_parquet"))
  receipt <- utils::read.delim(paste0(output, ".sources.tsv"), colClasses = "character")
  value <- function(field) receipt$value[receipt$field == field]
  expect_equal(value("source_parquet_sha256"),
               digest::digest(file = path_of("ancestry_reference_parquet"), algo = "sha256"))
  expect_equal(value("chain_sha256"),
               digest::digest(file = path_of("liftover_grch37_grch38_chain"), algo = "sha256"))
  expect_equal(value("source_fasta_sha256"),
               digest::digest(file = path_of("liftover_grch37_fasta"), algo = "sha256"))
  expect_equal(value("destination_fasta_sha256"),
               digest::digest(file = path_of("liftover_grch38_fasta"), algo = "sha256"))
  expect_equal(value("output_sha256"), digest::digest(file = output, algo = "sha256"))
  expect_equal(value("derivation_sha256"),
               digest::digest(plan$transform[plan$id == "ancestry_reference_grch38_parquet"],
                              algo = "sha256", serialize = FALSE))
  for (field in c("duckdb_version", "duckdb_r_version", "duckhts_htslib_version",
                  "rduckhts_version")) {
    expect_true(nzchar(value(field)))
  }
  expect_equal(as.integer(value("input_loci")), 10L)
  expect_equal(as.integer(value("output_loci")), 4L)
  expect_equal(as.integer(value("swapped")), 2L)
  expect_equal(as.integer(value("reverse_complemented")), 2L)
  expect_equal(as.integer(value("duplicate_destination_dropped")), 2L)
  expect_equal(as.integer(value("rejected")), 4L)
  reasons <- receipt[startsWith(receipt$field, "rejected_"), ]
  expect_equal(sum(as.integer(reasons$value)), 4L)
  expect_true(all(c("rejected_destination_allele_mismatch", "rejected_non_autosomal_destination",
                    "rejected_liftover_UnmappedAnchors") %in% reasons$field))
  expect_equal(as.integer(value("input_loci")),
               as.integer(value("rejected")) + as.integer(value("duplicate_destination_dropped")) +
                 as.integer(value("output_loci")))

  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  rows <- DBI::dbGetQuery(con, paste0("SELECT * FROM read_parquet(",
    as.character(DBI::dbQuoteString(con, output)), ")"))
  expect_equal(nrow(rows), 4L)
  expect_equal(ncol(rows), 41L)
  expect_equal(rows$chromosome, c(1L, 1L, 2L, 2L))
  # Destination positions: forward loci keep theirs, reverse-strand loci mirror.
  expect_equal(rows$position, c(10L, 20L, 61L, 71L))
  # The other allele is allele_a and the coded allele is allele_b, in destination
  # spelling; a role swap changes the key but never a frequency or loading.
  expect_equal(rows$allele_a, c("A", "A", "T", "T"))
  expect_equal(rows$allele_b, c("G", "G", "G", "G"))
  source_row <- c(1L, 2L, 9L, 8L)
  expect_equal(unname(as.matrix(rows[paste0("group", 1:21)])), unname(frequency[source_row, ]))
  expect_equal(unname(as.matrix(rows[paste0("PC", 1:16)])), unname(loading[source_row, ]))
  map_output <- sub("\\.parquet$", ".liftover.parquet", output)
  expect_equal(value("liftover_map_sha256"), digest::digest(file = map_output, algo = "sha256"))
  mapped <- DBI::dbGetQuery(con, paste0("SELECT * FROM read_parquet(",
    as.character(DBI::dbQuoteString(con, map_output)), ") ORDER BY 1, 2"))
  expect_equal(mapped$source_position, c(10L, 20L, 40L, 30L))
  expect_equal(mapped$source_chromosome, c(1L, 1L, 2L, 2L))
  expect_equal(mapped$swapped, c(FALSE, TRUE, TRUE, FALSE))
  expect_equal(mapped$reverse_complemented, c(FALSE, FALSE, TRUE, TRUE))
  expect_identical(duckhts_bench_stage_ancestry_grch38_parquet(), output)

  # A changed source binds nothing to the old receipt.
  fasta(path_of("liftover_grch37_fasta"), utils::modifyList(source_bases, list("4" = rep("C", 100L))))
  expect_error(duckhts_bench_stage_ancestry_grch38_parquet(), pattern = "identity")
  fasta(path_of("liftover_grch37_fasta"), source_bases)
  writeLines("corrupt", output)
  expect_error(duckhts_bench_stage_ancestry_grch38_parquet(), pattern = "identity")
  expect_true(file.exists(paste0(output, ".sources.tsv")))
})
