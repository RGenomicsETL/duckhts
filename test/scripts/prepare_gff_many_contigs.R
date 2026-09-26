#!/usr/bin/env Rscript
# Authored synthetic input; no external source. Reproduce with HTSlib 1.19.
# The benchmark registry records the committed BGZF and tabix identities.

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args) > 0L) args[[1L]] else "."
repo <- normalizePath(repo, winslash = "/", mustWork = TRUE)
output <- file.path(repo, "test/data/gff_many_contigs.gff.gz")
index <- paste0(output, ".tbi")

required_tools <- c("bgzip", "tabix", "sha256sum")
tool_paths <- Sys.which(required_tools)
if (any(!nzchar(tool_paths))) {
  stop("required tools are unavailable: ",
       paste(required_tools[!nzchar(tool_paths)], collapse = ", "),
       call. = FALSE)
}

contig_count <- 10000L
record_count <- 2L * contig_count
contigs <- sprintf("contig%05d", seq_len(contig_count))
gene_ids <- sprintf("gene%05d", seq_len(contig_count))
transcript_ids <- sprintf("transcript%05d", seq_len(contig_count))
records <- character(record_count)
records[seq.int(1L, record_count, by = 2L)] <- sprintf(
  "%s\tduckhts_fixture\tgene\t1\t100\t.\t+\t.\tID=%s;Parent=.;transcript_type=protein_coding;gene_type=protein_coding",
  contigs, gene_ids
)
records[seq.int(2L, record_count, by = 2L)] <- sprintf(
  "%s\tduckhts_fixture\tmRNA\t10\t90\t.\t+\t.\tID=%s;Parent=%s;transcript_type=protein_coding;gene_type=protein_coding",
  contigs, transcript_ids, gene_ids
)

plain <- tempfile("gff-many-contigs-", fileext = ".gff3")
writeLines(c(
  "##gff-version 3",
  "#!fixture-generator test/scripts/prepare_gff_many_contigs.R",
  "#!fixture-layout 10000-contigs;2-records-per-contig;keys=ID,Parent,transcript_type,gene_type",
  records
), plain, useBytes = TRUE)

unlink(c(output, index), force = TRUE)
if (system2(tool_paths[["bgzip"]], c("-@", "1", "-c", shQuote(plain)),
            stdout = output) != 0L) {
  stop("bgzip failed", call. = FALSE)
}
unlink(plain)
if (system2(tool_paths[["tabix"]], c("-f", "-p", "gff", shQuote(output))) != 0L) {
  stop("tabix failed", call. = FALSE)
}

lines <- readLines(gzfile(output), warn = FALSE)
data_lines <- lines[!startsWith(lines, "#")]
indexed_contigs <- system2(tool_paths[["tabix"]], c("-l", shQuote(output)), stdout = TRUE)
indexed_lines <- system2(
  tool_paths[["tabix"]], c(shQuote(output), contigs), stdout = TRUE
)
stopifnot(
  length(data_lines) == record_count,
  length(indexed_contigs) == contig_count,
  identical(indexed_contigs, contigs),
  identical(indexed_lines, data_lines),
  all(lengths(strsplit(data_lines, "\t", fixed = TRUE)) == 9L),
  all(vapply(c("ID=", "Parent=", "transcript_type=", "gene_type="),
             function(key) all(grepl(key, data_lines, fixed = TRUE)), logical(1L)))
)

sha256 <- function(path) {
  output <- system2(tool_paths[["sha256sum"]], shQuote(path), stdout = TRUE)
  strsplit(output, "[[:space:]]+")[[1L]][[1L]]
}
actual_bytes <- unname(file.info(c(output, index))$size)
actual_sha256 <- c(sha256(output), sha256(index))
provenance <- data.frame(
  artifact = c(output, index),
  bytes = actual_bytes,
  sha256 = actual_sha256,
  generator = "test/scripts/prepare_gff_many_contigs.R",
  R = R.version.string,
  bgzip = system2(tool_paths[["bgzip"]], "--version", stdout = TRUE)[[1L]],
  tabix = system2(tool_paths[["tabix"]], "--version", stdout = TRUE)[[1L]],
  stringsAsFactors = FALSE
)
print(provenance, row.names = FALSE)
