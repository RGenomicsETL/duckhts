#!/usr/bin/env Rscript
# Maintainer reconstruction: synthetic edge alignments, no downloaded inputs.
# samtools is the executable bundled by RBCFTools; Riker is a separately built
# v0.4.1 executable supplied as RIKER_BIN, never part of the extension build.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else ".", mustWork = TRUE)
riker <- Sys.getenv("RIKER_BIN")
stopifnot(file.exists(riker), requireNamespace("RBCFTools", quietly = TRUE))
source(file.path(root, "benchmarks/bam_mismatch_stage.R"))
output <- file.path(root, "test/data/insert_size")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
scratch <- tempfile("insert-size-fixtures-")
dir.create(scratch)
reference_plain <- file.path(scratch, "reference.fa")
reference <- file.path(output, "reference.fa.gz")
writeLines(c(">chr1", paste(rep("A", 2000100L), collapse = ""),
             ">chr2", paste(rep("A", 100000L), collapse = "")), reference_plain)
status <- system2(RBCFTools::bgzip_path(), shQuote(c("-c", reference_plain)), stdout = reference)
if (status != 0L) stop("reference compression failed", call. = FALSE)
bam_mismatch_run("samtools", c("faidx", reference))
header <- c("@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:2000100", "@SQ\tSN:chr2\tLN:100000",
            "@RG\tID:edges\tSM:edges")
row <- function(name, flag = 161L, pos = 100L, cigar = "10M", mate = "=", pnext = 190L,
                tlen = 100, seq = "AAAAAAAAAA", mapq = 0L) {
  paste(name, flag, "chr1", pos, mapq, cigar, mate, pnext,
        format(tlen, scientific = FALSE, trim = TRUE), seq,
        if (seq == "*") "*" else paste(rep("I", nchar(seq)), collapse = ""),
        "RG:Z:edges", sep = "\t")
}
# riker v0.4.1 commands/isize.rs:197-206 rejects first segments and QCFAIL.
# Improper pairs and MAPQ 0/255 remain eligible. No mate lookup is needed.
edges <- c(
  row("ordinary"), row("even_median", tlen = 200), row("tail", tlen = 1000),
  row("qcfail", flag = 161L + 512L), row("proper", flag = 163L),
  row("mapq255", mapq = 255L), row("no_segment_bit", flag = 33L),
  row("duplicate", flag = 161L + 1024L, tlen = 300),
  row("secondary", flag = 161L + 256L), row("supplementary", flag = 161L + 2048L),
  row("first", flag = 97L), row("unpaired", flag = 160L),
  row("mate_unmapped", flag = 169L), row("zero", tlen = 0),
  row("cross_reference", mate = "chr2"),
  row("forward_negative", tlen = -100), row("forward_equal", tlen = 1),
  row("reverse_cigar_end", flag = 145L, pnext = 104L, tlen = -200),
  row("reverse_equal", flag = 145L, cigar = "1M", pnext = 100L, tlen = -1, seq = "A"),
  row("reverse_skip", flag = 145L, cigar = "2S3M20N2M3S", pnext = 120L, tlen = -200),
  row("reverse_gap", flag = 145L, cigar = "3M5D7M", pnext = 111L, tlen = -100),
  row("missing_cigar", flag = 145L, cigar = "*", pnext = 100L, tlen = -2),
  row("tandem_forward", flag = 129L, tlen = 10),
  row("tandem_reverse", flag = 177L, tlen = -20),
  row("int32_min", flag = 177L, tlen = -2147483648),
  paste("unmapped", 165L, "*", 0L, 0L, "*", "chr1", 100L, 100L, "*", "*", "RG:Z:edges", sep = "\t")
)
fixtures <- list(edges = edges, empty = character(),
                 singleton = row("singleton"),
                 task_overlap = c(row("left", pos = 999995L, pnext = 1000004L),
                                  row("right", pos = 1000001L, pnext = 1000010L)),
                 long_cigar = row("long", flag = 145L,
                   cigar = paste(rep("1M1D", 32768L), collapse = ""),
                   seq = paste(rep("A", 32768L), collapse = ""),
                   pnext = 65634L, tlen = -65536),
                 even = c(row("a", tlen = 100), row("b", tlen = 200)),
                 mad_zero = c(row("a"), row("b"), row("tail", tlen = 1000)),
                 minority = c(vapply(seq_len(19L), function(i) row(paste0("fr", i)), ""),
                              row("rf", tlen = -100)))
# isize.rs:279-290 emits the complete raw histogram when no orientation
# passes min_frac. The raw cases are the independent integer-count oracle.
cases <- data.frame(name = c("default", "duplicates", "all_orientations", "median_only",
                             "filtered", "raw", "raw_duplicates"),
                    duplicates = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE),
                    deviations = c(10, 10, 10, 0, 10, 10, 10),
                    min_frac = c(.05, .05, 0, 0, .051, 1.1, 1.1))
write.table(cases, file.path(output, "cases.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
metrics <- list()
histograms <- list()
for (name in names(fixtures)) {
  sam <- file.path(scratch, paste0(name, ".sam"))
  writeLines(c(header, fixtures[[name]]), sam)
  bam <- file.path(output, paste0(name, ".bam"))
  bam_mismatch_run("samtools", c("sort", "--no-PG", "-o", bam, sam))
  bam_mismatch_run("samtools", c("index", bam))
  cram <- file.path(output, paste0(name, ".cram"))
  bam_mismatch_run("samtools", c("view", "--no-PG", "-C", "-T", reference, "-o", cram, bam))
  bam_mismatch_run("samtools", c("index", cram))
  for (format in c("bam", "cram")) {
    for (i in seq_len(nrow(cases))) {
      prefix <- file.path(scratch, paste(name, format, cases$name[[i]], sep = "."))
      arguments <- c("isize", "-i", paste0(file.path(output, name), ".", format),
                     "-r", reference, "-o", prefix, "--threads", "1", "--deviations",
                     cases$deviations[[i]], "--min-frac", cases$min_frac[[i]])
      if (cases$duplicates[[i]]) arguments <- c(arguments, "--include-duplicates")
      status <- system2(riker, shQuote(arguments))
      if (status != 0L) stop("Riker fixture oracle failed")
      identity <- data.frame(fixture = name, format = format, case_name = cases$name[[i]])
      metric_path <- paste0(prefix, ".isize-metrics.txt")
      histogram_path <- paste0(prefix, ".isize-histogram.txt")
      if (file.info(metric_path)$size > 0) {
        metrics[[length(metrics) + 1L]] <- cbind(identity,
          read.delim(metric_path, colClasses = "character", check.names = FALSE))
      }
      if (file.info(histogram_path)$size > 0) {
        histograms[[length(histograms) + 1L]] <- cbind(identity,
          read.delim(histogram_path, colClasses = "character", check.names = FALSE))
      }
    }
  }
}
write.table(do.call(rbind, metrics), file.path(output, "metrics.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(do.call(rbind, histograms), file.path(output, "histograms.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
provenance <- c(system2(riker, "--version", stdout = TRUE),
                "source: fulcrumgenomics/riker tag v0.4.1, 40d5e3260c5b1b9d210c28f7c90a75387013b4ce (MIT)",
                paste("riker_binary_sha256:", sub(" .*", "", system2("sha256sum", shQuote(riker), stdout = TRUE))),
                paste("samtools:", system2(bam_mismatch_tool("samtools"), "--version", stdout = TRUE)[[1L]]),
                "Input: synthetic SAM records from test/scripts/prepare_insert_size_fixtures.R",
                "BAM: samtools sort --no-PG; CRAM: samtools view --no-PG -C -T reference.fa.gz; both indexed.")
writeLines(provenance, file.path(output, "provenance.txt"))
stopifnot(file.copy(file.path(output, "edges.bam"), file.path(output, "unindexed.bam"), overwrite = TRUE))
stopifnot(file.copy(file.path(root, "test/sql/insert_size_reference.sql"),
                    output, overwrite = TRUE))
package_output <- file.path(root, "r/Rduckhts/inst/extdata/insert_size")
dir.create(package_output, recursive = TRUE, showWarnings = FALSE)
package_files <- c("edges.bam", "edges.bam.bai", "edges.cram", "edges.cram.crai",
                   "empty.bam", "empty.bam.bai", "empty.cram", "empty.cram.crai",
                   "unindexed.bam", "reference.fa.gz", "reference.fa.gz.fai", "reference.fa.gz.gzi",
                   "cases.tsv", "provenance.txt", "insert_size_reference.sql")
stopifnot(all(file.copy(file.path(output, package_files), package_output, overwrite = TRUE)))
for (oracle in c("metrics", "histograms")) {
  frame <- read.delim(file.path(output, paste0(oracle, ".tsv")), colClasses = "character")
  frame <- frame[frame$fixture %in% c("edges", "empty"), , drop = FALSE]
  write.table(frame, file.path(package_output, paste0(oracle, ".tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
}
unlink(scratch, recursive = TRUE)
