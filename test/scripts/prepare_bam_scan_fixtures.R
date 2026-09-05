# Maintainer-only, network-free reconstruction of the BAM scan regression corpus.
# Synthetic SAM records below are the source authority; samtools supplies BAM,
# BAI/CSI and reference-free CRAM/CRAI encoding. Run from the repository root.
# Originally generated with samtools 1.23 / HTSlib 1.23; no external reference.
stopifnot(file.exists("src/bam_reader.c"), nzchar(Sys.which("samtools")), nzchar(Sys.which("bgzip")))

prepare_bam_scan_fixtures <- function() {
  tmp <- tempfile("bam-scan-fixtures-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))
  header <- c("@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:1000")
  multi_header <- c(header, "@SQ\tSN:empty\tLN:1000", "@SQ\tSN:chr2\tLN:1000")
  mapped <- "mapped1\t0\tchr1\t10\t60\t1M\t*\t0\t0\tA\tI"
  placed <- "placed_unmapped\t4\tchr1\t20\t0\t*\t*\t0\t0\tC\tI"
  tail <- "unplaced\t4\t*\t0\t0\t*\t*\t0\t0\tG\tI"
  fixtures <- list(
    mixed = c(multi_header, mapped, placed,
              "mapped2\t0\tchr2\t30\t60\t1M\t*\t0\t0\tT\tI",
              tail, tail),
    single = c(header, mapped, placed, tail),
    all_unplaced = c(multi_header, rep(tail, 2053L)),
    empty = multi_header
  )
  for (name in names(fixtures)) {
    sam <- file.path(tmp, paste0(name, ".sam"))
    writeLines(fixtures[[name]], sam)
    base <- file.path("test/data", paste0("bam_scan_", name))
    bam <- paste0(base, ".bam")
    cram <- paste0(base, ".cram")
    run <- function(args) stopifnot(system2("samtools", shQuote(args)) == 0L)
    run(c("view", "--no-PG", "-b", "-o", bam, sam))
    run(c("index", "-b", bam))
    run(c("index", "-c", bam))
    # BAI/CSI allow omission of the final uint64 no-coordinate count.
    # Remove only those eight bytes, leaving all contig/chunk metadata intact.
    for (index in c("bai", "csi")) {
      indexed <- paste0(bam, ".", index)
      handle <- if (index == "csi") gzfile(indexed, "rb") else file(indexed, "rb")
      bytes <- readBin(handle, "raw", n = 100000L)
      close(handle)
      stopifnot(length(bytes) > 8L)
      plain <- file.path(tmp, paste0("legacy.", index))
      writeBin(head(bytes, -8L), plain)
      legacy <- paste0(bam, ".legacy.", index)
      if (index == "csi") {
        stopifnot(system2("bgzip", c("-c", shQuote(plain)), stdout = legacy) == 0L)
      } else {
        stopifnot(file.copy(plain, legacy, overwrite = TRUE))
      }
      stopifnot(file.copy(legacy, "r/Rduckhts/inst/extdata", overwrite = TRUE))
    }
    run(c("view", "--no-PG", "-C", "--output-fmt-option", "no_ref=1", "-o", cram, sam))
    run(c("index", cram))
    # The decoded SAM fields, including duplicate physical records, must survive.
    expected <- fixtures[[name]][!startsWith(fixtures[[name]], "@")]
    for (path in c(bam, cram)) {
      decoded <- system2("samtools", c("view", shQuote(path)), stdout = TRUE)
      stopifnot(identical(decoded, expected))
    }
    paths <- c(bam, paste0(bam, ".bai"), paste0(bam, ".csi"), cram, paste0(cram, ".crai"))
    stopifnot(all(file.copy(paths, "r/Rduckhts/inst/extdata", overwrite = TRUE)))
  }
}

prepare_bam_scan_fixtures()
stopifnot(file.copy("test/data/bam_scan_malformed.sam",
                   "r/Rduckhts/inst/extdata", overwrite = TRUE))
