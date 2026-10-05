# Maintainer-only, network-free reconstruction of the BAM scan regression corpus.
# Synthetic SAM records below are the source authority; samtools supplies BAM,
# BAI/CSI and reference-free CRAM/CRAI encoding. Run from the repository root.
# Originally generated with samtools 1.23 / HTSlib 1.23; no external reference.
stopifnot(file.exists("src/bam_reader.c"), nzchar(Sys.which("samtools")), nzchar(Sys.which("bgzip")))

# Copies a BAI without its per-reference statistics. The SAM specification
# makes the statistics optional: they are a pseudo-bin, number 37450, of each
# reference. The data bins, the linear index and the trailing no-coordinate
# count are copied unchanged.
write_bai_without_statistics <- function(source, output) {
  input <- file(source, "rb")
  on.exit(close(input), add = TRUE)
  copied <- rawConnection(raw(0L), "wb")
  on.exit(close(copied), add = TRUE)
  read_int32 <- function() readBin(input, "integer", n = 1L, size = 4L, endian = "little")
  write_int32 <- function(value) writeBin(as.integer(value), copied, size = 4L, endian = "little")
  magic <- readBin(input, "raw", n = 4L)
  stopifnot(identical(magic, as.raw(c(0x42, 0x41, 0x49, 0x01))))
  writeBin(magic, copied)
  reference_count <- read_int32()
  write_int32(reference_count)
  for (reference in seq_len(reference_count)) {
    bins <- lapply(seq_len(read_int32()), function(bin) {
      number <- read_int32()
      chunk_count <- read_int32()
      list(number = number, chunk_count = chunk_count,
           chunks = readBin(input, "raw", n = 16L * chunk_count))
    })
    kept <- Filter(function(bin) bin$number != 37450L, bins)
    write_int32(length(kept))
    for (bin in kept) {
      write_int32(bin$number)
      write_int32(bin$chunk_count)
      writeBin(bin$chunks, copied)
    }
    interval_count <- read_int32()
    write_int32(interval_count)
    writeBin(readBin(input, "raw", n = 8L * interval_count), copied)
  }
  writeBin(readBin(input, "raw", n = 8L), copied)
  stopifnot(length(readBin(input, "raw", n = 1L)) == 0L)
  writeBin(rawConnectionValue(copied), output)
}

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
    # Indexes with alignments of chr1 and chr2 and no statistics for them, with
    # and without the trailing no-coordinate count. htslib cannot locate the
    # reads without coordinates from either.
    if (name == "mixed") {
      without_statistics <- paste0(bam, ".nostats.bai")
      write_bai_without_statistics(paste0(bam, ".bai"), without_statistics)
      bytes <- readBin(without_statistics, "raw", n = file.size(without_statistics))
      without_count <- paste0(bam, ".nostats.legacy.bai")
      writeBin(head(bytes, -8L), without_count)
      stopifnot(all(file.copy(c(without_statistics, without_count),
                              "r/Rduckhts/inst/extdata", overwrite = TRUE)))
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
