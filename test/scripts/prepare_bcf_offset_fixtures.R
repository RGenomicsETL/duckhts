# Regenerates the read_bcf FILE_OFFSET fixtures of test/sql/bcf_offset.test in
# test/data and copies the ones the R tests use to r/Rduckhts/inst/extdata:
#   bcf_offset.vcf.gz (.tbi)        bgzip --binary: records straddle BGZF blocks
#   bcf_offset.htslib.vcf.gz (.tbi) HTSlib VCF writer: blocks end at record ends
#   bcf_offset.bcf (.csi)           HTSlib BCF writer: full blocks, so records
#                                   straddle them, and one record is lengthened
#                                   to end exactly at a block end
#   bcf_offset_tidy.vcf.gz          bgzip of tidy_chunk_boundary.vcf: one CSQ
#                                   record whose 2053 tidy rows cross a DuckDB
#                                   output chunk
#   bcf_record_ends.tsv             expected FILE_OFFSET of every record of
#                                   those four BGZF files
#   bcf_offset_plain.vcf            eight short records, uncompressed
#   bcf_offset_gzip.vcf.gz          the same, plain gzip (not BGZF)
#   bcf_offset_uncompressed.bcf     the same, as a BCF stream without BGZF
#                                   blocks (HTSlib writes uncompressed BCF, as
#                                   bcftools -Ou does, as level-0 BGZF instead)
# Run from the repository root: Rscript test/scripts/prepare_bcf_offset_fixtures.R
# It uses the bgzip, tabix and bcftools bundled with the RBCFTools package
# (HTSlib 1.24, bcftools 1.24). The output is deterministic.
#
# The expected positions do not come from HTSlib. They are read from the bytes
# of each file: the block size and uncompressed size in each BGZF block header
# and trailer, and the record ends in the stream that R's zlib decompresses
# (line ends for VCF, length prefixes for BCF). Each position is the BGZF
# virtual position just after the record; a record that ends its block counts
# as the start of the next block.

bgzip <- RBCFTools::bgzip_path()
tabix <- RBCFTools::tabix_path()
bcftools <- RBCFTools::bcftools_path()
out_dir <- "test/data"
package_dir <- "r/Rduckhts/inst/extdata"
stopifnot(dir.exists(out_dir), dir.exists(package_dir))

run <- function(command, args, stdout = "") {
  stopifnot(system2(command, shQuote(args), stdout = stdout) == 0L)
}
fixture <- function(name) file.path(out_dir, name)

little_endian <- function(bytes) {
  sum(as.numeric(bytes) * 256^(seq_along(bytes) - 1L))
}

# Compressed start and size, and uncompressed start and end, of each BGZF block.
bgzf_blocks <- function(bytes) {
  blocks <- NULL
  at <- 0
  stream_at <- 0
  while (at < length(bytes)) {
    # Gzip member with FEXTRA whose only subfield is BC: XLEN 6, 'B', 'C', SLEN 2.
    stopifnot(identical(as.integer(bytes[at + 1:4]), c(31L, 139L, 8L, 4L)),
              identical(as.integer(bytes[at + 11:16]), c(6L, 0L, 66L, 67L, 2L, 0L)))
    size <- little_endian(bytes[at + 17:18]) + 1
    stream_size <- little_endian(bytes[at + (size - 3):size])
    blocks <- rbind(blocks, data.frame(start = at, size = size, stream_start = stream_at,
                                       stream_end = stream_at + stream_size))
    at <- at + size
    stream_at <- stream_at + stream_size
  }
  stopifnot(at == length(bytes))
  blocks
}

# Stream positions (zero-based) at which each record starts and just after it
# ends, and its ID, in file order.
vcf_records <- function(stream) {
  line_ends <- which(stream == as.raw(10L)) # one-based index of '\n' = position after it
  line_starts <- c(0, line_ends[-length(line_ends)])
  found <- NULL
  for (i in seq_along(line_ends)) {
    line <- rawToChar(stream[(line_starts[i] + 1):(line_ends[i] - 1)])
    if (startsWith(line, "#")) next
    id <- strsplit(line, "\t", fixed = TRUE)[[1]][3]
    found <- rbind(found, data.frame(ID = id, stream_start = line_starts[i],
                                     stream_end = line_ends[i]))
  }
  found
}

bcf_records <- function(stream) {
  stopifnot(identical(rawToChar(stream[1:3]), "BCF"), identical(as.integer(stream[4:5]), c(2L, 2L)))
  at <- 9 + little_endian(stream[6:9])
  found <- NULL
  while (at < length(stream)) {
    shared_length <- little_endian(stream[at + 1:4])
    individual_length <- little_endian(stream[at + 5:8])
    # The ID follows 24 bytes of fixed fields: a typed string of under 15 chars.
    id_type <- as.integer(stream[at + 8 + 25])
    stopifnot(bitwAnd(id_type, 15L) == 7L, bitwShiftR(id_type, 4L) < 15L)
    id <- rawToChar(stream[at + 8 + 25 + seq_len(bitwShiftR(id_type, 4L))])
    end <- at + 8 + shared_length + individual_length
    found <- rbind(found, data.frame(ID = id, stream_start = at, stream_end = end))
    at <- end
  }
  stopifnot(at == length(stream))
  found
}

# The blocks and records of one BGZF fixture.
bgzf_layout <- function(name) {
  path <- fixture(name)
  bytes <- readBin(path, "raw", file.size(path))
  blocks <- bgzf_blocks(bytes)
  connection <- gzfile(path, "rb")
  stream <- readBin(connection, "raw", max(blocks$stream_end) + 1)
  close(connection)
  stopifnot(length(stream) == max(blocks$stream_end))
  records <- if (endsWith(name, ".bcf")) bcf_records(stream) else vcf_records(stream)
  list(blocks = blocks, records = records)
}

# The BGZF virtual position of the stream position just after a record.
virtual_position <- function(stream_end, blocks) {
  b <- which(blocks$stream_start < stream_end & stream_end <= blocks$stream_end)
  stopifnot(length(b) == 1L)
  if (stream_end == blocks$stream_end[b]) {
    return((blocks$start[b] + blocks$size[b]) * 65536)
  }
  blocks$start[b] * 65536 + (stream_end - blocks$stream_start[b])
}

# Index of the block holding a stream byte (zero-based position).
block_of_byte <- function(position, blocks) {
  which(blocks$stream_start <= position & position < blocks$stream_end)
}

# Three contigs, so that an indexed scan runs in parallel by contig. The PAD
# filler makes the stream several BGZF blocks long; its varying length moves
# the record ends within the blocks. Two records of chr2 share a position.
contigs <- data.frame(name = c("chr1", "chr2", "chr3"), n = c(40L, 40L, 40L))
bases <- c("A", "C", "G", "T")
header <- c(
  "##fileformat=VCFv4.2",
  sprintf("##contig=<ID=%s,length=1000000>", contigs$name),
  "##INFO=<ID=PAD,Number=1,Type=String,Description=\"Filler of varying length\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2"
)
make_records <- function(extra_pad) {
  records <- character()
  ordinal <- 0L
  for (k in seq_len(nrow(contigs))) {
    for (i in seq_len(contigs$n[k])) {
      ordinal <- ordinal + 1L
      pos <- if (contigs$name[k] == "chr2" && i == 21L) 20000L else i * 1000L
      pad_length <- 1200L + (ordinal * 37L) %% 800L + extra_pad[ordinal]
      stopifnot(pad_length <= 8000L)
      pad <- substr(strrep("ACGT", 2000L), 1L, pad_length)
      records <- c(records, paste(
        contigs$name[k], pos, sprintf("r%04d", ordinal),
        bases[ordinal %% 4L + 1L], bases[(ordinal + 1L) %% 4L + 1L], "50", "PASS",
        paste0("PAD=", pad), "GT",
        if (ordinal %% 2L == 0L) "0/1" else "0/0",
        if (ordinal %% 3L == 0L) "1/1" else "0/1",
        sep = "\t"
      ))
    }
  }
  records
}

work_dir <- tempfile("bcf-offset-")
dir.create(work_dir)
write_bgzf_fixtures <- function(records) {
  source <- file.path(work_dir, "bcf_offset.vcf")
  writeLines(c(header, records), source)
  run(bgzip, c("--binary", "-c", source), stdout = fixture("bcf_offset.vcf.gz"))
  run(tabix, c("-f", "-p", "vcf", fixture("bcf_offset.vcf.gz")))
  run(bcftools, c("view", "--no-version", "-Oz", "-o", fixture("bcf_offset.htslib.vcf.gz"), source))
  run(bcftools, c("index", "-f", "-t", fixture("bcf_offset.htslib.vcf.gz")))
  run(bcftools, c("view", "--no-version", "-Ob", "-o", fixture("bcf_offset.bcf"), source))
  run(bcftools, c("index", "-f", fixture("bcf_offset.bcf")))
}

# The BCF writer cuts full blocks, so a block end usually falls inside a record.
# Lengthening the record before the first such record by the bytes it lacks
# makes that record end exactly at the block end: one more PAD byte is one more
# BCF byte, since the PAD length keeps its two-byte encoding.
extra_pad <- integer(sum(contigs$n))
write_bgzf_fixtures(make_records(extra_pad))
layout <- bgzf_layout("bcf_offset.bcf")
inside <- vapply(layout$blocks$stream_end, function(end) {
  hit <- which(layout$records$stream_start < end & end < layout$records$stream_end)
  if (length(hit)) hit else NA_integer_
}, integer(1))
straddling <- inside[!is.na(inside)][1]
stopifnot(straddling >= 2L)
extra_pad[straddling - 1L] <- layout$blocks$stream_end[match(straddling, inside)] -
  layout$records$stream_start[straddling]
records <- make_records(extra_pad)
write_bgzf_fixtures(records)

run(bgzip, c("-c", fixture("tidy_chunk_boundary.vcf")), stdout = fixture("bcf_offset_tidy.vcf.gz"))

small_records <- sub("PAD=[ACGT]+", "PAD=ACGT", records[1:8])
writeLines(c(header, small_records), fixture("bcf_offset_plain.vcf"))
plain_gzip <- gzfile(fixture("bcf_offset_gzip.vcf.gz"), "wb")
writeLines(c(header, small_records), plain_gzip)
close(plain_gzip)
small_bcf <- file.path(work_dir, "small.bcf")
run(bcftools, c("view", "--no-version", "-Ob", "-o", small_bcf, fixture("bcf_offset_plain.vcf")))
connection <- gzfile(small_bcf, "rb")
naked_bcf <- readBin(connection, "raw", 1e6)
close(connection)
stopifnot(identical(rawToChar(naked_bcf[1:3]), "BCF"))
writeBin(naked_bcf, fixture("bcf_offset_uncompressed.bcf"))
unlink(work_dir, recursive = TRUE)

decoded <- function(name) {
  system2(bcftools, c("view", "-H", shQuote(fixture(name))), stdout = TRUE)
}
for (name in c("bcf_offset.vcf.gz", "bcf_offset.htslib.vcf.gz", "bcf_offset.bcf")) {
  stopifnot(identical(decoded(name), records))
}
for (name in c("bcf_offset_plain.vcf", "bcf_offset_gzip.vcf.gz", "bcf_offset_uncompressed.bcf")) {
  stopifnot(identical(decoded(name), small_records))
}
stopifnot(identical(decoded("bcf_offset_tidy.vcf.gz"), decoded("tidy_chunk_boundary.vcf")))

expected <- NULL
for (name in c("bcf_offset.vcf.gz", "bcf_offset.htslib.vcf.gz", "bcf_offset.bcf")) {
  layout <- bgzf_layout(name)
  blocks <- layout$blocks
  found <- layout$records
  stopifnot(identical(found$ID, sprintf("r%04d", seq_along(records))))
  # What the tests rely on: record ends in at least three blocks; records that
  # straddle a block end in the bgzip and BCF files and none in the HTSlib VCF;
  # a record other than the last ending exactly at a block end in the HTSlib
  # VCF and the BCF.
  first_block <- vapply(found$stream_start, block_of_byte, integer(1), blocks = blocks)
  last_block <- vapply(found$stream_end - 1, block_of_byte, integer(1), blocks = blocks)
  ends_block <- found$stream_end %in% blocks$stream_end
  stopifnot(length(unique(last_block)) >= 3L)
  if (name == "bcf_offset.htslib.vcf.gz") {
    stopifnot(all(first_block == last_block))
  } else {
    stopifnot(any(first_block != last_block))
  }
  if (name != "bcf_offset.vcf.gz") stopifnot(sum(ends_block[-nrow(found)]) >= 1L)
  offsets <- vapply(found$stream_end, virtual_position, numeric(1), blocks = blocks)
  stopifnot(all(diff(offsets) > 0))
  expected <- rbind(expected, data.frame(file = name, ID = found$ID,
                                         FILE_OFFSET = sprintf("%.0f", offsets)))
}
layout <- bgzf_layout("bcf_offset_tidy.vcf.gz")
stopifnot(identical(layout$records$ID, "tidy_chunk"))
expected <- rbind(expected, data.frame(
  file = "bcf_offset_tidy.vcf.gz", ID = "tidy_chunk",
  FILE_OFFSET = sprintf("%.0f", virtual_position(layout$records$stream_end, layout$blocks))))
write.table(expected, fixture("bcf_record_ends.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

bundled <- c("bcf_offset.vcf.gz", "bcf_offset.vcf.gz.tbi", "bcf_offset.bcf", "bcf_offset.bcf.csi",
             "bcf_offset_plain.vcf", "bcf_record_ends.tsv")
stopifnot(all(file.copy(fixture(bundled), package_dir, overwrite = TRUE)))
