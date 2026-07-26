# BigWig reader and wrapper tests
library(tinytest)
library(DBI)

test_bigwig <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  bigwig_path <- system.file(
    "extdata", "libbigwig_test.bw", package = "Rduckhts", mustWork = TRUE
  )

  read_little_endian <- function(bytes, offset0, width) {
    sum(as.integer(bytes[offset0 + seq_len(width)]) * 256^(0:(width - 1L)))
  }
  write_little_endian <- function(bytes, offset0, value, width) {
    bytes[offset0 + seq_len(width)] <- as.raw(
      floor(value / 256^(0:(width - 1L))) %% 256
    )
    bytes
  }
  write_bigwig_bytes <- function(bytes) {
    path <- tempfile("rduckhts_corrupt_", fileext = ".bw")
    writeBin(bytes, path)
    path
  }
  query_bigwig_count <- function(path) {
    dbGetQuery(con, paste0(
      "SELECT count(*) FROM read_bigwig(",
      as.character(dbQuoteString(con, path)),
      ")"
    ))
  }

  expect_silent(rduckhts_bigwig(
    con,
    "bigwig_intervals",
    bigwig_path,
    region = c("1:1-150", "1:100-250", "10:201-300"),
    blocks_per_iteration = 1L,
    overwrite = TRUE
  ))

  observed <- dbGetQuery(
    con,
    paste(
      "SELECT CHROM, START0, END0, round(VALUE::DOUBLE, 1) AS VALUE",
      "FROM bigwig_intervals ORDER BY CHROM, START0"
    )
  )
  expect_equal(nrow(observed), 6L)
  expect_equal(observed$CHROM, c("1", "1", "1", "1", "1", "10"))
  expect_equal(observed$START0, c(0L, 1L, 2L, 100L, 150L, 200L))
  expect_equal(observed$END0, c(1L, 2L, 3L, 150L, 151L, 300L))
  expect_equal(observed$VALUE, c(0.1, 0.2, 0.3, 1.4, 1.5, 2.0))

  expect_error(
    rduckhts_bigwig(con, "bigwig_intervals", bigwig_path),
    "already exists"
  )
  expect_error(
    rduckhts_bigwig(con, "bad_blocks", bigwig_path, blocks_per_iteration = 0),
    "between 1 and 1048576"
  )
  expect_error(
    rduckhts_bigwig(con, "bad_region", bigwig_path, region = character()),
    "region must contain"
  )

  fixture_bytes <- readBin(
    bigwig_path,
    what = "raw",
    n = file.info(bigwig_path)$size
  )
  corrupt_paths <- character()
  on.exit(unlink(corrupt_paths), add = TRUE)

  truncated_zoom_header <- write_little_endian(
    fixture_bytes,
    6L,
    1,
    2L
  )
  truncated_zoom_header <- truncated_zoom_header[seq_len(68L)]
  corrupt_paths <- c(
    corrupt_paths,
    write_bigwig_bytes(truncated_zoom_header)
  )
  expect_error(
    query_bigwig_count(corrupt_paths[[length(corrupt_paths)]]),
    "failed to open a valid BigWig file"
  )

  chrom_tree_offset <- read_little_endian(fixture_bytes, 8L, 8L)
  chrom_key_size <- read_little_endian(
    fixture_bytes,
    chrom_tree_offset + 8L,
    4L
  )
  chrom_count <- read_little_endian(
    fixture_bytes,
    chrom_tree_offset + 16L,
    8L
  )
  huge_chrom_count <- write_little_endian(
    fixture_bytes,
    chrom_tree_offset + 16L,
    1048577,
    8L
  )
  corrupt_paths <- c(corrupt_paths, write_bigwig_bytes(huge_chrom_count))
  expect_error(
    query_bigwig_count(corrupt_paths[[length(corrupt_paths)]]),
    "failed to open a valid BigWig file"
  )

  chrom_index_offset <- chrom_tree_offset + 36L + chrom_key_size
  bad_chrom_index <- write_little_endian(
    fixture_bytes,
    chrom_index_offset,
    chrom_count,
    4L
  )
  corrupt_paths <- c(corrupt_paths, write_bigwig_bytes(bad_chrom_index))
  expect_error(
    query_bigwig_count(corrupt_paths[[length(corrupt_paths)]]),
    "failed to open a valid BigWig file"
  )

  index_offset <- read_little_endian(fixture_bytes, 24L, 8L)
  index_root_offset <- index_offset + 48L
  expect_identical(
    read_little_endian(fixture_bytes, index_root_offset, 1L),
    1
  )
  bad_child_index <- write_little_endian(
    fixture_bytes,
    index_root_offset,
    0,
    1L
  )
  bad_child_index <- write_little_endian(
    bad_child_index,
    index_root_offset + 20L,
    length(fixture_bytes) + 1024,
    8L
  )
  corrupt_paths <- c(corrupt_paths, write_bigwig_bytes(bad_child_index))
  expect_error(
    query_bigwig_count(corrupt_paths[[length(corrupt_paths)]]),
    "failed to initialize an indexed range"
  )

  huge_block <- write_little_endian(
    fixture_bytes,
    index_root_offset + 28L,
    64 * 1024^2 + 1,
    8L
  )
  corrupt_paths <- c(corrupt_paths, write_bigwig_bytes(huge_block))
  expect_error(
    query_bigwig_count(corrupt_paths[[length(corrupt_paths)]]),
    "failed to initialize an indexed range"
  )

  huge_decode_buffer <- write_little_endian(
    fixture_bytes,
    52L,
    64 * 1024^2 + 1,
    4L
  )
  corrupt_paths <- c(corrupt_paths, write_bigwig_bytes(huge_decode_buffer))
  expect_error(
    query_bigwig_count(corrupt_paths[[length(corrupt_paths)]]),
    "failed to initialize an indexed range"
  )

  short_block <- write_little_endian(fixture_bytes, 52L, 0, 4L)
  short_block <- write_little_endian(
    short_block,
    index_root_offset + 28L,
    1,
    8L
  )
  corrupt_paths <- c(corrupt_paths, write_bigwig_bytes(short_block))
  expect_error(
    query_bigwig_count(corrupt_paths[[length(corrupt_paths)]]),
    "failed to initialize an indexed range|failed while reading an indexed range"
  )
}

test_bigwig()
