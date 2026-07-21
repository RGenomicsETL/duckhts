library(tinytest)
library(DBI)

test_variantkey_regionkey <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  expect_silent(rduckhts_load(con))

  contig_keys <- dbGetQuery(
    con,
    paste(
      "SELECT duckhts_contig_key(contig) AS key",
      "FROM (VALUES",
      "('chr1'), ('1'), ('chrM'), ('CHRMT'), ('M'), ('mt'), ('chrX'), ('chry'),",
      "('chrGL000220.1'), ('NC_000001.11'), ('23'), ('chrchr1'), ('chr'), (''), (NULL)",
      ") AS t(contig)"
    )
  )
  expect_equal(
    contig_keys$key,
    c("1", "1", "MT", "MT", "MT", "MT", "X", "Y", "GL000220.1", "NC_000001.11", "23", "chr1", "chr", "", NA_character_)
  )

  vkx <- dbGetQuery(
    con,
    "SELECT variantkey_hex(variantkey('chr1', 268435456, 'C', 'T')) AS vkx"
  )
  expect_equal(vkx$vkx[1], "0fffffff88b80000")

  raw_vk <- dbGetQuery(
    con,
    paste(
      "SELECT variantkey_hex(encode_variantkey(1, 324683, 145752064)) AS vkx,",
      "extract_variantkey_chrom(parse_variantkey_hex('08027a2588b00000')) AS chrom_code,",
      "extract_variantkey_pos(parse_variantkey_hex('08027a2588b00000')) AS pos0,",
      "extract_variantkey_refalt(parse_variantkey_hex('08027a2588b00000')) AS refalt_code"
    )
  )
  expect_equal(raw_vk$vkx[1], "08027a2588b00000")
  expect_equal(raw_vk$chrom_code[1], 1)
  expect_equal(raw_vk$pos0[1], 324683)
  expect_equal(raw_vk$refalt_code[1], 145752064)

  decode_vk <- dbGetQuery(
    con,
    paste(
      "SELECT dv.chrom_code, dv.pos0, dv.refalt_code",
      "FROM (SELECT decode_variantkey(parse_variantkey_hex('08027a2588b00000')) AS dv)"
    )
  )
  expect_equal(decode_vk$chrom_code[1], 1)
  expect_equal(decode_vk$pos0[1], 324683)
  expect_equal(decode_vk$refalt_code[1], 145752064)

  large_hash <- dbGetQuery(
    con,
    "SELECT variantkey_hex(variantkey('1', 324696, 'ATGCTTCTTTGCAGCCTCTCCAGGCCCAGCTCC', 'A')) AS vkx"
  )
  expect_equal(large_hash$vkx[1], "08027a2bed669117")

  symbolic <- dbGetQuery(
    con,
    paste(
      "SELECT rv.pos, rv.ref IS NULL AS ref_is_null, rv.reversible",
      "FROM (SELECT reverse_variantkey(variantkey('1', 101, 'A', '<DEL>')) AS rv)"
    )
  )
  expect_equal(symbolic$pos[1], 101)
  expect_true(isTRUE(symbolic$ref_is_null[1]))
  expect_false(isTRUE(symbolic$reversible[1]))

  reverse_vk <- dbGetQuery(
    con,
    paste(
      "SELECT rv.chrom, rv.chrom_code, rv.pos, rv.pos0, rv.ref, rv.alt, rv.reversible",
      "FROM (SELECT reverse_variantkey(parse_variantkey_hex('08027a2588b00000')) AS rv)"
    )
  )
  expect_equal(reverse_vk$chrom[1], "1")
  expect_equal(reverse_vk$chrom_code[1], 1)
  expect_equal(reverse_vk$pos[1], 324684)
  expect_equal(reverse_vk$pos0[1], 324683)
  expect_equal(reverse_vk$ref[1], "C")
  expect_equal(reverse_vk$alt[1], "G")
  expect_true(isTRUE(reverse_vk$reversible[1]))

  vk_range <- dbGetQuery(
    con,
    paste(
      "SELECT variantkey_hex(vkr.min) AS vk_min, variantkey_hex(vkr.max) AS vk_max",
      "FROM (SELECT variantkey_range('1', 100, 100) AS vkr)"
    )
  )
  expect_equal(vk_range$vk_min[1], "0800003180000000")
  expect_equal(vk_range$vk_max[1], "08000031ffffffff")

  invalid_vk <- dbGetQuery(
    con,
    paste(
      "SELECT",
      "variantkey('GL000220.1', 1, 'A', 'C') IS NULL AS bad_chrom_is_null,",
      "variantkey('1', 0, 'A', 'C') IS NULL AS bad_pos_is_null"
    )
  )
  expect_true(isTRUE(invalid_vk$bad_chrom_is_null[1]))
  expect_true(isTRUE(invalid_vk$bad_pos_is_null[1]))

  rkx <- dbGetQuery(
    con,
    "SELECT regionkey_hex(regionkey('X', 1007, 1807, 1)) AS rkx"
  )
  expect_equal(rkx$rkx[1], "b80001f78000387a")

  raw_rk <- dbGetQuery(
    con,
    paste(
      "SELECT regionkey_hex(encode_regionkey(23, 1007, 1807, 1)) AS rkx,",
      "extract_regionkey_chrom(parse_regionkey_hex('b80001f78000387a')) AS chrom_code,",
      "extract_regionkey_startpos(parse_regionkey_hex('b80001f78000387a')) AS start0,",
      "extract_regionkey_endpos(parse_regionkey_hex('b80001f78000387a')) AS end0,",
      "extract_regionkey_strand(parse_regionkey_hex('b80001f78000387a')) AS strand_code"
    )
  )
  expect_equal(raw_rk$rkx[1], "b80001f78000387a")
  expect_equal(raw_rk$chrom_code[1], 23)
  expect_equal(raw_rk$start0[1], 1007)
  expect_equal(raw_rk$end0[1], 1807)
  expect_equal(raw_rk$strand_code[1], 1)

  decode_rk <- dbGetQuery(
    con,
    paste(
      "SELECT dr.chrom_code, dr.start, dr.\"end\", dr.strand_code",
      "FROM (SELECT decode_regionkey(parse_regionkey_hex('b80001f78000387a')) AS dr)"
    )
  )
  expect_equal(decode_rk$chrom_code[1], 23)
  expect_equal(decode_rk$start[1], 1007)
  expect_equal(decode_rk$end[1], 1807)
  expect_equal(decode_rk$strand_code[1], 1)

  reverse_rk <- dbGetQuery(
    con,
    paste(
      "SELECT rr.chrom, rr.chrom_code, rr.start, rr.\"end\", rr.strand, rr.strand_code",
      "FROM (SELECT reverse_regionkey(parse_regionkey_hex('b80001f78000387a')) AS rr)"
    )
  )
  expect_equal(reverse_rk$chrom[1], "X")
  expect_equal(reverse_rk$chrom_code[1], 23)
  expect_equal(reverse_rk$start[1], 1007)
  expect_equal(reverse_rk$end[1], 1807)
  expect_equal(reverse_rk$strand[1], 1)
  expect_equal(reverse_rk$strand_code[1], 1)

  extend_rk <- dbGetQuery(
    con,
    paste(
      "SELECT rr.start, rr.\"end\", rr.strand",
      "FROM (SELECT reverse_regionkey(extend_regionkey(regionkey('X', 10000, 20000, -1), 1000)) AS rr)"
    )
  )
  expect_equal(extend_rk$start[1], 9000)
  expect_equal(extend_rk$end[1], 21000)
  expect_equal(extend_rk$strand[1], -1)

  overlaps <- dbGetQuery(
    con,
    paste(
      "SELECT",
      "are_overlapping_regions('1', 2, 4, '1', 3, 7) AS overlap_regions,",
      "are_overlapping_regionkeys(regionkey('X', 1007, 1807, 1), parse_regionkey_hex('b80001f78000387a')) AS overlap_keys,",
      "are_overlapping_region_regionkey('X', 2000, 2100, parse_regionkey_hex('b80001f78000387a')) AS overlap_miss"
    )
  )
  expect_true(isTRUE(overlaps$overlap_regions[1]))
  expect_true(isTRUE(overlaps$overlap_keys[1]))
  expect_false(isTRUE(overlaps$overlap_miss[1]))

  invalid_rk <- dbGetQuery(
    con,
    paste(
      "SELECT",
      "regionkey('GL000220.1', 0, 10) IS NULL AS bad_chrom_is_null,",
      "regionkey('1', 10, 0) IS NULL AS bad_span_is_null"
    )
  )
  expect_true(isTRUE(invalid_rk$bad_chrom_is_null[1]))
  expect_true(isTRUE(invalid_rk$bad_span_is_null[1]))
}

test_variantkey_regionkey()
