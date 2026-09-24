library(tinytest)
library(DBI)

test_cigar_utils <- function() {
  con <- rduckhts_connect()
  on.exit(
    {
      try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
    },
    add = TRUE
  )

  literal_metrics <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT",
      "cigar_left_soft_clip('5S90M5S') AS left_soft,",
      "cigar_right_soft_clip('5S90M5S') AS right_soft,",
      "cigar_query_length('5S90M5I') AS query_len,",
      "cigar_aligned_query_length('5S90M5I') AS aligned_query_len,",
      "cigar_reference_length('90M5D') AS ref_len,",
      "cigar_has_soft_clip('5S90M5S') AS has_soft,",
      "cigar_has_hard_clip('5H95M') AS has_hard,",
      "cigar_has_op('90M5D', 'D') AS has_d"
    )
  )

  expect_equal(literal_metrics$left_soft[[1]], 5)
  expect_equal(literal_metrics$right_soft[[1]], 5)
  expect_equal(literal_metrics$query_len[[1]], 100)
  expect_equal(literal_metrics$aligned_query_len[[1]], 90)
  expect_equal(literal_metrics$ref_len[[1]], 95)
  expect_true(isTRUE(literal_metrics$has_soft[[1]]))
  expect_true(isTRUE(literal_metrics$has_hard[[1]]))
  expect_true(isTRUE(literal_metrics$has_d[[1]]))

  invalid_cigars <- c(
    "10M0", "10M00", "18446744073709551617M", "9223372036854775807M1M",
    "9223372036854775807I1I", "9223372036854775807D1D", "9223372036854775807S1S",
    "1M0Q", "1M1Q", "1M1B", "0M", "M", "1M1", "-1M", "1M 1I", "1M1m", "1S0M"
  )
  invalid_rows <- paste0("(", DBI::dbQuoteString(con, invalid_cigars), ")", collapse = ", ")
  invalid <- DBI::dbGetQuery(con, paste0(
    "SELECT cigar_query_length(c) IS NULL AS query_null, ",
    "cigar_reference_length(c) IS NULL AS ref_null, ",
    "cigar_aligned_query_length(c) IS NULL AS aligned_null, ",
    "cigar_left_soft_clip(c) IS NULL AS left_null, ",
    "cigar_right_soft_clip(c) IS NULL AS right_null, ",
    "cigar_has_soft_clip(c) IS NULL AS soft_null, ",
    "cigar_has_hard_clip(c) IS NULL AS hard_null, ",
    "cigar_has_op(c, 'M') IS NULL AS match_null, ",
    "cigar_has_op(c, 'S') IS NULL AS clip_null ",
    "FROM (VALUES ", invalid_rows, ") AS invalid(c)"
  ))
  expect_equal(nrow(invalid), length(invalid_cigars))
  for (column in names(invalid)) {
    expect_true(all(invalid[[column]]), info = column)
  }

  validation_edges <- DBI::dbGetQuery(con, paste(
    "SELECT",
    "cigar_query_length('9223372036854775807M')::VARCHAR AS query_limit,",
    "cigar_reference_length('9223372036854775807M')::VARCHAR AS ref_limit,",
    "cigar_aligned_query_length('9223372036854775807M')::VARCHAR AS aligned_limit,",
    "cigar_query_length('9223372036854775807H9223372036854775807P') AS unconsumed,",
    "cigar_query_length('000000000000000000000000000000000000000001M') AS leading_zeroes,",
    "cigar_left_soft_clip('5H7S3M2S4H') AS hard_first,",
    "cigar_right_soft_clip('5H7S3M2S4H') AS hard_last,",
    "cigar_has_op('1M1D', 'd') AS lower_text,",
    "cigar_has_op([16, 18]::UINTEGER[], 'd') AS lower_packed,",
    "cigar_has_op('1M', 'Q') IS NULL AS invalid_requested_op,",
    "cigar_has_op('*', 'M') AS star_presence,",
    "cigar_has_op('', 'M') AS empty_text_presence,",
    "cigar_has_op([]::UINTEGER[], 'M') AS empty_packed_presence,",
    "cigar_has_op(NULL::VARCHAR, 'M') AS null_presence"
  ))
  expect_equal(validation_edges$query_limit[[1]], "9223372036854775807")
  expect_equal(validation_edges$ref_limit[[1]], "9223372036854775807")
  expect_equal(validation_edges$aligned_limit[[1]], "9223372036854775807")
  expect_equal(validation_edges$unconsumed[[1]], 0)
  expect_equal(validation_edges$leading_zeroes[[1]], 1)
  expect_equal(validation_edges$hard_first[[1]], 0)
  expect_equal(validation_edges$hard_last[[1]], 0)
  expect_true(validation_edges$lower_text[[1]])
  expect_true(validation_edges$lower_packed[[1]])
  expect_true(validation_edges$invalid_requested_op[[1]])
  expect_false(validation_edges$star_presence[[1]])
  expect_false(validation_edges$empty_text_presence[[1]])
  expect_false(validation_edges$empty_packed_presence[[1]])
  expect_true(is.na(validation_edges$null_presence[[1]]))

  invalid_packed <- DBI::dbGetQuery(con, paste(
    "SELECT cigar_query_length(c) IS NULL AS query_null,",
    "cigar_has_op(c, 'M') IS NULL AS presence_null",
    "FROM (VALUES ([16, 0]::UINTEGER[]), ([16, 25]::UINTEGER[]),",
    "([16, 31]::UINTEGER[]), ([16, NULL]::UINTEGER[])) invalid(c)"
  ))
  expect_equal(nrow(invalid_packed), 4L)
  expect_true(all(invalid_packed$query_null))
  expect_true(all(invalid_packed$presence_null))
  operators <- dbGetQuery(con, paste(
    "SELECT cigar_has_op('1M1I1D1N1S1H1P1=1X', op) AS text,",
    "cigar_has_op('1M1I1D1N1S1H1P1=1X', op, TRUE) AS strict_text,",
    "cigar_has_op([16,17,18,19,20,21,22,23,24]::UINTEGER[], op) AS packed,",
    "cigar_has_op([16,17,18,19,20,21,22,23,24]::UINTEGER[], op, TRUE) AS strict_packed",
    "FROM (VALUES ('m'), ('i'), ('d'), ('n'), ('s'), ('h'), ('p'), ('='), ('x')) t(op)"
  ))
  expect_equal(nrow(operators), 9L)
  expect_true(all(unlist(operators, use.names = FALSE)))
  test_cigar_strict_policy(con, invalid_cigars)
}

test_cigar_strict_policy <- function(con, invalid_cigars) {
  metric_expectations <- list(
    cigar_query_length = 26, cigar_reference_length = 22,
    cigar_aligned_query_length = 17, cigar_left_soft_clip = 0,
    cigar_right_soft_clip = 0, cigar_has_soft_clip = TRUE,
    cigar_has_hard_clip = TRUE
  )
  valid_inputs <- c(
    "'2H3S5M2I1P4=1X2D3N7M4S1H'",
    "[37,52,80,33,22,71,24,34,51,112,68,21]::UINTEGER[]"
  )
  for (fun in names(metric_expectations)) {
    for (input in valid_inputs) {
      observed <- dbGetQuery(con, sprintf(
        "SELECT %s(%s) AS omitted, %s(%s, FALSE) AS loose, %s(%s, TRUE) AS strict",
        fun, input, fun, input, fun, input
      ))
      expected <- metric_expectations[[fun]]
      expect_equal(observed$omitted[[1]], expected, info = paste(fun, input, "omitted"))
      expect_equal(observed$loose[[1]], expected, info = paste(fun, input, "FALSE"))
      expect_equal(observed$strict[[1]], expected, info = paste(fun, input, "TRUE"))
    }
    missing <- dbGetQuery(con, sprintf(paste(
      "SELECT %s(c, TRUE) IS NULL AS missing",
      "FROM (VALUES (NULL::VARCHAR), (''), ('*')) t(c)"
    ), fun))
    expect_equal(missing$missing, rep(TRUE, 3L), info = fun)
    missing <- dbGetQuery(con, sprintf(paste(
      "SELECT %s(c, TRUE) IS NULL AS missing",
      "FROM (VALUES (NULL::UINTEGER[]), ([]::UINTEGER[])) t(c)"
    ), fun))
    expect_equal(missing$missing, rep(TRUE, 2L), info = fun)
  }

  calls <- c(
    setNames(paste0(names(metric_expectations), "(%s%s)"), names(metric_expectations)),
    cigar_has_op = "cigar_has_op(%s, 'M'%s)",
    cigar_aligned_blocks = "cigar_aligned_blocks(%s, 1%s)"
  )
  invalid_inputs <- c(
    as.character(dbQuoteString(con, invalid_cigars)),
    "'1M' || chr(0) || '1M'",
    "[16,0]::UINTEGER[]", "[16,25]::UINTEGER[]", "[16,31]::UINTEGER[]",
    "[16,NULL]::UINTEGER[]", "[0]::UINTEGER[]", "[4294967295]::UINTEGER[]"
  )
  for (fun in names(calls)) {
    call <- calls[[fun]]
    for (input in invalid_inputs) {
      observed <- dbGetQuery(con, paste0(
        "SELECT ", sprintf(call, input, ""), " IS NULL AS omitted, ",
        sprintf(call, input, ", FALSE"), " IS NULL AS loose"
      ))
      expect_true(all(unlist(observed, use.names = FALSE)), info = paste(fun, input))
      expect_error(dbGetQuery(con, paste("SELECT", sprintf(call, input, ", TRUE"))),
                   paste0(fun, ":"), info = input)
    }
    for (input in c("'0M'", "[NULL,16]::UINTEGER[]")) {
      observed <- dbGetQuery(con, paste(
        "SELECT", sprintf(call, input, ", NULL::BOOLEAN"), "IS NULL AS null_flag"
      ))
      expect_true(observed$null_flag[[1]], info = fun)
    }
    expect_error(dbGetQuery(con, paste(
      "SELECT", sprintf(call, "[16,0]::UINTEGER[]", ", TRUE")
    )), "op index 2", info = fun)
    expect_error(dbGetQuery(con, paste(
      "SELECT", sprintf(call, "[16,16,NULL]::UINTEGER[]", ", TRUE")
    )), "op index 3", info = fun)
  }

  expect_error(dbGetQuery(con, "SELECT cigar_query_length('10M2Q', TRUE)"),
               "operation-start byte offset 4($|[^0-9])")
  expect_error(dbGetQuery(con,
    "SELECT cigar_query_length('9223372036854775807I1I', TRUE)"
  ), "operation-start byte offset 21($|[^0-9])")
  expect_error(dbGetQuery(con,
    "SELECT cigar_aligned_blocks('10M2M', 9223372036854775796, TRUE)"
  ), "operation-start byte offset 4($|[^0-9])")

  # Strict errors identify invalid arguments; top-level SQL NULL takes precedence.
  expect_error(dbGetQuery(con, "SELECT cigar_has_op('1M', 'Q', TRUE)"), "cigar_has_op:")
  expect_error(dbGetQuery(con, "SELECT cigar_has_op([16]::UINTEGER[], '', TRUE)"),
               "cigar_has_op:")
  expect_error(dbGetQuery(con, "SELECT cigar_has_op('1M', 'MM', TRUE)"), "cigar_has_op:")
  expect_error(dbGetQuery(con,
    "SELECT cigar_aligned_blocks('1M', 9223372036854775807, TRUE)"
  ), "cigar_aligned_blocks:")
  null_args <- dbGetQuery(con, paste(
    "SELECT cigar_has_op(NULL::VARCHAR, 'Q', TRUE) IS NULL AS null_cigar,",
    "cigar_has_op('0M', NULL::VARCHAR, TRUE) IS NULL AS null_op,",
    "cigar_aligned_blocks('0M', NULL::BIGINT, TRUE) IS NULL AS null_pos,",
    "cigar_has_op('*', 'M', TRUE) AS missing_presence"
  ))
  expect_true(null_args$null_cigar[[1]])
  expect_true(null_args$null_op[[1]])
  expect_true(null_args$null_pos[[1]])
  expect_false(null_args$missing_presence[[1]])
  for (input in valid_inputs) {
    geometry <- dbGetQuery(con, sprintf(paste(
      "SELECT (b).ref_start::VARCHAR AS ref_start,",
      "(b).query_start::VARCHAR AS query_start, (b).width::VARCHAR AS width",
      "FROM (SELECT cigar_aligned_blocks(%s, -5, TRUE) AS b)"
    ), input))
    expect_equal(geometry$ref_start[[1]], "[-5, 0, 4, 10]")
    expect_equal(geometry$query_start[[1]], "[3, 10, 14, 15]")
    expect_equal(geometry$width[[1]], "[5, 4, 1, 7]")
  }
  mixed <- dbGetQuery(con, paste(
    "SELECT count(*) AS rows,",
    "count(*) FILTER (WHERE i % 4 = 0 AND q IS NULL AND b IS NULL) AS invalid,",
    "count(*) FILTER (WHERE i % 4 = 2 AND q IS NULL AND b IS NULL) AS null_flag,",
    "count(*) FILTER (WHERE i % 2 = 1 AND q = 10 AND",
    "b.ref_start = [1,8] AND b.query_start = [0,5] AND b.width = [5,5]) AS valid",
    "FROM (SELECT i, cigar_query_length(c, strict) AS q,",
    "cigar_aligned_blocks(c, 1, strict) AS b FROM (",
    "SELECT i, CASE WHEN i % 4 = 0 THEN [16,0]::UINTEGER[]",
    "WHEN i % 4 = 2 THEN [NULL,16]::UINTEGER[] ELSE [80,34,80]::UINTEGER[] END AS c,",
    "CASE WHEN i % 4 = 2 THEN NULL ELSE i % 4 = 1 END AS strict",
    "FROM range(5000) t(i)))"
  ))
  expect_equal(mixed$rows[[1]], 5000)
  expect_equal(mixed$invalid[[1]], 1250)
  expect_equal(mixed$null_flag[[1]], 1250)
  expect_equal(mixed$valid[[1]], 2500)
}

test_cigar_geometry_and_flags <- function() {
  con <- rduckhts_connect()
  on.exit(try(dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)

  # binary CIGAR overload (UINTEGER[], oplen<<4|op): '5S90M5I' = [84, 1440, 81],
  # '90M5D' = [1440, 82]; bit-identical to the text path, empty -> NULL like '*'
  bin_metrics <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT",
      "cigar_query_length([84, 1440, 81]::UINTEGER[]) AS q,",
      "cigar_query_length([84, 1440, 81]::UINTEGER[]) = cigar_query_length('5S90M5I') AS q_eq,",
      "cigar_left_soft_clip([84, 1440, 81]::UINTEGER[]) AS lsc,",
      "cigar_reference_length([1440, 82]::UINTEGER[]) = cigar_reference_length('90M5D') AS ref_eq,",
      "cigar_has_op([84, 1440, 81]::UINTEGER[], 'S') = cigar_has_op('5S90M5I', 'S') AS hop_eq,",
      "cigar_query_length([]::UINTEGER[]) IS NULL AS empty_null"
    )
  )
  expect_equal(bin_metrics$q[[1]], 100)
  expect_true(isTRUE(bin_metrics$q_eq[[1]]))
  expect_equal(bin_metrics$lsc[[1]], 5)
  expect_true(isTRUE(bin_metrics$ref_eq[[1]]))
  expect_true(isTRUE(bin_metrics$hop_eq[[1]]))
  expect_true(isTRUE(bin_metrics$empty_null[[1]]))

  # cigar_aligned_blocks: one block per M/=/X op; ref_start carries pos's base,
  # query_start is the 0-based offset into the stored SEQ. Binary overload is
  # bit-identical; invalid input is NULL, a CIGAR with no aligned op is empty lists.
  blocks <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT",
      "(cigar_aligned_blocks('5S90M5S', 100)).ref_start::VARCHAR AS ref_start,",
      "(cigar_aligned_blocks('5S90M5S', 100)).query_start::VARCHAR AS query_start,",
      "(cigar_aligned_blocks('5S90M5S', 100)).width::VARCHAR AS width,",
      "(cigar_aligned_blocks('10M2I10M', 1)).ref_start::VARCHAR AS ins_ref,",
      "(cigar_aligned_blocks('10M2I10M', 1)).query_start::VARCHAR AS ins_query,",
      "cigar_aligned_blocks([84, 1440, 84]::UINTEGER[], 100) = cigar_aligned_blocks('5S90M5S', 100) AS bin_eq,",
      "cigar_aligned_blocks('*', 1) IS NULL AS star_null,",
      "cigar_aligned_blocks([]::UINTEGER[], 1) IS NULL AS empty_null,",
      "len((cigar_aligned_blocks('5S', 1)).width) AS clip_only_blocks"
    )
  )
  expect_equal(blocks$ref_start[[1]], "[100]")
  expect_equal(blocks$query_start[[1]], "[5]")
  expect_equal(blocks$width[[1]], "[90]")
  expect_equal(blocks$ins_ref[[1]], "[1, 11]")
  expect_equal(blocks$ins_query[[1]], "[0, 12]")
  expect_true(isTRUE(blocks$bin_eq[[1]]))
  expect_true(isTRUE(blocks$star_null[[1]]))
  expect_true(isTRUE(blocks$empty_null[[1]]))
  expect_equal(blocks$clip_only_blocks[[1]], 0)

  # Running spans are checked as they accumulate, trailing digits are NULL even
  # when zero-valued, and a rejected row does not disturb the next row's blocks.
  edge <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT",
      "cigar_aligned_blocks('9223372036854775807I1I1M', 1) IS NULL AS q_wrap_null,",
      "cigar_aligned_blocks('9223372036854775807D9223372036854775807D2D1M', 1) IS NULL AS r_wrap_null,",
      "(cigar_aligned_blocks('9223372036854775807M', 0)).width::VARCHAR AS limit_width,",
      "cigar_aligned_blocks('10M0', 1) IS NULL AS trailing_zero_null,",
      "cigar_aligned_blocks('10M00', 1) IS NULL AS trailing_zeros_null"
    )
  )
  expect_true(isTRUE(edge$q_wrap_null[[1]]))
  expect_true(isTRUE(edge$r_wrap_null[[1]]))
  expect_equal(edge$limit_width[[1]], "[9223372036854775807]")
  expect_true(isTRUE(edge$trailing_zero_null[[1]]))
  expect_true(isTRUE(edge$trailing_zeros_null[[1]]))
  rollback <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT (b IS NULL) AS is_null, (b).ref_start::VARCHAR AS ref_start",
      "FROM (SELECT cigar_aligned_blocks(c, 1) AS b, i",
      "      FROM (VALUES ('10M0', 1), ('5M2D5M', 2)) AS t(c, i) ORDER BY i)"
    )
  )
  expect_equal(rollback$is_null, c(TRUE, FALSE))
  expect_equal(rollback$ref_start[[2]], "[1, 8]")

  # On the bundled alignments the block count equals the number of M/=/X ops.
  bam_path <- system.file("extdata", "nanopore.bam", package = "Rduckhts")
  if (nzchar(bam_path)) {
    block_counts <- DBI::dbGetQuery(
      con,
      paste0(
        "SELECT count(*) AS reads, ",
        "count(*) FILTER (WHERE len((cigar_aligned_blocks(CIGAR, POS)).width) = ",
        "len(list_filter(CIGAR, lambda c: (c & 15) IN (0, 7, 8)))) AS agree ",
        "FROM read_bam(", as.character(DBI::dbQuoteString(con, bam_path)),
        ", cigar_representation := 'binary')"
      )
    )
    expect_equal(block_counts$agree[[1]], block_counts$reads[[1]])
  }

  # seq_hash_2bit nt16 overload (UTINYINT[]) is bit-identical; non-ACGT -> NULL
  hash_nt16 <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT",
      "seq_hash_2bit([1, 2, 4, 8]::UTINYINT[]) = seq_hash_2bit('ACGT') AS eq,",
      "seq_hash_2bit([1, 2, 15, 8]::UTINYINT[]) IS NULL AS n_null"
    )
  )
  expect_true(isTRUE(hash_nt16$eq[[1]]))
  expect_true(isTRUE(hash_nt16$n_null[[1]]))

  flag_bits <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT",
      "(sam_flag_bits(99)).is_paired AS is_paired,",
      "(sam_flag_bits(99)).is_proper_pair AS is_proper_pair,",
      "sam_flag_has(99, 2) AS has_proper_pair"
    )
  )

  expect_true(isTRUE(flag_bits$is_paired[[1]]))
  expect_true(isTRUE(flag_bits$is_proper_pair[[1]]))
  expect_true(isTRUE(flag_bits$has_proper_pair[[1]]))

  forward_bits <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT",
      "is_forward_aligned(0) AS forward_zero,",
      "is_forward_aligned(16) AS reverse_sixteen,",
      "is_forward_aligned(4) AS unmapped_four"
    )
  )

  expect_true(isTRUE(forward_bits$forward_zero[[1]]))
  expect_false(isTRUE(forward_bits$reverse_sixteen[[1]]))
  expect_true(is.na(forward_bits$unmapped_four[[1]]))
}

test_cigar_utils()
test_cigar_geometry_and_flags()
