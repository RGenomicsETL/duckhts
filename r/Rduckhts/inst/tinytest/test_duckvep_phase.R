library(tinytest)
library(DBI)

test_phase_preparation <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  dbExecute(con, "SET threads=4")
  path <- system.file("extdata", "geno_vcf44.vcf", package = "Rduckhts")
  stopifnot(nzchar(path))
  source <- paste0("(SELECT record_index, unnest(calls) c FROM read_geno(",
                   dbQuoteString(con, path), "))")
  actual <- dbGetQuery(con, paste0(
    "SELECT record_index::INTEGER r, c.sample_index::INTEGER s, ",
    "a.input_slot::INTEGER slot, a.haplotype_lane::INTEGER lane, a.phase_scope, a.status ",
    "FROM (SELECT record_index, c, unnest(duckvep_phase_call(c.alleles, c.phase_before)) a FROM ",
    source, ") ORDER BY r,s,slot"))
  expect_equal(actual, data.frame(
    r = c(rep(0L, 4), rep(1L, 2), rep(2L, 4), rep(3L, 6)),
    s = c(0L,0L,1L,1L,0L,1L,0L,0L,1L,1L,0L,0L,0L,1L,1L,1L),
    slot = c(1L,2L,1L,2L,1L,1L,1L,2L,1L,2L,1L,2L,3L,1L,2L,3L),
    lane = c(1L,2L,1L,2L,1L,1L,1L,2L,1L,2L,1L,2L,3L,1L,2L,3L),
    phase_scope = c(rep("phase_set", 4), rep("all_phase_sets", 2), rep("phase_set", 10)),
    status = c(rep("called", 6), "missing", "called", "missing", rep("called", 3),
               "missing", rep("called", 3))
  ))
  partial_path <- system.file("extdata", "geno_phase_partial.vcf", package = "Rduckhts")
  stopifnot(nzchar(partial_path))
  partial <- dbGetQuery(con, paste0(
    "SELECT record_index::INTEGER r,c,duckvep_phase_call(c.alleles,c.phase_before,",
    "phase_set:=c.phase_set) assignments FROM (SELECT record_index,calls[1] c FROM read_geno(",
    dbQuoteString(con, partial_path), ")) ORDER BY r"))
  expect_equal(partial$r, 0:11)
  expected_gt <- list(0:2,0:2,0:2,c(1L,0L,2L),c(1L,0L,2L),c(1L,0L,2L),
    c(NA_integer_,1L,2L),c(NA_integer_,1L,2L),c(NA_integer_,1L,2L),0:2,0:2,0:2)
  expected_flags <- list(c(FALSE,TRUE,FALSE),c(FALSE,TRUE,FALSE),c(TRUE,TRUE,FALSE),
    c(FALSE,TRUE,FALSE),c(FALSE,TRUE,FALSE),c(TRUE,TRUE,FALSE),c(FALSE,TRUE,FALSE),
    c(FALSE,TRUE,FALSE),c(TRUE,TRUE,FALSE),c(TRUE,FALSE,FALSE),c(FALSE,FALSE,TRUE),c(FALSE,FALSE,TRUE))
  expected_lanes <- list(c(NA,2L,NA),c(NA,2L,NA),1:3,c(NA,2L,NA),c(NA,2L,NA),1:3,
    c(NA,2L,NA),c(NA,2L,NA),1:3,c(1L,NA,NA),c(NA,NA,3L),c(NA,NA,3L))
  for (i in seq_len(nrow(partial))) {
    expect_equal(partial$c$alleles[[i]], expected_gt[[i]])
    expect_equal(partial$c$phase_before[[i]], expected_flags[[i]])
    slots <- partial$assignments[[i]]
    expect_equal(slots$input_slot, 1:3)
    expect_equal(slots$haplotype_lane, expected_lanes[[i]])
    expect_equal(slots$phase_set, ifelse(is.na(expected_lanes[[i]]), NA_real_, 10))
    expect_equal(slots$phase_scope, ifelse(is.na(expected_lanes[[i]]), "unresolved", "phase_set"))
    expect_equal(slots$status, ifelse(is.na(expected_gt[[i]]), "missing",
      ifelse(is.na(expected_lanes[[i]]), "unphased", "called")))
  }
  prepared <- function(gt, phase, policy = "strict", ps = "NULL") {
    dbGetQuery(con, paste0("SELECT a.input_slot::INTEGER slot, a.allele_index allele, ",
      "a.haplotype_lane::INTEGER lane, a.phase_set::VARCHAR ps, a.phase_scope, a.status ",
      "FROM (SELECT unnest(duckvep_phase_call(", gt, ",", phase, ",phase_set := ", ps,
      ",phase_policy := ", dbQuoteString(con, policy), ")) a) ORDER BY slot"))
  }
  ambiguous <- prepared("[0,1,NULL]", "[true,false,false]", ps = "10")
  expect_equal(ambiguous$allele, c(0L, 1L, NA_integer_))
  expect_equal(ambiguous$lane, c(1L, NA_integer_, NA_integer_))
  expect_equal(ambiguous$status, c("called", "unphased", "missing"))
  expect_equal(ambiguous$ps, c("10", NA, NA))
  compatible <- prepared("[0,1,NULL]", "[false,false,false]", "vep116_compat", "42")
  expect_equal(compatible$lane, c(1L, 2L, NA_integer_))
  expect_equal(compatible$phase_scope, rep("allele_slot", 3))
  expect_equal(compatible$status, c("called", "called", "missing"))
  expect_true(all(is.na(compatible$ps)))
  expect_equal(prepared("[NULL,1,NULL,2]", "NULL", "vep116_compat")$lane,
               c(NA_integer_, 1L, NA_integer_, 2L))
  expect_equal(prepared("[2,2,2]", "NULL", ps = "99")$phase_scope, rep("all_phase_sets", 3))
  expect_equal(prepared("[2,1,1]", "[true,false,false]", ps = "-1")$ps, rep("-1", 3))
  expect_true(all(is.na(prepared("[0,1]", "[NULL,NULL]")$lane)))
  expect_equal(prepared("[0,1]", "[NULL,NULL]"), data.frame(slot = 1:2, allele = 0:1,
    lane = NA_integer_, ps = NA_character_, phase_scope = "unresolved", status = "unphased"))
  expect_equal(prepared("[NULL,NULL]", "[true,true]"), data.frame(slot = 1:2,
    allele = NA_integer_, lane = 1:2, ps = NA_character_, phase_scope = "phase_set", status = "missing"))
  null_slots <- prepared("list_resize([NULL], 65535)", "NULL")
  expect_equal(nrow(null_slots), 65535L)
  expect_true(all(is.na(null_slots$allele)) && all(is.na(null_slots$lane)) &&
    all(null_slots$status == "missing"))
  expect_true(dbGetQuery(con, "SELECT duckvep_phase_call(NULL,NULL) IS NULL ok")$ok)
  expect_equal(dbGetQuery(con, paste(
    "SELECT count(*) n, count(*) FILTER(WHERE a.status='unphased') ambiguous",
    "FROM (SELECT unnest(duckvep_phase_call(CASE WHEN i%2=0 THEN [0,1] ELSE [1,1,1] END, NULL)) a",
    "FROM range(5000) r(i))")), data.frame(n = 12500, ambiguous = 5000))
  # Verify every row and field across changing list offsets and policy, not
  # only output cardinality. PS=0 must remain distinct from the missing set.
  many <- dbGetQuery(con, paste(
    "SELECT i::INTEGER i, a.input_slot::INTEGER slot, a.allele_index allele,",
    "a.haplotype_lane::INTEGER lane, a.phase_set::VARCHAR ps, a.phase_scope, a.status",
    "FROM (SELECT i, unnest(duckvep_phase_call(",
    "CASE WHEN i%3=2 THEN [NULL,2,0] ELSE [0,1] END,",
    "CASE WHEN i%3=2 THEN NULL WHEN i%3=0 THEN [true,true] ELSE [false,false] END,",
    "phase_set := i, phase_policy := CASE WHEN i%3=2 THEN 'vep116_compat' ELSE 'strict' END)) a",
    "FROM range(5000) r(i)) ORDER BY i,slot"))
  expected_many <- do.call(rbind, lapply(0:4999, function(i) {
    if (i %% 3L == 0L) {
      data.frame(i, slot = 1:2, allele = 0:1, lane = 1:2, ps = as.character(i),
                 phase_scope = "phase_set", status = "called")
    } else if (i %% 3L == 1L) {
      data.frame(i, slot = 1:2, allele = 0:1, lane = NA_integer_, ps = NA_character_,
                 phase_scope = "unresolved", status = "unphased")
    } else {
      data.frame(i, slot = 1:3, allele = c(NA_integer_,2L,0L), lane = c(NA_integer_,1L,2L),
                 ps = NA_character_, phase_scope = "allele_slot", status = c("missing","called","called"))
    }
  }))
  expect_equal(many, expected_many)
  expect_equal(nrow(prepared("list_transform(range(65535), x -> 1)", "NULL")), 65535L)
  expect_error(prepared("[]", "[]"), pattern = "ploidy")
  expect_error(prepared("[0,1]", "[true]"), pattern = "equal length")
  expect_error(prepared("[0,-1]", "NULL"), pattern = "non-negative")
  expect_error(prepared("[0,1]", "NULL", "guess"), pattern = "phase_policy")
  expect_error(prepared("list_transform(range(65536), x -> 1)", "NULL"), pattern = "ploidy")
  expect_equal(dbGetQuery(con, "SELECT 42 n")$n, 42L)
}

test_phase_preparation()
