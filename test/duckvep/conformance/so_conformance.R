#!/usr/bin/env Rscript
# so_conformance.R — FORMAL class-model tier of the conformance framework.
#
# The conformance framework's oracle is Ensembl VEP's authoritative consequence
# spec (the %OVERLAP_CONSEQUENCES table in Constants.pm), extracted to
# data/so_consequences.tsv by extract_so_spec.pl. This script proves the C
# engine's emitted SO vocabulary (duckvep_so.h/.c, dumped via so_dump.c) conforms
# to that spec: every term the engine can emit must exist in VEP's spec with the
# SAME IMPACT, severity rank, and evaluator tier. It is the part of the framework
# that runs against the C engine TODAY (the structural SO slice); the
# differential-fuzzer and stratified-CI tiers (generate_witnesses / fuzz_witnesses
# / stratified_conformance) prove per-variant behavior.
#
# Raw agree/total is reported per the project rule; a single mismatch fails and is
# never rounded away.
#
# Usage: Rscript test/duckvep/conformance/so_conformance.R [repository]

suppressPackageStartupMessages(library(glue))

args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(if (length(args) >= 1) args[[1]] else ".")

`%||%` <- function(x, y) if (is.null(x)) y else x
die <- function(...) stop(glue(..., .envir = parent.frame()), call. = FALSE)

run <- function(cmd, args = character()) {
  out <- suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = TRUE))
  list(status = attr(out, "status") %||% 0L, lines = out)
}

base <- file.path(repo, "test", "duckvep", "conformance")
spec_path <- file.path(base, "data", "so_consequences.tsv")
if (!file.exists(spec_path)) {
  die(
    "VEP SO spec missing: {spec_path}. Regenerate with ",
    "`make so-spec` (extract_so_spec.pl over VEP's Constants.pm)."
  )
}
spec <- read.delim(
  spec_path,
  stringsAsFactors = FALSE,
  na.strings = "\\N"
)
# Guard against a stale/corrupt oracle silently passing the check: require the
# spec's structure and an anchor term VEP has always carried.
if (
  !all(c("SO_term", "impact", "rank", "tier") %in% names(spec)) ||
    nrow(spec) < 30
) {
  die(
    "SO spec {spec_path} looks malformed (cols/rows). Regenerate with `make so-spec`."
  )
}
if (!identical(spec$impact[spec$SO_term == "missense_variant"], "MODERATE")) {
  die(
    "SO spec anchor failed: missense_variant is not MODERATE — wrong/stale oracle?"
  )
}

# --- the C engine's emitted SO vocabulary, via the real ABI -------------------
cc <- Sys.getenv("CC", "cc")
dump_bin <- file.path(tempdir(), "duckvep_so_dump")
build <- run(
  cc,
  c(
    "-std=c99",
    "-Wall",
    "-Wextra",
    glue("-I{repo}/src/duckvep/kernel/include"),
    glue("-I{repo}/src/duckvep/kernel/src"),
    glue("{base}/so_dump.c"),
    glue("{repo}/src/duckvep/kernel/src/duckvep_so.c"),
    "-o",
    dump_bin
  )
)
if (build$status != 0) {
  die("so_dump build failed:\n{paste(build$lines, collapse = '\n')}")
}
dump <- run(dump_bin)
if (dump$status != 0) {
  die("so_dump run failed:\n{paste(dump$lines, collapse = '\n')}")
}
eng <- read.delim(
  text = paste(dump$lines, collapse = "\n"),
  stringsAsFactors = FALSE
)

# --- conformance: term exists in the VEP spec with matching metadata ----------
spec_impact <- setNames(spec$impact, spec$SO_term)
spec_rank <- setNames(spec$rank, spec$SO_term)
spec_tier <- setNames(spec$tier, spec$SO_term)
eng$in_spec <- eng$so_term %in% spec$SO_term
eng$spec_impact <- unname(spec_impact[eng$so_term])
eng$spec_rank <- as.integer(unname(spec_rank[eng$so_term]))
eng$spec_tier <- as.integer(unname(spec_tier[eng$so_term]))
eng$impact_ok <- eng$in_spec & (eng$impact == eng$spec_impact)
eng$rank_ok <- eng$in_spec & (eng$rank == eng$spec_rank)
eng$tier_ok <- eng$in_spec & (eng$tier == eng$spec_tier)
eng$metadata_ok <- eng$impact_ok & eng$rank_ok & eng$tier_ok

n <- nrow(eng)
n_in_spec <- sum(eng$in_spec)
n_impact_ok <- sum(eng$impact_ok)
n_rank_ok <- sum(eng$rank_ok)
n_tier_ok <- sum(eng$tier_ok)
n_metadata_ok <- sum(eng$metadata_ok)

out_dir <- file.path(base, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(eng, file.path(out_dir, "so_conformance.csv"), row.names = FALSE)

# Honest scope: VOCABULARY / class-model conformance (does each declared term
# carry VEP's IMPACT, rank, and tier), NOT behavioral conformance (whether the
# annotator emits the right term per variant — the differential-fuzzer tier).
cat(
  glue(
    "
SO vocabulary conformance — C engine's declared SO table vs VEP {nrow(spec)}-term spec
  (class-model tier: term + IMPACT + rank + tier; NOT behavioral correctness)

  declared SO vocabulary terms    : {n}
  present in VEP spec             : {n_in_spec}/{n}
  IMPACT matches VEP spec         : {n_impact_ok}/{n}
  rank matches VEP spec           : {n_rank_ok}/{n}
  tier matches VEP spec           : {n_tier_ok}/{n}
  all metadata matches VEP spec   : {n_metadata_ok}/{n}
  verdict                         : {if (n_metadata_ok == n) 'VOCABULARY-CONFORMANT' else 'DIVERGENT'}
"
  ),
  "\n"
)

bad <- eng[
  !eng$metadata_ok,
  c(
    "so_term",
    "impact",
    "spec_impact",
    "rank",
    "spec_rank",
    "tier",
    "spec_tier",
    "in_spec"
  )
]
if (nrow(bad) > 0) {
  cat("divergences:\n")
  print(bad, row.names = FALSE)
  die(
    "SO class-model conformance FAILED: {n_metadata_ok}/{n} terms match VEP metadata"
  )
}
