#!/usr/bin/env Rscript
# Additional alignment lane over every complete sequence pair in an existing
# Haplosaurus receipt. No generator, selection rule or prior oracle is changed.
main <- function() {
  opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--artifacts", type = "character"),
    optparse::make_option("--public-artifacts", dest = "public_artifacts", type = "character"),
    optparse::make_option("--vep-prefix", dest = "vep_prefix",
      default = Sys.getenv("VEP_PREFIX", "/root/miniconda3/envs/vep"))
  )))
  source("scripts/duckvep_evidence.R", local = TRUE)
  root <- normalizePath(".")
  revision <- duckvep_evidence_revision(root)
  artifact <- normalizePath(opt$artifacts, mustWork = TRUE)
  public <- normalizePath(opt$public_artifacts, mustWork = TRUE)
  original <- jsonlite::read_json(file.path(artifact, "receipt.json"), simplifyVector = TRUE)
  public_receipt <- jsonlite::read_json(file.path(public, "receipt.json"), simplifyVector = TRUE)
  stopifnot(identical(public_receipt$source_receipt, file.path(artifact, "receipt.json")))
  original_paths <- file.path(artifact, c("inputs.rds", "native.rds", "summary.csv"))
  public_paths <- file.path(public, c("strict.rds", "vep116_compat.rds", "summary.csv"))
  for (p in original_paths) stopifnot(identical(duckvep_evidence_sha256(p), original$sha256[[p]]))
  for (p in public_paths) stopifnot(identical(duckvep_evidence_sha256(p), public_receipt$sha256[[p]]))
  stopifnot(read.csv(original_paths[3L])$failures == 0L,
    all(read.csv(public_paths[3L])$oracle_failures == 0L))
  mirrors <- setNames(normalizePath(c(".sync/ensembl-vep", ".sync/ensembl-variation")),
    c("vep", "variation"))
  revisions <- c(vep = "57ea5c52340acc1f156267f810ad162e26597082",
    variation = "2fb834b987ede3824e200197a838ce11e91aeb4b")
  for (name in names(revisions)) {
    stopifnot(identical(original$oracle_revisions[[name]], revisions[[name]]),
      identical(system2("git", c("-C", mirrors[[name]], "rev-parse", "HEAD"), stdout = TRUE), revisions[[name]]),
      !length(system2("git", c("-C", mirrors[[name]], "status", "--porcelain"), stdout = TRUE)))
  }
  out <- tempfile(paste0("sequence_diff_seed", original$seed, "_"), tmpdir = dirname(artifact))
  dir.create(out)
  message("Sequence-difference artifacts: ", out)
  sources <- c("test/duckvep/conformance/sequence_diff_probe.c",
    "src/duckvep/kernel/src/duckvep_sequence_diff.c")
  code <- c(sources, "src/duckvep/kernel/src/duckvep_sequence_diff.h",
    "test/duckvep/conformance/sequence_diff_oracle.pl", "test/duckvep/conformance/sequence_diff_differential.R")
  code_hashes <- vapply(code, duckvep_evidence_sha256, "")
  binding <- if (length(duckvep_evidence_tracked_changes(root)) ||
    length(duckvep_evidence_untracked_build_inputs(root))) "diagnostic_unbound" else "clean_checkout"
  prefix <- normalizePath(opt$vep_prefix, mustWork = TRUE)
  conda <- function(...) do.call(blit::conda,
    c(as.list(duckvep_blit_quote(c(...))), list(conda = duckvep_blit_quote("micromamba"))))
  environment <- file.path(out, "environment.txt")
  stopifnot(blit::cmd_run(conda("list", "-p", prefix, "--explicit"),
    stdout = environment, stderr = "2>&1", verbose = FALSE) == 0L,
    identical(duckvep_evidence_explicit_packages(readLines(environment)),
      duckvep_evidence_explicit_packages(readLines(
        "test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt"))))
  cases <- readRDS(original_paths[1L])$cases
  native <- readRDS(original_paths[2L])
  pairs <- vector("list", length(cases) * 12L)
  at <- 0L
  for (case in cases) {
    paths <- native[[case$transcript]]
    stopifnot(length(paths) == 6L)
    reference_paths <- which(!lengths(lapply(paths, `[[`, "contributors")))
    stopifnot(length(reference_paths) > 0L)
    reference <- paths[[reference_paths[[1L]]]]
    stopifnot(identical(reference$cds, case$cds))
    for (lane in seq_along(paths)) for (axis in c("cds", "protein")) {
      at <- at + 1L
      pairs[[at]] <- list(id = paste(case$transcript, lane, axis, sep = "/"),
        reference = reference[[axis]], alternate = paths[[lane]][[axis]],
        align = as.integer("indel" %in% paths[[lane]]$flags))
    }
  }
  stopifnot(at == length(pairs))
  pair_file <- file.path(out, "pairs.jsonl")
  writeLines(vapply(pairs, jsonlite::toJSON, "", auto_unbox = TRUE), pair_file)
  perl_lib <- paste(c(file.path(mirrors, "modules"),
    file.path(prefix, "share/ensembl-vep-116.0-0")), collapse = ":")
  oracle_file <- file.path(out, "oracle.jsonl")
  run <- conda("run", "--clean-env", "--env", paste0("PERL5LIB=", perl_lib), "-p", prefix,
    "perl", file.path(root, "test/duckvep/conformance/sequence_diff_oracle.pl"), pair_file)
  stopifnot(blit::cmd_run(run, stdout = oracle_file, stderr = file.path(out, "oracle.log"),
    stdin = NULL, verbose = FALSE) == 0L)
  oracle <- lapply(readLines(oracle_file), jsonlite::fromJSON)
  pair_ids <- vapply(pairs, `[[`, "", "id")
  stopifnot(identical(vapply(oracle, `[[`, "", "id"), pair_ids), !anyDuplicated(pair_ids))
  shared <- file.path(out, paste0("sequence_diff_probe", .Platform$dynlib.ext))
  stopifnot(system2(Sys.getenv("CC", "cc"), shQuote(c("-std=c99", "-O1", "-g", "-Wall", "-Wextra",
    "-fPIC", "-shared", "-I", "src/duckvep/kernel/src", sources, "-o", shared)),
    stdout = file.path(out, "compiler.log"), stderr = file.path(out, "compiler.log")) == 0L)
  dll <- dyn.load(shared)
  on.exit(dyn.unload(shared), add = TRUE)
  symbol <- getNativeSymbolInfo("duckhts_test_sequence_differences", PACKAGE = dll)
  canonical <- function(x) data.frame(ref_start0 = as.numeric(x$ref_start0),
    alt_start0 = as.numeric(x$alt_start0), reference = as.character(x$reference),
    alternate = as.character(x$alternate), alignment_start0 = as.numeric(x$alignment_start0))
  expected <- setNames(lapply(oracle, function(x) canonical(x$differences)), pair_ids)
  actual <- setNames(lapply(pairs, function(p) {
    capacity <- as.integer(nchar(p$reference) + nchar(p$alternate) + 1L)
    x <- .C(symbol, p$reference, p$alternate, as.integer(p$align), spans = double(5L * capacity),
      capacity, status = integer(1L), facts = double(3L))
    stopifnot(x$status == 0L)
    spans <- matrix(x$spans, ncol = 5L, byrow = TRUE)[seq_len(x$facts[1L]), , drop = FALSE]
    slice <- function(text, start, length) ifelse(length == 0, "", substring(text, start + 1, start + length))
    canonical(data.frame(ref_start0 = spans[, 1L], alt_start0 = spans[, 2L],
      reference = slice(p$reference, spans[, 1L], spans[, 3L]),
      alternate = slice(p$alternate, spans[, 2L], spans[, 4L]), alignment_start0 = spans[, 5L]))
  }), pair_ids)
  equal <- function(x) identical(x, expected)
  failures <- !vapply(seq_along(actual), function(i) identical(actual[[i]], expected[[i]]), TRUE)
  saveRDS(list(actual = actual, expected = expected, failed_pairs = pair_ids[failures]), file.path(out, "pairs.rds"))
  controls <- c(missing_pair = !equal(actual[-1L]), duplicate_pair = !equal(c(actual, actual[1L])))
  first <- which(lengths(lapply(actual, row.names)) > 0L)[1L]
  for (field in names(actual[[first]])) {
    changed <- actual
    value <- changed[[first]][[field]][1L]
    changed[[first]][[field]][1L] <- if (is.numeric(value)) value + 1 else paste0(value, "X")
    controls[field] <- !equal(changed)
  }
  changed <- actual; changed[[first]] <- changed[[first]][-1L, , drop = FALSE]
  controls["missing_run"] <- !equal(changed)
  changed <- actual; changed[[first]] <- rbind(changed[[first]], changed[[first]][1L, ])
  controls["extra_run"] <- !equal(changed)
  public_checks <- list()
  for (policy in c("strict", "vep116_compat")) {
    leaves <- readRDS(file.path(public, paste0(policy, ".rds")))$leaves
    for (axis in c("cds", "protein")) {
      checks <- vapply(seq_len(nrow(leaves)), function(i) {
        case <- cases[[as.integer(leaves$transcript_index[i]) + 1L]]
        paths <- native[[case$transcript]]
        same <- which(vapply(paths, function(p) p$cds == leaves$cds[i] &&
          ("indel" %in% p$flags) == (bitwAnd(as.integer(leaves$sequence_flags[i]), 1L) != 0L), TRUE))
        if (!length(same)) return(FALSE)
        wanted <- expected[[paste(case$transcript, same[1L], axis, sep = "/")]]
        identical(canonical(leaves[[paste0(axis, "_differences")]][[i]]), wanted)
      }, TRUE)
      public_checks[[paste(policy, axis)]] <- data.frame(policy, axis,
        leaves = nrow(leaves), failures = sum(!checks))
    }
  }
  summary <- data.frame(seed = original$seed, transcripts = length(cases), sequence_pairs = length(pairs),
    difference_runs = sum(vapply(expected, nrow, 1L)), native_failures = sum(failures),
    controls_rejected = sum(controls))
  write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
  write.csv(do.call(rbind, public_checks), file.path(out, "public.csv"), row.names = FALSE)
  write.csv(data.frame(control = names(controls), rejected = controls), file.path(out, "controls.csv"), row.names = FALSE)
  identities <- c(original_paths, public_paths, file.path(artifact, "receipt.json"),
    file.path(public, "receipt.json"), code,
    file.path(out, c("pairs.jsonl", "oracle.jsonl", "pairs.rds", "summary.csv", "public.csv", "controls.csv", "environment.txt")))
  jsonlite::write_json(list(source_revision = revision, source_binding = binding, oracle_revisions = revisions,
    scope = "exact_sequence_alignment_given_receipted_pairs_and_public_cds_and_protein_differences_not_hgvs_or_combined_so",
    source_receipt = file.path(artifact, "receipt.json"), public_receipt = file.path(public, "receipt.json"),
    sha256 = as.list(vapply(identities, duckvep_evidence_sha256, ""))), file.path(out, "receipt.json"),
    pretty = TRUE, auto_unbox = TRUE)
  print(summary); print(do.call(rbind, public_checks))
  stopifnot(!any(failures), all(controls), all(vapply(public_checks, function(x) x$failures == 0L, TRUE)),
    identical(code_hashes, vapply(code, duckvep_evidence_sha256, "")),
    identical(revision, duckvep_evidence_revision(root)))
}
main()
