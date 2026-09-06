#!/usr/bin/env Rscript
# Exact presentation comparison with pinned VEP JSON. The independent-event
# generator and existing consequence/HGVS campaigns are deliberately unchanged.
duckvep_projection_expected <- function(path, gff) {
  records <- lapply(readLines(path), jsonlite::fromJSON, simplifyVector = FALSE)
  features <- utils::read.delim(gff, header = FALSE, comment.char = "#", quote = "")
  exon_total <- sum(features$V3 == "exon")
  scalar <- function(value) if (is.null(value)) NA_character_ else as.character(value)
  ordinals <- function(value) {
    if (is.null(value)) return(c(NA_character_, NA_character_))
    span <- strsplit(strsplit(value, "/", fixed = TRUE)[[1L]][[1L]], "-", fixed = TRUE)[[1L]]
    c(span[[1L]], tail(span, 1L))
  }
  alleles <- function(value) {
    if (is.null(value)) return(c(NA_character_, NA_character_))
    pair <- strsplit(value, "/", fixed = TRUE)[[1L]]
    stopifnot(length(pair) %in% 1:2)
    c(pair[[1L]], tail(pair, 1L))
  }
  rows <- lapply(records, function(record) {
    # Each fixture has one transcript; never discard unexpected oracle rows.
    stopifnot(length(record$transcript_consequences) == 1L)
    tx <- record$transcript_consequences[[1L]]
    exon <- ordinals(tx$exon)
    intron <- ordinals(tx$intron)
    aa <- alleles(tx$amino_acids)
    codon <- alleles(tx$codons)
    coordinates <- c("cdna_start", "cdna_end", "cds_start", "cds_end", "protein_start", "protein_end")
    c(ID = strsplit(record$input, "\t", fixed = TRUE)[[1L]][[3L]],
      transcript_id = scalar(tx$transcript_id), output_allele = scalar(tx$variant_allele),
      interbase = as.character(record$start > record$end),
      vapply(stats::setNames(coordinates, coordinates),
        function(field) scalar(tx[[field]]), character(1L)),
      exon_first = exon[[1L]], exon_last = exon[[2L]], exon_total = as.character(exon_total),
      intron_first = intron[[1L]], intron_last = intron[[2L]], intron_total = as.character(exon_total - 1L),
      transcript_distance = scalar(tx$distance),
      cds_start_nf = as.character("cds_start_NF" %in% tx$flags),
      cds_end_nf = as.character("cds_end_NF" %in% tx$flags),
      reference_amino_acids = aa[[1L]], alternate_amino_acids = aa[[2L]],
      reference_codons = codon[[1L]], alternate_codons = codon[[2L]])
  })
  as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE, check.names = FALSE)
}

duckvep_projection_consequences_expected <- function(path) {
  records <- lapply(readLines(path), jsonlite::fromJSON, simplifyVector = FALSE)
  rows <- lapply(records, function(record) {
    stopifnot(length(record$transcript_consequences) == 1L)
    tx <- record$transcript_consequences[[1L]]
    data.frame(ID = strsplit(record$input, "\t", fixed = TRUE)[[1L]][[3L]],
      transcript_id = tx$transcript_id,
      consequence_terms = paste(sort(unlist(tx$consequence_terms)), collapse = "&"),
      status = "supported", reason = NA_character_)
  })
  do.call(rbind, rows)
}

duckvep_projection_compare <- function(actual, expected) {
  keys <- c("ID", "transcript_id")
  stopifnot(identical(names(actual), names(expected)))
  key <- function(x) do.call(paste, c(x[keys], sep = "\t"))
  a <- key(actual)
  e <- key(expected)
  problems <- list()
  add <- function(field, id, actual, expected) {
    problems[[length(problems) + 1L]] <<- data.frame(field, ID = id, actual, expected)
  }
  for (side in c("actual", "expected")) {
    x <- if (side == "actual") actual else expected
    k <- if (side == "actual") a else e
    bad <- is.na(x$ID) | is.na(x$transcript_id) | duplicated(k)
    if (any(bad)) add(paste0(side, "_key"), k[bad], "missing or duplicate", "unique non-null key")
  }
  if (length(setdiff(a, e))) add("unexpected_pair", setdiff(a, e), "present", "absent")
  if (length(setdiff(e, a))) add("missing_pair", setdiff(e, a), "absent", "present")
  selected <- match(e, a)
  present <- !is.na(selected)
  for (field in setdiff(names(expected), keys)) {
    observed <- actual[[field]][selected[present]]
    oracle <- expected[[field]][present]
    bad <- is.na(observed) != is.na(oracle) |
      (!is.na(observed) & !is.na(oracle) & observed != oracle)
    if (any(bad)) add(field, e[present][bad], observed[bad], oracle[bad])
  }
  if (!length(problems)) return(data.frame(field = character(), ID = character(),
    actual = character(), expected = character()))
  do.call(rbind, problems)
}

main <- function() {
  opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--random-cases", dest = "random_cases", type = "integer", default = 1000L),
    optparse::make_option("--max-random-length", dest = "max_random_length", type = "integer", default = 10L),
    optparse::make_option("--seed", type = "integer", default = 173L),
    optparse::make_option("--consequences", action = "store_true", default = FALSE,
      help = "also require native SO/status agreement on every retained record/transcript pair"),
    optparse::make_option("--extension", default = "build/release/duckhts.duckdb_extension"),
    optparse::make_option("--extension-receipt", dest = "extension_receipt", default = NULL),
    optparse::make_option("--summary-output", dest = "summary_output", default = NULL),
    optparse::make_option("--vep-prefix", dest = "vep_prefix", default = Sys.getenv("VEP_PREFIX")),
    optparse::make_option("--vep-git", dest = "vep_git", default = ".sync/ensembl-vep"),
    optparse::make_option("--variation-git", dest = "variation_git", default = ".sync/ensembl-variation")
  )))
  stopifnot(!is.na(opt$random_cases), opt$random_cases >= 0L, !is.na(opt$seed),
    !is.na(opt$max_random_length), opt$max_random_length >= 1L, opt$max_random_length <= 65534L)
  root <- normalizePath(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
  source(file.path(root, "scripts/duckvep_evidence.R"), local = TRUE)
  source(file.path(root, "test/duckvep/conformance/projection_fixtures.R"), local = TRUE)
  revision <- duckvep_evidence_revision(root)
  extension <- normalizePath(opt$extension, mustWork = TRUE)
  extension_hash <- duckvep_evidence_sha256(extension)
  binding <- "diagnostic_unbound"
  allowed_outputs <- if (is.null(opt$summary_output)) character() else opt$summary_output
  if (length(allowed_outputs) && is.null(opt$extension_receipt))
    stop("--summary-output requires a source-bound extension receipt")
  if (!is.null(opt$extension_receipt)) {
    duckvep_evidence_assert_checkout(root, revision, allowed_outputs)
    binding <- duckvep_evidence_read_extension_receipt(opt$extension_receipt, root,
      extension, revision)$binding
  }
  results <- file.path(root, "test/duckvep/conformance/results")
  dir.create(results, showWarnings = FALSE)
  directory <- tempfile(paste0("projection_seed", opt$seed, "_"), tmpdir = results)
  dir.create(directory)
  message("Artifacts: ", directory)
  inputs <- duckhtsbench::duckhts_bench_stage_repository_fixtures(root, "duckvep-projection")
  mirrors <- c(vep = normalizePath(opt$vep_git, mustWork = TRUE),
    variation = normalizePath(opt$variation_git, mustWork = TRUE))
  revisions <- c(vep = "57ea5c52340acc1f156267f810ad162e26597082",
    variation = "2fb834b987ede3824e200197a838ce11e91aeb4b")
  command <- function(exe, args) duckvep_evidence_command(exe, args, paste("failed:", exe))
  check_oracle <- function() {
    for (name in names(mirrors)) {
      stopifnot(identical(command("git", c("-C", mirrors[[name]], "rev-parse", "HEAD")),
        revisions[[name]]), !length(command("git", c("-C", mirrors[[name]], "status", "--porcelain"))))
    }
  }
  check_oracle()
  prefix <- normalizePath(opt$vep_prefix, mustWork = TRUE)
  environment <- command("micromamba", c("list", "-p", prefix, "--explicit"))
  lock <- file.path(root, "test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt")
  stopifnot(identical(duckvep_evidence_explicit_packages(environment),
    duckvep_evidence_explicit_packages(readLines(lock))))
  writeLines(environment, file.path(directory, "environment.txt"))
  vep_directory <- file.path(directory, "vep")
  dir.create(vep_directory)
  perl_library <- paste(c(file.path(mirrors, "modules"),
    file.path(prefix, "share/ensembl-vep-116.0-0")), collapse = .Platform$path.sep)
  con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(con, paste("LOAD", DBI::dbQuoteString(con, extension)))
  DBI::dbExecute(con, "SET threads=1")
  summaries <- controls <- consequence_summaries <- list()
  for (case in duckvep_projection_cases) {
    gff <- duckvep_projection_fixture(con, root, inputs, case, directory)
    vcf <- file.path(directory, paste0(case, ".vcf"))
    generated <- file.path(directory, paste0(case, ".generated.vcf"))
    json <- file.path(directory, paste0(case, ".json"))
    command("Rscript", c(file.path(root, "test/duckvep/conformance/generate_witnesses.R"),
      "--gff", gff, "--fasta", inputs[["projection_reference"]], "--ext", extension,
      "--out", generated, "--random-cases", opt$random_cases, "--seed", opt$seed,
      "--max-random-length", opt$max_random_length))
    duckvep_projection_label_records(generated, vcf)
    status <- system2("bgzip", c("-c", shQuote(gff)), stdout = paste0(gff, ".gz"))
    stopifnot(status == 0L)
    command("tabix", c("-p", "gff", paste0(gff, ".gz")))
    log <- command("micromamba", c("run", "--clean-env", "--env", paste0("HOME=", vep_directory),
      "--env", paste0("PERL5LIB=", perl_library), "-p", prefix, "perl", file.path(mirrors[["vep"]], "vep"),
      "-i", vcf, "--fasta", inputs[["projection_reference"]], "--gff", paste0(gff, ".gz"),
      "--dir", vep_directory, "--json", "--numbers", "--distance", "5000", "--no_stats",
      "--force_overwrite", "-o", json))
    writeLines(log, file.path(directory, paste0(case, ".vep.log")))
    duckvep_projection_input(con, vcf)
    expected <- duckvep_projection_expected(json, gff)
    actual <- DBI::dbGetQuery(con, "SELECT * EXCLUDE(event_index, transcript_index) FROM projection_actual")
    actual <- as.data.frame(lapply(actual[names(expected)], as.character), stringsAsFactors = FALSE)
    utils::write.csv(actual, file.path(directory, paste0(case, ".actual.csv")), row.names = FALSE)
    utils::write.csv(expected, file.path(directory, paste0(case, ".expected.csv")), row.names = FALSE)
    input_count <- DBI::dbGetQuery(con, "SELECT count(*) n FROM projection_events")$n
    stopifnot(input_count == nrow(expected))
    mismatch <- duckvep_projection_compare(actual, expected)
    utils::write.csv(mismatch, file.path(directory, paste0(case, ".mismatch.csv")), row.names = FALSE)
    summaries[[case]] <- data.frame(case, seed = opt$seed, max_random_length = opt$max_random_length,
      input_alleles = input_count,
      expected_pairs = nrow(expected), actual_pairs = nrow(actual), mismatches = nrow(mismatch))
    print(summaries[[case]])
    if (opt$consequences) {
      expected_so <- duckvep_projection_consequences_expected(json)
      actual_so <- DBI::dbGetQuery(con, "SELECT e.ID, n.transcript_id,
        coalesce((SELECT string_agg(s.consequence, '&' ORDER BY s.consequence)
          FROM duckvep_so_terms() s
          WHERE (a.consequence_mask & s.consequence_mask) != 0), '') consequence_terms,
        a.duckvep_status status, a.duckvep_reason reason
        FROM duckvep_annotate('projection_events', 'projection', rich := true) a
        JOIN projection_events e USING(event_index)
        LEFT JOIN duckvep_transcript_names n USING(transcript_index)")
      so_mismatch <- duckvep_projection_compare(actual_so, expected_so)
      for (name in c("actual", "expected", "mismatch")) {
        value <- switch(name, actual = actual_so, expected = expected_so,
          mismatch = so_mismatch)
        utils::write.csv(value, file.path(directory, paste0(case, ".consequence.", name, ".csv")),
          row.names = FALSE)
      }
      consequence_summaries[[case]] <- data.frame(case, input_alleles = input_count,
        expected_pairs = nrow(expected_so), actual_pairs = nrow(actual_so),
        mismatched_pairs = length(unique(so_mismatch$ID)),
        mismatched_fields = nrow(so_mismatch),
        unresolved_pairs = sum(is.na(actual_so$status) | actual_so$status != "supported"))
      print(consequence_summaries[[case]])
    }
    # A valid comparison must fail on every changed field, lost/extra pair and
    # duplicate key, and on changing a value to or from NULL. Errors do not pass.
    mutations <- lapply(names(expected), function(field) {
      corrupt <- expected
      corrupt[[field]][[1L]] <- "DELIBERATE_CORRUPTION"
      corrupt
    })
    names(mutations) <- names(expected)
    mutations$missing_pair <- expected[-1L, ]
    mutations$duplicate_pair <- rbind(expected, expected[1L, ])
    extra <- expected[1L, ]
    extra$ID <- "DELIBERATE_EXTRA_KEY"
    mutations$extra_pair <- rbind(expected, extra)
    for (field in setdiff(names(expected), c("ID", "transcript_id"))) {
      for (missing in c(FALSE, TRUE)) {
        index <- which(is.na(expected[[field]]) == missing)[1L]
        if (!is.na(index)) {
          corrupt <- expected
          corrupt[[field]][[index]] <- if (missing) "DELIBERATE_VALUE" else NA_character_
          mutations[[paste(field, if (missing) "from_null" else "to_null", sep = "_")]] <- corrupt
        }
      }
    }
    detected <- vapply(mutations, function(corrupt)
      nrow(duckvep_projection_compare(corrupt, expected)) > 0L, logical(1L))
    if (case %in% c("noncoding_first_exon_phase1", "noncoding_first_exon_phase2")) {
      # Deliberately substitute the translation-start phase for VEP's first-
      # transcript-exon display phase. This is a wrong-policy control, not an
      # alteration of the oracle or the valid model loaded by the kernel.
      DBI::dbExecute(con, "CREATE OR REPLACE TABLE wrong_projection_phase AS
        SELECT * REPLACE(list_transform(exons, e -> struct_pack(
          exon_start := e.exon_start, exon_end := e.exon_end,
          exon_cdna_start := e.exon_cdna_start, exon_cdna_end := e.exon_cdna_end,
          phase := CASE WHEN e.exon_cdna_start=1 THEN exons[2].phase ELSE e.phase END,
          end_phase := e.end_phase)) AS exons) FROM projection_transcripts")
      wrong <- DBI::dbGetQuery(con, "SELECT p.* EXCLUDE(event_index, transcript_index), e.ID,
        n.transcript_id FROM duckvep_transcript_projection(
          'projection_events', 'projection_annotations', 'wrong_projection_phase') p
        JOIN projection_events e USING(event_index)
        LEFT JOIN duckvep_transcript_names n USING(transcript_index)")
      wrong <- as.data.frame(lapply(wrong[names(expected)], as.character), stringsAsFactors = FALSE)
      wrong_comparison <- duckvep_projection_compare(wrong, expected)
      utils::write.csv(wrong_comparison, file.path(directory, paste0(case, ".wrong-phase.csv")), row.names = FALSE)
      detected <- c(detected, coding_exon_phase_substitution = any(wrong_comparison$field == "cds_start"))
    }
    stopifnot(all(detected), nrow(duckvep_projection_compare(expected, expected)) == 0L)
    controls[[case]] <- data.frame(case, control = names(detected), detected)
    DBI::dbGetQuery(con, "SELECT duckvep_model_drop('projection')")
  }
  summary <- do.call(rbind, summaries)
  utils::write.csv(summary, file.path(directory, "summary.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, controls), file.path(directory, "controls.csv"), row.names = FALSE)
  consequence_summary <- do.call(rbind, consequence_summaries)
  if (opt$consequences) utils::write.csv(consequence_summary,
    file.path(directory, "consequence_summary.csv"), row.names = FALSE)
  artifacts <- c(inputs, lock, file.path(root, c("test/duckvep/conformance/projection_differential.R",
    "test/duckvep/conformance/projection_fixtures.R", "test/duckvep/conformance/minimal_model.sql",
    "test/duckvep/conformance/generate_witnesses.R", "scripts/duckvep_evidence.R")),
    list.files(directory, full.names = TRUE, pattern = "\\.(gff3|vcf|json|csv)$"))
  utils::write.table(data.frame(path = artifacts,
    sha256 = vapply(artifacts, duckvep_evidence_sha256, character(1L))),
    file.path(directory, "artifacts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  receipt <- c(source_revision = revision, build_binding = binding, extension_sha256 = extension_hash,
    revisions, seed = opt$seed, random_cases = opt$random_cases, max_random_length = opt$max_random_length,
    input_alleles = sum(summary$input_alleles),
    expected_pairs = sum(summary$expected_pairs), actual_pairs = sum(summary$actual_pairs),
    mismatches = sum(summary$mismatches))
  if (opt$consequences) receipt <- c(receipt, native_consequences = TRUE,
    consequence_mismatched_pairs = sum(consequence_summary$mismatched_pairs),
    consequence_mismatched_fields = sum(consequence_summary$mismatched_fields),
    consequence_unresolved_pairs = sum(consequence_summary$unresolved_pairs))
  utils::write.table(data.frame(field = names(receipt), value = unname(receipt)),
    file.path(directory, "receipt.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  check_oracle()
  stopifnot(identical(extension_hash, duckvep_evidence_sha256(extension)))
  if (!is.null(opt$extension_receipt)) duckvep_evidence_assert_checkout(root, revision, allowed_outputs)
  stopifnot(sum(summary$mismatches) == 0L)
  if (opt$consequences) stopifnot(sum(consequence_summary$mismatched_fields) == 0L,
    sum(consequence_summary$unresolved_pairs) == 0L)
  if (!is.null(opt$summary_output)) {
    evidence <- cbind(source_revision = revision, build_binding = binding,
      extension_sha256 = extension_hash, summary)
    evidence$detected_controls <- vapply(controls, nrow, integer(1L))
    for (artifact in c("gff3", "generated.vcf", "vcf", "json", "actual.csv", "expected.csv")) {
      evidence[[paste0(gsub("\\.", "_", artifact), "_sha256")]] <- vapply(summary$case,
        function(case) duckvep_evidence_sha256(file.path(directory, paste0(case, ".", artifact))),
        character(1L))
    }
    if (file.exists(opt$summary_output)) evidence <- rbind(
      utils::read.csv(opt$summary_output, stringsAsFactors = FALSE), evidence)
    dir.create(dirname(opt$summary_output), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(evidence, opt$summary_output, row.names = FALSE)
  }
}

if (sys.nframe() == 0L) main()
