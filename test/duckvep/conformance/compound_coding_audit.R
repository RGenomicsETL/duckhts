#!/usr/bin/env Rscript
# Additional support audit, not a compound-SO conformance certificate. Every
# original lane is retained, including reference lanes and unsupported results.
main <- function() {
  opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--artifacts", type = "character"),
    optparse::make_option("--bcftools-mirror", dest = "bcftools_mirror",
      default = ".sync/RBCFTools")
  )))
  source("scripts/duckvep_evidence.R", local = TRUE)
  root <- normalizePath(".")
  revision <- duckvep_evidence_revision(root)
  artifact <- normalizePath(opt$artifacts, mustWork = TRUE)
  original <- jsonlite::read_json(file.path(artifact, "receipt.json"), simplifyVector = TRUE)
  inputs <- file.path(artifact, c("inputs.rds", "native.rds", "summary.csv",
    "carriers.vcf", "reference.fa", "model.gff3.gz"))
  hashes <- vapply(inputs, duckvep_evidence_sha256, "")
  for (p in inputs) stopifnot(identical(hashes[[p]], original$sha256[[p]]))
  stopifnot(read.csv(inputs[3L])$failures == 0L)
  out <- tempfile(paste0("compound_coding_seed", original$seed, "_"), tmpdir = dirname(artifact))
  dir.create(out)
  message("Compound coding audit: ", out)

  # Export exactly the pinned source tree: never reuse an unbound local binary,
  # rebuild a user's mirror in place, or let a newer executable change semantics.
  mirror <- normalizePath(opt$bcftools_mirror, mustWork = TRUE)
  upstream <- "9adeaf4cfcc3bff40efca6237749fefb53391678"
  archive <- file.path(out, "bcftools-source.tar")
  stopifnot(system2("git", shQuote(c("-C", mirror, "archive", "--format=tar",
    paste0("--output=", archive), upstream, "src/bcftools-1.23"))) == 0L)
  utils::untar(archive, exdir = out)
  build <- file.path(out, "src/bcftools-1.23")
  stopifnot(system2("make", shQuote(c("-C", build, "-j2", "bcftools", "PLUGINS_ENABLED=no")),
    stdout = file.path(out, "bcftools-build.log"), stderr = file.path(out, "bcftools-build.log")) == 0L)
  bcftools <- file.path(build, "bcftools")
  version <- system2(bcftools, "--version", stdout = TRUE)
  stopifnot(identical(version[1:2], c("bcftools 1.23", "Using htslib 1.23")))
  csq_file <- file.path(out, "bcftools-csq.tsv")
  command <- c("csq", "-f", inputs[5L], "-g", inputs[6L], "-p", "a", "-Ot", "-o", csq_file, inputs[4L])
  bcftools_status <- system2(bcftools, shQuote(command), stdout = file.path(out, "bcftools-stdout.log"),
    stderr = file.path(out, "bcftools-stderr.log"))
  csq <- read.delim(csq_file, comment.char = "#", header = FALSE, quote = "",
    col.names = c("kind", "sample", "lane", "chrom", "pos", "annotation"))
  stopifnot(all(csq$kind == "CSQ"))
  annotations <- strsplit(csq$annotation, "|", fixed = TRUE)
  # A crash can leave a syntactically split but truncated final field.
  bytes <- file(csq_file, "rb")
  seek(bytes, file.info(csq_file)$size - 1, origin = "start")
  terminated <- identical(readBin(bytes, "raw", n = 1L), charToRaw("\n"))
  close(bytes)
  complete <- lengths(annotations) >= 4L
  if (length(complete) && !terminated) complete[length(complete)] <- FALSE
  field <- function(i) vapply(annotations,
    function(x) if (length(x) >= i) x[[i]] else NA_character_, "")
  csq$transcript <- field(3L)
  csq$terms <- field(1L)
  csq$key <- paste(csq$transcript, csq$sample, csq$lane, sep = "/")
  write.csv(csq, file.path(out, "bcftools-rows.csv"), row.names = FALSE)

  sources <- c("test/duckvep/conformance/compound_coding_probe.c",
    paste0("src/duckvep/kernel/src/duckvep_", c("delta", "haplotype", "codon", "coding", "projection"), ".c"))
  code <- c(sources, list.files("src/duckvep/kernel/src", "\\.h$", full.names = TRUE),
    "src/duckvep/kernel/include/duckvep_kernel.h", "test/duckvep/conformance/compound_coding_audit.R")
  code_hashes <- vapply(code, duckvep_evidence_sha256, "")
  shared <- file.path(out, paste0("compound_coding_probe", .Platform$dynlib.ext))
  stopifnot(system2(Sys.getenv("CC", "cc"), shQuote(c("-std=c99", "-O1", "-g", "-Wall", "-Wextra",
    "-fPIC", "-shared", "-I", "src/duckvep/kernel/src", "-I", "src/duckvep/kernel/include",
    sources, "-o", shared)), stdout = file.path(out, "compiler.log"),
    stderr = file.path(out, "compiler.log")) == 0L)
  dll <- dyn.load(shared)
  on.exit(dyn.unload(shared), add = TRUE)
  symbol <- getNativeSymbolInfo("duckhts_test_compound_coding", PACKAGE = dll)
  cases <- readRDS(inputs[1L])$cases
  native <- readRDS(inputs[2L])
  rows <- vector("list", 6L * length(cases))
  at <- 0L
  for (case in cases) {
    paths <- native[[case$transcript]]
    stopifnot(length(paths) == 6L)
    for (slot in seq_along(paths)) {
      path <- paths[[slot]]
      edits <- case$edits[match(path$contributors, case$edits$id), , drop = FALSE]
      stopifnot(!anyNA(edits))
      edits <- edits[order(edits$start, decreasing = TRUE), , drop = FALSE]
      capacity <- as.integer(nchar(case$cds) + sum(edits$alt_len) + 2L)
      x <- .C(symbol, case$cds, case$strand, as.integer(edits$start), edits$ref, edits$alt,
        as.integer(nrow(edits)), capacity, cds = raw(capacity), protein = raw(capacity),
        statuses = integer(2L), facts = integer(17L), lengths = integer(2L))
      actual_cds <- rawToChar(x$cds[seq_len(x$lengths[1L])])
      full_protein <- rawToChar(x$protein[seq_len(x$lengths[2L])])
      displayed <- sub("(\\*).*", "\\1", full_protein)
      names(x$facts) <- c("valid", "sequence_status", "synonymous", "missense", "stop_gained",
        "stop_lost", "stop_retained", "start_lost", "start_retained", "frameshift",
        "inframe_deletion", "inframe_insertion", "protein_altering", "coding_unknown",
        "partial_codon", "flags", "applied_edits")
      key <- paste(case$transcript, path$sample, (slot - 1L) %% 2L + 1L, sep = "/")
      matched <- which(csq$key == key)
      at <- at + 1L
      rows[[at]] <- cbind(data.frame(key, shape = case$shape, strand = case$strand,
        occupied = nrow(edits) > 0L, edit_count = nrow(edits),
        net_length = sum(edits$alt_len - edits$ref_len),
        compound_indel = nrow(edits) > 1L && any(edits$alt_len != edits$ref_len),
        build_status = x$statuses[1L], delta_status = x$statuses[2L],
        cds_match = identical(actual_cds, path$cds), protein_match = identical(displayed, path$protein),
        cds = actual_cds, protein = displayed, full_protein,
        contributors = paste(path$contributors, collapse = ";"),
        bcftools_rows = length(matched),
        bcftools_terms = paste(csq$terms[matched], collapse = ";")), as.data.frame(as.list(x$facts)))
    }
  }
  stopifnot(at == length(rows))
  rows <- do.call(rbind, rows)
  stopifnot(!anyDuplicated(rows$key))
  extra <- csq[!csq$key %in% rows$key, , drop = FALSE]
  write.csv(rows, file.path(out, "lanes.csv"), row.names = FALSE)
  write.csv(extra, file.path(out, "extra-bcftools-rows.csv"), row.names = FALSE)
  # This gate checks a declared support limit, NOT agreement between VEP and BCSQ.
  # Unsupported results stay in the denominator and are never counted as SO matches.
  false_support <- with(rows, compound_indel & (delta_status == 0L | valid != 0L))
  summary <- data.frame(seed = original$seed, transcripts = length(cases), lanes = nrow(rows),
    occupied = sum(rows$occupied), context_failures = sum(rows$build_status != 0L),
    sequence_failures = sum(!rows$cds_match | !rows$protein_match),
    unsupported = sum(rows$delta_status == 2L), compound_indels = sum(rows$compound_indel),
    falsely_supported_compound_indels = sum(false_support), bcftools_exit_status = bcftools_status,
    # Noncoding/splice consequences have four fields; strand and peptide are optional.
    bcftools_rows = nrow(csq), incomplete_bcftools_rows = sum(!complete),
    missing_occupied = sum(rows$occupied & rows$bcftools_rows == 0L), extra_bcftools_rows = nrow(extra),
    restored_with_stop_gained_frameshift = sum(rows$net_length == 0L & rows$compound_indel &
      grepl("stop_gained&frameshift", rows$bcftools_terms, fixed = TRUE)))
  write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
  identities <- c(inputs, file.path(artifact, "receipt.json"), code, archive, bcftools, shared,
    file.path(out, c("bcftools-csq.tsv", "bcftools-rows.csv", "lanes.csv", "summary.csv",
      "extra-bcftools-rows.csv", "bcftools-stderr.log", "bcftools-build.log")))
  jsonlite::write_json(list(source_revision = revision, source_binding = "diagnostic_unbound",
    scope = "coding_context_support_audit_and_complementary_bcftools_observation_not_compound_so_conformance",
    source_receipt = file.path(artifact, "receipt.json"), bcftools_source_revision = upstream,
    bcftools_source_tree = "src/bcftools-1.23", bcftools_version = version, bcftools_command = command,
    sha256 = as.list(vapply(identities, duckvep_evidence_sha256, ""))), file.path(out, "receipt.json"),
    pretty = TRUE, auto_unbox = TRUE)
  print(summary)
  stopifnot(identical(hashes, vapply(inputs, duckvep_evidence_sha256, "")),
    identical(code_hashes, vapply(code, duckvep_evidence_sha256, "")),
    identical(revision, duckvep_evidence_revision(root)),
    !any(false_support), all(rows$build_status == 0L), all(rows$cds_match & rows$protein_match),
    bcftools_status == 0L, summary$incomplete_bcftools_rows == 0L,
    summary$missing_occupied == 0L, nrow(extra) == 0L)
}
main()
