#!/usr/bin/env Rscript
# Finite raw-GT audit against unmodified Haplosaurus. Disagreements are retained,
# not accepted as conformance. The existing seeded replay corpus is untouched.

genotypes <- function(max_ploidy) {
  unlist(lapply(seq_len(max_ploidy), function(n) {
    alleles <- as.matrix(expand.grid(rep(list(c("0", "1", "2", ".")), n),
      stringsAsFactors = FALSE))
    separators <- if (n == 1L) matrix(character(), 1L, 0L) else
      as.matrix(expand.grid(rep(list(c("|", "/")), n - 1L), stringsAsFactors = FALSE))
    unlist(lapply(seq_len(nrow(alleles)), function(i) {
      vapply(seq_len(nrow(separators)), function(j) {
        paste0(paste0(alleles[i, ], c(separators[j, ], "")), collapse = "")
      }, "")
    }), use.names = FALSE)
  }), use.names = FALSE)
}

canonical <- function(rows) {
  if (!length(rows)) return(list())
  keys <- vapply(rows, function(x) jsonlite::toJSON(list(cds = x$cds, protein = x$protein),
    auto_unbox = TRUE, na = "null"), "")
  lapply(split(rows, keys), function(group) list(cds = group[[1L]]$cds,
    protein = group[[1L]]$protein, count = sum(vapply(group, function(x) as.numeric(x$count), 0)),
    contributors = sort(unique(unlist(lapply(group, `[[`, "contributors"), use.names = FALSE)))))
}

phase_equal <- function(expected, observed) {
  stopifnot(identical(names(expected), names(observed)), nrow(expected) == nrow(observed))
  Reduce(`&`, Map(`==`, expected, observed))
}

main <- function() {
  opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--max-ploidy", dest = "max_ploidy", type = "integer", default = 4L),
    optparse::make_option("--vep-prefix", dest = "vep_prefix",
      default = Sys.getenv("VEP_PREFIX", "/root/miniconda3/envs/vep")),
    optparse::make_option("--extension-receipt", dest = "extension_receipt", default = NULL)
  )))
  stopifnot(opt$max_ploidy >= 2L, opt$max_ploidy <= 4L)
  source("scripts/duckvep_evidence.R", local = TRUE)
  root <- normalizePath(".")
  revision <- duckvep_evidence_revision(root)
  extension <- normalizePath("build/release/duckhts.duckdb_extension")
  binding <- "diagnostic_unbound"
  if (!is.null(opt$extension_receipt)) {
    duckvep_evidence_assert_checkout(root, revision)
    binding <- duckvep_evidence_read_extension_receipt(opt$extension_receipt,
      root, extension, revision)$binding
  }
  mirrors <- c(vep = normalizePath(".sync/ensembl-vep"),
    variation = normalizePath(".sync/ensembl-variation"))
  pins <- c(vep = "57ea5c52340acc1f156267f810ad162e26597082",
    variation = "2fb834b987ede3824e200197a838ce11e91aeb4b")
  for (name in names(pins)) {
    stopifnot(identical(duckvep_evidence_command("git", c("-C", mirrors[[name]],
      "rev-parse", "HEAD"), "oracle revision"), pins[[name]]),
      !length(duckvep_evidence_command("git", c("-C", mirrors[[name]], "status",
        "--porcelain"), "oracle worktree")))
  }
  prefix <- normalizePath(opt$vep_prefix)
  out <- tempfile("haplotype_phase_", tmpdir = "test/duckvep/conformance/results")
  dir.create(out)
  message("Phase artifacts: ", out)
  environment <- duckvep_evidence_command("micromamba", c("list", "-p", prefix, "--explicit"),
    "oracle environment")
  writeLines(environment, file.path(out, "environment.txt"))
  stopifnot(identical(duckvep_evidence_explicit_packages(environment),
    duckvep_evidence_explicit_packages(readLines(
      "test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt"))))
  Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(root, "r/duckhtsbench/inst/benchmark_registry.tsv"))
  paths <- duckhtsbench::duckhts_bench_stage_repository_fixtures(root, "duckvep-haplotypes")
  cds <- readLines(paths[["haplotype_benchmark_reference"]])[2L]
  stopifnot(nchar(cds) == 180L, substr(cds, 31L, 31L) == "G", substr(cds, 34L, 34L) == "G")
  plain <- genotypes(opt$max_ploidy)
  gt <- c(plain, paste0("|", plain), paste0("/", plain))
  stopifnot(!anyDuplicated(gt), length(gt) == 3 * sum(4^(1:opt$max_ploidy) * 2^(0:(opt$max_ploidy - 1L))))
  ploidy <- lengths(strsplit(sub("^[|/]", "", gt), "[|/]"))
  cases <- data.frame(transcript_index = seq_along(gt) - 1L, seq_region = seq_along(gt) - 1L,
    chrom = sprintf("chrP%05d", seq_along(gt)), transcript = sprintf("HP%05d", seq_along(gt)),
    GT = gt, ploidy, prefix = grepl("^[|/]", gt), missing = grepl(".", gt, fixed = TRUE),
    mixed = grepl("|", gt, fixed = TRUE) & grepl("/", gt, fixed = TRUE), cds = cds)
  write.table(cases, file.path(out, "cases.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(as.vector(rbind(paste0(">", cases$chrom), paste0(strrep("A", 10), cds, strrep("A", 10)))),
    file.path(out, "reference.fa"))
  gff <- unlist(lapply(seq_len(nrow(cases)), function(i) {
    id <- cases$transcript[i]
    attrs <- c(paste0("ID=gene:", id, ";biotype=protein_coding"),
      paste0("ID=transcript:", id, ";Parent=gene:", id, ";biotype=protein_coding"),
      paste0("ID=exon:", id, ";Parent=transcript:", id), paste0("Parent=transcript:", id))
    paste(cases$chrom[i], "phase", c("gene", "mRNA", "exon", "CDS"), 11, 190, ".", "+",
      c(".", ".", ".", "0"), attrs, sep = "\t")
  }))
  writeLines(c("##gff-version 3", gff), file.path(out, "model.gff3"))
  vcf <- unlist(lapply(seq_len(nrow(cases)), function(i) c(
    paste(cases$chrom[i], 41, "a", "G", "A,T", ".", "PASS", ".", "GT:PS", paste0(gt[i], ":10"), sep = "\t"),
    paste(cases$chrom[i], 44, "b", "G", "C", ".", "PASS", ".", "GT:PS",
      paste0(paste(rep("1", ploidy[i]), collapse = "|"), ":20"), sep = "\t"))))
  writeLines(c("##fileformat=VCFv4.4", paste0("##contig=<ID=", cases$chrom, ",length=200>"),
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample", vcf), file.path(out, "calls.vcf"))
  run <- function(command, args, name) {
    status <- system2(command, shQuote(args), stdout = file.path(out, paste0(name, ".stdout")),
      stderr = file.path(out, paste0(name, ".stderr")))
    if (status != 0L) stop(name, " failed; retained artifacts: ", out)
  }
  run("samtools", c("faidx", file.path(out, "reference.fa")), "faidx")
  run("bcftools", c("norm", "-c", "e", "-f", file.path(out, "reference.fa"),
    "-o", file.path(out, "ref_checked.vcf"), file.path(out, "calls.vcf")), "reference_check")
  run("bgzip", file.path(out, "model.gff3"), "bgzip")
  run("tabix", c("-p", "gff", file.path(out, "model.gff3.gz")), "tabix")
  perl_lib <- paste(c(file.path(mirrors, "modules"), file.path(prefix, "share/ensembl-vep-116.0-0")), collapse = ":")
  run("micromamba", c("run", "--clean-env", "--env", paste0("PERL5LIB=", perl_lib), "-p", prefix,
    "perl", "test/duckvep/conformance/haplotype_oracle.pl", file.path(out, "calls.vcf"),
    file.path(out, "reference.fa"), file.path(out, "model.gff3.gz"),
    file.path(out, "phase.jsonl")), "oracle")
  oracle <- lapply(readLines(file.path(out, "oracle.stdout")), jsonlite::fromJSON, simplifyVector = FALSE)
  names(oracle) <- vapply(oracle, `[[`, "", "transcript")
  stopifnot(!anyDuplicated(names(oracle)), setequal(names(oracle), cases$transcript))
  for (o in oracle) {
    stopifnot(sum(vapply(o$haplotypes, function(h) as.numeric(h$count), 0)) == o$total_haplotype_count,
      all(vapply(o$haplotypes, function(h)
        identical(names(h$samples), "sample") && h$samples$sample == h$count, TRUE)))
  }

  # Observe actual upstream objects, not a second implementation of its parser.
  # Source ploidy/missingness is independently known from the generated grammar;
  # effective file ploidy and retained allele slots come from Haplosaurus.
  phase <- lapply(readLines(file.path(out, "phase.jsonl")), jsonlite::fromJSON, simplifyVector = FALSE)
  names(phase) <- vapply(phase, `[[`, "", "transcript")
  stopifnot(!anyDuplicated(names(phase)), setequal(names(phase), cases$transcript))
  expected_phase <- do.call(rbind, lapply(seq_len(nrow(cases)), function(i) {
    p <- phase[[cases$transcript[i]]]
    stopifnot(p$default_ploidy == 2L, identical(names(p$sample_ploidy), "sample"),
      p$sample_ploidy$sample == 2L)
    ids <- vapply(p$calls, `[[`, "", "source_id")
    stopifnot(!anyDuplicated(ids), all(ids %in% c("a", "b")),
      all(vapply(p$calls, function(x) identical(x$sample, "sample"), TRUE)))
    do.call(rbind, lapply(c("a", "b"), function(id) {
      call <- p$calls[ids == id]
      alleles <- if (id == "a") c("G", "A", "T") else c("G", "C")
      indices <- if (length(call)) match(unlist(call[[1L]]$genotype), alleles) - 1L else integer()
      stopifnot(!anyNA(indices))
      data.frame(status = 0L, retained = length(call) == 1L, ploidy = cases$ploidy[i],
        missing = as.integer(id == "a" && cases$missing[i]), slots = length(indices),
        first = if (length(indices)) indices[1L] else -1L,
        second = if (length(indices) > 1L) indices[2L] else -1L)
    }))
  }))
  raw_gt <- as.vector(rbind(cases$GT, vapply(cases$ploidy,
    function(n) paste(rep("1", n), collapse = "|"), "")))
  probe_sources <- c("test/duckvep/conformance/phase_probe.c", "src/duckvep/kernel/src/duckvep_phase.c")
  probe <- file.path(out, paste0("phase_probe", .Platform$dynlib.ext))
  compiler <- Sys.getenv("CC", "cc")
  run(compiler, c("-std=c99", "-O2", "-Wall", "-Wextra", "-Werror", "-fPIC", "-shared",
    "-Isrc/duckvep/kernel/src", probe_sources, "-o", probe), "phase_compile")
  writeLines(duckvep_evidence_command(compiler, "--version", "phase compiler"),
    file.path(out, "phase_compiler.txt"))
  dll <- dyn.load(probe)
  parsed <- .C(getNativeSymbolInfo("duckhts_test_vep116_raw_gt", dll), raw_gt,
    rep(c(2L, 1L), nrow(cases)), as.integer(length(raw_gt)), status = integer(length(raw_gt)),
    disposition = integer(length(raw_gt)), ploidy = integer(length(raw_gt)),
    missing = integer(length(raw_gt)), slots = integer(length(raw_gt)),
    first = integer(length(raw_gt)), second = integer(length(raw_gt)))
  dyn.unload(dll[["path"]])
  observed_phase <- as.data.frame(parsed[c("status", "ploidy", "missing", "slots", "first", "second")])
  observed_phase$retained <- parsed$disposition == 3L
  observed_phase <- observed_phase[names(expected_phase)]
  phase_matches <- phase_equal(expected_phase, observed_phase)
  phase_keys <- data.frame(transcript = rep(cases$transcript, each = 2L),
    source_id = rep(c("a", "b"), nrow(cases)), GT = raw_gt)
  saveRDS(list(keys = phase_keys, expected = expected_phase, observed = observed_phase,
    disposition = parsed$disposition), file.path(out, "phase_comparisons.rds"))
  write.csv(cbind(phase_keys, equal = phase_matches), file.path(out, "phase_summary.csv"), row.names = FALSE)
  write.csv(cbind(phase_keys[!phase_matches, ], expected_phase[!phase_matches, ],
    observed_phase[!phase_matches, ]), file.path(out, "phase_mismatches.csv"), row.names = FALSE)
  phase_rejected <- vapply(names(expected_phase), function(field) {
    corrupt <- expected_phase
    corrupt[[field]][1L] <- if (is.logical(corrupt[[field]])) !corrupt[[field]][1L] else
      corrupt[[field]][1L] + 1L
    !all(phase_equal(expected_phase, corrupt))
  }, TRUE)
  stopifnot(all(phase_rejected))
  write.csv(data.frame(control = names(phase_rejected), rejected = phase_rejected),
    file.path(out, "phase_controls.csv"), row.names = FALSE)

  con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  q <- function(x) as.character(DBI::dbQuoteString(con, x))
  DBI::dbExecute(con, paste("LOAD", q(extension)))
  DBI::dbExecute(con, "SET threads=4")
  DBI::dbWriteTable(con, "models", cases)
  queries <- c("SELECT seq_region::UINTEGER seq_region FROM models ORDER BY seq_region",
    "SELECT transcript_index::UINTEGER transcript_index,seq_region::UINTEGER seq_region,
     11::UBIGINT transcript_start,190::UBIGINT transcript_end,1::TINYINT strand,0::UINTEGER gene_index,
     3::UBIGINT transcript_flags,11::UBIGINT cds_start,190::UBIGINT cds_end,cds::BLOB cds_sequence,
     1::UTINYINT codon_table FROM models ORDER BY transcript_index",
    "SELECT transcript_index::UINTEGER transcript_index,11::UBIGINT exon_start,190::UBIGINT exon_end,
     1::UBIGINT exon_cdna_start,180::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase
     FROM models ORDER BY transcript_index")
  stopifnot(DBI::dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('phase',",
    paste(q(queries), collapse = ","), ")"))$loaded)
  DBI::dbExecute(con, paste0("CREATE TABLE records AS SELECT * FROM read_geno(",
    q(normalizePath(file.path(out, "calls.vcf"))), ")"))
  DBI::dbExecute(con, "CREATE TABLE calls AS SELECT r.record_index*3+a.i event_index,
    m.seq_region,r.POS AS position,r.REF reference,r.ALT[a.i] alternate,a.i alt_index,
    m.transcript_index,c.sample_index,c.alleles,c.phase_before,c.phase_set
    FROM records r JOIN models m ON m.chrom=r.CHROM,unnest(r.calls) u(c),range(1,len(r.ALT)+1) a(i)")
  records <- DBI::dbGetQuery(con, "SELECT * FROM records")
  calls <- DBI::dbGetQuery(con, "SELECT * FROM calls")
  stopifnot(nrow(records) == 2L * nrow(cases), nrow(calls) == 3L * nrow(cases),
    all(lengths(calls$alleles) == cases$ploidy[calls$transcript_index + 1L]))
  actual <- DBI::dbGetQuery(con, "SELECT * FROM duckvep_haplotypes('SELECT * FROM calls',
    'phase',phase_policy:='vep116_compat')")
  stopifnot(all(actual$carrier_count == vapply(actual$carriers, nrow, 1L)),
    all(vapply(actual$carriers, function(c) all(c$sample_index == 0L), TRUE)))
  saveRDS(list(records = records, calls = calls, actual = actual), file.path(out, "native.rds"))
  comparisons <- lapply(seq_len(nrow(cases)), function(i) {
    o <- oracle[[cases$transcript[i]]]
    a <- actual[actual$transcript_index == cases$transcript_index[i], ]
    observed <- lapply(seq_len(nrow(a)), function(j) list(cds = a$cds[j], protein = a$protein[j],
      count = a$carrier_count[j], contributors = records$ID[match(
        (a$contributors[[j]]$event_index - 1) %/% 3, records$record_index)]))
    expected <- canonical(o$haplotypes)
    observed <- canonical(observed)
    list(expected = expected, observed = observed, equal = identical(expected, observed),
      oracle_lanes = o$total_haplotype_count, native_lanes = sum(a$carrier_count),
      native_unknown = sum(is.na(a$cds)),
      native_unavailable_carriers = sum(a$carrier_count[is.na(a$cds)]))
  })
  saveRDS(comparisons, file.path(out, "comparisons.rds"))
  summary <- cases[setdiff(names(cases), "cds")]
  for (name in c("equal", "oracle_lanes", "native_lanes", "native_unknown", "native_unavailable_carriers"))
    summary[[name]] <- vapply(comparisons, `[[`, if (name == "equal") TRUE else 0, name)
  write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
  # A fully called ordinary diploid GT is a positive control, not an exclusion:
  # all other profiles and their failures remain in the complete comparison table.
  control <- with(summary, ploidy == 2L & !prefix & !mixed & !missing)
  stopifnot(sum(control) == 18L, all(summary$equal[control]))
  witness <- oracle[[cases$transcript[which(control)[1L]]]]$haplotypes
  controls <- list(duplicate = c(witness, witness[1L]), cds = witness,
    protein = witness, contributor = witness)
  controls$cds[[1L]]$cds <- paste0("C", substring(witness[[1L]]$cds, 2L))
  controls$protein[[1L]]$protein <- paste0("X", substring(witness[[1L]]$protein, 2L))
  controls$contributor[[1L]]$contributors <- c(witness[[1L]]$contributors, "deliberate_extra_event")
  rejected <- vapply(controls, function(x) !identical(canonical(x), canonical(witness)), TRUE)
  stopifnot(all(rejected))
  write.csv(data.frame(control = names(rejected), rejected),
    file.path(out, "controls.csv"), row.names = FALSE)
  decoded <- DBI::dbGetQuery(con, "SELECT m.transcript_index,to_json(struct_pack(
    alleles:=c.alleles,phase_before:=c.phase_before))::VARCHAR decoded_gt
    FROM records r JOIN models m ON m.chrom=r.CHROM,unnest(r.calls) u(c) WHERE r.ID='a'
    ORDER BY m.transcript_index")
  stopifnot(identical(as.integer(decoded$transcript_index), cases$transcript_index))
  summary$decoded_gt <- decoded$decoded_gt
  summary$oracle_signature <- vapply(comparisons, function(x)
    jsonlite::toJSON(x$expected, auto_unbox = TRUE, na = "null"), "")
  collisions <- split(summary, summary$decoded_gt)
  collisions <- collisions[vapply(collisions, function(x) length(unique(x$oracle_signature)) > 1L, TRUE)]
  saveRDS(collisions, file.path(out, "decoded_collisions.rds"))
  inputs <- c(paths, extension, probe_sources, "src/duckvep/kernel/src/duckvep_phase.h",
    "test/duckvep/conformance/haplotype_phase_differential.R",
    "test/duckvep/conformance/haplotype_oracle.pl", file.path(prefix,
      "share/ensembl-vep-116.0-0/Bio/EnsEMBL/IO/Parser/BaseVCF4.pm"))
  jsonlite::write_json(list(source_revision = revision, extension_build_binding = binding,
    oracle_revisions = as.list(pins), scope = "raw_GT_finite_phase_audit_not_conformance",
    max_ploidy = opt$max_ploidy, cases = nrow(summary), disagreements = sum(!summary$equal),
    decoded_collision_groups = length(collisions), controls_rejected = sum(rejected),
    raw_parser_calls = length(raw_gt), raw_parser_disagreements = sum(!phase_matches),
    raw_parser_controls_rejected = sum(phase_rejected),
    input_records = nrow(records), source_alt_events = 3L * nrow(cases),
    input_genotype_calls = nrow(records), input_allele_slots = 2L * sum(cases$ploidy),
    candidate_alt_calls = nrow(calls), native_leaves = nrow(actual),
    oracle_lanes = sum(summary$oracle_lanes), native_lanes = sum(summary$native_lanes),
    native_unavailable_carriers = sum(summary$native_unavailable_carriers),
    sha256 = as.list(vapply(unique(c(inputs, list.files(out, full.names = TRUE))),
      duckvep_evidence_sha256, ""))), file.path(out, "receipt.json"), pretty = TRUE, auto_unbox = TRUE)
  if (!is.null(opt$extension_receipt)) duckvep_evidence_assert_checkout(root, revision)
  print(aggregate(cbind(cases = rep(1L, nrow(summary)), disagreements = as.integer(!summary$equal)) ~
    ploidy + prefix + missing + mixed, data = summary, FUN = sum), row.names = FALSE)
  message("Decoded-equivalence collision groups: ", length(collisions))
  message("Raw parser: ", sum(!phase_matches), " disagreements / ", length(raw_gt), " calls")
  if (any(!phase_matches)) stop("Raw parser disagreements retained: ", out, call. = FALSE)
  if (any(!summary$equal)) stop("Raw-GT compatibility disagreements retained: ", out, call. = FALSE)
}
main()
