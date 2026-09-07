#!/usr/bin/env Rscript
# Additional public-SQL lane over unchanged, receipted Haplosaurus artifacts.
# No generation, sampling, oracle modification or replacement of prior results.
main <- function() {
  opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--artifacts", type = "character"),
    optparse::make_option("--extension", default = "build/release/duckhts.duckdb_extension"),
    optparse::make_option("--extension-receipt", dest = "extension_receipt", default = NULL)
  )))
  source("scripts/duckvep_evidence.R", local = TRUE)
  root <- normalizePath(".")
  revision <- duckvep_evidence_revision(root)
  artifact <- normalizePath(opt$artifacts, mustWork = TRUE)
  extension <- normalizePath(opt$extension, mustWork = TRUE)
  binding <- "diagnostic_unbound"
  if (!is.null(opt$extension_receipt)) {
    duckvep_evidence_assert_checkout(root, revision)
    binding <- duckvep_evidence_read_extension_receipt(opt$extension_receipt, root, extension, revision)$binding
  }
  original <- jsonlite::read_json(file.path(artifact, "receipt.json"), simplifyVector = TRUE)
  stopifnot(identical(original$oracle_revisions$vep, "57ea5c52340acc1f156267f810ad162e26597082"),
            identical(original$oracle_revisions$variation, "2fb834b987ede3824e200197a838ce11e91aeb4b"))
  files <- file.path(artifact, c("inputs.rds", "genotypes.rds", "native.rds", "oracle.jsonl", "summary.csv"))
  hashes <- vapply(files, duckvep_evidence_sha256, "")
  stopifnot(identical(unname(hashes), unname(unlist(original$sha256[files]))))
  baseline <- read.csv(files[5L])
  stopifnot(baseline$failures == 0L, baseline$verifier_controls_rejected == 7L,
            baseline$genotype_routing_controls_rejected == 6L)
  cases <- readRDS(files[1L])$cases
  inputs <- readRDS(files[2L])
  native <- readRDS(files[3L])
  oracle <- lapply(readLines(files[4L]), jsonlite::fromJSON, simplifyVector = FALSE)
  names(oracle) <- vapply(oracle, `[[`, "", "transcript")
  con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(con, paste("LOAD", DBI::dbQuoteString(con, extension)))
  DBI::dbExecute(con, "SET threads=4")
  regions <- data.frame(seq_region = seq_along(cases) - 1L, chrom = vapply(cases, `[[`, "", "chrom"))
  transcripts <- data.frame(transcript_index = seq_along(cases) - 1L, seq_region = regions$seq_region,
    strand = vapply(cases, `[[`, 1L, "strand"), cds = vapply(cases, `[[`, "", "cds"))
  exons <- do.call(rbind, lapply(seq_along(cases), function(i) {
    e <- cases[[i]]$exons[order(cases[[i]]$exons$start, decreasing = cases[[i]]$strand < 0L), ]
    ends <- cumsum(e$end - e$start + 1L)
    data.frame(transcript_index = i - 1L, start = e$start, end = e$end,
               cdna_start = c(1L, head(ends, -1L) + 1L), cdna_end = ends)
  }))
  DBI::dbWriteTable(con, "hap_regions", regions)
  DBI::dbWriteTable(con, "hap_transcripts", transcripts)
  DBI::dbWriteTable(con, "hap_exons", exons)
  DBI::dbWriteTable(con, "hap_genotypes", inputs$calls)
  queries <- c("SELECT seq_region::UINTEGER seq_region FROM hap_regions ORDER BY seq_region",
    paste("SELECT transcript_index::UINTEGER transcript_index,seq_region::UINTEGER seq_region,",
      "11::UBIGINT transcript_start,220::UBIGINT transcript_end,strand::TINYINT strand,",
      "transcript_index::UINTEGER gene_index,3::UBIGINT transcript_flags,11::UBIGINT cds_start,",
      "220::UBIGINT cds_end,cds::BLOB cds_sequence,1::UTINYINT codon_table FROM hap_transcripts ORDER BY seq_region"),
    paste('SELECT transcript_index::UINTEGER transcript_index,start::UBIGINT exon_start,"end"::UBIGINT exon_end,',
      "cdna_start::UBIGINT exon_cdna_start,cdna_end::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase",
      "FROM hap_exons ORDER BY transcript_index,exon_cdna_start"))
  stopifnot(DBI::dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('public_hap',",
    paste(DBI::dbQuoteString(con, queries), collapse = ","), ")"))$loaded)
  calls_query <- paste('SELECT record_index event_index,r.seq_region,g.POS AS position,g.REF AS reference,',
    'g.ALT[1] AS alternate,1 alt_index,t.transcript_index,g.sample_index,g.alleles,g.phase_before,g.phase_set',
    'FROM hap_genotypes g JOIN hap_regions r ON g.CHROM=r.chrom JOIN hap_transcripts t USING(seq_region)')
  stopifnot(DBI::dbGetQuery(con, paste0("SELECT count(*) n FROM (", calls_query, ")"))$n == nrow(inputs$calls))
  events <- unique(inputs$calls[c("record_index", "ID")])
  stopifnot(!anyDuplicated(events$record_index), !anyDuplicated(events$ID))
  stopifnot(length(cases) == baseline$transcripts, nrow(events) == baseline$input_records,
            sum(lengths(inputs$calls$alleles)) == baseline$input_allele_slots,
            6L * length(cases) == baseline$haplotype_lanes)
  event_names <- setNames(events$ID, as.character(events$record_index))
  out <- tempfile(paste0("haplotype_sql_seed", original$seed, "_"), tmpdir = dirname(artifact))
  dir.create(out)
  message("Public SQL artifacts: ", out)
  summaries <- list()
  controls <- list()
  for (policy in c("strict", "vep116_compat")) {
    DBI::dbExecute(con, paste0("CREATE OR REPLACE TABLE hap_output AS SELECT * FROM duckvep_haplotypes(",
      DBI::dbQuoteString(con, calls_query), ",'public_hap',phase_policy:=", DBI::dbQuoteString(con, policy), ")"))
    leaves <- DBI::dbGetQuery(con, "SELECT * FROM hap_output")
    rows <- DBI::dbGetQuery(con, paste("SELECT transcript_index,cds,protein,sequence_flags,projection_status,sequence_status,",
      "c.sample_index,c.phase_set,c.haplotype_lane,c.ploidy,",
      "list_transform(contributors,x -> x.event_index) event_ids FROM hap_output,unnest(carriers) u(c)"))
    saveRDS(list(leaves = leaves, carriers = rows), file.path(out, paste0(policy, ".rds")))
    expected <- do.call(rbind, lapply(seq_along(cases), function(i) {
      paths <- native[[cases[[i]]$transcript]]
      stopifnot(length(paths) == 6L)
      data.frame(transcript_index = i - 1L, sample_index = rep(0:2, each = 2L), lane = rep(1:2, 3L),
        cds = vapply(paths, `[[`, "", "cds"), protein = vapply(paths, `[[`, "", "protein"),
        flags = vapply(paths, function(x) paste(x$flags, collapse = ","), ""),
        events = vapply(paths, function(x) paste(sort(x$contributors), collapse = ","), ""))
    }))
    # Reference lanes are intentionally implicit. Required occupied keys come
    # from independent complete paths, not from whichever SQL rows survived.
    occupied <- expected[nzchar(expected$events), ]
    key <- function(x) paste(x$transcript_index, x$sample_index, x$lane, sep = "/")
    convert <- function(actual) {
      data.frame(transcript_index = as.integer(actual$transcript_index),
        sample_index = as.integer(actual$sample_index), lane = as.integer(actual$haplotype_lane),
        cds = actual$cds, protein = actual$protein,
        flags = vapply(as.integer(actual$sequence_flags), function(x)
          paste(sort(c("indel", "frameshift", "resolved_frameshift")[bitwAnd(x,c(1L,2L,4L)) != 0L]), collapse = ","), ""),
        events = vapply(actual$event_ids, function(ids) {
          names <- unname(event_names[as.character(ids)])
          if (anyNA(names)) return("UNKNOWN_EVENT")
          paste(sort(names), collapse = ",")
        }, ""))
    }
    equal <- function(actual) {
      x <- convert(actual)
      if (anyDuplicated(key(x)) || !setequal(key(x),key(occupied))) return(FALSE)
      x <- x[match(key(occupied),key(x)), ]
      isTRUE(all.equal(unname(as.list(x)),unname(as.list(occupied)),check.attributes = FALSE)) &&
        all(actual$ploidy == 2L) && all(actual$projection_status == "ok") && all(actual$sequence_status == "ok") &&
        if (policy == "strict") all(actual$phase_set == 10) else all(is.na(actual$phase_set))
    }
    mechanics_equal <- equal(rows)
    checks <- logical()
    reconstructed <- rbind(expected[!nzchar(expected$events), ], convert(rows))
    for (i in seq_along(cases)) {
      paths <- reconstructed[reconstructed$transcript_index == i - 1L, ]
      o <- oracle[[cases[[i]]$transcript]]
      checks <- c(checks, o$total_haplotype_count == nrow(paths),
        setequal(paths$cds, vapply(o$haplotypes, `[[`, "", "cds")))
      for (h in o$haplotypes) {
        p <- paths[paths$cds == h$cds, ]
        counts <- table(inputs$samples$sample_name[p$sample_index + 1L])
        observed_counts <- unlist(h$samples)
        contributors <- unique(unlist(strsplit(p$events[nzchar(p$events)], ",", fixed = TRUE)))
        checks <- c(checks, nrow(p) == h$count, all(p$protein == h$protein),
          all(p$flags == paste(sort(as.character(unlist(h$flags))), collapse = ",")),
          setequal(contributors, unlist(h$contributors)),
          setequal(names(counts), names(observed_counts)) &&
            identical(as.integer(counts[sort(names(counts))]), as.integer(observed_counts[sort(names(counts))])))
      }
    }
    mutations <- list(missing = rows[-1L, ], duplicate = rbind(rows, rows[1L, ]),
      protein = rows, contributor = rows, ploidy = rows, phase_set = rows)
    mutations$protein$protein[1L] <- paste0(rows$protein[1L], "X")
    mutations$contributor$event_ids[[1L]] <- -1
    mutations$ploidy$ploidy[1L] <- 3L
    mutations$phase_set$phase_set[1L] <- 99
    rejected <- vapply(mutations, function(x) !equal(x), TRUE)
    controls[[policy]] <- data.frame(policy, control = names(rejected), rejected)
    summaries[[policy]] <- data.frame(policy, transcripts = length(cases), input_records = nrow(events),
      input_call_rows = nrow(inputs$calls), input_allele_slots = sum(lengths(inputs$calls$alleles)),
      complete_lanes = nrow(expected), occupied_carriers = nrow(rows), output_leaves = nrow(leaves),
      native_paths_equal = mechanics_equal, oracle_comparisons = length(checks),
      oracle_failures = sum(!checks), controls_rejected = sum(rejected))
  }
  summary <- do.call(rbind, summaries)
  write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
  write.csv(do.call(rbind, controls), file.path(out, "controls.csv"), row.names = FALSE)
  identities <- c(files, file.path(artifact,"receipt.json"), extension,
    "test/duckvep/conformance/haplotype_sql_differential.R", file.path(out,c("strict.rds","vep116_compat.rds","summary.csv","controls.csv")))
  jsonlite::write_json(list(source_revision = revision, extension_build_binding = binding,
    scope = "public_literal_sequence_replay_not_combined_consequence_hgvs_or_broad_phase_compatibility",
    source_receipt = file.path(artifact,"receipt.json"), oracle_revisions = original$oracle_revisions,
    seed = original$seed, sha256 = as.list(vapply(identities,duckvep_evidence_sha256,""))),
    file.path(out,"receipt.json"), pretty = TRUE, auto_unbox = TRUE)
  print(summary, row.names = FALSE)
  stopifnot(all(summary$native_paths_equal), all(summary$oracle_failures == 0L), all(summary$controls_rejected == 6L),
    identical(hashes, vapply(files,duckvep_evidence_sha256,"")), identical(revision,duckvep_evidence_revision(root)))
  if (!is.null(opt$extension_receipt)) duckvep_evidence_assert_checkout(root,revision)
}
main()
