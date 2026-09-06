#!/usr/bin/env Rscript
# Executable Haplosaurus comparison for the existing edit-set mechanics.
# The declared model supplies CDS coordinates to C; Haplosaurus independently
# parses the VCF and projects it through GFF. This does not certify SQL carrier
# grouping, strict phase, compound SO/HGVS, or arbitrary-ploidy compatibility.

main <- function() {
  opt <- optparse::parse_args(optparse::OptionParser(
    option_list = list(
      optparse::make_option("--cases", type = "integer", default = 128L),
      optparse::make_option("--seed", type = "integer", default = 173L),
      optparse::make_option(
        "--vep-prefix",
        dest = "vep_prefix",
        default = Sys.getenv(
          "VEP_PREFIX",
          "/root/miniconda3/envs/vep"
        )
      ),
      optparse::make_option(
        "--vep-git",
        dest = "vep_git",
        default = ".sync/ensembl-vep"
      ),
      optparse::make_option(
        "--variation-git",
        dest = "variation_git",
        default = ".sync/ensembl-variation"
      )
    )
  ))
  stopifnot(opt$cases >= 10L, opt$cases <= 100000L, !is.na(opt$seed))
  root <- normalizePath(system2(
    "git",
    c("rev-parse", "--show-toplevel"),
    stdout = TRUE
  ))
  source(file.path(root, "scripts/duckvep_evidence.R"), local = TRUE)
  results <- file.path(root, "test/duckvep/conformance/results")
  dir.create(results, showWarnings = FALSE)
  out <- tempfile(paste0("haplotype_seed", opt$seed, "_"), tmpdir = results)
  dir.create(out)
  message("Artifacts: ", out)

  # Fail closed on a different oracle or a modified source checkout. Environment
  # matching uses the same exact-package authority as independent-event campaigns.
  revisions <- c(
    vep = "57ea5c52340acc1f156267f810ad162e26597082",
    variation = "2fb834b987ede3824e200197a838ce11e91aeb4b"
  )
  mirrors <- c(
    vep = normalizePath(opt$vep_git, mustWork = TRUE),
    variation = normalizePath(opt$variation_git, mustWork = TRUE)
  )
  for (name in names(mirrors)) {
    revision <- system2(
      "git",
      c("-C", shQuote(mirrors[[name]]), "rev-parse", "HEAD"),
      stdout = TRUE
    )
    dirty <- system2(
      "git",
      c("-C", shQuote(mirrors[[name]]), "status", "--porcelain"),
      stdout = TRUE
    )
    stopifnot(identical(revision, revisions[[name]]), !length(dirty))
  }
  prefix <- normalizePath(opt$vep_prefix, mustWork = TRUE)
  conda <- function(...) {
    do.call(
      blit::conda,
      c(
        as.list(duckvep_blit_quote(c(...))),
        list(conda = duckvep_blit_quote("micromamba"))
      )
    )
  }
  environment_path <- file.path(out, "environment.txt")
  stopifnot(
    blit::cmd_run(
      conda("list", "-p", prefix, "--explicit"),
      stdout = environment_path,
      stderr = "2>&1",
      verbose = FALSE
    ) ==
      0L
  )
  lock <- file.path(
    root,
    "test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt"
  )
  stopifnot(identical(
    duckvep_evidence_explicit_packages(readLines(environment_path)),
    duckvep_evidence_explicit_packages(readLines(lock))
  ))

  # A small native test bridge, not another production annotation implementation.
  shared <- file.path(out, paste0("haplotype_probe", .Platform$dynlib.ext))
  sources <- file.path(
    root,
    c(
      "test/duckvep/conformance/haplotype_probe.c",
      "src/duckvep/kernel/src/duckvep_haplotype.c",
      "src/duckvep/kernel/src/duckvep_carriers.c",
      "src/duckvep/kernel/src/duckvep_codon.c"
    )
  )
  compiler <- Sys.getenv("CC", "cc")
  stopifnot(
    system2(
      compiler,
      c(
        "-std=c99",
        "-O1",
        "-g",
        "-Wall",
        "-Wextra",
        "-fPIC",
        "-shared",
        "-I",
        shQuote(file.path(root, "src/duckvep/kernel/src")),
        "-I",
        shQuote(file.path(root, "src/duckvep/kernel/include")),
        shQuote(sources),
        "-o",
        shQuote(shared)
      ),
      stdout = file.path(out, "compiler.log"),
      stderr = file.path(out, "compiler.log")
    ) ==
      0L
  )
  dll <- dyn.load(shared)
  on.exit(dyn.unload(shared), add = TRUE)
  symbol <- getNativeSymbolInfo("duckhts_test_haplotype", PACKAGE = dll)
  carrier_symbol <- getNativeSymbolInfo(
    "duckhts_test_carrier_haplotypes",
    PACKAGE = dll
  )
  reverse_complement <- function(x) {
    paste(rev(strsplit(chartr("ACGT", "TGCA", x), "")[[1L]]), collapse = "")
  }
  dna <- function(n) {
    paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
  }
  codons <- apply(
    expand.grid(rep(list(c("A", "C", "G", "T")), 3L)),
    1L,
    paste,
    collapse = ""
  )
  codons <- setdiff(codons, c("TAA", "TAG", "TGA"))
  shapes <- c(
    "same_codon",
    "restored_frame",
    "open_frame",
    "in_frame",
    "random_set"
  )
  set.seed(opt$seed)
  cases <- vector("list", opt$cases)
  fasta <- gff <- vcf <- vector("list", opt$cases)
  for (i in seq_len(opt$cases)) {
    chrom <- sprintf("chrH%06d", i)
    transcript <- sprintf("DHT%06d", i)
    strand <- if (i %% 2L) 1L else -1L
    shape <- shapes[(i - 1L) %% length(shapes) + 1L]
    cds <- paste0(
      "ATG",
      paste(sample(codons, 58L, replace = TRUE), collapse = ""),
      "TAA"
    )
    genome_cds <- if (strand == 1L) cds else reverse_complement(cds)
    intron <- if (strand == 1L) {
      paste0("GT", strrep("A", 26L), "AG")
    } else {
      paste0("CT", strrep("T", 26L), "AC")
    }
    genome <- paste0(
      strrep("A", 10L),
      substr(genome_cds, 1L, 90L),
      intron,
      substr(genome_cds, 91L, 180L),
      strrep("A", 10L)
    )
    stopifnot(nchar(cds) == 180L, nchar(genome) == 230L)
    fasta[[i]] <- c(paste0(">", chrom), genome)
    features <- data.frame(
      feature = c("gene", "mRNA", "exon", "CDS", "exon", "CDS"),
      start = c(11L, 11L, 11L, 11L, 131L, 131L),
      end = c(220L, 220L, 100L, 100L, 220L, 220L),
      phase = c(".", ".", ".", "0", ".", "0"),
      attributes = c(
        paste0("ID=gene:", transcript, ";biotype=protein_coding"),
        paste0(
          "ID=transcript:",
          transcript,
          ";Parent=gene:",
          transcript,
          ";biotype=protein_coding"
        ),
        paste0("ID=exon:", transcript, "a;Parent=transcript:", transcript),
        paste0("Parent=transcript:", transcript),
        paste0("ID=exon:", transcript, "b;Parent=transcript:", transcript),
        paste0("Parent=transcript:", transcript)
      )
    )
    gff[[i]] <- apply(features, 1L, function(row) {
      paste(
        c(
          chrom,
          "duckvep",
          row["feature"],
          trimws(row["start"]),
          trimws(row["end"]),
          ".",
          if (strand == 1L) "+" else "-",
          row["phase"],
          row["attributes"]
        ),
        collapse = "\t"
      )
    })
    edits <- switch(
      shape,
      same_codon = data.frame(start = c(4L, 5L), ref_len = 1L, alt_len = 1L),
      restored_frame = data.frame(
        start = c(10L, 40L),
        ref_len = c(0L, 1L),
        alt_len = c(1L, 0L)
      ),
      open_frame = data.frame(
        start = c(10L, 40L),
        ref_len = c(0L, 1L),
        alt_len = c(2L, 1L)
      ),
      in_frame = data.frame(
        start = c(10L, 40L),
        ref_len = c(0L, 3L),
        alt_len = c(3L, 0L)
      ),
      random_set = {
        starts <- sort(sample(
          c(seq(4L, 82L, 6L), seq(94L, 172L, 6L)),
          sample(2:10, 1L)
        ))
        data.frame(
          start = starts,
          ref_len = sample(0:3, length(starts), TRUE),
          alt_len = sample(0:4, length(starts), TRUE)
        )
      }
    )
    edits$ref <- substring(cds, edits$start, edits$start + edits$ref_len - 1L)
    edits$ref[edits$ref_len == 0L] <- ""
    edits$alt <- vapply(edits$alt_len, dna, "")
    for (j in seq_len(nrow(edits))) {
      while (edits$ref[j] == edits$alt[j]) {
        if (edits$alt_len[j] == 0L) {
          edits$alt_len[j] <- 1L
        }
        edits$alt[j] <- dna(edits$alt_len[j])
      }
    }
    # C sees unanchored edits in original-CDS coordinates, with genomic-strand
    # alleles. Haplosaurus sees ordinary left-anchored genomic VCF records.
    if (strand == -1L) {
      edits$ref <- vapply(edits$ref, reverse_complement, "")
      edits$alt <- vapply(edits$alt, reverse_complement, "")
    }
    edits$id <- paste0(transcript, "_e", seq_len(nrow(edits)))
    genomic_position <- function(p) {
      if (strand == 1L) {
        ifelse(p <= 90L, p + 10L, p + 40L)
      } else {
        ifelse(p <= 90L, 221L - p, 191L - p)
      }
    }
    variants <- vector("list", nrow(edits))
    for (j in seq_len(nrow(edits))) {
      edit <- edits[j, ]
      pos <- if (strand == 1L || edit$ref_len == 0L) {
        genomic_position(edit$start)
      } else {
        genomic_position(edit$start + edit$ref_len - 1L)
      }
      ref <- edit$ref
      alt <- edit$alt
      if (edit$ref_len == 0L || edit$alt_len == 0L) {
        pos <- pos - as.integer(edit$ref_len != 0L || strand == 1L)
        anchor <- substring(genome, pos, pos)
        ref <- paste0(anchor, ref)
        alt <- paste0(anchor, alt)
      }
      stopifnot(substring(genome, pos, pos + nchar(ref) - 1L) == ref)
      variants[[j]] <- data.frame(
        pos = pos,
        row = paste(
          c(
            chrom,
            pos,
            edit$id,
            ref,
            alt,
            ".",
            "PASS",
            ".",
            "GT:PS",
            "1|0:10",
            if (j %% 2L) "1|0:10" else "0|1:10",
            "1|0:10"
          ),
          collapse = "\t"
        )
      )
    }
    variants <- do.call(rbind, variants)
    vcf[[i]] <- variants$row[order(variants$pos)]
    cases[[i]] <- list(
      transcript = transcript,
      shape = shape,
      strand = strand,
      cds = cds,
      positions = variants$pos,
      edits = edits
    )
  }
  fasta <- unlist(fasta, use.names = FALSE)
  gff <- unlist(gff, use.names = FALSE)
  vcf <- unlist(vcf, use.names = FALSE)
  writeLines(fasta, file.path(out, "reference.fa"))
  writeLines(c("##gff-version 3", gff), file.path(out, "model.gff3"))
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      sprintf("##contig=<ID=chrH%06d,length=230>", seq_len(opt$cases)),
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcis\ttrans\tshared",
      vcf
    ),
    file.path(out, "carriers.vcf")
  )
  saveRDS(list(seed = opt$seed, cases = cases), file.path(out, "inputs.rds"))
  stopifnot(
    system2("samtools", c("faidx", shQuote(file.path(out, "reference.fa")))) ==
      0L,
    system2("bgzip", shQuote(file.path(out, "model.gff3"))) == 0L,
    system2(
      "tabix",
      c("-p", "gff", shQuote(file.path(out, "model.gff3.gz")))
    ) ==
      0L
  )
  perl_lib <- paste(
    c(
      file.path(mirrors, "modules"),
      file.path(prefix, "share/ensembl-vep-116.0-0")
    ),
    collapse = ":"
  )
  run <- conda(
    "run",
    "--clean-env",
    "--env",
    paste0("PERL5LIB=", perl_lib),
    "-p",
    prefix,
    "perl",
    file.path(root, "test/duckvep/conformance/haplotype_oracle.pl"),
    file.path(out, "carriers.vcf"),
    file.path(out, "reference.fa"),
    file.path(out, "model.gff3.gz")
  )
  stopifnot(
    blit::cmd_run(
      run,
      stdout = file.path(out, "oracle.jsonl"),
      stderr = file.path(out, "oracle.log"),
      stdin = NULL,
      verbose = FALSE
    ) ==
      0L
  )
  observed <- lapply(
    readLines(file.path(out, "oracle.jsonl")),
    jsonlite::fromJSON,
    simplifyVector = FALSE
  )
  ids <- vapply(observed, `[[`, "", "transcript")
  stopifnot(
    !anyDuplicated(ids),
    setequal(ids, vapply(cases, `[[`, "", "transcript"))
  )
  names(observed) <- ids

  apply_edits <- function(case, selected) {
    edits <- case$edits[selected, , drop = FALSE]
    edits <- edits[order(edits$start, decreasing = TRUE), , drop = FALSE]
    capacity <- nchar(case$cds) + sum(nchar(edits$alt)) + 1L
    answer <- .C(
      symbol,
      charToRaw(case$cds),
      as.integer(nchar(case$cds)),
      case$strand,
      as.integer(edits$start),
      edits$ref,
      edits$alt,
      as.integer(nrow(edits)),
      cds = raw(capacity),
      as.integer(capacity),
      protein = raw(capacity),
      as.integer(capacity),
      status = integer(1L),
      cds_length = integer(1L),
      protein_length = integer(1L),
      flags = integer(1L)
    )
    stopifnot(answer$status == 0L)
    list(
      cds = rawToChar(answer$cds[seq_len(answer$cds_length)]),
      protein = rawToChar(answer$protein[seq_len(answer$protein_length)]),
      flags = sort(c("indel", "frameshift", "resolved_frameshift")[
        bitwAnd(answer$flags, c(1L, 2L, 4L)) != 0L
      ]),
      contributors = sort(edits$id)
    )
  }
  stream_edits <- function(case) {
    edits <- case$edits
    n <- nrow(edits)
    capacity <- nchar(case$cds) + sum(nchar(edits$alt)) + 1L
    lanes <- rbind(
      rep(1L, n),
      1L + as.integer(seq_len(n) %% 2L == 0L),
      rep(1L, n)
    )
    answer <- .C(
      carrier_symbol,
      charToRaw(case$cds),
      as.integer(nchar(case$cds)),
      case$strand,
      as.integer(case$positions),
      as.integer(edits$start),
      edits$ref,
      edits$alt,
      as.integer(n),
      as.integer(order(case$positions, seq_len(n)) - 1L),
      as.integer(lanes),
      as.integer(capacity),
      cds = raw(6L * capacity),
      protein = raw(6L * capacity),
      cds_lengths = integer(6L),
      protein_lengths = integer(6L),
      flags = integer(6L),
      contributors = integer(6L * n),
      contributor_counts = integer(6L),
      metrics = double(6L),
      status = integer(1L)
    )
    if (answer$status != 0L) {
      saveRDS(
        list(input = case, output = answer),
        file.path(out, "carrier_failure.rds")
      )
      stop("carrier index failed for ", case$transcript, ": ", answer$status)
    }
    paths <- lapply(seq_len(6L), function(slot) {
      offset <- (slot - 1L) * capacity
      contributors <- answer$contributors[
        (slot - 1L) * n + seq_len(answer$contributor_counts[slot])
      ]
      list(
        cds = rawToChar(answer$cds[offset + seq_len(answer$cds_lengths[slot])]),
        protein = rawToChar(answer$protein[
          offset + seq_len(answer$protein_lengths[slot])
        ]),
        flags = sort(c("indel", "frameshift", "resolved_frameshift")[
          bitwAnd(answer$flags[slot], c(1L, 2L, 4L)) != 0L
        ]),
        contributors = sort(edits$id[contributors + 1L]),
        sample = rep(c("cis", "trans", "shared"), each = 2L)[slot]
      )
    })
    list(
      paths = paths,
      metrics = data.frame(
        transcript = case$transcript,
        input_events = n,
        input_carriers = 3L * n,
        peak_transcripts = answer$metrics[1L],
        peak_calls = answer$metrics[2L],
        peak_prefixes = answer$metrics[3L],
        completed_event_paths = answer$metrics[4L],
        translated_bases = answer$metrics[5L],
        completed_carriers = answer$metrics[6L]
      )
    )
  }
  compare_haplotypes <- function(expected, oracle) {
    groups <- split(expected, vapply(expected, `[[`, "", "cds"))
    checks <- logical()
    compare <- function(condition, field) {
      checks <<- c(checks, setNames(isTRUE(condition), field))
    }
    compare(oracle$total_haplotype_count == 6L, "total_haplotype_count")
    compare(
      setequal(names(groups), vapply(oracle$haplotypes, `[[`, "", "cds")),
      "complete_CDS_set"
    )
    for (haplotype in oracle$haplotypes) {
      carriers <- groups[[haplotype$cds]]
      if (is.null(carriers)) {
        next
      }
      compare(
        all(vapply(
          carriers,
          function(x) identical(x$protein, haplotype$protein),
          TRUE
        )),
        "complete_protein"
      )
      compare(
        all(vapply(
          carriers,
          function(x) {
            identical(x$flags, sort(as.character(unlist(haplotype$flags))))
          },
          TRUE
        )),
        "frame_flags"
      )
      compare(
        setequal(
          unique(unlist(lapply(carriers, `[[`, "contributors"))),
          unlist(haplotype$contributors)
        ),
        "complete_contributors"
      )
      counts <- table(vapply(carriers, `[[`, "", "sample"))
      actual_counts <- unlist(haplotype$samples)
      compare(
        identical(sort(names(counts)), sort(names(actual_counts))) &&
          identical(
            as.integer(counts[sort(names(counts))]),
            as.integer(actual_counts[sort(names(actual_counts))])
          ),
        "sample_counts"
      )
      compare(haplotype$count == length(carriers), "carrier_count")
    }
    checks
  }
  failures <- list()
  comparisons <- 0L
  native <- vector("list", length(cases))
  carrier_metrics <- vector("list", length(cases))
  names(native) <- vapply(cases, `[[`, "", "transcript")
  for (case in cases) {
    expected <- list()
    for (sample in c("cis", "trans", "shared")) {
      for (lane in 1:2) {
        selected <- if (sample == "trans") {
          seq_len(nrow(case$edits)) %% 2L == (lane %% 2L)
        } else {
          rep(lane == 1L, nrow(case$edits))
        }
        haplotype <- apply_edits(case, selected)
        haplotype$sample <- sample
        expected[[length(expected) + 1L]] <- haplotype
      }
    }
    streamed <- stream_edits(case)
    if (!identical(expected, streamed$paths)) {
      saveRDS(
        list(input = case, direct = expected, streamed = streamed),
        file.path(out, "carrier_failure.rds")
      )
      stop(
        "carrier index disagrees with direct edit sets for ",
        case$transcript
      )
    }
    native[[case$transcript]] <- streamed$paths
    carrier_metrics[[match(case$transcript, names(native))]] <- streamed$metrics
    checks <- compare_haplotypes(streamed$paths, observed[[case$transcript]])
    comparisons <- comparisons + length(checks)
    if (any(!checks)) {
      failures[[length(failures) + 1L]] <- data.frame(
        transcript = case$transcript,
        shape = case$shape,
        strand = case$strand,
        field = names(checks)[!checks]
      )
    }
  }
  failures <- if (length(failures)) {
    do.call(rbind, failures)
  } else {
    data.frame(
      transcript = character(),
      shape = character(),
      strand = integer(),
      field = character()
    )
  }
  saveRDS(native, file.path(out, "native.rds"))
  write.csv(
    do.call(rbind, carrier_metrics),
    file.path(out, "carrier_metrics.csv"),
    row.names = FALSE
  )
  write.csv(failures, file.path(out, "failures.csv"), row.names = FALSE)
  # Exercise every comparison axis with a deliberate corruption. These seven
  # verifier controls are separate from the engine-comparison denominators.
  controls <- data.frame(
    field = c(
      "total_haplotype_count",
      "complete_CDS_set",
      "complete_protein",
      "frame_flags",
      "complete_contributors",
      "sample_counts",
      "carrier_count"
    ),
    rejected = FALSE
  )
  first <- observed[[cases[[1L]]$transcript]]
  if (all(compare_haplotypes(native[[1L]], first))) {
    for (i in seq_len(nrow(controls))) {
      mutant <- first
      switch(
        controls$field[i],
        total_haplotype_count = {
          mutant$total_haplotype_count <- first$total_haplotype_count + 1L
        },
        complete_CDS_set = {
          mutant$haplotypes[[1L]]$cds <- "N"
        },
        complete_protein = {
          mutant$haplotypes[[1L]]$protein <- "WRONG"
        },
        frame_flags = {
          mutant$haplotypes[[1L]]$flags <- list("WRONG")
        },
        complete_contributors = {
          mutant$haplotypes[[1L]]$contributors <- list("WRONG")
        },
        sample_counts = {
          mutant$haplotypes[[1L]]$samples <- list(WRONG = 1L)
        },
        carrier_count = {
          mutant$haplotypes[[1L]]$count <- first$haplotypes[[1L]]$count + 1L
        }
      )
      checks <- compare_haplotypes(native[[1L]], mutant)
      controls$rejected[i] <- any(!checks[names(checks) == controls$field[i]])
    }
  }
  write.csv(
    controls,
    file.path(out, "verifier_controls.csv"),
    row.names = FALSE
  )
  summary <- data.frame(
    seed = opt$seed,
    transcripts = length(cases),
    input_records = length(vcf),
    input_allele_slots = 6L * length(vcf),
    haplotype_lanes = 6L * length(cases),
    comparisons = comparisons,
    failures = nrow(failures),
    verifier_controls_rejected = sum(controls$rejected)
  )
  write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
  identities <- c(
    sources,
    file.path(root, "test/duckvep/conformance/haplotype_differential.R"),
    file.path(root, "test/duckvep/conformance/haplotype_oracle.pl"),
    file.path(
      root,
      c(
        "scripts/duckvep_evidence.R",
        "src/duckvep/kernel/src/duckvep_haplotype.h",
        "src/duckvep/kernel/src/duckvep_carriers.h",
        "src/duckvep/kernel/include/duckvep_kernel.h",
        "src/duckvep/kernel/src/duckvep_codon.h",
        "src/duckvep/kernel/src/duckvep_dna.h"
      )
    ),
    lock,
    file.path(
      out,
      c(
        "reference.fa",
        "reference.fa.fai",
        "model.gff3.gz",
        "model.gff3.gz.tbi",
        "carriers.vcf",
        "inputs.rds",
        "native.rds",
        "carrier_metrics.csv",
        "oracle.jsonl",
        "summary.csv",
        "failures.csv",
        "verifier_controls.csv",
        "environment.txt"
      )
    ),
    shared
  )
  hashes <- vapply(identities, duckvep_evidence_sha256, "")
  jsonlite::write_json(
    list(
      scope = "diploid_carrier_prefix_and_edit_mechanics_not_public_phased_execution",
      source_revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
      dirty_worktree = length(system2(
        "git",
        c("status", "--porcelain"),
        stdout = TRUE
      )) !=
        0L,
      oracle_revisions = as.list(revisions),
      seed = opt$seed,
      compiler = system2(compiler, "--version", stdout = TRUE),
      sha256 = as.list(hashes)
    ),
    file.path(out, "receipt.json"),
    auto_unbox = TRUE,
    pretty = TRUE
  )
  print(summary, row.names = FALSE)
  if (!all(controls$rejected)) {
    stop("Haplosaurus comparison controls failed; inspect ", out, call. = FALSE)
  }
  if (nrow(failures)) {
    stop(
      "Haplosaurus mechanics differential failed; inspect ",
      out,
      call. = FALSE
    )
  }
}
main()
