# Generate transcript-oriented CDS splits and genomic source records from the
# registered reference. GFF phase is the complement of Ensembl exon phase.
haplotype_geometry_cases <- function(cds, quota, region_offset, transcript_offset,
    interaction = FALSE, prefix = if (interaction) "P" else "G") {
  stopifnot(is.character(cds), length(cds) == 1L, !is.na(cds), nchar(cds) >= 36L,
    grepl("^[ACGT]+$", cds), length(quota) == 1L, !is.na(quota),
    quota >= 1L, quota == as.integer(quota),
    length(prefix) == 1L, !is.na(prefix), grepl("^[A-Za-z0-9_]+$", prefix))
  cds_length <- nchar(cds)
  rc <- function(x) paste(rev(strsplit(chartr("ACGT", "TGCA", x), "", fixed = TRUE)[[1L]]), collapse = "")
  strata <- expand.grid(
    exon_count = c(2L, 3L, 5L, 7L), split_phase = 0:2,
    utr = c("none", "inside", "separate"), strand = c(1L, -1L),
    shape = c(
      "coding", "cds_start", "cds_end", "exon_exit", "exon_entry",
      "whole_exon", "junction_insertion"
    ), stringsAsFactors = FALSE
  )
  if (interaction) strata <- merge(strata, expand.grid(
    second_location = c("overlap", "same_exon", "other_exon"),
    second_edit = c("snv", "insertion", "deletion", "replacement"),
    stringsAsFactors = FALSE
  ), by = NULL)
  cases <- strata[rep(seq_len(nrow(strata)), each = quota), ]
  cases$draw <- rep(seq_len(quota), nrow(strata))
  models <- exons <- records <- fasta <- gff <- vector("list", nrow(cases))
  for (i in seq_len(nrow(cases))) {
    profile <- cases[i, ]
    tx_index <- transcript_offset + i - 1L
    region <- region_offset + i - 1L
    tx <- paste0("T", prefix, i)
    chr <- sprintf(paste0("chr", prefix, "%04d"), i)
    first <- sample(seq(12L + profile$split_phase,
      min(32L, cds_length - profile$exon_count + 1L), 3L), 1L)
    cuts <- c(0L, first, sort(sample(seq.int(first + 1L, cds_length - 1L),
      profile$exon_count - 2L)), cds_length)
    u5 <- if (profile$utr == "inside") sample(1:9, 1L) else 0L
    u3 <- if (profile$utr == "inside") sample(1:9, 1L) else 0L
    genome <- strrep("A", 10L)
    span <- data.frame(start = integer(), end = integer(), coding_start = integer(), coding_end = integer())
    if (profile$utr == "separate") {
      genome <- paste0(genome, "ACACACA", strrep("A", sample(2:15, 1L)))
      span <- rbind(span, data.frame(start = 11L, end = 17L, coding_start = NA_integer_, coding_end = NA_integer_))
    }
    for (j in seq_len(profile$exon_count)) {
      s <- nchar(genome) + 1L
      part <- paste0(
        if (j == 1L) strrep("C", u5) else "",
        substr(cds, cuts[j] + 1L, cuts[j + 1L]), if (j == profile$exon_count) strrep("G", u3) else ""
      )
      genome <- paste0(genome, part)
      e <- nchar(genome)
      span <- rbind(span, data.frame(
        start = s, end = e, coding_start = s + if (j == 1L) u5 else 0L,
        coding_end = e - if (j == profile$exon_count) u3 else 0L
      ))
      if (j < profile$exon_count) genome <- paste0(genome, strrep("A", sample(1:24, 1L)))
    }
    if (profile$utr == "separate") {
      genome <- paste0(genome, strrep("A", sample(2:15, 1L)))
      s <- nchar(genome) + 1L
      genome <- paste0(genome, "GAGAGAG")
      span <- rbind(span, data.frame(start = s, end = nchar(genome), coding_start = NA_integer_, coding_end = NA_integer_))
    }
    genome <- paste0(genome, strrep("A", 10L))
    coding <- span[!is.na(span$coding_start), ]
    j <- sample(seq_len(nrow(coding)), 1L)
    pos <- switch(profile$shape,
      coding = coding$coding_start[j],
      cds_start = coding$coding_start[1L] - 1L,
      cds_end = tail(coding$coding_end, 1L) - 1L,
      exon_exit = coding$coding_end[j],
      exon_entry = coding$coding_start[j] - 1L,
      whole_exon = coding$coding_start[j] - 1L,
      junction_insertion = coding$coding_end[j]
    )
    ref_length <- switch(profile$shape,
      coding = 1L,
      cds_start = 3L,
      cds_end = 3L,
      exon_exit = 2L,
      exon_entry = 2L,
      whole_exon = coding$coding_end[j] - coding$coding_start[j] + 3L,
      junction_insertion = 1L
    )
    ref <- substr(genome, pos, pos + ref_length - 1L)
    alt <- if (profile$shape == "junction_insertion") paste0(ref, "GC") else chartr("ACGT", "CGTA", ref)
    anchor <- tail(coding$coding_start, 1L)
    anchor_ref <- substr(genome, anchor, anchor)
    r <- data.frame(
      position = c(pos, anchor), reference = c(ref, anchor_ref),
      alt = c(alt, chartr("ACGT", "CGTA", anchor_ref)), source_id = c("edge", "anchor"),
      s0 = c("1|0", "1|1"), s1 = c("0|1", "1|1"), s2 = c(".|1", "1|1")
    )
    if (interaction) {
      edge_exon <- switch(profile$shape, cds_start = 1L, cds_end = nrow(coding), j)
      other_exons <- setdiff(seq_len(nrow(coding)), edge_exon)
      partner_exon <- if (profile$second_location == "other_exon")
        other_exons[sample.int(length(other_exons), 1L)] else edge_exon
      partner_positions <- if (profile$second_location == "overlap")
        seq.int(max(pos, min(span$start)), min(pos + ref_length - 1L, max(span$end))) else
        seq.int(coding$coding_start[partner_exon], coding$coding_end[partner_exon])
      partner_pos <- partner_positions[sample.int(length(partner_positions), 1L)]
      partner_length <- switch(profile$second_edit,
        snv = 1L, insertion = 1L, deletion = sample(2:5, 1L), replacement = sample(2:6, 1L))
      partner_ref <- substr(genome, partner_pos, partner_pos + partner_length - 1L)
      partner_alt <- switch(profile$second_edit,
        insertion = paste0(partner_ref, paste0(sample(c("A", "C", "G", "T"),
          sample(1:3, 1L), replace = TRUE), collapse = "")),
        deletion = substr(partner_ref, 1L, 1L), chartr("ACGT", "CGTA", partner_ref))
      stopifnot(nchar(partner_ref) == partner_length, partner_ref != partner_alt,
        profile$second_location != "other_exon" || partner_exon != edge_exon)
      partner <- data.frame(position = partner_pos, reference = partner_ref,
        alt = partner_alt, source_id = "partner", s0 = "1|0", s1 = "1|0", s2 = "1|1")
      r <- rbind(r[1L, ], partner, r[2L, ])
      stopifnot(all(r$position <= max(span$end) &
        r$position + nchar(r$reference) - 1L >= min(span$start)))
    }
    span$ce <- cumsum(span$end - span$start + 1L)
    span$cs <- c(1L, head(span$ce, -1L) + 1L)
    coding_before <- cumsum(c(0L, head(ifelse(is.na(span$coding_start), 0L,
      span$coding_end - span$coding_start + 1L
    ), -1L)))
    span$phase <- ifelse(is.na(span$coding_start), -1L, coding_before %% 3L)
    span$end_phase <- ifelse(is.na(span$coding_start), -1L,
      (coding_before + span$coding_end - span$coding_start + 1L) %% 3L
    )
    if (profile$strand < 0L) {
      n <- nchar(genome)
      span[c("start", "end")] <- data.frame(start = n - span$end + 1L, end = n - span$start + 1L)
      span[c("coding_start", "coding_end")] <- data.frame(
        coding_start = n - span$coding_end + 1L,
        coding_end = n - span$coding_start + 1L
      )
      r$position <- n - r$position - nchar(r$reference) + 2L
      r$reference <- vapply(r$reference, rc, "")
      r$alt <- vapply(r$alt, rc, "")
      genome <- rc(genome)
    }
    r$chrom <- chr
    r$seq_region <- region
    records[[i]] <- r
    span$transcript_index <- tx_index
    exons[[i]] <- span
    models[[i]] <- data.frame(
      transcript_index = tx_index, seq_region = region, transcript = tx,
      strand = profile$strand, transcript_start = min(span$start), transcript_end = max(span$end),
      cds_start = min(span$coding_start, na.rm = TRUE), cds_end = max(span$coding_end, na.rm = TRUE), cds = cds
    )
    gf <- data.frame(
      feature = c("gene", "mRNA"), start = min(span$start), end = max(span$end), phase = ".",
      attr = c(
        paste0("ID=gene:", tx, ";biotype=protein_coding"),
        paste0("ID=transcript:", tx, ";Parent=gene:", tx, ";biotype=protein_coding")
      )
    )
    for (j in seq_len(nrow(span))) {
      gf <- rbind(gf, data.frame(
        feature = "exon", start = span$start[j], end = span$end[j], phase = ".",
        attr = paste0("ID=exon:", tx, "_", j, ";Parent=transcript:", tx)
      ))
      if (!is.na(span$coding_start[j])) {
        gf <- rbind(gf, data.frame(
          feature = "CDS",
          start = span$coding_start[j], end = span$coding_end[j], phase = (3L - span$phase[j]) %% 3L,
          attr = paste0("Parent=transcript:", tx)
        ))
      }
    }
    gff[[i]] <- with(gf, paste(chr, "geometry", feature, start, end, ".", if (profile$strand > 0) "+" else "-", phase, attr, sep = "\t"))
    fasta[[i]] <- c(paste0(">", chr), genome)
  }

  models <- do.call(rbind, models)
  exons <- do.call(rbind, exons)
  records <- do.call(rbind, records)
  layouts <- models[c("transcript_index", "transcript_start", "transcript_end", "cds_start", "cds_end")]
  phases <- exons[c("transcript_index", "cs", "phase", "end_phase")]
  models$kind <- "geometry"
  models <- models[c("transcript_index", "seq_region", "transcript", "strand", "kind", "cds")]
  exons <- exons[c("transcript_index", "start", "end", "ce", "cs")]
  key <- function(x) do.call(paste, c(x[names(strata)], sep = ":"))
  coverage <- strata
  coverage$required <- quota
  coverage$observed <- tabulate(match(key(cases), key(strata)), nrow(strata))
  stopifnot(nrow(strata) == if (interaction) 6048L else 504L,
    all(coverage$observed == quota))
  list(
    cases = cases, models = models, exons = exons, records = records,
    layouts = layouts, phases = phases, fasta = fasta, gff = gff, coverage = coverage
  )
}

# Repeat the registered CDS's internal codons, retaining ATG and a raw TAA suffix.
# Partial terminal codons exercise the pinned reference/alternate translation rules.
haplotype_length_cases <- function(cds, quota, region_offset, transcript_offset) {
  stopifnot(length(cds) == 1L, !is.na(cds), nchar(cds) == 180L,
    startsWith(cds, "ATG"), endsWith(cds, "TAA"), grepl("^[ACGT]+$", cds))
  lengths <- c(36L, 37L, 38L, 2047L, 2048L, 2049L, 6143L, 6144L, 6145L)
  internal <- substr(cds, 4L, nchar(cds) - 3L)
  cohorts <- vector("list", length(lengths))
  for (i in seq_along(lengths)) {
    n <- lengths[i]
    sequence <- paste0("ATG", substr(strrep(internal, ceiling((n - 6L) / nchar(internal))),
      1L, n - 6L), "TAA")
    added <- haplotype_geometry_cases(sequence, quota, region_offset, transcript_offset,
      prefix = paste0("L", n, "_"))
    added$cases$cds_length <- added$coverage$cds_length <- n
    cohorts[[i]] <- added
    region_offset <- region_offset + nrow(added$models)
    transcript_offset <- transcript_offset + nrow(added$models)
  }
  result <- lapply(names(cohorts[[1L]]), function(name)
    do.call(if (name %in% c("fasta", "gff")) c else rbind, lapply(cohorts, `[[`, name)))
  names(result) <- names(cohorts[[1L]])
  stopifnot(nrow(result$coverage) == 4536L,
    all(result$coverage$required == quota), all(result$coverage$observed == quota),
    identical(nchar(result$models$cds), result$cases$cds_length))
  result
}
