# Generate transcript-oriented CDS splits and genomic source records from the
# registered reference. GFF phase is the complement of Ensembl exon phase.
haplotype_geometry_cases <- function(cds, quota, region_offset, transcript_offset) {
  stopifnot(nchar(cds) == 180L, quota >= 1L)
  rc <- function(x) paste(rev(strsplit(chartr("ACGT", "TGCA", x), "", fixed = TRUE)[[1L]]), collapse = "")
  strata <- expand.grid(
    exon_count = c(2L, 3L, 5L, 7L), split_phase = 0:2,
    utr = c("none", "inside", "separate"), strand = c(1L, -1L),
    shape = c(
      "coding", "cds_start", "cds_end", "exon_exit", "exon_entry",
      "whole_exon", "junction_insertion"
    ), stringsAsFactors = FALSE
  )
  cases <- strata[rep(seq_len(nrow(strata)), each = quota), ]
  cases$draw <- rep(seq_len(quota), nrow(strata))
  models <- exons <- records <- fasta <- gff <- vector("list", nrow(cases))
  for (i in seq_len(nrow(cases))) {
    profile <- cases[i, ]
    tx_index <- transcript_offset + i - 1L
    region <- region_offset + i - 1L
    tx <- paste0("TG", i)
    chr <- sprintf("chrG%04d", i)
    first <- sample(seq(12L + profile$split_phase, 32L, 3L), 1L)
    cuts <- c(0L, first, sort(sample(seq.int(first + 1L, 179L), profile$exon_count - 2L)), 180L)
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
  stopifnot(nrow(strata) == 504L, all(coverage$observed == quota))
  list(
    cases = cases, models = models, exons = exons, records = records,
    layouts = layouts, phases = phases, fasta = fasta, gff = gff, coverage = coverage
  )
}
