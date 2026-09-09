#!/usr/bin/env Rscript
# Generator contracts: genomic REF, spliced CDS, ranked phases and cell quotas.
source("test/duckvep/conformance/haplotype_model_geometry.R")

check_geometry <- function(generated, cds) {
  n <- nrow(generated$models)
  stopifnot(nrow(generated$cases) == n, sum(generated$coverage$observed) == n,
            all(generated$coverage$required == 1L), all(generated$coverage$observed == 1L),
            identical(generated$models$transcript_index, 10L + seq_len(n)),
            identical(generated$models$seq_region, 6L + seq_len(n)),
            !anyDuplicated(generated$models$transcript),
            !anyDuplicated(vapply(generated$fasta, `[`, "", 1L)))
  exons <- split(generated$exons, generated$exons$transcript_index)
  phases <- split(generated$phases, generated$phases$transcript_index)
  records <- split(generated$records, generated$records$seq_region)
  reverse_complement <- function(x)
    paste0(rev(strsplit(chartr("ACGT", "TGCA", x), "", fixed = TRUE)[[1L]]), collapse = "")
  for (i in seq_len(nrow(generated$models))) {
    model <- generated$models[i, ]
    layout <- generated$layouts[i, ]
    exon <- exons[[as.character(model$transcript_index)]]
    phase <- phases[[as.character(model$transcript_index)]]
    genome <- generated$fasta[[i]][2L]
    record <- records[[as.character(model$seq_region)]]
    stopifnot(identical(substring(genome, record$position,
      record$position + nchar(record$reference) - 1L), record$reference))
    first <- pmax(exon$start, layout$cds_start)
    last <- pmin(exon$end, layout$cds_end)
    lengths <- pmax(0L, last - first + 1L)
    parts <- substring(genome, first[lengths > 0L], last[lengths > 0L])
    if (model$strand < 0L) parts <- vapply(parts, reverse_complement, "", USE.NAMES = FALSE)
    stopifnot(identical(paste0(parts, collapse = ""), cds[i]), identical(model$cds, cds[i]),
              sum(lengths) == nchar(cds[i]), all(exon$cs == c(1L, head(exon$ce, -1L) + 1L)),
              all(exon$ce - exon$cs == exon$end - exon$start))
    before <- c(0L, head(cumsum(lengths), -1L))
    stopifnot(identical(as.integer(phase$phase), ifelse(lengths > 0L, before %% 3L, -1L)),
              identical(as.integer(phase$end_phase), ifelse(lengths > 0L, (before + lengths) %% 3L, -1L)))
    gff <- do.call(rbind, strsplit(generated$gff[[i]], "\t", fixed = TRUE))
    coding <- gff[gff[, 3L] == "CDS", , drop = FALSE]
    stopifnot(identical(as.integer(coding[, 4L]), first[lengths > 0L]),
              identical(as.integer(coding[, 5L]), last[lengths > 0L]),
              identical(as.integer(coding[, 8L]), (3L - before[lengths > 0L] %% 3L) %% 3L),
              all(gff[, 1L] == substring(generated$fasta[[i]][1L], 2L)),
              all(gff[, 7L] == if (model$strand > 0L) "+" else "-"))
  }
  invisible(generated)
}

for (cds_length in c(180L, 36L, 37L, 38L, 2047L, 2048L, 2049L, 6143L, 6144L, 6145L)) {
  cds <- substr(paste0("ATG", strrep("GCCATC", ceiling(cds_length / 6))), 1L, cds_length)
  set.seed(173L)
  generated <- haplotype_geometry_cases(cds, 1L, 7L, 11L)
  stopifnot(nrow(generated$models) == 504L)
  check_geometry(generated, rep(cds, 504L))
}
for (cds in list(character(), NA_character_, "", strrep("A", 35L),
                rep(strrep("A", 36L), 2L), strrep("?", 36L))) {
  error <- tryCatch({haplotype_geometry_cases(cds, 1L, 0L, 0L); NULL}, error = identity)
  stopifnot(inherits(error, "error"))
}
cds <- readLines("test/data/duckvep/haplotype_benchmark.fa")[2L]
set.seed(20260909L)
generated <- haplotype_length_cases(cds, 1L, 7L, 11L)
stopifnot(nrow(generated$models) == 4536L, nrow(generated$coverage) == 4536L,
          identical(sort(unique(generated$cases$cds_length)),
            c(36L, 37L, 38L, 2047L, 2048L, 2049L, 6143L, 6144L, 6145L)),
          identical(nchar(generated$models$cds), generated$cases$cds_length),
          all(startsWith(generated$models$cds, "ATG")),
          all(endsWith(generated$models$cds, "TAA")))
check_geometry(generated, generated$models$cds)
corrupt <- list(ref = generated, model = generated, phase = generated,
                gff = generated, quota = generated)
corrupt$ref$records$reference[1L] <- "INVALID"
corrupt$model$models$cds[1L] <- paste0(generated$models$cds[1L], "A")
corrupt$phase$phases$phase[1L] <- 9L
corrupt$gff$gff[[1L]] <- sub("\tCDS\t", "\tintron\t", generated$gff[[1L]], fixed = TRUE)
corrupt$quota$coverage$observed[1L] <- 0L
for (bad in corrupt) {
  error <- tryCatch({check_geometry(bad, generated$models$cds); NULL}, error = identity)
  stopifnot(inherits(error, "error"))
}
cat("Haplotype geometry: 9,576 models retain genomic REF, spliced CDS, phases and quotas; 5 corruption controls: OK\n")
