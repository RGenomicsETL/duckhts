source("benchmarks/count_mixture_lowpass/stage.R")

run_staging_test <- function(extension) {
  cache <- tempfile("lowpass-cache-")
  dir.create(cache)
  old <- Sys.getenv("DUCKHTS_CACHE_DIR", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("DUCKHTS_CACHE_DIR") else Sys.setenv(DUCKHTS_CACHE_DIR = old)
    unlink(cache, recursive = TRUE)
  })
  Sys.setenv(DUCKHTS_CACHE_DIR = cache)
  ids <- c(lowpass_declaration$source_artifacts, lowpass_declaration$count_artifacts)
  paths <- setNames(vapply(ids, duckhts_bench_artifact_path, ""), ids)
  stopifnot(all(startsWith(paths, paste0(cache, "/"))), !any(file.exists(paths)),
            startsWith(duckhts_bench_artifact_path("count_mixture_lowpass_chr20_haplotypes"), paste0(cache, "/")))
  for (path in paths) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  pedigree <- data.frame(FamilyID = paste0("F", 1:120), SampleID = paste0("S", 1:120),
                         FatherID = "0", MotherID = "0", Population = rep(c("CEU", "YRI"), each = 60))
  write.table(pedigree, paths[["roh_ancestry_pedigree"]], quote = FALSE, row.names = FALSE, sep = "\t")
  genotypes <- rbind(c("0|0", "0|1", "1|1", "0|1"),
                     c("1|0", "1|1", "0|0", "1|0"),
                     c("0|1", "1|0", "0|1", "1|1"))
  genotypes <- genotypes[, rep(1:4, each = 30)]
  header <- c("##fileformat=VCFv4.2", "##contig=<ID=20,length=4000000>",
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
              paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", pedigree$SampleID), collapse = "\t"))
  positions <- c(1L, 1000001L, 2000001L)
  records <- vapply(seq_along(positions), function(i) paste(
    c("20", positions[[i]], ".", "A", "C", ".", "PASS", ".", "GT", genotypes[i, ]), collapse = "\t"), "")
  compressed <- gzfile(paths[["roh_ancestry_chr20_source"]], "wt")
  writeLines(c(header, records), compressed)
  close(compressed)
  con <- lowpass_connect(extension)
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  sites <- data.frame(chrom = "20", pos = positions, ref = "A", alt = "C")
  dbWriteTable(con, "sites", sites)
  dbExecute(con, paste("COPY sites TO", dbQuoteString(con, paths[["roh_ancestry_chr20_af_sites"]]), "(FORMAT PARQUET)"))
  for (id in lowpass_declaration$count_artifacts) {
    dbExecute(con, paste("COPY sites TO", dbQuoteString(con, paths[[id]]), "(FORMAT PARQUET)"))
  }
  destination <- lowpass_stage_haplotypes(extension)
  observed <- readRDS(destination)
  stopifnot(identical(as.integer(observed$sites$pos), positions), length(observed$genotypes) == 120L,
            identical(observed$sites$f, c(.25, .75, .5)), identical(observed$sites$g, c(.75, .25, .75)),
            identical(as.integer(observed$genotypes$S1), c(0L, 2L, 1L)))
  stopifnot(identical(lowpass_stage_haplotypes(extension), destination))
  unlink(paths[["roh_ancestry_chr20_source"]])
  unlink(destination)
  missing <- try(lowpass_stage_haplotypes(extension), silent = TRUE)
  stopifnot(inherits(missing, "try-error"), grepl("Missing registered inputs", missing))
  cat("low-pass staging: 11 checks, no network, registered derivation and missing-input failure OK\n")
}
run_staging_test(commandArgs(TRUE)[[1L]])
