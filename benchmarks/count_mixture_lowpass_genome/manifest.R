source("benchmarks/count_mixture_lowpass_genome/common.R")

genome_manifest <- function(extension, destination) {
  source_identity <- genome_input_identity()
  phase <- list()
  for (id in source_identity$id) {
    path <- duckhts_bench_artifact_path(id)
    duckhts_bench_validate_identity(id, path)
    if (grepl("_children_bcf$", id)) {
      temporary <- tempfile()
      status <- system2(RBCFTools::bcftools_path(), c(
        "query", "-i", shQuote('GT~"/"'), "-f", shQuote("%CHROM\\t%POS\\n"), shQuote(path)
      ), stdout = temporary)
      counterexamples <- readLines(temporary, n = 10L)
      unlink(temporary)
      if (status != 0L || length(counterexamples)) stop("Children phase audit failed: ", id)
      phase[[length(phase) + 1L]] <- data.frame(
        chrom = as.integer(sub(".*chr([0-9]+)_children_bcf", "\\1", id)),
        artifact = id, unphased_records = 0L
      )
    }
  }
  artifacts <- do.call(rbind, lapply(genome_settings$chromosomes, function(chrom) {
    path <- duckhts_bench_artifact_path(genome_id(chrom))
    panel <- readRDS(path)
    expected <- source_identity[source_identity$id %in% c(
      "roh_ancestry_pedigree",
      paste0("roh_ancestry_chr", chrom, c("_children_bcf", "_af_sites"))
    ), ]
    stopifnot(identical(panel$input_identity, expected), identical(panel$settings, genome_settings))
    data.frame(
      chrom = chrom, artifact = genome_id(chrom), sha256 = digest::digest(file = path, algo = "sha256"),
      sites = nrow(panel$sites), bytes = file.info(path)$size
    )
  }))
  sources <- c(
    "benchmarks/count_mixture_lowpass/common.R",
    list.files("benchmarks/count_mixture_lowpass_genome", "\\.R$", full.names = TRUE)
  )
  manifest <- list(
    revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    extension_sha256 = digest::digest(file = extension, algo = "sha256"), settings = genome_settings,
    input_identity = source_identity, artifacts = artifacts, phase = do.call(rbind, phase),
    driver_sha256 = setNames(vapply(sources, function(path) digest::digest(file = path, algo = "sha256"), ""), sources),
    runtime = list(R = R.version.string, duckdb = as.character(packageVersion("duckdb")))
  )
  saveRDS(manifest, destination)
  invisible(manifest)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  genome_manifest(args[[1L]], args[[2L]])
}
