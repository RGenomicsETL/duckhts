Sys.setenv(DUCKHTSBENCH_REGISTRY = normalizePath("r/duckhtsbench/inst/benchmark_registry.tsv"))
source("benchmarks/count_mixture_lowpass/common.R")
source("r/duckhtsbench/R/registry.R")
source("r/duckhtsbench/R/stage.R")

genome_settings <- list(
  seed = 332202611L, chromosomes = 1:22,
  populations = c(sample = "CEU", other = "YRI"), block_bases = 10000000L,
  min_blocks = 30L, max_depth = 8L, artefact_weight = .002
)

genome_id <- function(chrom) paste0("count_mixture_lowpass_chr", chrom, "_children_haplotypes")

genome_connect <- function(extension, threads = 1L) {
  con <- lowpass_connect(extension, threads)
  dbExecute(con, paste("SET memory_limit =", dbQuoteString(con, "6GB")))
  dbExecute(con, "SET preserve_insertion_order = false")
  con
}

genome_fit <- function(con, error) {
  dbGetQuery(con, sprintf(
    paste(
      "SELECT * FROM duckhts_count_mixture_lowpass('counts', %.17g, %.17g, 0,",
      "block_bases := %d, min_blocks := %d) ORDER BY candidate, relationship_model"
    ),
    error, genome_settings$artefact_weight, genome_settings$block_bases,
    genome_settings$min_blocks
  ))
}

genome_input_identity <- function() {
  ids <- c("roh_ancestry_pedigree", unlist(lapply(
    genome_settings$chromosomes,
    function(chrom) paste0("roh_ancestry_chr", chrom, c("_children_bcf", "_af_sites"))
  )))
  registry <- duckhts_bench_registry()
  rows <- registry[match(ids, registry$id), ]
  stopifnot(!anyNA(rows$id))
  rows[c("id", "release", "supplier_identity")]
}

genome_select_children <- function(bcf) {
  ids <- system2(RBCFTools::bcftools_path(), c("query", "-l", shQuote(bcf)), stdout = TRUE)
  ped <- read.table(duckhts_bench_artifact_path("roh_ancestry_pedigree"),
    header = TRUE,
    colClasses = "character", stringsAsFactors = FALSE
  )
  ped <- ped[match(ids, ped$SampleID), ]
  stopifnot(!anyNA(ped), all(ped$FatherID != "0" & ped$MotherID != "0"))
  ped <- ped[order(ped$SampleID), ]
  ped <- ped[!duplicated(ped$FamilyID) & ped$Population %in% genome_settings$populations, ]
  stopifnot(!anyDuplicated(ped$FamilyID))
  selected <- lapply(genome_settings$populations, function(pop) ped$SampleID[ped$Population == pop])
  stopifnot(all(lengths(selected) >= 2L))
  list(selected = selected, pedigree = ped)
}

genome_simulate <- function(panel, setting, identity) {
  sites <- panel$sites
  h <- panel$genotypes[[identity$sample]]
  dosage <- function(code) code %/% 2L + code %% 2L
  a <- dosage(as.integer(h)) / 2
  donor <- as.integer(panel$genotypes[[identity$donor]])
  if (setting$relationship == "ibd1") {
    transmitted <- if (identity$sample_haplotype == 1L) as.integer(h) %/% 2L else as.integer(h) %% 2L
    independent <- if (identity$donor_haplotype == 1L) donor %/% 2L else donor %% 2L
    b <- (transmitted + independent) / 2
  } else {
    b <- dosage(donor) / 2
  }
  p <- (1 - setting$share) * a + setting$share * b
  p <- setting$error + (1 - 2 * setting$error) * p
  n <- rpois(nrow(sites), setting$depth)
  alt <- rbinom(length(n), n, p)
  artefact <- runif(length(n)) < genome_settings$artefact_weight
  alt[artefact] <- floor(runif(sum(artefact)) * (n[artefact] + 1L))
  data.frame(
    chrom = sites$chrom, pos = sites$pos, ref_count = n - alt,
    alt_count = alt, other_count = 0L, sample_af = sites$f,
    candidate_af = if (setting$relationship == "other") sites$g else sites$f
  )
}

genome_draw_identity <- function(selected, relationship) {
  sample_id <- sample(selected$sample, 1L)
  pool <- if (relationship == "other") selected$other else setdiff(selected$sample, sample_id)
  list(
    sample = sample_id, donor = sample(pool, 1L), sample_haplotype = sample(1:2, 1L),
    donor_haplotype = sample(1:2, 1L)
  )
}
