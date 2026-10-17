source("benchmarks/count_mixture_lowpass/common.R")
args <- commandArgs(TRUE)
kind <- args[[1L]]
id <- as.integer(args[[2L]])
output <- args[[3L]]
extension <- args[[4L]]
con <- lowpass_connect(extension)
lowpass_moment_macro(con)
set.seed(lowpass_declaration$seed + id + if (kind == "thinning") 10000L else 0L)

hap <- function(x, copy) if (copy == 1L) x %/% 2L else x %% 2L

if (kind == "synthetic") {
  stratum <- lowpass_strata()[id, ]
  genomes <- readRDS(duckhts_bench_artifact_path("count_mixture_lowpass_chr20_haplotypes"))
  sites <- genomes$sites
  error <- stratum$seq_error
  replicates <- lowpass_declaration$replicates
  simulate <- function() {
    founders <- sample(genomes$selected$sample, 4L)
    copies <- sample.int(2L, 4L, replace = TRUE)
    first <- hap(genomes$genotypes[[founders[[1L]]]], copies[[1L]])
    second <- hap(genomes$genotypes[[founders[[2L]]]], copies[[2L]])
    receiver <- first + second
    if (stratum$relationship == "ibd1") {
      transmitted <- if (sample.int(2L, 1L) == 1L) first else second
      donor <- transmitted + hap(genomes$genotypes[[founders[[3L]]]], copies[[3L]])
    } else {
      donors <- if (stratum$relationship == "other") sample(genomes$selected$other, 2L) else founders[3:4]
      donor <- hap(genomes$genotypes[[donors[[1L]]]], copies[[3L]]) +
        hap(genomes$genotypes[[donors[[2L]]]], copies[[4L]])
    }
    n <- rpois(nrow(sites), stratum$depth)
    probability <- error + (1-2*error)*((1-stratum$share)*receiver + stratum$share*donor)/2
    alt <- rbinom(length(n), n, probability)
    artefact <- runif(length(n)) < lowpass_declaration$artefact_weight
    alt[artefact] <- floor(runif(sum(artefact))*(n[artefact]+1))
    data.frame(sample_id = "s", candidate = stratum$relationship, chrom = "20", pos = sites$pos,
               ref_count = n-alt, alt_count = alt, other_count = 0L, sample_af = sites$f,
               candidate_af = if (stratum$relationship == "other") sites$g else sites$f)
  }
} else if (kind == "thinning") {
  strata <- lowpass_thinning_strata()
  stratum <- strata[id, ]
  paths <- setNames(vapply(lowpass_declaration$count_artifacts, duckhts_bench_artifact_path, ""),
                    lowpass_declaration$count_samples)
  receiver <- dbGetQuery(con, paste("SELECT * FROM read_parquet(", dbQuoteString(con, paths[[stratum$receiver]]), ") ORDER BY pos"))
  donor <- dbGetQuery(con, paste("SELECT * FROM read_parquet(", dbQuoteString(con, paths[[stratum$donor]]), ") ORDER BY pos"))
  stopifnot(identical(receiver$pos, donor$pos), identical(receiver$ref, donor$ref), identical(receiver$alt, donor$alt))
  error <- lowpass_declaration$thinning_error
  replicates <- lowpass_declaration$thinning_replicates
  simulate <- function() {
    keep_receiver <- stratum$probability*(1-stratum$share)
    keep_donor <- stratum$probability*stratum$share
    thin <- function(column) rbinom(nrow(receiver), receiver[[column]], keep_receiver) +
      rbinom(nrow(donor), donor[[column]], keep_donor)
    data.frame(sample_id = "s", candidate = "donor", chrom = receiver$chrom, pos = receiver$pos,
               ref_count = thin("ref_count"), alt_count = thin("alt_count"), other_count = thin("other_count"),
               sample_af = receiver$af, candidate_af = donor$af)
  }
} else stop("Unknown worker kind")

rows <- vector("list", replicates)
checks <- NULL
for (replicate in seq_len(replicates)) {
  counts <- simulate()
  dbWriteTable(con, "counts", counts, overwrite = TRUE)
  fit <- lowpass_fit(con, error)
  if (replicate == 1L) {
    moments <- lowpass_moments(con)
    oracle <- lowpass_oracle(counts, moments)
    dbExecute(con, "SET threads=4")
    four <- lowpass_fit(con, error)
    four_moments <- lowpass_moments(con)
    thread <- lowpass_thread_check(fit, four, moments, four_moments)
    dbExecute(con, "SET threads=1")
    checks <- cbind(data.frame(kind = kind, id = id), as.data.frame(c(oracle, thread)))
  }
  fit$replicate <- replicate
  fit$id <- id
  fit$primary <- if (kind == "synthetic") {
    fit$relationship_model == if (stratum$relationship == "ibd1") "ibd1" else "unrelated"
  } else fit$relationship_model == "unrelated"
  rows[[replicate]] <- fit
}
saveRDS(list(stratum = stratum, fits = do.call(rbind, rows), checks = checks), output, compress = "gzip")
dbDisconnect(con, shutdown = TRUE)
