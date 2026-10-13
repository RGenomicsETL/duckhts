source("benchmarks/count_mixture_lowpass/common.R")

lowpass_stage_paths <- function() {
  ids <- c(lowpass_declaration$source_artifacts, lowpass_declaration$count_artifacts)
  paths <- setNames(vapply(ids, duckhts_bench_artifact_path, character(1L)), ids)
  if (!all(file.exists(paths))) stop("Missing registered inputs: ", paste(names(paths)[!file.exists(paths)], collapse = ", "))
  paths
}

lowpass_stage_haplotypes <- function(extension, destination = duckhts_bench_artifact_path("count_mixture_lowpass_chr20_haplotypes")) {
  if (file.exists(destination)) return(invisible(destination))
  paths <- lowpass_stage_paths()
  pedigree <- read.table(paths[["roh_ancestry_pedigree"]], header = TRUE,
                         colClasses = "character", stringsAsFactors = FALSE)
  founders <- pedigree[pedigree$FatherID == "0" & pedigree$MotherID == "0" &
                        pedigree$Population %in% lowpass_declaration$populations, ]
  set.seed(lowpass_declaration$seed)
  selected <- lapply(lowpass_declaration$populations, function(pop) {
    pool <- sort(founders$SampleID[founders$Population == pop])
    stopifnot(length(pool) >= lowpass_declaration$founders_per_population)
    sample(pool, lowpass_declaration$founders_per_population)
  })
  con <- lowpass_connect(extension)
  on.exit(dbDisconnect(con, shutdown = TRUE))
  q <- function(x) as.character(dbQuoteString(con, x))
  column <- function(sample) paste0('"FORMAT_GT_', sample, '"')
  allele_sum <- function(sample) sprintf(
    "((substr(%1$s,1,1)='1')::INTEGER+(substr(%1$s,3,1)='1')::INTEGER)", column(sample))
  frequencies <- vapply(names(selected), function(name) {
    members <- founders$SampleID[founders$Population == lowpass_declaration$populations[[name]]]
    sprintf("greatest(.001,least(.999,(%s)/%d.0)) AS %s",
            paste(vapply(members, allele_sum, ""), collapse = "+"), 2L*length(members),
            if (name == "sample") "f" else "g")
  }, "")
  codes <- vapply(unlist(selected), function(sample) sprintf(
    "(2*(substr(%1$s,1,1)='1')::INTEGER+(substr(%1$s,3,1)='1')::INTEGER)::TINYINT AS %2$s",
    column(sample), column(sample)), "")
  invalid <- paste(sprintf("(NOT regexp_full_match(%s,'[01][|][01]'))::INTEGER",
                           vapply(founders$SampleID, column, "")), collapse = "+")
  panel <- dbGetQuery(con, paste0("SELECT v.POS AS pos, ", paste(frequencies, collapse = ","),
    ", (", invalid, ") AS invalid, ", paste(codes, collapse = ","),
    " FROM read_bcf(", q(paths[["roh_ancestry_chr20_source"]]),
    ", samples := ", q(paste(founders$SampleID, collapse = ",")), ") v ",
    "JOIN read_parquet(", q(paths[["roh_ancestry_chr20_af_sites"]]),
    ") p ON v.POS=p.pos AND v.REF=p.ref AND len(v.ALT)=1 AND v.ALT[1]=p.alt ORDER BY v.POS"))
  sites <- panel[c("pos", "f", "g")]
  stopifnot(all(panel$invalid == 0), !anyNA(panel), !anyDuplicated(sites$pos))
  genotypes <- setNames(lapply(unlist(selected), function(sample) panel[[paste0("FORMAT_GT_", sample)]]),
                        unlist(selected))
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(sites = sites, genotypes = genotypes, selected = selected,
               seed = lowpass_declaration$seed, source_ids = names(paths)), destination, compress = "gzip")
  invisible(destination)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  lowpass_stage_haplotypes(args[[1L]])
}
