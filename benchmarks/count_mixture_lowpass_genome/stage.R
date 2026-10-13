source("benchmarks/count_mixture_lowpass_genome/common.R")

genome_stage_chromosome <- function(chrom, extension, destination = duckhts_bench_artifact_path(genome_id(chrom))) {
  source_ids <- c("roh_ancestry_pedigree", paste0("roh_ancestry_chr", chrom, c("_children_bcf", "_af_sites")))
  paths <- setNames(vapply(source_ids, duckhts_bench_artifact_path, ""), source_ids)
  for (id in source_ids) duckhts_bench_validate_identity(id, paths[[id]])
  identity <- genome_input_identity()
  identity <- identity[identity$id %in% source_ids, ]
  if (file.exists(destination)) {
    existing <- readRDS(destination)
    stopifnot(
      identical(existing$input_identity, identity),
      identical(existing$settings, genome_settings)
    )
    return(invisible(destination))
  }
  bcf <- paths[[paste0("roh_ancestry_chr", chrom, "_children_bcf")]]
  children <- genome_select_children(bcf)
  con <- genome_connect(extension)
  on.exit(dbDisconnect(con, shutdown = TRUE))
  q <- function(x) as.character(dbQuoteString(con, x))
  members <- unlist(children$selected, use.names = FALSE)
  fields <- paste0('v."FORMAT_GT_', members, '"')
  query <- paste0(
    "SELECT p.chrom, p.pos, p.INFO_AF_EUR[1] AS f, p.INFO_AF_AFR[1] AS g, ",
    paste(fields, collapse = ", "), " FROM read_bcf(", q(bcf), ", samples := ",
    q(paste(members, collapse = ",")), ") v JOIN read_parquet(",
    q(paths[[paste0("roh_ancestry_chr", chrom, "_af_sites")]]),
    ") p ON v.POS=p.pos AND v.REF=p.ref AND len(v.ALT)=1 AND v.ALT[1]=p.alt ORDER BY p.pos"
  )
  panel <- dbGetQuery(con, query)
  stopifnot(nrow(panel) > 0L, !anyNA(panel), !anyDuplicated(panel$pos))
  gt <- panel[paste0("FORMAT_GT_", members)]
  phased <- vapply(gt, function(x) all(grepl("^[01]\\|[01]$", x)), logical(1L))
  if (!all(phased)) stop("Non-phased or unavailable GTs: ", paste(members[!phased], collapse = ", "))
  genotypes <- setNames(lapply(gt, function(x) {
    as.raw(2L * (substr(x, 1L, 1L) == "1") +
      (substr(x, 3L, 3L) == "1"))
  }), members)
  sites <- panel[c("chrom", "pos", "f", "g")]
  stopifnot(all(sites$f >= 0 & sites$f <= 1), all(sites$g >= 0 & sites$g <= 1))
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(destination, ".tmp")
  saveRDS(
    list(
      sites = sites, genotypes = genotypes, selected = children$selected,
      pedigree = children$pedigree, input_identity = identity, settings = genome_settings
    ),
    temporary,
    compress = "gzip"
  )
  stopifnot(file.rename(temporary, destination))
  cat(
    "chr", chrom, nrow(sites), "phased sites;", lengths(children$selected), "children;",
    file.info(destination)$size, "bytes\n"
  )
  invisible(destination)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  genome_stage_chromosome(as.integer(args[[1L]]), args[[2L]])
}
