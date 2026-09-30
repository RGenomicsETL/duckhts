#' Stage the GRCh38 ancestry acceptance inputs.
#'
#' Builds the deterministic genome-wide ancestry panel of the GRCh38 product
#' with `rduckhts_ancestry_panel()` and stages, for a fixed set of 1000 Genomes
#' phase-3 individuals, both sides of the acceptance comparison at exactly those
#' sites:
#'
#' * the GRCh37 phase-3 genotypes at each panel site's GRCh37 source locus
#'   (from the product's liftover map), read by indexed region access to the
#'   public per-chromosome phase-3 VCFs with `read_bcf(region := ...)`; no
#'   chromosome is downloaded, and
#' * a small indexed CRAM per individual holding the reads that overlap a panel
#'   site, cut from the registered 30x GRCh38 CRAM with `samtools view -M -L`
#'   against the registered GRCh38 FASTA.
#'
#' The receipt (`<genotypes>.sources.tsv`, columns `field` and `value`) binds the
#' panel SHA-256, the GRCh38 product and liftover map hashes, the identity
#' (bytes and ETag) of every remote VCF and CRAM read, the panel parameters, the
#' row counts, and the SHA-256 of each staged file.
#'
#' @param samples Phase-3 sample identifiers with a registered 30x CRAM.
#' @param max_sites,spacing_bp Panel size and window width.
#' @param vcf_url,cram_url Optional overrides for tests: a function from a
#'   chromosome number to a VCF path or URL, and a named character vector of
#'   CRAM paths or URLs. The default sources are the registered public files.
#' @return The staged genotype Parquet path.
#' @export
duckhts_bench_stage_ancestry_grch38_acceptance <- function(
    samples = c("HG00188", "HG00403", "NA18507"), max_sites = 1000L,
    spacing_bp = 5000L, vcf_url = NULL, cram_url = NULL) {
  for (package in c("DBI", "duckdb", "Rduckhts")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(package, " is required for GRCh38 acceptance staging", call. = FALSE)
    }
  }
  samtools <- Sys.which("samtools")
  if (!nzchar(samtools)) stop("samtools is required for CRAM staging", call. = FALSE)
  registry <- duckhts_bench_registry()
  cram_ids <- c(NA18507 = "ancestry_30x_na18507", HG00188 = "riker_hg00188_cram",
                HG00403 = "ancestry_30x_hg00403")
  registered_cram <- is.null(cram_url)
  if (registered_cram) {
    if (!all(samples %in% names(cram_ids))) stop("unregistered 30x CRAM sample", call. = FALSE)
    cram_url <- stats::setNames(registry$locator[match(cram_ids[samples], registry$id)], samples)
  }
  if (is.null(vcf_url)) {
    template <- registry$locator[registry$id == "ancestry_1000g_chr22"]
    vcf_url <- function(chrom) sub("chr22", paste0("chr", chrom), template, fixed = TRUE)
  }
  product <- duckhts_bench_stage_ancestry_grch38_parquet()
  product_receipt <- utils::read.delim(paste0(product, ".sources.tsv"), colClasses = "character")
  product_value <- function(field) product_receipt$value[product_receipt$field == field]
  map_path <- sub("\\.parquet$", ".liftover.parquet", product)
  destination_fasta <- duckhts_bench_fetch("liftover_grch38_fasta")
  output <- duckhts_bench_artifact_path("ancestry_grch38_acceptance_genotypes")
  directory <- dirname(output)
  receipt <- paste0(output, ".sources.tsv")
  cram_paths <- stats::setNames(file.path(directory, paste0(samples, ".panel-sites.cram")), samples)

  quote_str <- function(con, x) as.character(DBI::dbQuoteString(con, x))
  con <- Rduckhts::rduckhts_connect()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(con, "SET threads=4")
  columns <- DBI::dbGetQuery(con, paste0("DESCRIBE SELECT * FROM read_parquet(",
                                         quote_str(con, product), ")"))$column_name
  pcs <- paste0("PC", seq_len(16L))
  groups <- setdiff(columns, c("chromosome", "position", "allele_a", "allele_b", pcs))
  unpivot <- function(names, name, value) paste0(
    "(UNPIVOT read_parquet(", quote_str(con, product), ") ON ",
    paste0('"', names, '"', collapse = ", "), " INTO NAME ", name, " VALUE ", value, ")")
  DBI::dbExecute(con, paste0(
    "CREATE VIEW long_reference AS SELECT 'chr' || chromosome::VARCHAR AS chromosome, ",
    "position, allele_a, allele_b, group_id, frequency FROM ",
    unpivot(groups, "group_id", "frequency")))
  DBI::dbExecute(con, paste0(
    "CREATE VIEW long_loadings AS SELECT 'chr' || chromosome::VARCHAR AS chromosome, ",
    "position, allele_a, allele_b, replace(pc, 'PC', '')::INTEGER AS pc, loading FROM ",
    unpivot(pcs, "pc", "loading")))
  panel_summary <- Rduckhts::rduckhts_ancestry_panel(
    con, "long_reference", "long_loadings", "acceptance_panel", "GRCh38",
    spacing_bp = spacing_bp, max_sites = max_sites)
  DBI::dbExecute(con, paste0(
    "CREATE TEMP TABLE panel_sites AS SELECT p.site_index, p.region AS chromosome, ",
    "p.position, p.allele_a, p.allele_b, m.source_chromosome::VARCHAR AS source_chromosome, ",
    "m.source_position::BIGINT AS source_position, m.source_allele_a, m.source_allele_b, ",
    "m.swapped, m.reverse_complemented FROM acceptance_panel p JOIN read_parquet(",
    quote_str(con, map_path), ") m ON 'chr' || m.chromosome::VARCHAR = p.region ",
    "AND m.position = p.position ORDER BY p.site_index"))
  sites <- DBI::dbGetQuery(con, "SELECT count(*) AS n FROM panel_sites")$n
  if (sites != panel_summary$sites) stop("panel loci are missing from the liftover map", call. = FALSE)

  expected <- data.frame(
    field = c("panel_sha256", "panel_sites", "max_sites", "spacing_bp", "samples",
              "ancestry_reference_grch38_output_sha256",
              "ancestry_reference_grch38_liftover_map_sha256"),
    value = c(panel_summary$panel_sha256, as.character(sites), as.character(max_sites),
              as.character(spacing_bp), paste(samples, collapse = ","),
              product_value("output_sha256"), product_value("liftover_map_sha256")),
    stringsAsFactors = FALSE)
  sha256 <- function(path) unname(digest::digest(file = path, algo = "sha256"))
  staged <- c(output, cram_paths, paste0(cram_paths, ".crai"))
  if (any(file.exists(c(staged, receipt)))) {
    stored <- tryCatch(utils::read.delim(receipt, colClasses = "character"),
                       error = function(e) NULL)
    ok <- !is.null(stored) && all(file.exists(c(staged, receipt))) &&
      identical(stored$field[seq_len(nrow(expected))], expected$field) &&
      identical(stored$value[seq_len(nrow(expected))], expected$value)
    if (ok) {
      files <- stats::setNames(staged, c("genotypes", paste0("cram_", samples),
                                         paste0("crai_", samples)))
      recorded <- stored$value[match(paste0(names(files), "_sha256"), stored$field)]
      ok <- !anyNA(recorded) && identical(unname(vapply(files, sha256, "")), recorded)
    }
    if (ok) return(output)
    stop("cached GRCh38 acceptance inputs do not match their sources: ", output, call. = FALSE)
  }
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)

  remote_identity <- function(source) {
    if (!grepl("^https?://", source)) {
      return(c(bytes = as.character(file.info(source)$size), etag = ""))
    }
    header <- system2(Sys.which("curl"), c("-fsSIL", "--max-time", "60", shQuote(source)),
                      stdout = TRUE)
    field <- function(name) {
      hit <- utils::tail(header[grepl(paste0("^", name, ":"), header, ignore.case = TRUE)], 1L)
      if (!length(hit)) return("")
      gsub('["\r]', "", trimws(sub("^[^:]+:", "", hit)))
    }
    c(bytes = field("content-length"), etag = field("etag"))
  }
  provenance <- list()
  chromosomes <- DBI::dbGetQuery(con,
    "SELECT DISTINCT source_chromosome FROM panel_sites")$source_chromosome
  chromosomes <- chromosomes[order(as.integer(chromosomes))]
  temporary <- paste0(output, ".partial-", Sys.getpid())
  dir.create(temporary)
  on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
  genotype_parts <- character()
  for (chrom in chromosomes) {
    source <- vcf_url(chrom)
    identity <- remote_identity(source)
    provenance[[paste0("vcf_chr", chrom)]] <- paste0(
      source, ";bytes=", identity[["bytes"]], ";etag=", identity[["etag"]])
    positions <- DBI::dbGetQuery(con, paste0(
      "SELECT DISTINCT source_position FROM panel_sites WHERE source_chromosome = ",
      quote_str(con, chrom), " ORDER BY 1"))$source_position
    regions <- paste0(chrom, ":", positions, "-", positions, collapse = ",")
    part <- file.path(temporary, paste0("genotypes-", chrom, ".parquet"))
    DBI::dbExecute(con, paste0(
      "COPY (SELECT s.site_index, s.chromosome, s.position, s.allele_a, s.allele_b, ",
      "s.source_chromosome, s.source_position, s.source_allele_a, s.source_allele_b, ",
      "s.swapped, v.SAMPLE_ID AS sample_id, v.REF AS ref, v.ALT[1] AS alt, ",
      "v.FORMAT_GT AS gt FROM panel_sites s JOIN (SELECT CHROM, POS, REF, ALT, SAMPLE_ID, ",
      "FORMAT_GT FROM read_bcf(", quote_str(con, source), ", region := ",
      quote_str(con, regions), ", samples := ", quote_str(con, paste(samples, collapse = ",")),
      ", tidy_format := true) WHERE len(ALT) = 1) v ON v.CHROM = s.source_chromosome ",
      "AND v.POS = s.source_position AND ((v.REF = s.source_allele_a AND ",
      "v.ALT[1] = s.source_allele_b) OR (v.REF = s.source_allele_b AND ",
      "v.ALT[1] = s.source_allele_a)) ORDER BY s.site_index, v.SAMPLE_ID) TO ",
      quote_str(con, part), " (FORMAT PARQUET, COMPRESSION ZSTD)"))
    genotype_parts <- c(genotype_parts, part)
  }
  temporary_output <- paste0(output, ".partial-", Sys.getpid(), ".parquet")
  on.exit(unlink(temporary_output), add = TRUE)
  DBI::dbExecute(con, paste0(
    "COPY (SELECT * FROM read_parquet([", paste(quote_str(con, genotype_parts), collapse = ", "),
    "]) ORDER BY site_index, sample_id) TO ", quote_str(con, temporary_output),
    " (FORMAT PARQUET, COMPRESSION ZSTD)"))
  genotype_rows <- DBI::dbGetQuery(con, paste0("SELECT count(*) AS n FROM read_parquet(",
                                               quote_str(con, temporary_output), ")"))$n

  bed <- file.path(temporary, "panel-sites.bed")
  panel_sites <- DBI::dbGetQuery(con,
    "SELECT chromosome, position FROM panel_sites ORDER BY site_index")
  utils::write.table(data.frame(panel_sites$chromosome, panel_sites$position - 1L,
                                panel_sites$position), bed, sep = "\t",
                     row.names = FALSE, col.names = FALSE, quote = FALSE)
  identities <- lapply(cram_url[samples], remote_identity)
  if (registered_cram) {
    for (sample in samples) {
      pinned <- duckhts_bench_identity_fields(
        registry$supplier_identity[registry$id == cram_ids[[sample]]])
      if (!identical(identities[[sample]][["bytes"]], pinned[["bytes"]]) ||
          ("etag" %in% names(pinned) && !identical(identities[[sample]][["etag"]], pinned[["etag"]]))) {
        stop("registered 30x CRAM identity does not match the remote file: ", sample, call. = FALSE)
      }
    }
  }
  partial_cram <- stats::setNames(paste0(cram_paths, ".partial-", Sys.getpid()), samples)
  on.exit(unlink(c(partial_cram, paste0(partial_cram, ".crai"))), add = TRUE)
  cut_cram <- function(sample) {
    status <- system2(samtools, c("view", "-C", "-M", "-L", shQuote(bed), "-T",
                                  shQuote(destination_fasta), "-o", shQuote(partial_cram[[sample]]),
                                  shQuote(cram_url[[sample]])))
    if (status == 0L) status <- system2(samtools, c("index", shQuote(partial_cram[[sample]])))
    status
  }
  status <- parallel::mclapply(samples, cut_cram, mc.cores = length(samples))
  if (any(vapply(status, function(x) !identical(x, 0L), NA))) {
    stop("could not cut the panel-site CRAMs", call. = FALSE)
  }
  for (sample in samples) {
    provenance[[paste0("cram_", sample)]] <- paste0(
      cram_url[[sample]], ";bytes=", identities[[sample]][["bytes"]],
      ";etag=", identities[[sample]][["etag"]])
  }
  hashes <- c(genotypes = sha256(temporary_output),
              stats::setNames(vapply(partial_cram, sha256, ""), paste0("cram_", samples)),
              stats::setNames(vapply(paste0(partial_cram, ".crai"), sha256, ""),
                              paste0("crai_", samples)))
  full <- rbind(
    expected,
    data.frame(field = c(paste0("source_", names(provenance)), "genotype_rows",
                         paste0(names(hashes), "_sha256"),
                         "duckdb_version", "duckhts_htslib_version", "rduckhts_version",
                         "samtools_version"),
               value = c(unlist(provenance), as.character(genotype_rows), hashes,
                         DBI::dbGetQuery(con, "SELECT version() AS v")$v,
                         DBI::dbGetQuery(con, "SELECT duckhts_htslib_version() AS v")$v,
                         as.character(utils::packageVersion("Rduckhts")),
                         system2(samtools, "--version", stdout = TRUE)[[1L]]),
               stringsAsFactors = FALSE))
  temporary_receipt <- paste0(receipt, ".partial-", Sys.getpid())
  on.exit(unlink(temporary_receipt), add = TRUE)
  utils::write.table(full, temporary_receipt, sep = "\t", row.names = FALSE, quote = FALSE)
  moves <- c(temporary_output, partial_cram, paste0(partial_cram, ".crai"), temporary_receipt)
  targets <- c(output, cram_paths, paste0(cram_paths, ".crai"), receipt)
  for (i in seq_along(moves)) {
    if (!file.rename(moves[[i]], targets[[i]])) {
      unlink(targets)
      stop("cannot publish GRCh38 acceptance inputs", call. = FALSE)
    }
  }
  output
}
