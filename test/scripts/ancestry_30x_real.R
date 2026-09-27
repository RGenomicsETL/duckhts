# Indexed 30x CRAM chr22 versus phase-3 GRCh37 VCF at lifted, allele-checked sites.
# Run from the repository root; never fetch a whole public CRAM.
source("test/scripts/ancestry_1000g_real.R", local = TRUE)

bundle <- c(source_fasta = duckhts_bench_artifact_path("liftover_grch37_fasta"),
  destination_fasta = duckhts_bench_fetch("liftover_grch38_fasta"),
  chain = duckhts_bench_fetch("liftover_grch37_grch38_chain"))
stopifnot(all(file.exists(bundle)))
con <- rduckhts_connect()
dbExecute(con, "SET threads=4")
source_sites <- data.frame(chrom = "22", pos = as.integer(site$pos),
                           ref = site$a, alt = site$b)
dbWriteTable(con, "source_sites", source_sites, overwrite = TRUE)
lifted <- as.data.table(rduckhts_liftover(
  con, "source_sites", bundle[["chain"]], bundle[["destination_fasta"]],
  chrom_col = "chrom", pos_col = "pos", ref_col = "ref", alt_col = "alt",
  src_fasta_ref = bundle[["source_fasta"]]))
lifted <- merge(lifted, site, by.x = "src_pos", by.y = "pos", sort = FALSE)
message("Lifted ", sum(lifted$mapped), "/", nrow(site), " loci; ",
        sum(lifted$swap != 0L, na.rm = TRUE), " source/destination reference swaps")
mapped <- lifted[mapped & dest_chrom == "chr22"]
aligned <- mapped[swap == 0L & src_ref == a & src_alt == b]
snvs <- aligned[nchar(dest_ref) == 1L & nchar(dest_alt) == 1L]
snvs[, copies := .N, by = dest_pos]
usable <- snvs[copies == 1L][order(dest_pos)]
message("Mapped on chr22: ", nrow(mapped), "; unswapped source alleles: ",
        nrow(aligned), "; SNVs: ", nrow(snvs),
        "; unique destination loci: ", nrow(usable), "/", nrow(site))
stopifnot(nrow(usable) >= 17000L)
usable <- usable[unique(as.integer(round(seq(1, nrow(usable), length.out = 17000L))))]
usable[, allele_a := fifelse(a0 == a, dest_ref, dest_alt)]
usable[, allele_b := fifelse(a1 == b, dest_alt, dest_ref)]
stopifnot(all(usable$allele_a != usable$allele_b))
cache <- dirname(duckhts_bench_artifact_path("ancestry_ref_freqs"))
region <- sprintf("chr22:%d-%d", max(1L, min(usable$dest_pos) - 1000L),
                  max(usable$dest_pos) + 1000L)
message("Indexed remote CRAM region: ", region)

ids <- c(NA18507 = "ancestry_30x_na18507", HG00188 = "riker_hg00188_cram",
         HG00403 = "ancestry_30x_hg00403")
requested <- Sys.getenv("ANCESTRY_30X_SAMPLES", "")
if (nzchar(requested)) {
  ids <- ids[strsplit(requested, ",", fixed = TRUE)[[1L]]]
  stopifnot(length(ids) > 0L, !anyNA(ids))
}
registry <- duckhts_bench_registry()
cram_sources <- setNames(registry$locator[match(ids, registry$id)], names(ids))
stopifnot(all(startsWith(cram_sources, "https://")))
for (sample in names(cram_sources)) {
  response <- system2("curl", c("-fIs", "--max-time", "30",
    shQuote(cram_sources[[sample]])), stdout = TRUE, stderr = TRUE)
  stopifnot(is.null(attr(response, "status")))
  metadata <- registry$supplier_identity[[match(ids[[sample]], registry$id)]]
  fields <- strsplit(metadata, ";", fixed = TRUE)[[1L]]
  values <- strsplit(fields, "=", fixed = TRUE)
  names(values) <- vapply(values, `[[`, character(1L), 1L)
  content_length <- tail(response[grepl("^content-length:", response,
    ignore.case = TRUE)], 1L)
  stopifnot(length(content_length) == 1L,
    identical(trimws(sub("^[^:]+:", "", content_length)), values$bytes[[2L]]))
  etag <- tail(response[grepl("^etag:", response, ignore.case = TRUE)], 1L)
  if ("etag" %in% names(values)) {
    stopifnot(length(etag) == 1L,
      identical(gsub('"', "", trimws(sub("^[^:]+:", "", etag)), fixed = TRUE),
                values$etag[[2L]]))
  }
}

# Local indexed regional CRAMs keep repeated panel/depth calls bounded and
# retain a single physical source for native and downsampled measurements.
for (sample in names(cram_sources)) {
  path <- file.path(cache, paste0(sample, ".chr22-region.cram"))
  start <- proc.time()[["elapsed"]]
  status <- system2("samtools", c("view", "-C", "-T", shQuote(bundle[["destination_fasta"]]),
    "-o", shQuote(path), shQuote(cram_sources[[sample]]), shQuote(region)))
  stopifnot(status == 0L, file.info(path)$size > 0L)
  unlink(paste0(path, ".crai"))
  stopifnot(system2("samtools", c("index", shQuote(path))) == 0L)
  message(sample, ": staged chr22 region in ",
          round(proc.time()[["elapsed"]] - start, 2), " s; ",
          file.info(path)$size, " bytes")
  downsampled <- file.path(cache, paste0(sample, ".chr22-quarter.cram"))
  status <- system2("samtools", c("view", "-C", "-s", "42.25", "-T",
    shQuote(bundle[["destination_fasta"]]), "-o", shQuote(downsampled), shQuote(path)))
  stopifnot(status == 0L)
  unlink(paste0(downsampled, ".crai"))
  stopifnot(system2("samtools", c("index", shQuote(downsampled))) == 0L)
}

dbWriteTable(con, "ancestry_correction", data.frame(pc = seq_along(correction),
  coefficient = correction), overwrite = TRUE)
results <- list()
for (n in c(1000L, 5000L, 17000L)) {
  selected <- usable[unique(as.integer(round(seq(1, nrow(usable), length.out = n))))]
  panel <- data.frame(assembly = "GRCh38", site_index = seq_len(n) - 1L,
    region = selected$dest_chrom, position = selected$dest_pos,
    allele_a = pmin(selected$dest_ref, selected$dest_alt),
    allele_b = pmax(selected$dest_ref, selected$dest_alt))
  key <- data.frame(chromosome = selected$dest_chrom, position = selected$dest_pos,
    allele_a = selected$allele_a, allele_b = selected$allele_b)
  index <- selected$ref_index
  frequencies <- cbind(key[rep(seq_len(n), times = ncol(ref) - 5L), ],
    group_id = rep(names(ref)[-(1:5)], each = n),
    frequency = unlist(ref[index, -(1:5), with = FALSE], use.names = FALSE))
  loadings <- cbind(key[rep(seq_len(n), times = ncol(pc) - 5L), ],
    pc = rep(seq_len(ncol(pc) - 5L), each = n),
    loading = unlist(pc[index, -(1:5), with = FALSE], use.names = FALSE))
  dbWriteTable(con, "ancestry_panel", panel, overwrite = TRUE)
  panel_sha256 <- dbGetQuery(con,
    "SELECT duckhts_somalier_panel_sha256('ancestry_panel') AS panel_sha256")$panel_sha256[[1L]]
  dbWriteTable(con, "ancestry_reference", frequencies, overwrite = TRUE)
  dbWriteTable(con, "ancestry_loadings", loadings, overwrite = TRUE)
  for (sample in names(cram_sources)) {
    gt <- selected[[sample]]
    dosage <- fifelse(gt %in% c("0|0", "0/0"), 0,
                      fifelse(gt %in% c("0|1", "1|0", "0/1", "1/0"), 1,
                              fifelse(gt %in% c("1|1", "1/1"), 2, NA_real_)))
    genotype <- data.frame(sample_id = sample, chromosome = selected$dest_chrom,
      position = selected$dest_pos, allele_a = selected$dest_ref,
      allele_b = selected$dest_alt, frequency = dosage / 2)
    dbWriteTable(con, "ancestry_genotype", genotype, overwrite = TRUE)
    vcf <- rduckhts_ancestry_proportions(con, "ancestry_genotype",
      "ancestry_reference", "ancestry_loadings", "ancestry_correction")
    for (depth in c("native", "quarter")) {
      path <- file.path(cache, paste0(sample, ".chr22-region.cram"))
      if (depth == "quarter") {
        path <- file.path(cache, paste0(sample, ".chr22-quarter.cram"))
      }
      rduckhts_somalier_bam_counts(con, path, sample,
        bundle[["destination_fasta"]], panel_table = "ancestry_panel",
        table_name = "ancestry_depth", worker_count = 2L)
      observed_depth <- dbGetQuery(con, paste(
        "SELECT median(a+b) AS median_depth, sum(a+b) AS allele_reads,",
        "count(*) FILTER (WHERE a+b >= 7) AS depth_eligible FROM ancestry_depth"
      ))
      dbExecute(con, "DROP TABLE ancestry_depth")
      for (method in c("allele_fraction", "called_genotype")) {
        start <- proc.time()[["elapsed"]]
        result <- rduckhts_ancestry_bam(con, path, sample,
          bundle[["destination_fasta"]], "ancestry_panel", "ancestry_reference",
          "ancestry_loadings", "ancestry_correction", frequency_method = method,
          min_depth = 7L, worker_count = 2L)
        elapsed <- proc.time()[["elapsed"]] - start
        status <- unique(result$status)
        if (length(status) != 1L) stop("mixed status across groups")
        delta <- if (identical(status, "ok") && all(vcf$status == "ok")) {
          max(abs(result$proportion[match(vcf$group_id, result$group_id)] - vcf$proportion))
        } else NA_real_
        rss <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
        peak_rss_mib <- as.numeric(strsplit(trimws(rss), "[[:space:]]+")[[1L]][[2L]]) / 1024
        row <- data.frame(sample_id = sample, superpopulation = superpopulation[[sample]],
          sites = n, panel_sha256 = panel_sha256, depth = depth, mode = method,
          vcf_used = vcf$used_variants[[1L]],
          bam_used = result$used_variants[[1L]],
          depth_eligible = observed_depth$depth_eligible[[1L]],
          median_depth = observed_depth$median_depth[[1L]],
          allele_reads = observed_depth$allele_reads[[1L]],
          cor_pred = result$cor_pred[[1L]],
          vcf_cor_pred = vcf$cor_pred[[1L]], status = status,
          max_group_difference = delta, seconds = elapsed,
          process_peak_rss_mib = peak_rss_mib)
        results[[length(results) + 1L]] <- row
        print(row)
      }
    }
  }
}
sensitivity <- do.call(rbind, results)
# For every simplex coefficient vector, the smallest singular value of the
# centered frequency matrix bounds its prediction norm from below. The
# seven-decimal QP rounding perturbation is at most 21 * 5e-8 per site.
rounding_bounds <- setNames(vapply(c(1000L, 5000L, 17000L), function(n) {
  selected <- usable[unique(as.integer(round(seq(1, nrow(usable), length.out = n))))]
  frequencies <- as.matrix(ref[selected$ref_index, -(1:5), with = FALSE])
  centered <- sweep(frequencies, 2L, colMeans(frequencies), "-")
  minimum_norm <- min(svd(centered, nu = 0L, nv = 0L)$d) / sqrt(ncol(centered))
  perturbation <- sqrt(n) * ncol(centered) * 5e-8
  2 * perturbation / (minimum_norm - perturbation)
}, numeric(1L)), c("1000", "5000", "17000"))
sensitivity$rounding_bound <- unname(rounding_bounds[as.character(sensitivity$sites)])
sensitivity$gate_distance <- abs(sensitivity$cor_pred - 0.4)
measurable <- is.finite(sensitivity$gate_distance)
message("Gate failures: ", sum(sensitivity$status != "ok"), "/", nrow(sensitivity),
        "; maximum rounding bound: ", max(sensitivity$rounding_bound),
        "; minimum measured gate distance: ",
        min(sensitivity$gate_distance[measurable], na.rm = TRUE))
stopifnot(all(sensitivity$gate_distance[measurable] >
              sensitivity$rounding_bound[measurable]))
if (length(ids) == 3L) {
  utils::write.table(sensitivity, "benchmarks/ancestry_30x_sensitivity.tsv",
    sep = "\t", row.names = FALSE, quote = FALSE)
}
dbDisconnect(con, shutdown = TRUE)
