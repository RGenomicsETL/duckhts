# GRCh37 phase-3 chr22 genotype parity against bigsnpr 1.12.21.
# Run from the repository root with bcftools, Rduckhts, bigsnpr and duckhtsbench.
library(data.table)
library(DBI)
library(Rduckhts)
library(duckhtsbench)

Sys.setenv(DUCKHTSBENCH_REGISTRY = normalizePath("r/duckhtsbench/inst/benchmark_registry.tsv"))
ref_path <- duckhts_bench_fetch("ancestry_ref_freqs")
pc_path <- duckhts_bench_fetch("ancestry_projection")
vcf_path <- duckhts_bench_fetch("ancestry_1000g_chr22")
load_start <- proc.time()[["elapsed"]]
ref <- as.data.table(bigreadr::fread2(ref_path))
pc <- as.data.table(bigreadr::fread2(pc_path))
load_seconds <- proc.time()[["elapsed"]] - load_start
message("Loaded ", nrow(ref), " reference rows and ", nrow(pc),
        " loading rows in ", round(load_seconds, 2), " s")
stopifnot(nrow(ref) == nrow(pc), identical(ref$pos, pc$pos),
          identical(ref$a0, pc$a0), identical(ref$a1, pc$a1))
correction <- data.table::fread("test/data/ancestry_correction_bigsnpr_1.12.21.tsv")$coefficient
samples <- c("HG01879", "HG02561", "NA18507", "HG00188", "NA12878",
             "HG00403", "HG02057", "HG03742", "HG03866", "HG00731", "HG01112")
sample_panel <- fread(duckhts_bench_fetch("ancestry_1000g_panel"), header = FALSE,
  col.names = c("sample_id", "population", "superpopulation", "sex"))
superpopulation <- setNames(sample_panel$superpopulation[
  match(samples, sample_panel$sample_id)], samples)
stopifnot(identical(unname(superpopulation),
  c(rep("AFR", 3L), rep("EUR", 2L), rep("EAS", 2L), rep("SAS", 2L), rep("AMR", 2L))))
selected_path <- file.path(dirname(ref_path), "chr22.ten-samples.vcf.gz")
sample_file <- tempfile("ancestry-samples-")
writeLines(samples, sample_file)
status <- system2("bcftools", c("view", "--threads", "2", "-S",
  shQuote(sample_file), "-v", "snps", "-m2", "-M2", "-Oz", "-o",
  shQuote(selected_path), shQuote(vcf_path)))
unlink(sample_file)
stopifnot(status == 0L)
stopifnot(system2("bcftools", c("index", "-f", "-t", shQuote(selected_path))) == 0L)
vcf_samples <- system2("bcftools", c("query", "-l", shQuote(selected_path)), stdout = TRUE)
stopifnot(setequal(samples, vcf_samples))
format <- "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n"
query <- sprintf("bcftools query -f %s %s", shQuote(format), shQuote(selected_path))
vcf <- fread(cmd = query, header = FALSE, sep = "\t",
             col.names = c("chromosome", "pos", "a", "b", vcf_samples))
stopifnot(nrow(vcf) > 100000L)
vcf[, copies := .N, by = pos]
vcf <- vcf[copies == 1L & nchar(a) == 1L & nchar(b) == 1L]
# One unambiguous chr22 SNV admission predicate owns both reference rows and
# their original product indices; the latter index the paired PC loadings.
reference_index <- which(ref$chr == 22L & nchar(ref$a0) == 1L &
  nchar(ref$a1) == 1L & !paste0(ref$a0, ref$a1) %in% c("AT", "TA", "CG", "GC"))
reference <- copy(ref[reference_index])
reference[, ref_index := reference_index]
reference[, copies := .N, by = pos]
reference <- reference[copies == 1L]
joined <- merge(reference[, .(pos, a0, a1, ref_index)],
                vcf, by = "pos", sort = TRUE)
# Accept either allele orientation at one physical position; exclude
# strand-ambiguous SNVs before this match.
joined <- joined[(a0 == a & a1 == b) | (a0 == b & a1 == a)]
joined[, reversed := a0 == b]
message("chr22 VCF biallelic SNVs: ", nrow(vcf), "; reference overlap: ", nrow(joined))
stopifnot(nrow(joined) >= 20000L)
site <- joined[unique(as.integer(round(seq(1, nrow(joined), length.out = 20000L))))]
index <- site$ref_index
keys <- data.frame(chromosome = "22", position = site$pos,
                   allele_a = site$a0, allele_b = site$a1)
ref_long <- cbind(keys[rep(seq_len(nrow(keys)), times = ncol(ref) - 5L), ],
  group_id = rep(names(ref)[-(1:5)], each = nrow(keys)),
  frequency = unlist(ref[index, -(1:5), with = FALSE], use.names = FALSE))
pc_long <- cbind(keys[rep(seq_len(nrow(keys)), times = ncol(pc) - 5L), ],
  pc = rep(seq_len(ncol(pc) - 5L), each = nrow(keys)),
  loading = unlist(pc[index, -(1:5), with = FALSE], use.names = FALSE))
con <- rduckhts_connect()
dbExecute(con, "SET threads=4")
dbWriteTable(con, "real_reference", ref_long, overwrite = TRUE)
dbWriteTable(con, "real_loadings", pc_long, overwrite = TRUE)
dbWriteTable(con, "real_correction", data.frame(pc = seq_along(correction),
  coefficient = correction), overwrite = TRUE)
input <- rbindlist(lapply(samples, function(sample) {
  gt <- site[[sample]]
  dosage <- fifelse(gt %in% c("0|0", "0/0"), 0,
                    fifelse(gt %in% c("0|1", "1|0", "0/1", "1/0"), 1,
                            fifelse(gt %in% c("1|1", "1/1"), 2, NA_real_)))
  data.table(sample_id = sample, chromosome = "22", position = site$pos,
             allele_a = site$a, allele_b = site$b, frequency = dosage / 2)
}))
dbWriteTable(con, "real_input", input, overwrite = TRUE)
start <- proc.time()[["elapsed"]]
result <- rduckhts_ancestry_proportions(con, "real_input", "real_reference",
                                        "real_loadings", "real_correction", min_cor = 0.4)
elapsed <- proc.time()[["elapsed"]] - start
comparisons <- rbindlist(lapply(samples, function(sample) {
  rows <- input[sample_id == sample]
  valid <- which(!is.na(rows$frequency))
  observed <- ifelse(site$reversed[valid], 1 - rows$frequency[valid],
                     rows$frequency[valid])
  oracle <- tryCatch(bigsnpr::snp_ancestry_summary(
    observed, ref[index[valid], -(1:5), with = FALSE],
    as.matrix(pc[index[valid], -(1:5), with = FALSE]), correction),
    error = function(e) e)
  if (inherits(oracle, "error")) {
    message("bigsnpr gate: ", sample, ": ", conditionMessage(oracle))
    return(data.table(sample_id = sample, group_id = names(ref)[-(1:5)],
                      oracle = NA_real_, oracle_cor_pred = NA_real_,
                      oracle_error = conditionMessage(oracle), oracle_sites = length(valid)))
  }
  data.table(sample_id = sample, group_id = names(oracle),
             oracle = as.numeric(oracle), oracle_cor_pred = attr(oracle, "cor_pred"),
             oracle_error = "", oracle_sites = length(valid))
}))
parity <- merge(as.data.table(result), comparisons, by = c("sample_id", "group_id"))
parity[, `:=`(superpopulation = superpopulation[sample_id],
             absolute_difference = abs(proportion - oracle),
             correlation_difference = abs(cor_pred - oracle_cor_pred))]
summary <- parity[, .(superpopulation = first(superpopulation),
  input_rows = first(input_variants), duckhts_sites = first(used_variants),
  bigsnpr_sites = first(oracle_sites), max_group_difference = max(absolute_difference),
  cor_pred = first(cor_pred), oracle_cor_pred = first(oracle_cor_pred),
  status = first(status), oracle_error = first(oracle_error)), by = sample_id]
summary[, rounding_bound := vapply(sample_id, function(sample) {
  rows <- input[sample_id == sample]
  valid <- which(!is.na(rows$frequency))
  frequencies <- as.matrix(ref[index[valid], -(1:5), with = FALSE])
  coefficients <- parity[sample_id == sample][match(names(ref)[-(1:5)], group_id)]$oracle
  prediction <- as.vector(frequencies %*% coefficients)
  centered_norm <- sqrt(sum((prediction - mean(prediction))^2))
  perturbation <- sqrt(length(valid)) * ncol(frequencies) * 5e-8
  2 * perturbation / (centered_norm - perturbation)
}, numeric(1L))]
summary[, gate_distance := abs(cor_pred - 0.4)]
print(summary)
message("DuckHTS chr22 query time: ", elapsed, " s; ", nrow(parity),
        " output group rows; ", nrow(input), " genotype input rows; ",
        "maximum rounding bound: ", max(summary$rounding_bound),
        "; minimum gate distance: ", min(summary$gate_distance))
stopifnot(nrow(parity) == length(samples) * 21L,
          all(summary$status == "ok"), all(summary$oracle_error == ""),
          all(summary$max_group_difference < 2e-4),
          all(summary$gate_distance > summary$rounding_bound))
dbDisconnect(con, shutdown = TRUE)
