# Public epilepsy summary-statistics comparison with the bigsnpr 1.12.21 oracle.
# Run from the repository root after installing Rduckhts, duckhtsbench and bigsnpr.
library(data.table)
library(DBI)
library(Rduckhts)
library(duckhtsbench)

Sys.setenv(DUCKHTSBENCH_REGISTRY = normalizePath("r/duckhtsbench/inst/benchmark_registry.tsv"))
paths <- setNames(vapply(duckhts_bench_stage_plan("ancestry-reference")$id,
                         duckhts_bench_artifact_path, character(1L)),
                  duckhts_bench_stage_plan("ancestry-reference")$id)
for (id in names(paths)) duckhts_bench_validate_identity(id)
correction <- scan("test/data/ancestry_correction_bigsnpr_1.12.21.tsv",
                   what = numeric(), sep = "\t", skip = 1L, quiet = TRUE)
# The committed correction has two columns (PC ordinal and coefficient).
correction <- correction[seq(2L, length(correction), 2L)]

message("Reading public reference frequencies and loadings")
start <- proc.time()[["elapsed"]]
ref <- bigreadr::fread2(paths[["ancestry_ref_freqs"]])
pc <- bigreadr::fread2(paths[["ancestry_projection"]])
setDT(ref)
setDT(pc)
stopifnot(nrow(ref) == 5816590L, nrow(pc) == nrow(ref),
          identical(ref[["chr"]], pc[["chr"]]),
          identical(ref[["pos"]], pc[["pos"]]),
          identical(ref[["a0"]], pc[["a0"]]),
          identical(ref[["a1"]], pc[["a1"]]),
          identical(ref[["rsid"]], pc[["rsid"]]))
load_seconds <- proc.time()[["elapsed"]] - start
message("Loaded ", nrow(ref), " sites in ", round(load_seconds, 2), " s")
sumstats <- bigreadr::fread2(paths[["ancestry_epilepsy"]],
  select = c("CHR", "BP", "Allele2", "Allele1", "Freq1"),
  col.names = c("chr", "pos", "a0", "a1", "freq"))
setDT(sumstats)
sumstats[, `:=`(chr = as.integer(chr), a0 = toupper(a0),
                a1 = toupper(a1), beta = 1)]
message("Summary input rows: ", nrow(sumstats))
matched <- bigsnpr::snp_match(sumstats, ref[, 1:5, with = FALSE],
                               return_flip_and_rev = TRUE)
setDT(matched)
matched[, freq := fifelse(`_REV_`, 1 - freq, freq)]
message("bigsnpr matched: ", nrow(matched), "; reversed: ",
        sum(matched$`_REV_`), "; flipped: ", sum(matched$`_FLIP_`))
stopifnot(nrow(matched) > 1000L)
# Compare on unique, unambiguous SNVs with finite observed frequency.
start <- proc.time()[["elapsed"]]
upstream <- bigsnpr::snp_ancestry_summary(
  matched$freq, ref[matched$`_NUM_ID_`, -(1:5), with = FALSE],
  as.matrix(pc[matched$`_NUM_ID_`, -(1:5), with = FALSE]), correction)
message("Published bigsnpr solution: ", proc.time()[["elapsed"]] - start, " s")
print(data.frame(group_id = names(upstream), proportion = as.numeric(upstream),
                 cor_pred = attr(upstream, "cor_pred")))

# Use the same physical variants and orientation for a solver parity comparison.
# A separate audit of the unfiltered input is required for matching-policy parity.
site <- ref[matched$`_NUM_ID_`, 1:4, with = FALSE]
# SNV admission, physical-locus uniqueness, and observed-frequency validity
# are separate reasons for exclusion from the shared solver denominator.
snv <- nchar(site$a0) == 1L & nchar(site$a1) == 1L &
  !paste0(site$a0, site$a1) %in% c("AT", "TA", "CG", "GC")
locus <- paste(site$chr, site$pos)
unique_locus <- !duplicated(locus) & !duplicated(locus, fromLast = TRUE)
keep <- snv & unique_locus & is.finite(matched$freq)
chr22_only <- Sys.getenv("ANCESTRY_EPILEPSY_CHR22_ONLY", "true")
stopifnot(chr22_only %in% c("true", "false"))
if (chr22_only == "true") keep <- keep & site$chr == 22L
eligible <- sum(keep)
limit <- if (chr22_only == "true") min(17000L, eligible) else eligible
keep <- which(keep)[seq_len(limit)]
message("Shared unambiguous SNV sites: ", length(keep),
        "/", eligible, " eligible; ", nrow(matched), " genome-wide bigsnpr matches")
index <- matched$`_NUM_ID_`[keep]
shared <- bigsnpr::snp_ancestry_summary(
  matched$freq[keep], ref[index, -(1:5), with = FALSE],
  as.matrix(pc[index, -(1:5), with = FALSE]), correction)
con <- rduckhts_connect()
dbExecute(con, "SET threads=4")
input <- data.frame(sample_id = "epilepsy", chromosome = paste0("chr", site$chr[keep]),
                    position = site$pos[keep], allele_a = site$a0[keep],
                    allele_b = site$a1[keep], frequency = matched$freq[keep])
reference_parquet <- duckhts_bench_stage_ancestry_parquet()
dbExecute(con, paste0("CREATE VIEW real_reference AS SELECT * FROM read_parquet(",
                      as.character(DBI::dbQuoteString(con, reference_parquet)), ")"))
dbWriteTable(con, "real_input", input, overwrite = TRUE)
input_parquet <- Sys.getenv("DUCKHTS_ANCESTRY_INPUT_PARQUET")
if (nzchar(input_parquet)) {
  dbExecute(con, paste0("COPY real_input TO ",
                        as.character(DBI::dbQuoteString(con, input_parquet)),
                        " (FORMAT PARQUET)"))
}
dbWriteTable(con, "real_correction", data.frame(pc = seq_along(correction),
                                               coefficient = correction), overwrite = TRUE)
start <- proc.time()[["elapsed"]]
result <- rduckhts_ancestry_proportions(con, "real_input", "real_reference",
                                             "real_correction")
elapsed <- proc.time()[["elapsed"]] - start
parity <- merge(result, data.frame(group_id = names(shared), oracle = as.numeric(shared)),
                by = "group_id", sort = FALSE)
parity$absolute_difference <- abs(parity$proportion - parity$oracle)
print(parity[, c("group_id", "proportion", "oracle", "absolute_difference",
                 "cor_pred", "status", "used_variants", "input_variants")])
reference_matrix <- as.matrix(ref[index, -(1:5), with = FALSE])
predicted <- as.vector(reference_matrix %*% as.numeric(shared))
predicted <- predicted - mean(predicted)
rounding_norm_bound <- sqrt(length(predicted)) * ncol(reference_matrix) * 5e-8
correlation_bound <- 2 * rounding_norm_bound /
  (sqrt(sum(predicted^2)) - rounding_norm_bound)
message("shared bigsnpr cor_pred: ", attr(shared, "cor_pred"),
        "; DuckHTS cor_pred: ", result$cor_pred[[1L]],
        "; maximum seven-decimal QP-rounding correlation change: ",
        correlation_bound, "; gate 0.4 distance: ",
        abs(result$cor_pred[[1L]] - 0.4),
        "; DuckHTS time: ", elapsed, " s")
stopifnot(all(parity$status == "ok"),
          all(round(parity$proportion, 7L) == round(parity$oracle, 7L)),
          abs(result$cor_pred[[1L]] - 0.4) > correlation_bound)
dbDisconnect(con, shutdown = TRUE)
