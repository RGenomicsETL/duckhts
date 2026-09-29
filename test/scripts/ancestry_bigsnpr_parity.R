#!/usr/bin/env Rscript
# Fixed synthetic matching audit and pinned bigsnpr 1.12.21 oracle.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: ancestry_bigsnpr_parity.R private_library")
.libPaths(c(args[[1L]], .libPaths()))
stopifnot(as.character(packageVersion("bigsnpr")) == "1.12.21")
library(DBI)
library(Rduckhts)
con <- rduckhts_connect()
on.exit(dbDisconnect(con, shutdown = TRUE))
dbExecute(con, paste(
  "CREATE TABLE ancestry_parity_ref AS SELECT 'chr1' AS chromosome, (100+i)::INTEGER AS position,",
  "CASE WHEN i % 2 = 0 THEN 'A' ELSE 'C' END AS allele_a,",
  "CASE WHEN i % 2 = 0 THEN 'G' ELSE 'T' END AS allele_b, g.group_id,",
  "CASE WHEN g.group_id = 'A' THEN 0.15 + i*0.06 ELSE 0.8 - i*0.04 END AS frequency",
  "FROM range(8) t(i), (VALUES ('A'), ('B')) g(group_id)"
))
dbExecute(con, paste(
  "CREATE TABLE ancestry_parity_pc AS SELECT DISTINCT chromosome, position, allele_a, allele_b,",
  "pc, CASE WHEN pc=1 THEN position - 103.0 ELSE sin(position) END AS loading",
  "FROM ancestry_parity_ref, (VALUES (1), (2)) p(pc)"
))
dbExecute(con, "CREATE TABLE ancestry_parity_correction AS SELECT 1 AS pc, 1.0 AS coefficient UNION ALL SELECT 2, 1.0")
dbExecute(con, paste(
  "CREATE TABLE ancestry_parity_input AS SELECT 'S' AS sample_id, chromosome, position, allele_a, allele_b,",
  "sum(CASE WHEN group_id='A' THEN 0.7 ELSE 0.3 END * frequency) AS frequency",
  "FROM ancestry_parity_ref GROUP BY chromosome, position, allele_a, allele_b"
))
dbExecute(con, "UPDATE ancestry_parity_input SET allele_a='T', allele_b='C' WHERE position=100")
dbExecute(con, "UPDATE ancestry_parity_input SET allele_a='T', allele_b='C', frequency=1-frequency WHERE position=101")
dbExecute(con, "UPDATE ancestry_parity_input SET frequency=frequency+0.04 WHERE position=102")
dbExecute(con, "UPDATE ancestry_parity_input SET frequency=NULL WHERE position=107")
dbExecute(con, "INSERT INTO ancestry_parity_input SELECT * FROM ancestry_parity_input WHERE position=106")
dbExecute(con, "INSERT INTO ancestry_parity_input VALUES ('S','chr1',999,'A','G',0.5)")
observed <- rduckhts_ancestry_proportions(
  con, "ancestry_parity_input", "ancestry_parity_ref", "ancestry_parity_pc",
  "ancestry_parity_correction", min_cor = 0
)
input <- dbGetQuery(con, "SELECT chromosome,position,allele_a,allele_b,frequency FROM ancestry_parity_input")
ref <- dbGetQuery(con, "SELECT chromosome,position,allele_a,allele_b,group_id,frequency FROM ancestry_parity_ref ORDER BY position,group_id")
reference_sites <- ref[!duplicated(ref$position), ]
sumstats <- data.frame(chr = input$chromosome, pos = input$position,
                       a0 = input$allele_a, a1 = input$allele_b,
                       beta = input$frequency)
sumstats <- sumstats[!is.na(sumstats$beta), ]
matched <- bigsnpr::snp_match(
  sumstats, data.frame(chr = reference_sites$chromosome,
                       pos = reference_sites$position,
                       a0 = reference_sites$allele_a,
                       a1 = reference_sites$allele_b),
  return_flip_and_rev = TRUE
)
stopifnot(nrow(matched) == 6L, sum(matched$`_FLIP_`) == 1L,
          sum(matched$`_REV_`) == 1L,
          unique(observed$used_variants) == nrow(matched),
          unique(observed$flipped_variants) == sum(matched$`_FLIP_`),
          unique(observed$reversed_variants) == sum(matched$`_REV_`))
frequencies <- sumstats$beta[match(matched$pos, sumstats$pos)]
frequencies[matched$`_REV_`] <- 1 - frequencies[matched$`_REV_`]
reference_matrix <- matrix(ref$frequency[ref$position %in% matched$pos],
                           ncol = 2L, byrow = TRUE)
projection <- cbind(matched$pos - 103, sin(matched$pos))
oracle <- bigsnpr::snp_ancestry_summary(frequencies, reference_matrix,
                                         projection, correction = c(1, 1), min_cor = 0)
error <- max(abs(unname(oracle) - observed$proportion))
cat(sprintf("source=bigsnpr_1.12.21 matched=%d input=%d dropped=%d flip=%d reverse=%d max_proportion_error=%.9g cor_pred_error=%.9g\n",
            nrow(matched), unique(observed$input_variants),
            unique(observed$dropped_variants), sum(matched$`_FLIP_`),
            sum(matched$`_REV_`), error,
            abs(attr(oracle, "cor_pred") - observed$cor_pred[[1L]])))
for (i in seq_len(nrow(observed))) {
  cat(sprintf("group=%s duckhts=%.7f bigsnpr=%.7f difference=%.9g cor_each_error=%.9g\n",
              observed$group_id[[i]], observed$proportion[[i]], oracle[[i]],
              abs(observed$proportion[[i]] - oracle[[i]]),
              abs(observed$cor_each[[i]] - attr(oracle, "cor_each")[[i]])))
}
stopifnot(error < 2e-4,
          max(abs(observed$cor_each - attr(oracle, "cor_each"))) < 1e-8,
          abs(attr(oracle, "cor_pred") - observed$cor_pred[[1L]]) < 1e-6)
