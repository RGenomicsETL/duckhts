#!/usr/bin/env Rscript
# Stage the declared VariantKey provider inputs. The registry is the authority
# for every remote locator, cache destination, derivation, and consumer.

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("must be invoked as an R script", call. = FALSE)
package_dir <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1L]])), ".."), mustWork = TRUE)
if (!nzchar(Sys.getenv("DUCKHTSBENCH_REGISTRY", unset = ""))) {
  Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(package_dir, "inst", "benchmark_registry.tsv"))
}
source(file.path(package_dir, "R", "registry.R"))
source(file.path(package_dir, "R", "stage.R"))

if (length(args) && args[[1L]] == "--plan") {
  print(duckhts_bench_stage_plan("variantkey-providers"), row.names = FALSE)
  quit(status = 0L)
}
if (length(args)) stop("usage: stage_variantkey_providers.R [--plan]", call. = FALSE)
if (!requireNamespace("DBI", quietly = TRUE) || !requireNamespace("duckdb", quietly = TRUE)) {
  stop("staging requires the installed R packages DBI and duckdb", call. = FALSE)
}

plan <- duckhts_bench_stage_plan("variantkey-providers")
artifact <- duckhts_bench_artifact_path
# Direct source downloads are the only generic operation. All other transforms
# below name their input/output registry IDs explicitly.
for (id in plan$id[plan$transform == "direct_download"]) {
  duckhts_bench_fetch(id)
}

fasta_gz <- artifact("ensembl116_grch38_fasta")
fasta <- artifact("ensembl116_grch38_fasta_fa")
if (!file.exists(fasta) || !file.info(fasta)$size) {
  dir.create(dirname(fasta), recursive = TRUE, showWarnings = FALSE)
  gzip <- Sys.which("gzip")
  if (!nzchar(gzip)) stop("gzip is required to derive the Ensembl FASTA", call. = FALSE)
  status <- system2(gzip, c("--decompress", "--keep", "--stdout", fasta_gz), stdout = fasta)
  if (status != 0L) stop("could not decompress registered Ensembl FASTA", call. = FALSE)
}
samtools <- Sys.which("samtools")
if (!nzchar(samtools)) stop("samtools is required to index the Ensembl FASTA", call. = FALSE)
if (!file.exists(paste0(fasta, ".fai")) || !file.info(paste0(fasta, ".fai"))$size) {
  status <- system2(samtools, c("faidx", fasta))
  if (status != 0L) stop("could not index registered Ensembl FASTA", call. = FALSE)
}
duckhts_bench_write_provenance("ensembl116_grch38_fasta_fa", fasta)

revel_zip <- artifact("revel_v13_source_zip")
revel_parquet <- artifact("revel_v13_grch37")
if (!file.exists(revel_parquet) || !file.info(revel_parquet)$size) {
  dir.create(dirname(revel_parquet), recursive = TRUE, showWarnings = FALSE)
  extracted_dir <- tempfile("duckhtsbench-revel-")
  on.exit(unlink(extracted_dir, recursive = TRUE, force = TRUE), add = TRUE)
  utils::unzip(revel_zip, exdir = extracted_dir)
  candidates <- list.files(extracted_dir, pattern = "^revel_with_transcript_ids$", recursive = TRUE, full.names = TRUE)
  if (length(candidates) != 1L) stop("REVEL v1.3 archive lacks a unique revel_with_transcript_ids member", call. = FALSE)
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(con, sprintf(
    paste(
      "COPY (SELECT chr AS chrom, CAST(hg19_pos AS BIGINT) AS pos, ref, alt,",
      "CAST(REVEL AS DOUBLE) AS revel",
      "FROM read_csv_auto('%s', header = TRUE) WHERE TRY_CAST(hg19_pos AS BIGINT) IS NOT NULL)",
      "TO '%s' (FORMAT PARQUET, COMPRESSION ZSTD)"
    ),
    gsub("'", "''", candidates[[1L]], fixed = TRUE),
    gsub("'", "''", revel_parquet, fixed = TRUE)
  ))
}
duckhts_bench_write_provenance("revel_v13_grch37", revel_parquet)

clinvarbitration_archive <- artifact("clinvarbitration_202508_source_archive")
clinvarbitration_tsv <- artifact("clinvarbitration_202508")
if (!file.exists(clinvarbitration_tsv) || !file.info(clinvarbitration_tsv)$size) {
  dir.create(dirname(clinvarbitration_tsv), recursive = TRUE, showWarnings = FALSE)
  members <- utils::untar(clinvarbitration_archive, list = TRUE)
  member <- members[grepl("(^|/)clinvar_decisions\\.tsv$", members)]
  if (length(member) != 1L) stop("ClinvArbitration archive lacks a unique clinvar_decisions.tsv member", call. = FALSE)
  work <- tempfile("duckhtsbench-clinvarbitration-")
  on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)
  utils::untar(clinvarbitration_archive, files = member, exdir = work)
  source <- file.path(work, member)
  if (!file.exists(source)) stop("could not extract ClinvArbitration decision table", call. = FALSE)
  if (!file.copy(source, clinvarbitration_tsv, overwrite = TRUE)) stop("could not stage ClinvArbitration decision table", call. = FALSE)
}
duckhts_bench_write_provenance("clinvarbitration_202508", clinvarbitration_tsv)

model <- artifact("duckvep_ensembl116_model")
if (!file.exists(model) || !file.info(model)$size) {
  stop(
    "the registered DuckVEP Ensembl 116 model is absent: ", model,
    ". Its public producer is tracked at https://github.com/RGenomicsETL/duckhts/issues/179",
    call. = FALSE
  )
}
duckhts_bench_write_provenance("duckvep_ensembl116_model", model)

export_relation <- function(id, query) {
  output <- artifact(id)
  if (!file.exists(output) || !file.info(output)$size) {
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
    DBI::dbExecute(con, sprintf("ATTACH '%s' AS model (READ_ONLY)", gsub("'", "''", model, fixed = TRUE)))
    DBI::dbExecute(con, sprintf(
      "COPY (%s) TO '%s' (FORMAT PARQUET, COMPRESSION ZSTD)",
      query, gsub("'", "''", output, fixed = TRUE)
    ))
  }
  duckhts_bench_write_provenance(id, output)
}

export_relation(
  "duckvep_ensembl116_regulatory_parquet",
  paste(
    "SELECT f.stable_id, r.seq_region_name AS chrom, f.feature_start AS start, f.feature_end AS end,",
    "f.feature_class AS feature_type, f.feature_so_term AS so_term, f.feature_so_accession AS so_accession",
    "FROM model.duckvep_regulation_features f JOIN model.model_regions r USING (seq_region)"
  )
)
export_relation(
  "duckvep_ensembl116_transcripts_parquet",
  "SELECT gene_stable_id AS gene_id FROM model.model_transcripts WHERE gene_stable_id IS NOT NULL"
)

message("VariantKey provider inputs staged under ", duckhts_bench_cache_dir())
