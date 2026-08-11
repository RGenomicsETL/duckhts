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
source(file.path(package_dir, "R", "duckvep.R"))

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
temporary_path <- function(output) paste0(output, ".partial-", Sys.getpid())
publish <- function(temporary, output) {
  if (file.exists(output) && unlink(output, force = TRUE) != 0L) {
    stop("could not replace staged artifact: ", output, call. = FALSE)
  }
  if (!file.rename(temporary, output)) stop("could not publish staged artifact: ", output, call. = FALSE)
}

# Direct source downloads are the only generic operation. All other transforms
# below name their input/output registry IDs explicitly.
for (id in plan$id[plan$transform == "direct_download"]) {
  duckhts_bench_fetch(id)
}

fasta_gz <- artifact("ensembl116_grch38_fasta")
fasta <- artifact("ensembl116_grch38_fasta_fa")
dir.create(dirname(fasta), recursive = TRUE, showWarnings = FALSE)
gzip <- Sys.which("gzip")
if (!nzchar(gzip)) stop("gzip is required to derive the Ensembl FASTA", call. = FALSE)
temporary <- temporary_path(fasta)
unlink(temporary, force = TRUE)
status <- system2(gzip, c("--decompress", "--keep", "--stdout", shQuote(fasta_gz)), stdout = temporary)
if (status != 0L) stop("could not decompress registered Ensembl FASTA", call. = FALSE)
publish(temporary, fasta)
samtools <- Sys.which("samtools")
if (!nzchar(samtools)) stop("samtools is required to index the Ensembl FASTA", call. = FALSE)
index <- paste0(fasta, ".fai")
if (file.exists(index) && unlink(index, force = TRUE) != 0L) stop("could not replace Ensembl FASTA index", call. = FALSE)
status <- system2(samtools, c("faidx", shQuote(fasta)))
if (status != 0L || !file.exists(index) || !file.info(index)$size) stop("could not index registered Ensembl FASTA", call. = FALSE)
duckhts_bench_write_provenance("ensembl116_grch38_fasta_fa", fasta)

revel_zip <- artifact("revel_v13_source_zip")
revel_parquet <- artifact("revel_v13_grch37")
dir.create(dirname(revel_parquet), recursive = TRUE, showWarnings = FALSE)
extracted_dir <- tempfile("duckhtsbench-revel-")
utils::unzip(revel_zip, exdir = extracted_dir)
candidates <- list.files(extracted_dir, pattern = "^revel_with_transcript_ids$", recursive = TRUE, full.names = TRUE)
if (length(candidates) != 1L) stop("REVEL v1.3 archive lacks a unique revel_with_transcript_ids member", call. = FALSE)
temporary <- temporary_path(revel_parquet)
unlink(temporary, force = TRUE)
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
DBI::dbExecute(con, sprintf(
  paste(
    "COPY (SELECT chr AS chrom, CAST(hg19_pos AS BIGINT) AS pos, ref, alt,",
    "CAST(REVEL AS DOUBLE) AS revel",
    "FROM read_csv_auto('%s', header = TRUE, types = {'chr': 'VARCHAR'})",
    "WHERE TRY_CAST(hg19_pos AS BIGINT) IS NOT NULL)",
    "TO '%s' (FORMAT PARQUET, COMPRESSION ZSTD)"
  ),
  gsub("'", "''", candidates[[1L]], fixed = TRUE),
  gsub("'", "''", temporary, fixed = TRUE)
))
DBI::dbDisconnect(con, shutdown = TRUE)
unlink(extracted_dir, recursive = TRUE, force = TRUE)
publish(temporary, revel_parquet)
duckhts_bench_write_provenance("revel_v13_grch37", revel_parquet)

clinvarbitration_archive <- artifact("clinvarbitration_202508_source_archive")
clinvarbitration_tsv <- artifact("clinvarbitration_202508")
dir.create(dirname(clinvarbitration_tsv), recursive = TRUE, showWarnings = FALSE)
members <- utils::untar(clinvarbitration_archive, list = TRUE)
member <- members[grepl("(^|/)clinvar_decisions\\.tsv$", members)]
if (length(member) != 1L) stop("ClinvArbitration archive lacks a unique clinvar_decisions.tsv member", call. = FALSE)
work <- tempfile("duckhtsbench-clinvarbitration-")
utils::untar(clinvarbitration_archive, files = member, exdir = work)
source <- file.path(work, member)
if (!file.exists(source)) stop("could not extract ClinvArbitration decision table", call. = FALSE)
temporary <- temporary_path(clinvarbitration_tsv)
unlink(temporary, force = TRUE)
if (!file.copy(source, temporary, overwrite = TRUE)) stop("could not stage ClinvArbitration decision table", call. = FALSE)
unlink(work, recursive = TRUE, force = TRUE)
publish(temporary, clinvarbitration_tsv)
duckhts_bench_write_provenance("clinvarbitration_202508", clinvarbitration_tsv)

extension <- Sys.getenv(
  "DUCKHTS_EXTENSION",
  unset = file.path(package_dir, "..", "..", "build", "release", "duckhts.duckdb_extension")
)
model <- duckhts_bench_stage_duckvep_ensembl116_model(extension)

export_relation <- function(id, query) {
  output <- artifact(id)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  temporary <- temporary_path(output)
  unlink(temporary, force = TRUE)
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  DBI::dbExecute(con, sprintf("ATTACH '%s' AS model (READ_ONLY)", gsub("'", "''", model, fixed = TRUE)))
  DBI::dbExecute(con, sprintf(
    "COPY (%s) TO '%s' (FORMAT PARQUET, COMPRESSION ZSTD)",
    query, gsub("'", "''", temporary, fixed = TRUE)
  ))
  DBI::dbDisconnect(con, shutdown = TRUE)
  publish(temporary, output)
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
