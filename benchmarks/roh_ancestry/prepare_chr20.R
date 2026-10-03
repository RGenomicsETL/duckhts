# The installed Rduckhts; load_all() would compile r/Rduckhts in place.
library(Rduckhts)
cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-ancestry")
con <- rduckhts_connect(extension_path = normalizePath("build/release/duckhts.duckdb_extension"),
                        config = list(memory_limit = "12GB", threads = "4",
                                      temp_directory = file.path(cache, "duckdb-tmp")))
q <- function(x) as.character(DBI::dbQuoteString(con, x))
reference_path <- duckhtsbench::duckhts_bench_artifact_path("ancestry_reference_grch38_parquet")
bcf_path <- file.path(cache, "chr20.children.bcf")
# The children decoded by every arm: the BCF header samples, in header order.
bcftools <- Sys.getenv("BCFTOOLS", "bcftools")
children <- system2(bcftools, c("query", "-l", shQuote(bcf_path)), stdout = TRUE)
if (!is.null(attr(children, "status")) || length(children) != 377L) {
  stop("chr20 children BCF must list the 377 eligible children", call. = FALSE)
}
writeLines(children, file.path(cache, "chr20.children.present.txt"))
long_output <- file.path(cache, "chr20.reference_long.parquet")
af_output <- file.path(cache, "chr20.af_reference_sites.parquet")
DBI::dbExecute(con, sprintf(
  "CREATE TEMP VIEW ref20 AS SELECT * FROM read_parquet(%s) WHERE try_cast(regexp_replace(CAST(chromosome AS VARCHAR), '^chr', '') AS INTEGER) = 20",
  q(reference_path)))
columns <- DBI::dbGetQuery(con, "DESCRIBE ref20")$column_name
pc_columns <- paste0("PC", seq_len(16L))
group_columns <- setdiff(columns, c("chromosome", "position", "allele_a", "allele_b", pc_columns))
if (length(group_columns) != 21L) stop("expected 21 ancestry groups", call. = FALSE)
unpivot_columns <- paste(as.character(DBI::dbQuoteIdentifier(con, group_columns)), collapse = ", ")
long_query <- sprintf(paste0(
  "UNPIVOT ref20 ON %s INTO NAME group_id VALUE frequency"), unpivot_columns)
DBI::dbExecute(con, sprintf(
  "COPY (SELECT chromosome::VARCHAR AS chromosome, position::BIGINT AS position, allele_a, allele_b, group_id, frequency::DOUBLE AS frequency FROM (%s)) TO %s (FORMAT PARQUET, COMPRESSION ZSTD)",
  long_query, q(long_output)))
DBI::dbExecute(con, sprintf(paste0(
  "COPY (SELECT v.CHROM AS chrom, v.POS AS pos, v.REF AS ref, v.ALT[1] AS alt, ",
  "v.INFO_AF, v.INFO_AF_AFR, v.INFO_AF_AMR, v.INFO_AF_EAS, v.INFO_AF_EUR ",
  "FROM read_bcf(%s) v JOIN (SELECT DISTINCT position, allele_a, allele_b FROM ref20) r ",
  "ON v.POS = r.position AND ((v.REF = r.allele_a AND v.ALT[1] = r.allele_b) OR ",
  "(v.REF = r.allele_b AND v.ALT[1] = r.allele_a)) ",
  "WHERE len(v.ALT) = 1 AND NOT ((r.allele_a IN ('A','T') AND r.allele_b IN ('A','T')) OR ",
  "(r.allele_a IN ('C','G') AND r.allele_b IN ('C','G'))) ",
  "AND v.INFO_AF IS NOT NULL AND v.INFO_AF_AFR IS NOT NULL AND v.INFO_AF_AMR IS NOT NULL ",
  "AND v.INFO_AF_EAS IS NOT NULL AND v.INFO_AF_EUR IS NOT NULL) TO %s ",
  "(FORMAT PARQUET, COMPRESSION ZSTD)"), q(bcf_path), q(af_output)))
cat("reference rows:", DBI::dbGetQuery(con, sprintf("SELECT count(*) n FROM read_parquet(%s)", q(long_output)))$n, "\n")
cat("eligible reference sites:", DBI::dbGetQuery(con, sprintf("SELECT count(*) n FROM read_parquet(%s)", q(af_output)))$n, "\n")
cat("INFO fields and site rows:\n")
print(DBI::dbGetQuery(con, sprintf("SELECT count(*) n, count(DISTINCT (chrom,pos,ref,alt)) sites FROM read_parquet(%s)", q(af_output))))
DBI::dbDisconnect(con, shutdown = TRUE)
