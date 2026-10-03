cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-ancestry")
driver <- duckdb::duckdb(dbdir = ":memory:", shared_home = FALSE,
  allow_extensions = TRUE, config = list(memory_limit = "12GB", threads = "4",
    allow_unsigned_extensions = "true", autoload_known_extensions = "false",
    autoinstall_known_extensions = "false"))
connection <- DBI::dbConnect(driver)
DBI::dbExecute(connection, sprintf("LOAD %s", as.character(DBI::dbQuoteString(
  connection, normalizePath("build/release/duckhts.duckdb_extension")))))
q <- function(value) as.character(DBI::dbQuoteString(connection, value))
sql <- sprintf(paste0(
  "WITH bcf_sites AS (",
  "SELECT DISTINCT regexp_replace(CHROM, '^chr', '')::VARCHAR AS chrom, POS::BIGINT AS pos, REF::VARCHAR AS ref, ",
  "ALT[1]::VARCHAR AS alt, INFO_AF, INFO_AF_AFR, INFO_AF_AMR, INFO_AF_EAS, INFO_AF_EUR ",
  "FROM read_bcf(%s) WHERE len(ALT) = 1), ",
  "ref_pairs AS (SELECT DISTINCT try_cast(regexp_replace(chromosome, '^chr', '') AS VARCHAR) AS chrom, ",
  "position::BIGINT AS pos, least(allele_a, allele_b) AS lo, greatest(allele_a, allele_b) AS hi ",
  "FROM read_parquet(%s) WHERE NOT ((allele_a='A' AND allele_b='T') OR ",
  "(allele_a='T' AND allele_b='A') OR (allele_a='C' AND allele_b='G') OR ",
  "(allele_a='G' AND allele_b='C'))), ",
  "ancestry_keys AS (SELECT DISTINCT c.chrom, c.pos, least(c.ref,c.alt) AS lo, ",
  "greatest(c.ref,c.alt) AS hi FROM bcf_sites c JOIN ref_pairs r ",
  "ON c.chrom=r.chrom AND c.pos=r.pos AND least(c.ref,c.alt)=r.lo AND greatest(c.ref,c.alt)=r.hi), ",
  "af_keys AS (SELECT DISTINCT c.chrom, c.pos, least(c.ref,c.alt) AS lo, ",
  "greatest(c.ref,c.alt) AS hi FROM bcf_sites c JOIN ref_pairs r ",
  "ON c.chrom=r.chrom AND c.pos=r.pos AND least(c.ref,c.alt)=r.lo AND greatest(c.ref,c.alt)=r.hi ",
  "WHERE c.INFO_AF IS NOT NULL AND c.INFO_AF_AFR IS NOT NULL AND c.INFO_AF_AMR IS NOT NULL ",
  "AND c.INFO_AF_EAS IS NOT NULL AND c.INFO_AF_EUR IS NOT NULL), ",
  "af_only AS (SELECT * FROM af_keys EXCEPT SELECT * FROM ancestry_keys), ",
  "ancestry_only AS (SELECT * FROM ancestry_keys EXCEPT SELECT * FROM af_keys) ",
  "SELECT (SELECT count(*) FROM af_keys) af_site_count, ",
  "(SELECT count(*) FROM ancestry_keys) ancestry_site_count, ",
  "(SELECT count(*) FROM af_only) af_only_count, ",
  "(SELECT count(*) FROM ancestry_only) ancestry_only_count"),
  q(file.path(cache, "chr20.children.bcf")),
  q(file.path(cache, "chr20.reference_long.parquet")))
result <- DBI::dbGetQuery(connection, sql)
print(result)
utils::write.csv(result,
  "benchmarks/results/roh-ancestry/chr20/site_set_equivalence.csv", row.names = FALSE)
stopifnot(result$af_only_count == 0L, result$ancestry_only_count == 0L,
          result$af_site_count == result$ancestry_site_count)
DBI::dbDisconnect(connection, shutdown = TRUE)
