cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-ancestry")
connection <- Rduckhts::rduckhts_connect(
  extension_path = normalizePath("build/release/duckhts.duckdb_extension"),
  config = list(memory_limit = "12GB", threads = "4",
                temp_directory = file.path(cache, "duckdb-tmp")))
quote_string <- function(value) as.character(DBI::dbQuoteString(connection, value))
results <- lapply(1:22, function(chromosome) {
  bcf_path <- file.path(cache, sprintf("chr%d.children.bcf", chromosome))
  reference_path <- file.path(cache, sprintf("chr%d.reference_long.parquet", chromosome))
  sql <- sprintf(paste0(
    "WITH bcf_sites AS (SELECT DISTINCT regexp_replace(CHROM, '^chr', '')::VARCHAR AS chrom, ",
    "POS::BIGINT AS pos, REF::VARCHAR AS ref, ALT[1]::VARCHAR AS alt, ",
    "INFO_AF, INFO_AF_AFR, INFO_AF_AMR, INFO_AF_EAS, INFO_AF_EUR ",
    "FROM read_bcf(%s) WHERE len(ALT)=1), ",
    "ref_pairs AS (SELECT DISTINCT regexp_replace(chromosome, '^chr', '')::VARCHAR AS chrom, ",
    "position::BIGINT AS pos, least(allele_a,allele_b) AS lo, ",
    "greatest(allele_a,allele_b) AS hi FROM read_parquet(%s) ",
    "WHERE NOT ((allele_a='A' AND allele_b='T') OR (allele_a='T' AND allele_b='A') ",
    "OR (allele_a='C' AND allele_b='G') OR (allele_a='G' AND allele_b='C'))), ",
    "ancestry_keys AS (SELECT DISTINCT c.chrom,c.pos,least(c.ref,c.alt) AS lo, ",
    "greatest(c.ref,c.alt) AS hi FROM bcf_sites c JOIN ref_pairs r ",
    "ON c.chrom=r.chrom AND c.pos=r.pos AND least(c.ref,c.alt)=r.lo ",
    "AND greatest(c.ref,c.alt)=r.hi), ",
    "af_keys AS (SELECT DISTINCT c.chrom,c.pos,least(c.ref,c.alt) AS lo, ",
    "greatest(c.ref,c.alt) AS hi FROM bcf_sites c JOIN ref_pairs r ",
    "ON c.chrom=r.chrom AND c.pos=r.pos AND least(c.ref,c.alt)=r.lo ",
    "AND greatest(c.ref,c.alt)=r.hi WHERE c.INFO_AF IS NOT NULL ",
    "AND c.INFO_AF_AFR IS NOT NULL AND c.INFO_AF_AMR IS NOT NULL ",
    "AND c.INFO_AF_EAS IS NOT NULL AND c.INFO_AF_EUR IS NOT NULL), ",
    "af_only AS (SELECT * FROM af_keys EXCEPT SELECT * FROM ancestry_keys), ",
    "ancestry_only AS (SELECT * FROM ancestry_keys EXCEPT SELECT * FROM af_keys) ",
    "SELECT %d AS chromosome, (SELECT count(*) FROM af_keys) AS af_site_count, ",
    "(SELECT count(*) FROM ancestry_keys) AS ancestry_site_count, ",
    "(SELECT count(*) FROM af_only) AS af_only_count, ",
    "(SELECT count(*) FROM ancestry_only) AS ancestry_only_count"),
    quote_string(bcf_path), quote_string(reference_path), chromosome)
  DBI::dbGetQuery(connection, sql)
})
result <- do.call(rbind, results)
if (any(result$af_only_count != 0L) || any(result$ancestry_only_count != 0L) ||
    any(result$af_site_count != result$ancestry_site_count)) {
  stop("AF and ancestry site keys differ on one or more autosomes", call. = FALSE)
}
utils::write.csv(result,
  "benchmarks/results/roh-ancestry/autosome_site_set_equivalence.csv",
  row.names = FALSE)
cat("all autosomes have identical AF and ancestry site keys;",
    sum(result$af_site_count), "sites total\n")
print(result)
DBI::dbDisconnect(connection, shutdown = TRUE)
