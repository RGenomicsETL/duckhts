# The installed Rduckhts; load_all() would compile r/Rduckhts in place.
library(Rduckhts)
cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-ancestry")
con <- rduckhts_connect(extension_path = normalizePath("build/release/duckhts.duckdb_extension"),
                        config = list(memory_limit = "12GB", threads = "4",
                                      temp_directory = file.path(cache, "duckdb-tmp")))
q <- function(x) as.character(DBI::dbQuoteString(con, x))
bcf <- file.path(cache, "chr20.children.bcf")
af <- file.path(cache, "chr20.af_reference_sites.parquet")
ref <- file.path(cache, "chr20.reference_long.parquet")
query <- sprintf(paste0(
  "WITH a AS (SELECT SAMPLE_ID AS sample_id, count(*) AS sites FROM read_bcf(%s, tidy_format := true) b ",
  "JOIN read_parquet(%s) f ON b.CHROM = f.chrom AND b.POS = f.pos AND b.REF = f.ref AND b.ALT[1] = f.alt ",
  "GROUP BY SAMPLE_ID), r AS (SELECT DISTINCT position, allele_a, allele_b FROM read_parquet(%s)), ",
  "z AS (SELECT SAMPLE_ID AS sample_id, count(*) AS sites FROM read_bcf(%s, tidy_format := true) b ",
  "JOIN r ON b.POS = r.position AND ((b.REF = r.allele_a AND b.ALT[1] = r.allele_b) OR ",
  "(b.REF = r.allele_b AND b.ALT[1] = r.allele_a)) ",
  "AND NOT ((r.allele_a IN ('A','T') AND r.allele_b IN ('A','T')) OR ",
  "(r.allele_a IN ('C','G') AND r.allele_b IN ('C','G'))) GROUP BY SAMPLE_ID) ",
  "SELECT a.sample_id, a.sites AS pooled, a.sites AS single_population, z.sites AS ancestry_tuned, ",
  "a.sites = z.sites AS parity FROM a JOIN z USING(sample_id) ORDER BY a.sample_id"),
  q(bcf), q(af), q(ref), q(bcf))
site_sets <- DBI::dbGetQuery(con, query)
utils::write.csv(site_sets, file.path(cache, "chr20.site_sets.csv"), row.names = FALSE)
print(table(site_sets$pooled, site_sets$ancestry_tuned))
print(table(site_sets$parity))
cat("Rows per arm (unique):\n")
print(sapply(site_sets[c("pooled", "single_population", "ancestry_tuned")], unique))
DBI::dbDisconnect(con, shutdown = TRUE)
