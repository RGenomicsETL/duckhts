# Check the staged synthetic children against their source through read_bcf,
# apart from the staging code, which edits bcftools text.
#
# Each BCF is scanned once: 849,143 records by 377 samples, 320,126,911 genotype
# rows. The 1,131 truth rows are the build side of a join on the sample, and the
# rows aggregate to one group per sample and planted segment plus one group per
# sample for the genotypes outside every segment: 1,508 groups.
#
# Inside a segment the synthetic child has no heterozygous genotype, and is
# homozygous for the alternate allele exactly where haplotype 1 of the source
# carries it. Outside, the genotypes are the source's: the counts and an
# order-free digest of (POS, GT) agree.
results <- "benchmarks/results/roh-synthetic"
dir.create(results, recursive = TRUE, showWarnings = FALSE)
cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-synthetic")
dir.create(file.path(cache, "duckdb-tmp"), recursive = TRUE, showWarnings = FALSE)
driver <- duckdb::duckdb(dbdir = ":memory:", shared_home = FALSE,
  allow_extensions = TRUE, config = list(memory_limit = "12GB", threads = "4",
    temp_directory = file.path(cache, "duckdb-tmp"),
    allow_unsigned_extensions = "true", autoinstall_known_extensions = "false",
    autoload_known_extensions = "false"))
con <- DBI::dbConnect(driver)
q <- function(x) as.character(DBI::dbQuoteString(con, x))
DBI::dbExecute(con, sprintf("LOAD %s", q(normalizePath("build/release/duckhts.duckdb_extension"))))
DBI::dbExecute(con, sprintf(
  "CREATE TEMP TABLE truth AS SELECT * FROM read_csv(%s)",
  q(duckhtsbench::duckhts_bench_artifact_path("roh_synthetic_chr20_truth"))))

# Segment 0 is the group of genotypes outside every planted segment.
genotype_groups <- function(artifact) {
  DBI::dbGetQuery(con, sprintf(paste0(
    "SELECT g.SAMPLE_ID AS sample_id, coalesce(t.segment, 0)::INTEGER AS segment, ",
    "count(*)::DOUBLE AS records, ",
    "count(*) FILTER (WHERE g.FORMAT_GT IN ('0|1', '1|0'))::DOUBLE AS heterozygous, ",
    "count(*) FILTER (WHERE g.FORMAT_GT IN ('1|0', '1|1'))::DOUBLE AS first_allele_alt, ",
    "count(*) FILTER (WHERE g.FORMAT_GT = '1|1')::DOUBLE AS homozygous_alt, ",
    "sum(hash(g.POS, g.FORMAT_GT)::HUGEINT)::VARCHAR AS digest ",
    "FROM read_bcf(%s, tidy_format := true) AS g LEFT JOIN truth AS t ",
    "ON t.sample_id = g.SAMPLE_ID AND g.POS BETWEEN t.start AND t.\"end\" ",
    "GROUP BY ALL ORDER BY sample_id, segment"),
    q(duckhtsbench::duckhts_bench_artifact_path(artifact))))
}
source_groups <- genotype_groups("roh_ancestry_chr20_children_bcf")
synthetic_groups <- genotype_groups("roh_synthetic_chr20_children_bcf")
truth <- DBI::dbGetQuery(con, paste(
  "SELECT sample_id, segment::INTEGER AS segment, records::DOUBLE AS records,",
  "source_heterozygous::DOUBLE AS source_heterozygous FROM truth ORDER BY sample_id, segment"))
DBI::dbDisconnect(con, shutdown = TRUE)

same_groups <- function(left, right) {
  nrow(left) == nrow(right) && all(left$sample_id == right$sample_id) &&
    all(left$segment == right$segment)
}
inside <- source_groups$segment > 0
if (!same_groups(source_groups, synthetic_groups) ||
    !same_groups(source_groups[inside, , drop = FALSE], truth)) {
  stop("source, synthetic and truth groups differ", call. = FALSE)
}
checks <- data.frame(
  check = c("planted segments hold the truth's record counts",
            "planted segments replaced the truth's heterozygous counts",
            "synthetic planted segments hold no heterozygous genotype",
            "synthetic planted segments are homozygous for haplotype 1's allele",
            "genotypes outside planted segments are unchanged"),
  groups = c(rep(sum(inside), 4L), sum(!inside)),
  genotypes = c(rep(sum(source_groups$records[inside]), 4L),
                sum(source_groups$records[!inside])),
  failures = c(
    sum(source_groups$records[inside] != truth$records |
          synthetic_groups$records[inside] != truth$records),
    sum(source_groups$heterozygous[inside] != truth$source_heterozygous),
    sum(synthetic_groups$heterozygous[inside] != 0),
    sum(synthetic_groups$homozygous_alt[inside] != source_groups$first_allele_alt[inside] |
          synthetic_groups$first_allele_alt[inside] != source_groups$first_allele_alt[inside]),
    sum(source_groups$records[!inside] != synthetic_groups$records[!inside] |
          source_groups$heterozygous[!inside] != synthetic_groups$heterozygous[!inside] |
          source_groups$digest[!inside] != synthetic_groups$digest[!inside])))
utils::write.csv(checks, file.path(results, "planting_verification.csv"), row.names = FALSE)
print(checks)
if (any(checks$failures != 0)) stop("the planted children differ from the plan", call. = FALSE)
