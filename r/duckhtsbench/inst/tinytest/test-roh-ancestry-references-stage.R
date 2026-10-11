library(tinytest)

registry <- duckhts_bench_registry()
transform <- "duckhts_bench_stage_roh_ancestry_references"
# Every autosome but chr20 is staged by the function, from the registered
# reference and its children BCF; chr20 keeps the rows of prepare_chr20.R.
for (chromosome in setdiff(1:22, 20L)) {
  long <- registry[registry$id == sprintf("roh_ancestry_chr%d_reference_long", chromosome), ]
  sites <- registry[registry$id == sprintf("roh_ancestry_chr%d_af_sites", chromosome), ]
  expect_equal(nrow(long), 1L)
  expect_equal(nrow(sites), 1L)
  expect_equal(c(long$transform, sites$transform), c(transform, transform))
  expect_equal(long$locator, "artifact:ancestry_reference_grch38_parquet")
  expect_equal(sites$locator, sprintf(
    "artifact:roh_ancestry_chr%d_children_bcf;artifact:ancestry_reference_grch38_parquet",
    chromosome))
  expect_equal(long$cache_relpath,
               sprintf("benchmarks/roh-ancestry/chr%d.reference_long.parquet", chromosome))
  expect_equal(sites$cache_relpath,
               sprintf("benchmarks/roh-ancestry/chr%d.af_reference_sites.parquet", chromosome))
}
expect_false(any(registry$transform[registry$id %in% c(
  "roh_ancestry_chr20_reference_long", "roh_ancestry_chr20_af_sites")] == transform))
# The genome-wide drivers validate roh-ancestry files against these
# identities, so each must parse as key=value fields.
identity_fields <- getFromNamespace("duckhts_bench_identity_fields", "duckhtsbench")
for (identity in registry$supplier_identity[registry$workload == "roh-ancestry"]) {
  expect_silent(identity_fields(identity))
}

if (requireNamespace("DBI", quietly = TRUE) && requireNamespace("duckdb", quietly = TRUE) &&
    requireNamespace("Rduckhts", quietly = TRUE)) local({
  root <- tempfile("roh-ancestry-references-")
  dir.create(file.path(root, "cache"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE))
  ids <- c("ancestry_reference_grch38_parquet",
           paste0("roh_ancestry_chr", rep(c(20L, 21L, 22L), each = 3L), "_",
                  c("children_bcf", "reference_long", "af_sites")))
  plan <- registry[match(ids, registry$id), , drop = FALSE]
  plan$cache_relpath <- file.path("cache", basename(plan$cache_relpath))
  path_of <- function(id) file.path(root, plan$cache_relpath[plan$id == id])
  sha256 <- function(path) unname(digest::digest(file = path, algo = "sha256"))
  registry_path <- file.path(root, "registry.tsv")
  write_registry <- function() {
    utils::write.table(plan, registry_path, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  previous <- Sys.getenv(c("DUCKHTSBENCH_REGISTRY", "DUCKHTS_CACHE_DIR"), unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) {
      Sys.unsetenv(name)
    } else {
      do.call(Sys.setenv, as.list(previous[name]))
    }
  }, add = TRUE)
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry_path, DUCKHTS_CACHE_DIR = root)

  # A keyed wide reference with the GRCh38 product's schema. Each chr21 site
  # meets one rule: kept forward, kept reversed, palindromic, allele mismatch,
  # a missing INFO field, two ALT alleles in the BCF, and no BCF record. The
  # chr22 site at position 100 repeats chr21's first site.
  sites <- data.frame(
    chromosome = c(rep(21L, 7L), 22L, 22L),
    position = c(100L, 200L, 300L, 400L, 500L, 600L, 700L, 100L, 150L),
    allele_a = c("A", "C", "A", "A", "A", "C", "G", "A", "G"),
    allele_b = c("G", "T", "T", "C", "G", "A", "T", "G", "A"),
    stringsAsFactors = FALSE)
  groups <- sprintf("Group %02d (test)", seq_len(21L))
  frequency <- matrix(seq_len(nrow(sites) * 21L) / 1000, nrow(sites))
  loading <- matrix(seq_len(nrow(sites) * 16L) / 100, nrow(sites))
  wide <- cbind(sites, stats::setNames(as.data.frame(loading), paste0("PC", 1:16)),
                stats::setNames(as.data.frame(frequency), groups))
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  quote_str <- function(x) as.character(DBI::dbQuoteString(con, x))
  DBI::dbWriteTable(con, "wide", wide)
  reference <- path_of("ancestry_reference_grch38_parquet")
  DBI::dbExecute(con, paste0(
    "COPY (SELECT chromosome::UTINYINT AS chromosome, position::UINTEGER AS position, ",
    "* EXCLUDE (chromosome, position) FROM wide ORDER BY chromosome, position) TO ",
    quote_str(reference), " (FORMAT PARQUET)"))
  write_reference_receipt <- function(output_sha256) {
    utils::write.table(data.frame(field = "output_sha256", value = output_sha256),
                       paste0(reference, ".sources.tsv"), sep = "\t", row.names = FALSE,
                       quote = FALSE)
  }
  write_reference_receipt(sha256(reference))

  # The children "BCFs" are text VCFs, which read_bcf detects by content.
  info <- function(af) sprintf("AF=%s;AF_AFR=0.5;AF_AMR=0.6;AF_EAS=0.7;AF_EUR=0.8", af)
  record <- function(chrom, pos, ref, alt, info_field) {
    paste(chrom, pos, ".", ref, alt, ".", "PASS", info_field, "GT", "0/1", sep = "\t")
  }
  header <- c("##fileformat=VCFv4.2", "##contig=<ID=chr21,length=1000>",
              "##contig=<ID=chr22,length=1000>",
              sprintf('##INFO=<ID=%s,Number=A,Type=Float,Description="%s">',
                      c("AF", "AF_AFR", "AF_AMR", "AF_EAS", "AF_EUR"), "frequency"),
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1")
  chr21_records <- function(first_af) c(
    record("chr21", 100, "A", "G", info(first_af)),
    record("chr21", 200, "T", "C", info("0.2")),
    record("chr21", 300, "A", "T", info("0.3")),
    record("chr21", 400, "A", "G", info("0.4")),
    record("chr21", 500, "A", "G", "AF=0.5;AF_AFR=0.5;AF_AMR=0.6;AF_EUR=0.8"),
    record("chr21", 600, "C", "A,T", paste0("AF=0.1,0.2;AF_AFR=0.1,0.2;AF_AMR=0.1,0.2;",
                                            "AF_EAS=0.1,0.2;AF_EUR=0.1,0.2")),
    # A record of another contig at a chr21 reference position is not used.
    record("chr22", 100, "A", "G", info("0.9")))
  chr22_records <- c(record("chr22", 100, "A", "G", info("0.11")),
                     record("chr22", 150, "G", "A", info("0.15")))
  write_bcf <- function(chromosome, records) {
    id <- sprintf("roh_ancestry_chr%d_children_bcf", chromosome)
    writeLines(c(header, records), path_of(id))
    plan$supplier_identity[plan$id == id] <<- paste0("records=", length(records),
                                                     ";sha256=", sha256(path_of(id)))
    write_registry()
  }
  write_bcf(21L, chr21_records("0.1"))
  write_bcf(22L, chr22_records)
  read_output <- function(path) {
    DBI::dbGetQuery(con, paste0("SELECT * FROM read_parquet(", quote_str(path), ")"))
  }
  receipt_of <- function(path) {
    receipt <- utils::read.delim(paste0(path, ".sources.tsv"), colClasses = "character")
    stats::setNames(receipt$value, receipt$field)
  }

  staged <- duckhts_bench_stage_roh_ancestry_references(c(21L, 22L))
  expect_equal(staged$id, c("roh_ancestry_chr21_reference_long", "roh_ancestry_chr21_af_sites",
                            "roh_ancestry_chr22_reference_long", "roh_ancestry_chr22_af_sites"))
  expect_equal(staged$path, vapply(staged$id, path_of, character(1L), USE.NAMES = FALSE))
  expect_equal(staged$rows, c(7 * 21, 2, 2 * 21, 2))
  expect_false(any(staged$reused))

  long21 <- read_output(path_of("roh_ancestry_chr21_reference_long"))
  expect_equal(names(long21), c("chromosome", "position", "allele_a", "allele_b", "group_id",
                                "frequency"))
  expect_equal(unique(long21$chromosome), "21")
  expect_equal(long21$position, rep(c(100, 200, 300, 400, 500, 600, 700), each = 21L))
  expect_equal(long21$group_id, rep(sort(groups), 7L))
  # Each row keeps its site's group frequency unchanged.
  expected_frequency <- frequency[cbind(rep(1:7, each = 21L),
                                        rep(match(sort(groups), groups), 7L))]
  expect_equal(long21$frequency, expected_frequency)

  sites21 <- read_output(path_of("roh_ancestry_chr21_af_sites"))
  expect_equal(names(sites21), c("chrom", "pos", "ref", "alt", "INFO_AF", "INFO_AF_AFR",
                                 "INFO_AF_AMR", "INFO_AF_EAS", "INFO_AF_EUR"))
  expect_equal(sites21$chrom, c("chr21", "chr21"))
  expect_equal(sites21$pos, c(100, 200))
  expect_equal(sites21$ref, c("A", "T"))
  expect_equal(sites21$alt, c("G", "C"))
  expect_equal(unlist(sites21$INFO_AF), c(0.1, 0.2), tolerance = 1e-6)
  expect_equal(unlist(sites21$INFO_AF_EUR), c(0.8, 0.8), tolerance = 1e-6)
  sites22 <- read_output(path_of("roh_ancestry_chr22_af_sites"))
  expect_equal(paste(sites22$chrom, sites22$pos, sites22$ref, sites22$alt),
               c("chr22 100 A G", "chr22 150 G A"))

  extension <- system.file("duckhts_extension", "build", "duckhts.duckdb_extension",
                           package = "Rduckhts")
  receipt <- receipt_of(path_of("roh_ancestry_chr21_af_sites"))
  expect_equal(receipt[["artifact_id"]], "roh_ancestry_chr21_af_sites")
  expect_equal(receipt[["chromosome"]], "21")
  expect_equal(receipt[["reference_sha256"]], sha256(reference))
  expect_equal(receipt[["children_bcf_sha256"]],
               sha256(path_of("roh_ancestry_chr21_children_bcf")))
  expect_equal(receipt[["extension_sha256"]], sha256(extension))
  expect_equal(receipt[["output_sha256"]], sha256(path_of("roh_ancestry_chr21_af_sites")))
  expect_true(grepl("^[0-9a-f]{64}$", receipt[["derivation_sha256"]]))
  expect_true(nzchar(receipt[["duckdb_version"]]))
  expect_equal(receipt[["rows"]], "2")
  expect_equal(receipt[["reference_sites"]], "6")
  receipt <- receipt_of(path_of("roh_ancestry_chr21_reference_long"))
  expect_equal(receipt[["reference_sha256"]], sha256(reference))
  expect_equal(receipt[["output_sha256"]], sha256(path_of("roh_ancestry_chr21_reference_long")))
  expect_equal(c(receipt[["rows"]], receipt[["sites"]]), c("147", "7"))
  expect_false("children_bcf_sha256" %in% names(receipt))

  # Matching receipts are reused.
  expect_true(all(duckhts_bench_stage_roh_ancestry_references(c(21L, 22L))$reused))

  # A changed output no longer matches its receipt and is rebuilt.
  writeLines("corrupt", path_of("roh_ancestry_chr21_af_sites"))
  rebuilt <- duckhts_bench_stage_roh_ancestry_references(21L)
  expect_equal(rebuilt$reused, c(TRUE, FALSE))
  expect_equal(read_output(path_of("roh_ancestry_chr21_af_sites")), sites21)

  # A BCF that differs from its registry identity is an error; a newly
  # registered BCF rebuilds only the AF sites.
  writeLines(c(header, chr21_records("0.15")), path_of("roh_ancestry_chr21_children_bcf"))
  expect_error(duckhts_bench_stage_roh_ancestry_references(21L), "SHA-256")
  write_bcf(21L, chr21_records("0.15"))
  changed <- duckhts_bench_stage_roh_ancestry_references(21L)
  expect_equal(changed$reused, c(TRUE, FALSE))
  expect_equal(unlist(read_output(path_of("roh_ancestry_chr21_af_sites"))$INFO_AF),
               c(0.15, 0.2), tolerance = 1e-6)

  # A failed derivation publishes nothing: two records of one site leave the
  # previous output and receipt in place and no partial file behind.
  before <- sha256(path_of("roh_ancestry_chr22_af_sites"))
  write_bcf(22L, c(chr22_records, chr22_records[[2L]]))
  expect_error(duckhts_bench_stage_roh_ancestry_references(22L), "one row per site")
  expect_equal(sha256(path_of("roh_ancestry_chr22_af_sites")), before)
  expect_equal(receipt_of(path_of("roh_ancestry_chr22_af_sites"))[["output_sha256"]], before)
  expect_false(any(grepl("partial|spill", list.files(file.path(root, "cache")))))

  # chr20 is not staged here; requests are checked before any work.
  expect_error(duckhts_bench_stage_roh_ancestry_references(20L), "not staged by")
  expect_error(duckhts_bench_stage_roh_ancestry_references(c(21L, 21L)), "distinct")
  expect_error(duckhts_bench_stage_roh_ancestry_references(23L), "distinct")
  expect_error(duckhts_bench_stage_roh_ancestry_references(1.5), "distinct")

  # The reference must match its own receipt.
  write_reference_receipt(strrep("0", 64L))
  expect_error(duckhts_bench_stage_roh_ancestry_references(21L), "receipt")
})
