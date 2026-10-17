suppressPackageStartupMessages(library(tinytest))
suppressPackageStartupMessages(library(RBCFTools))
setwd(tinytest::get_call_wd())
source("benchmarks/count_mixture_lowpass_genome/stage.R")
test_genome_staging <- function() {
  original_registry <- Sys.getenv("DUCKHTSBENCH_REGISTRY")
  original_cache <- Sys.getenv("DUCKHTS_CACHE_DIR")
  original_settings <- genome_settings
  root <- tempfile("lowpass-genome-stage-")
  dir.create(root)
  on.exit({
    Sys.setenv(DUCKHTSBENCH_REGISTRY = original_registry, DUCKHTS_CACHE_DIR = original_cache)
    unlink(root, recursive = TRUE)
  })
  registry <- duckhts_bench_registry()
  ids <- c("roh_ancestry_pedigree", "roh_ancestry_chr20_children_bcf", "roh_ancestry_chr20_af_sites")
  registry <- registry[match(ids, registry$id), ]
  registry$cache_relpath <- c("pedigree.txt", "children.bcf", "sites.parquet")
  pedigree <- data.frame(
    FamilyID = paste0("family", 1:4), SampleID = paste0("child", 1:4),
    FatherID = paste0("father", 1:4), MotherID = paste0("mother", 1:4), Sex = "1",
    Population = c("CEU", "CEU", "YRI", "YRI"), Superpopulation = c("EUR", "EUR", "AFR", "AFR")
  )
  write.table(pedigree, file.path(root, "pedigree.txt"), quote = FALSE, row.names = FALSE)
  vcf <- file.path(root, "children.vcf")
  writeLines(c(
    "##fileformat=VCFv4.2", "##contig=<ID=chr20,length=100>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tchild1\tchild2\tchild3\tchild4",
    "chr20\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0|0\t0|1\t1|0\t1|1",
    "chr20\t2\t.\tG\tA\t.\tPASS\t.\tGT\t1|1\t1|0\t0|1\t0|0"
  ), vcf)
  bcf <- file.path(root, "children.bcf")
  status <- system2(RBCFTools::bcftools_path(), c("view", "-Ob", "-o", shQuote(bcf), shQuote(vcf)))
  expect_equal(status, 0L)
  con <- dbConnect(duckdb(shared_home = FALSE))
  parquet <- file.path(root, "sites.parquet")
  dbExecute(con, sprintf(paste(
    "COPY (SELECT 'chr20' AS chrom, 1::BIGINT AS pos, 'A' AS ref, 'C' AS alt,",
    "[0.0]::FLOAT[] AS INFO_AF_EUR, [1.0]::FLOAT[] AS INFO_AF_AFR UNION ALL",
    "SELECT 'chr20', 2, 'G', 'A', [1.0]::FLOAT[], [0.0]::FLOAT[]) TO '%s' (FORMAT PARQUET)"
  ), parquet))
  dbDisconnect(con, shutdown = TRUE)
  registry$supplier_identity <- paste0("sha256=", vapply(file.path(root, registry$cache_relpath), function(path) digest::digest(file = path, algo = "sha256"), ""))
  registry_path <- file.path(root, "registry.tsv")
  write.table(registry, registry_path, sep = "\t", quote = FALSE, row.names = FALSE)
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry_path, DUCKHTS_CACHE_DIR = root)
  genome_settings$chromosomes <<- 20L
  on.exit({
    Sys.setenv(DUCKHTSBENCH_REGISTRY = original_registry, DUCKHTS_CACHE_DIR = original_cache)
    genome_settings <<- original_settings
    unlink(root, recursive = TRUE)
  })
  extension <- Sys.getenv("DUCKHTS_EXTENSION", "build/release/duckhts.duckdb_extension")
  destination <- file.path(root, "haplotypes.rds")
  genome_stage_chromosome(20L, extension, destination)
  panel <- readRDS(destination)
  expect_equal(panel$sites$f, c(0, 1))
  expect_equal(panel$sites$g, c(1, 0))
  expect_identical(panel$genotypes$child2, as.raw(c(1, 2)))
  expect_identical(panel$genotypes$child3, as.raw(c(2, 1)))
  expect_equal(lengths(panel$selected), c(sample = 2L, other = 2L))
  expect_true(!anyDuplicated(panel$pedigree$FamilyID))
  genome_stage_chromosome(20L, extension, destination)
  expect_identical(readRDS(destination), panel)
  expect_identical(panel$genotypes$child1, as.raw(c(0, 3)))
  expect_identical(panel$genotypes$child4, as.raw(c(3, 0)))
  for (relationship in c("same", "other", "ibd1")) {
    set.seed(332L)
    drawn <- genome_draw_identity(panel$selected, relationship)
    expect_true(drawn$sample != drawn$donor)
    expect_true(drawn$donor %in% if (relationship == "other") panel$selected$other else panel$selected$sample)
  }
  lines <- readLines(vcf)
  lines <- gsub("0|1", "0/1", lines, fixed = TRUE)
  writeLines(lines, vcf)
  status <- system2(RBCFTools::bcftools_path(), c("view", "-Ob", "-o", shQuote(bcf), shQuote(vcf)))
  expect_equal(status, 0L)
  registry$supplier_identity[registry$id == "roh_ancestry_chr20_children_bcf"] <-
    paste0("sha256=", digest::digest(file = bcf, algo = "sha256"))
  write.table(registry, registry_path, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_error(genome_stage_chromosome(20L, extension, paste0(destination, ".unphased")), "Non-phased")
  cat("changed\n", file = file.path(root, "pedigree.txt"), append = TRUE)
  expect_error(genome_stage_chromosome(20L, extension, destination), "supplier SHA-256 identity")
}
test_genome_staging()
