library(tinytest)

registry_path <- Sys.getenv("DUCKHTSBENCH_REGISTRY", unset = "")
if (!nzchar(registry_path)) {
  registry_path <- system.file("benchmark_registry.tsv", package = "duckhtsbench")
}
expect_true(file.exists(registry_path))
Sys.setenv(DUCKHTSBENCH_REGISTRY = registry_path)

registry <- duckhts_bench_registry()
expect_true(all(c("id", "workload", "release", "locator", "access", "cache_relpath", "transform", "consumer") %in% names(registry)))
expect_true(all(nzchar(registry$id)))
expect_equal(length(unique(registry$id)), nrow(registry))
expect_true(all(!grepl("^/", registry$cache_relpath)))
expect_true(all(!grepl("\\.\\.", registry$cache_relpath)))
expect_true(any(registry$id == "revel_v13_grch37"))
expect_match(duckhts_bench_artifact_path("revel_v13_grch37"), "revel_grch37\\.parquet$")
expect_equal(nrow(duckhts_bench_stage_plan("variantkey-providers")), sum(registry$workload == "variantkey-providers"))

old_registry <- Sys.getenv("DUCKHTSBENCH_REGISTRY", unset = NA_character_)
old_cache <- Sys.getenv("DUCKHTS_CACHE_DIR", unset = NA_character_)
tmp <- tempfile("duckhtsbench-fetch-")
dir.create(tmp)
source_path <- file.path(tmp, "source.txt")
writeLines("registry fetch", source_path)
mini_registry <- file.path(tmp, "registry.tsv")
writeLines(c(
  paste(names(registry), collapse = "\t"),
  paste(c("fixture", "fixture", "fixture", "fixture-1", paste0("file://", source_path), "public", "fixture/output.txt", "direct_download", "tinytest", "1", ""), collapse = "\t")
), mini_registry)
Sys.setenv(DUCKHTSBENCH_REGISTRY = mini_registry, DUCKHTS_CACHE_DIR = file.path(tmp, "cache"))
output <- duckhts_bench_fetch("fixture")
expect_equal(readLines(output), "registry fetch")
expect_true(file.exists(paste0(output, ".provenance.tsv")))
if (is.na(old_registry)) Sys.unsetenv("DUCKHTSBENCH_REGISTRY") else Sys.setenv(DUCKHTSBENCH_REGISTRY = old_registry)
if (is.na(old_cache)) Sys.unsetenv("DUCKHTS_CACHE_DIR") else Sys.setenv(DUCKHTS_CACHE_DIR = old_cache)

source_vcf <- file.path(tmp, "HG001.vcf.gz")
writeLines("##fileformat=VCFv4.2", source_vcf)
giab_registry <- file.path(tmp, "giab.tsv")
writeLines(c(
  paste(names(registry), collapse = "\t"),
  paste(c("giab_hg001_grch37_v421", "giab-v4.2.1", "benchmark_vcf", "test", paste0("file://", source_vcf), "public", "datasets/giab/HG001.vcf.gz", "direct_download", "tinytest", "1", ""), collapse = "\t")
), giab_registry)
bcftools <- file.path(tmp, "bcftools")
writeLines(c("#!/usr/bin/env sh", "touch \"$4.tbi\""), bcftools)
Sys.chmod(bcftools, "0755")
Sys.setenv(DUCKHTSBENCH_REGISTRY = giab_registry, DUCKHTS_CACHE_DIR = file.path(tmp, "giab-cache"))
giab <- duckhts_bench_stage_giab("HG001", bcftools)
expect_true(file.exists(giab[["HG001"]]))
expect_true(file.exists(paste0(giab[["HG001"]], ".tbi")))
if (is.na(old_registry)) Sys.unsetenv("DUCKHTSBENCH_REGISTRY") else Sys.setenv(DUCKHTSBENCH_REGISTRY = old_registry)
if (is.na(old_cache)) Sys.unsetenv("DUCKHTS_CACHE_DIR") else Sys.setenv(DUCKHTS_CACHE_DIR = old_cache)

riker_reference <- file.path(tmp, "reference.fa")
riker_cram <- file.path(tmp, "input.cram")
writeLines(">chr1\nA", riker_reference)
writeLines("CRAM", riker_cram)
riker_registry <- file.path(tmp, "riker.tsv")
writeLines(c(
  paste(names(registry), collapse = "\t"),
  paste(c("riker_hg00188_reference", "riker-wgs", "reference", "test", paste0("file://", riker_reference), "public", "unused.fa", "direct_download", "tinytest", "1", ""), collapse = "\t"),
  paste(c("riker_hg00188_cram", "riker-wgs", "input_cram", "test", paste0("file://", riker_cram), "public", "unused.cram", "direct_download", "tinytest", "2", ""), collapse = "\t")
), riker_registry)
samtools <- file.path(tmp, "samtools")
writeLines(c("#!/usr/bin/env sh", "case \"$1\" in", "faidx) touch \"$2.fai\" ;;", "quickcheck) exit 0 ;;", "view) while [ \"$1\" != \"-o\" ]; do shift; done; touch \"$2\" ;;", "index) shift 3; touch \"$1.bai\" ;;", "esac"), samtools)
Sys.chmod(samtools, "0755")
Sys.setenv(DUCKHTSBENCH_REGISTRY = riker_registry, DUCKHTS_CACHE_DIR = file.path(tmp, "riker-cache"))
riker_bam <- duckhts_bench_stage_riker(file.path(tmp, "riker"), 1L, samtools)
expect_true(file.exists(riker_bam))
expect_true(file.exists(paste0(riker_bam, ".bai")))
expect_true(file.exists(file.path(tmp, "riker", "provenance.tsv")))
if (is.na(old_registry)) Sys.unsetenv("DUCKHTSBENCH_REGISTRY") else Sys.setenv(DUCKHTSBENCH_REGISTRY = old_registry)
if (is.na(old_cache)) Sys.unsetenv("DUCKHTS_CACHE_DIR") else Sys.setenv(DUCKHTS_CACHE_DIR = old_cache)
