library(tinytest)

registry_path <- Sys.getenv("DUCKHTSBENCH_REGISTRY", unset = "")
if (!nzchar(registry_path)) {
  registry_path <- system.file("benchmark_registry.tsv", package = "duckhtsbench")
}
expect_true(file.exists(registry_path))
Sys.setenv(DUCKHTSBENCH_REGISTRY = registry_path)
missing_registry <- tempfile("duckhtsbench-missing-registry-")
Sys.setenv(DUCKHTSBENCH_REGISTRY = missing_registry)
expect_error(duckhts_bench_registry(), "does not exist")
Sys.setenv(DUCKHTSBENCH_REGISTRY = registry_path)

registry <- duckhts_bench_registry()
expect_true(all(c("id", "workload", "release", "locator", "access", "cache_relpath", "transform", "consumer") %in% names(registry)))
expect_true(all(nzchar(registry$id)))
expect_equal(length(unique(registry$id)), nrow(registry))
expect_true(all(!grepl("^/", registry$cache_relpath)))
expect_true(all(!grepl("\\.\\.", registry$cache_relpath)))
expect_true(any(registry$id == "variantkey_giab_hg002_v421"))
expect_match(duckhts_bench_artifact_path("variantkey_giab_hg002_v421"), "HG002_GRCh38_1_22_v4\\.2\\.1_benchmark\\.vcf\\.gz$")
expect_equal(registry$id[registry$workload == "giab-benchmark-vcf"], "variantkey_giab_hg002_v421")
expect_equal(nrow(duckhts_bench_stage_plan("giab-benchmark-vcf")), 1L)
expect_match(registry$consumer[registry$id == "variantkey_giab_hg002_v421"],
  "benchmark_bcf_record_cache.Rmd", fixed = TRUE)
expect_match(registry$consumer[registry$id == "variantkey_giab_hg002_v421"],
  "benchmark_genotype_format.Rmd", fixed = TRUE)
expect_equal(registry$id[registry$workload == "genotype-reader"],
             c("geno_hprc_vcfgz", "geno_hprc_bcf"))
expect_equal(registry$id[registry$workload == "genotype-source"],
             c("geno_hprc_source", "geno_hprc_source_tbi"))
expect_equal(registry$transform[registry$id == "geno_hprc_source_tbi"],
             "direct_download")

old_registry <- Sys.getenv("DUCKHTSBENCH_REGISTRY", unset = NA_character_)
old_cache <- Sys.getenv("DUCKHTS_CACHE_DIR", unset = NA_character_)
tmp <- tempfile("duckhtsbench-fetch-")
dir.create(tmp)
source_path <- file.path(tmp, "source.txt")
writeLines("registry fetch", source_path)
source_bytes <- file.info(source_path)$size
mini_registry <- file.path(tmp, "registry.tsv")
writeLines(c(
  paste(names(registry), collapse = "\t"),
  paste(c("fixture", "fixture", "fixture", "fixture-1", paste0("file://", source_path), "public", "fixture/output.txt", "direct_download", "tinytest", "1", paste0("bytes=", source_bytes)), collapse = "\t")
), mini_registry)
Sys.setenv(DUCKHTSBENCH_REGISTRY = mini_registry, DUCKHTS_CACHE_DIR = file.path(tmp, "cache"))
expect_equal(duckhts_bench_cache_path("nested/artifact"), file.path(tmp, "cache", "nested", "artifact"))
expect_error(duckhts_bench_cache_path("../outside"), "cache-contained")
expect_error(duckhts_bench_cache_path("nested\\..\\outside"), "cache-contained")
expect_error(duckhts_bench_cache_path("C:outside"), "cache-contained")
expect_error(duckhts_bench_cache_path("\\\\server\\share"), "cache-contained")
output <- duckhts_bench_fetch("fixture")
expect_equal(readLines(output), "registry fetch")
expect_true(file.exists(paste0(output, ".provenance.tsv")))
writeLines("poisoned cache entry", output)
expect_equal(readLines(duckhts_bench_fetch("fixture")), "registry fetch")
writeLines(c(
  paste(names(registry), collapse = "\t"),
  paste(c("escape", "fixture", "fixture", "fixture-1", paste0("file://", source_path), "public", "../outside", "direct_download", "tinytest", "1", ""), collapse = "\t")
), mini_registry)
expect_error(duckhts_bench_artifact_path("escape"), "cache-contained")
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
writeLines(c("#!/usr/bin/env sh", "printf \"idx\\n\" >\"$4.tbi\""), bcftools)
Sys.chmod(bcftools, "0755")
Sys.setenv(DUCKHTSBENCH_REGISTRY = giab_registry, DUCKHTS_CACHE_DIR = file.path(tmp, "giab-cache"))
giab <- duckhts_bench_stage_giab("HG001", bcftools)
expect_true(file.exists(giab[["HG001"]]))
expect_true(file.exists(paste0(giab[["HG001"]], ".tbi")))
writeLines("stale", paste0(giab[["HG001"]], ".tbi"))
giab <- duckhts_bench_stage_giab("HG001", bcftools)
expect_equal(readLines(paste0(giab[["HG001"]], ".tbi")), "idx")
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
  paste(c("riker_hg00188_cram", "riker-wgs", "input_cram", "test", paste0("file://", riker_cram), "public", "unused.cram", "direct_download", "tinytest", "2", ""), collapse = "\t"),
  paste(c("riker_hg00188_bam", "riker-wgs", "derived_bam", "test", "artifact:riker_hg00188_reference;artifact:riker_hg00188_cram", "local_derived", "relocated/derived.bam", "samtools_view_reference_cram;samtools_index", "tinytest", "3", ""), collapse = "\t")
), riker_registry)
samtools <- file.path(tmp, "samtools")
writeLines(c("#!/usr/bin/env sh", "case \"$1\" in", "faidx) printf \"idx\\n\" >\"$2.fai\" ;;", "quickcheck) exit 0 ;;", "view) while [ \"$1\" != \"-o\" ]; do shift; done; printf \"BAM\\n\" >\"$2\" ;;", "index) shift 3; printf \"BAI\\n\" >\"$1.bai\" ;;", "esac"), samtools)
Sys.chmod(samtools, "0755")
Sys.setenv(DUCKHTSBENCH_REGISTRY = riker_registry, DUCKHTS_CACHE_DIR = file.path(tmp, "riker-cache"))
riker_bam <- duckhts_bench_stage_riker(file.path(tmp, "riker"), 1L, samtools)
expect_true(file.exists(riker_bam))
expect_true(file.exists(paste0(riker_bam, ".bai")))
expect_true(file.exists(file.path(tmp, "riker", "provenance.tsv")))
writeLines("stale", riker_bam)
riker_bam <- duckhts_bench_stage_riker(file.path(tmp, "riker"), 1L, samtools)
expect_equal(readLines(riker_bam), "BAM")
Sys.setenv(DUCKHTS_CACHE_DIR = file.path(tmp, "riker-registry-cache"))
riker_registry_bam <- duckhts_bench_stage_riker(NULL, 1L, samtools)
expect_equal(riker_registry_bam, file.path(tmp, "riker-registry-cache", "relocated", "derived.bam"))
if (is.na(old_registry)) Sys.unsetenv("DUCKHTSBENCH_REGISTRY") else Sys.setenv(DUCKHTSBENCH_REGISTRY = old_registry)
if (is.na(old_cache)) Sys.unsetenv("DUCKHTS_CACHE_DIR") else Sys.setenv(DUCKHTS_CACHE_DIR = old_cache)

lift_gz <- file.path(tmp, "grch37.fa.gz")
connection <- gzfile(lift_gz, "wb")
writeLines(c(">1", "ACGT"), connection)
close(connection)
lift_dst <- file.path(tmp, "grch38.fa")
lift_chain <- file.path(tmp, "chain.gz")
writeLines(c(">chr1", "ACGT"), lift_dst)
writeLines("chain", lift_chain)
lift_registry <- file.path(tmp, "liftover.tsv")
writeLines(c(
  paste(names(registry), collapse = "\t"),
  paste(c("liftover_grch37_fasta_gz", "liftover", "source", "test", paste0("file://", lift_gz), "public", "source.fa.gz", "direct_download", "tinytest", "1", ""), collapse = "\t"),
  paste(c("liftover_grch37_fasta", "liftover", "source", "test", "artifact:liftover_grch37_fasta_gz", "local_derived", "source.fa", "decompress", "tinytest", "2", ""), collapse = "\t"),
  paste(c("liftover_grch38_fasta", "liftover", "destination", "test", paste0("file://", lift_dst), "public", "destination.fa", "direct_download", "tinytest", "3", ""), collapse = "\t"),
  paste(c("liftover_grch37_grch38_chain", "liftover", "chain", "test", paste0("file://", lift_chain), "public", "chain.gz", "direct_download", "tinytest", "3", ""), collapse = "\t")
), lift_registry)
Sys.setenv(DUCKHTSBENCH_REGISTRY = lift_registry, DUCKHTS_CACHE_DIR = file.path(tmp, "lift-cache"))
lift_paths <- duckhts_bench_stage_liftover(file.path(tmp, "lift"), samtools, Sys.which("gzip"))
expect_true(all(file.exists(lift_paths)))
expect_true(all(file.exists(paste0(lift_paths[1:2], ".fai"))))
expect_true(file.exists(file.path(tmp, "lift", "provenance.tsv")))
writeLines("stale", lift_paths[["source_fasta"]])
lift_paths <- duckhts_bench_stage_liftover(file.path(tmp, "lift"), samtools, Sys.which("gzip"))
expect_equal(readLines(lift_paths[["source_fasta"]])[[1L]], ">1")
Sys.setenv(DUCKHTS_CACHE_DIR = file.path(tmp, "lift-registry-cache"))
lift_registry_paths <- duckhts_bench_stage_liftover(NULL, samtools, Sys.which("gzip"))
expect_equal(lift_registry_paths[["source_fasta"]], file.path(tmp, "lift-registry-cache", "source.fa"))
expect_equal(lift_registry_paths[["destination_fasta"]], file.path(tmp, "lift-registry-cache", "destination.fa"))
expect_equal(lift_registry_paths[["chain"]], file.path(tmp, "lift-registry-cache", "chain.gz"))
if (is.na(old_registry)) Sys.unsetenv("DUCKHTSBENCH_REGISTRY") else Sys.setenv(DUCKHTSBENCH_REGISTRY = old_registry)
if (is.na(old_cache)) Sys.unsetenv("DUCKHTS_CACHE_DIR") else Sys.setenv(DUCKHTS_CACHE_DIR = old_cache)

norm_source <- file.path(tmp, "HG00096.g.vcf.gz")
writeLines("##fileformat=VCFv4.2", norm_source)
norm_registry <- file.path(tmp, "norm.tsv")
writeLines(c(
  paste(names(registry), collapse = "\t"),
  paste(c("norm_hg00096_raw_gvcf", "norm", "raw_gvcf", "test", paste0("file://", norm_source), "public", "norm/HG00096.raw.g.vcf.gz", "direct_download", "tinytest", "1", ""), collapse = "\t"),
  paste(c("norm_hg00096_chr22_20m_30m", "norm", "gvcf_slice", "test", "artifact:norm_hg00096_raw_gvcf", "local_derived", "norm/HG00096.g.vcf.gz", "bcftools_view_region_sample;tabix_index", "tinytest", "2", ""), collapse = "\t")
), norm_registry)
norm_bcftools <- file.path(tmp, "norm-bcftools")
writeLines(c("#!/usr/bin/env sh", "while [ \"$1\" != \"-o\" ]; do shift; done", "printf \"VCF\\n\" >\"$2\""), norm_bcftools)
norm_tabix <- file.path(tmp, "norm-tabix")
writeLines(c("#!/usr/bin/env sh", "for arg in \"$@\"; do target=$arg; done", "printf TBI >\"$target.tbi\""), norm_tabix)
Sys.chmod(c(norm_bcftools, norm_tabix), "0755")
Sys.setenv(DUCKHTSBENCH_REGISTRY = norm_registry, DUCKHTS_CACHE_DIR = file.path(tmp, "norm-cache"))
norm_output <- duckhts_bench_stage_norm(file.path(tmp, "norm-stage", "HG00096.g.vcf.gz"), norm_bcftools, norm_tabix, 1L)
expect_true(file.exists(norm_output))
expect_true(file.exists(paste0(norm_output, ".tbi")))
expect_true(file.exists(paste0(norm_output, ".provenance.tsv")))
writeLines("stale", norm_output)
norm_output <- duckhts_bench_stage_norm(norm_output, norm_bcftools, norm_tabix, 1L)
expect_equal(readLines(norm_output), "VCF")
if (is.na(old_registry)) Sys.unsetenv("DUCKHTSBENCH_REGISTRY") else Sys.setenv(DUCKHTSBENCH_REGISTRY = old_registry)
if (is.na(old_cache)) Sys.unsetenv("DUCKHTS_CACHE_DIR") else Sys.setenv(DUCKHTS_CACHE_DIR = old_cache)

gff_wheel <- file.path(tmp, "gffbase.whl")
writeLines("verified wheel", gff_wheel)
gff_registry <- file.path(tmp, "gffbase.tsv")
writeLines(c(
  paste(names(registry), collapse = "\t"),
  paste(c("gffbase_010_linux_x86_64_wheel", "gffbase", "python_wheel", "test", paste0("file://", gff_wheel), "public", "wheel.whl", "direct_download", "tinytest", "1", ""), collapse = "\t"),
  paste(c("gffbase_010", "gffbase", "python_package", "test", "artifact:gffbase_010_linux_x86_64_wheel", "local_derived", "site", "pip_install_verified_wheel", "tinytest", "2", ""), collapse = "\t")
), gff_registry)
python <- file.path(tmp, "python3")
writeLines(c("#!/usr/bin/env sh", "if [ \"$1\" = \"-c\" ]; then exit 0; fi", "while [ \"$1\" != \"--target\" ]; do shift; done", "mkdir -p \"$2/gffbase\""), python)
Sys.chmod(python, "0755")
Sys.setenv(DUCKHTSBENCH_REGISTRY = gff_registry, DUCKHTS_CACHE_DIR = file.path(tmp, "gff-cache"))
gff_site <- duckhts_bench_stage_gffbase(file.path(tmp, "gff-site"), python)
expect_true(dir.exists(file.path(gff_site, "gffbase")))
expect_true(file.exists(file.path(gff_site, "provenance.tsv")))

featuredb_registry <- file.path(tmp, "featuredb.tsv")
featuredb_input <- file.path(tmp, "featuredb.gff3.gz")
writeLines("##gff-version 3", featuredb_input)
checksum <- function(path) {
  sha256_bin <- Sys.which("sha256sum")
  if (!nzchar(sha256_bin)) sha256_bin <- Sys.which("shasum")
  args <- if (basename(sha256_bin) == "shasum") c("-a", "256", shQuote(path)) else shQuote(path)
  fields <- strsplit(trimws(system2(sha256_bin, args, stdout = TRUE)),
                     "[[:space:]]+")[[1L]]
  paste0("sha256=", fields[[1L]], ";bytes=", file.info(path)$size)
}
featuredb_rows <- list(
  c("mane_v15_ensembl_gff3", "gffbase-featuredb", "gff3", "test",
    paste0("file://", featuredb_input), "public", "inputs/mane.gff3.gz",
    "direct_download", "tinytest", "1", checksum(featuredb_input)),
  c("gencode_v49_basic_gff3", "gffbase-featuredb", "gff3", "test",
    paste0("file://", featuredb_input), "public", "inputs/gencode.gff3.gz",
    "direct_download", "tinytest", "2", checksum(featuredb_input)),
  c("gffbase_021_linux_x86_64_wheel", "gffbase-featuredb", "python_wheel", "test",
    paste0("file://", gff_wheel), "public", "wheel.whl", "direct_download",
    "tinytest", "3", checksum(gff_wheel)),
  c("gffbase_021", "gffbase-featuredb", "python_package", "test",
    "artifact:gffbase_021_linux_x86_64_wheel", "local_derived", "site-021",
    "pip_install_verified_wheel", "tinytest", "4", "")
)
writeLines(c(paste(names(registry), collapse = "\t"),
             vapply(featuredb_rows, paste, character(1), collapse = "\t")), featuredb_registry)
Sys.setenv(DUCKHTSBENCH_REGISTRY = featuredb_registry,
           DUCKHTS_CACHE_DIR = file.path(tmp, "featuredb-cache"))
for (id in c("mane_v15_ensembl_gff3", "gencode_v49_basic_gff3")) {
  staged <- duckhts_bench_fetch(id)
  expect_equal(readLines(staged), "##gff-version 3")
  expect_true(file.exists(paste0(staged, ".provenance.tsv")))
}
featuredb_site <- duckhts_bench_stage_gffbase(
  python = python, artifact_id = "gffbase_021")
expect_equal(featuredb_site, duckhts_bench_artifact_path("gffbase_021"))
expect_true(dir.exists(file.path(featuredb_site, "gffbase")))
expect_true(any(grepl("version\\t0.2.1", readLines(file.path(featuredb_site, "provenance.tsv")))))
if (is.na(old_registry)) Sys.unsetenv("DUCKHTSBENCH_REGISTRY") else Sys.setenv(DUCKHTSBENCH_REGISTRY = old_registry)
if (is.na(old_cache)) Sys.unsetenv("DUCKHTS_CACHE_DIR") else Sys.setenv(DUCKHTS_CACHE_DIR = old_cache)

duckbed_registry <- file.path(tmp, "duckbedqc.tsv")
commit <- "118fc21c6cde9d680989dd4d1b613789539469f3"
writeLines(c(
  paste(names(registry), collapse = "\t"),
  paste(c("duckbedqc_118fc21", "cgranges", "beds", "test", "https://example.test/DuckBedQC.git", "public", "checkout", "git_clone", "tinytest", "1", paste0("git_commit=", commit)), collapse = "\t")
), duckbed_registry)
git <- file.path(tmp, "git")
writeLines(c("#!/usr/bin/env sh", "case \"$1\" in", "clone) mkdir -p \"$3/.git\" \"$3/data\"; touch \"$3/data/GRCh38_exons.bed\" \"$3/data/GRCh38_illumina_clinical_regions_v100.39.0.bed\" ;;", "-C) if [ \"$3\" = \"rev-parse\" ]; then printf '%s\\n' '118fc21c6cde9d680989dd4d1b613789539469f3'; fi ;;", "esac"), git)
Sys.chmod(git, "0755")
Sys.setenv(DUCKHTSBENCH_REGISTRY = duckbed_registry, DUCKHTS_CACHE_DIR = file.path(tmp, "duckbed-cache"))
duckbed_dir <- duckhts_bench_stage_duckbedqc(file.path(tmp, "duckbed"), git)
expect_true(file.exists(file.path(duckbed_dir, "provenance.tsv")))
if (is.na(old_registry)) Sys.unsetenv("DUCKHTSBENCH_REGISTRY") else Sys.setenv(DUCKHTSBENCH_REGISTRY = old_registry)
if (is.na(old_cache)) Sys.unsetenv("DUCKHTS_CACHE_DIR") else Sys.setenv(DUCKHTS_CACHE_DIR = old_cache)

identity_registry <- file.path(tmp, "identity.tsv")
writeLines(c(
  paste(names(registry), collapse = "\t"),
  paste(c("identity_fixture", "fixture", "fixture", "test", "file:///dev/null", "public", "unused", "direct_download", "tinytest", "1", "bytes=999"), collapse = "\t")
), identity_registry)
identity_file <- file.path(tmp, "identity.txt")
writeLines("not 999 bytes", identity_file)
Sys.setenv(DUCKHTSBENCH_REGISTRY = identity_registry)
expect_error(duckhts_bench_validate_identity("identity_fixture", identity_file), "supplier byte identity")
sha256_bin <- Sys.which("sha256sum")
if (!nzchar(sha256_bin)) sha256_bin <- Sys.which("shasum")
sha256_args <- if (basename(sha256_bin) == "shasum") c("-a", "256", shQuote(identity_file)) else shQuote(identity_file)
identity_sha256 <- strsplit(trimws(system2(sha256_bin, sha256_args, stdout = TRUE)), "[[:space:]]+")[[1L]][[1L]]
writeLines(c(
  paste(names(registry), collapse = "\t"),
  paste(c("sha256_fixture", "fixture", "fixture", "test", "file:///dev/null", "public", "unused", "direct_download", "tinytest", "1", paste0("sha256=", identity_sha256)), collapse = "\t")
), identity_registry)
Sys.setenv(DUCKHTSBENCH_REGISTRY = identity_registry)
expect_true(duckhts_bench_validate_identity("sha256_fixture", identity_file))
local({
  checksum_bin <- file.path(tmp, "checksum commands")
  dir.create(checksum_bin)
  checksum_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = checksum_path), add = TRUE)
  Sys.setenv(PATH = paste(checksum_bin, checksum_path, sep = .Platform$path.sep))
  for (algorithm in c("sha256", "sum")) {
    identity <- if (algorithm == "sha256") paste0("sha256=", identity_sha256) else "sum=123;blocks=1"
    expected_output <- if (algorithm == "sha256") paste(identity_sha256, "file") else "123 1 file"
    writeLines(c(paste(names(registry), collapse = "\t"),
      paste(c("checksum_fixture", "fixture", "fixture", "test", "file:///dev/null",
        "public", "unused", "direct_download", "tinytest", "1", identity), collapse = "\t")
    ), identity_registry)
    checksum_command <- file.path(checksum_bin, if (algorithm == "sha256") "sha256sum" else "sum")
    for (command_status in c(0L, 23L)) {
      writeLines(c("#!/bin/sh", paste("printf '%s\\n'", shQuote(expected_output)),
        paste("exit", command_status)), checksum_command)
      Sys.chmod(checksum_command, "0755")
      if (command_status == 0L) {
        expect_true(duckhts_bench_validate_identity("checksum_fixture", identity_file))
      } else {
        expect_error(duckhts_bench_validate_identity("checksum_fixture", identity_file),
          "supplier checksum command failed")
      }
    }
    for (output_lines in list(character(), c(expected_output, expected_output))) {
      writeLines(c("#!/bin/sh", if (length(output_lines))
        paste("printf '%s\\n'", shQuote(output_lines)), "exit 0"), checksum_command)
      Sys.chmod(checksum_command, "0755")
      expect_error(duckhts_bench_validate_identity("checksum_fixture", identity_file),
        "supplier checksum command failed")
    }
  }
})
if (is.na(old_registry)) Sys.unsetenv("DUCKHTSBENCH_REGISTRY") else Sys.setenv(DUCKHTSBENCH_REGISTRY = old_registry)
