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
expect_equal(duckhts_bench_cache_path("nested/artifact"), file.path(tmp, "cache", "nested", "artifact"))
expect_error(duckhts_bench_cache_path("../outside"), "cache-contained")
expect_error(duckhts_bench_cache_path("nested\\..\\outside"), "cache-contained")
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
  paste(c("liftover_grch38_fasta", "liftover", "destination", "test", paste0("file://", lift_dst), "public", "destination.fa", "direct_download", "tinytest", "2", ""), collapse = "\t"),
  paste(c("liftover_grch37_grch38_chain", "liftover", "chain", "test", paste0("file://", lift_chain), "public", "chain.gz", "direct_download", "tinytest", "3", ""), collapse = "\t")
), lift_registry)
Sys.setenv(DUCKHTSBENCH_REGISTRY = lift_registry, DUCKHTS_CACHE_DIR = file.path(tmp, "lift-cache"))
lift_paths <- duckhts_bench_stage_liftover(file.path(tmp, "lift"), samtools, Sys.which("gzip"))
expect_true(all(file.exists(lift_paths)))
expect_true(all(file.exists(paste0(lift_paths[1:2], ".fai"))))
expect_true(file.exists(file.path(tmp, "lift", "provenance.tsv")))
if (is.na(old_registry)) Sys.unsetenv("DUCKHTSBENCH_REGISTRY") else Sys.setenv(DUCKHTSBENCH_REGISTRY = old_registry)
if (is.na(old_cache)) Sys.unsetenv("DUCKHTS_CACHE_DIR") else Sys.setenv(DUCKHTS_CACHE_DIR = old_cache)

gff_registry <- file.path(tmp, "gffbase.tsv")
writeLines(c(
  paste(names(registry), collapse = "\t"),
  paste(c("gffbase_010", "gffbase", "python_package", "test", "https://example.test/gffbase", "public", "site", "pip_install", "tinytest", "1", ""), collapse = "\t")
), gff_registry)
python <- file.path(tmp, "python3")
writeLines(c("#!/usr/bin/env sh", "while [ \"$1\" != \"--target\" ]; do shift; done", "mkdir -p \"$2/gffbase\""), python)
Sys.chmod(python, "0755")
Sys.setenv(DUCKHTSBENCH_REGISTRY = gff_registry, DUCKHTS_CACHE_DIR = file.path(tmp, "gff-cache"))
gff_site <- duckhts_bench_stage_gffbase(file.path(tmp, "gff-site"), python)
expect_true(dir.exists(file.path(gff_site, "gffbase")))
expect_true(file.exists(file.path(gff_site, "provenance.tsv")))
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
if (is.na(old_registry)) Sys.unsetenv("DUCKHTSBENCH_REGISTRY") else Sys.setenv(DUCKHTSBENCH_REGISTRY = old_registry)
