#!/usr/bin/env Rscript

script <- grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]
package_dir <- normalizePath(file.path(dirname(sub("^--file=", "", script)), ".."))
repo <- normalizePath(file.path(package_dir, "..", ".."))
if (!nzchar(Sys.getenv("DUCKHTSBENCH_REGISTRY", unset = ""))) {
  Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(package_dir, "inst/benchmark_registry.tsv"))
}
for (file in c("registry.R", "stage.R", "vep_cache.R")) source(file.path(package_dir, "R", file))
parser <- optparse::OptionParser(option_list = list(
  optparse::make_option("--artifact", default = "vep116_grch38_cache_chr21"),
  optparse::make_option("--archive", default = NULL,
    help = "local archive, only for a registry-pinned content checksum")
))
options <- optparse::parse_args(parser)
print(duckhts_bench_stage_vep_cache(repo, options$artifact, options$archive))
