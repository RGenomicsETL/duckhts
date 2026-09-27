#!/usr/bin/env Rscript
root <- normalizePath(Sys.getenv("DUCKHTS_ROOT", unset = "."), mustWork = TRUE)
package_root <- file.path(root, "r", "duckhtsbench")
Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(package_root, "inst", "benchmark_registry.tsv"))
for (file in c("registry.R", "stage.R", "gffbase.R")) {
  source(file.path(package_root, "R", file))
}
# The duplicate-ID test writes its own fixture and needs only the pinned wheel.
if (!identical(Sys.getenv("GFFBASE_FEATUREDB_WHEEL_ONLY"), "1")) {
  for (id in c("mane_v15_ensembl_gff3", "gencode_v49_basic_gff3")) {
    duckhts_bench_fetch(id)
  }
}
site <- duckhts_bench_stage_gffbase(
  site_dir = duckhts_bench_artifact_path("gffbase_021"),
  python = Sys.getenv("GFFBASE_FEATUREDB_PYTHON", unset = Sys.which("python3")),
  artifact_id = "gffbase_021"
)
cat("GFFBASE_FEATUREDB_SITE=", site, "\n", sep = "")
