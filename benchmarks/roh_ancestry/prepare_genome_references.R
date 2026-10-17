# Stage the per-chromosome reference relations of the genome-wide ROH ancestry
# evaluation: registry artifacts roh_ancestry_chr<N>_reference_long and
# roh_ancestry_chr<N>_af_sites, through
# duckhtsbench::duckhts_bench_stage_roh_ancestry_references(). It reuses an
# output only when its receipt matches the current sources. chr20's two
# artifacts are made by prepare_chr20.R. Every chromosome's two artifacts are
# then checked against the SHA-256 pinned in the registry.
source("benchmarks/roh_ancestry/genome_common.R")
chromosomes <- genome_chromosomes()
staged <- setdiff(chromosomes, 20L)
if (length(staged)) {
  result <- duckhtsbench::duckhts_bench_stage_roh_ancestry_references(
    staged, extension_path = normalizePath("build/release/duckhts.duckdb_extension",
                                           mustWork = TRUE))
  print(result[c("chromosome", "id", "rows", "reused")], row.names = FALSE)
}
for (chromosome in chromosomes) {
  for (id in sprintf(c("roh_ancestry_chr%d_reference_long", "roh_ancestry_chr%d_af_sites"),
                     chromosome)) {
    duckhtsbench::duckhts_bench_validate_identity(id)
    cat(id, "matches its registered SHA-256\n")
  }
}
