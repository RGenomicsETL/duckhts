#!/usr/bin/env Rscript

# Stage the GENCODE mouse vM25 inputs and a deterministic BED-like 16-column sample.
if (!nzchar(Sys.getenv("DUCKHTSBENCH_REGISTRY", unset = ""))) {
  Sys.setenv(DUCKHTSBENCH_REGISTRY = normalizePath(
    "r/duckhtsbench/inst/benchmark_registry.tsv", mustWork = TRUE
  ))
}
ids <- c("tabix_split_gff3", "tabix_split_gtf")
offline <- identical(commandArgs(trailingOnly = TRUE), "--offline")
for (id in ids) {
  path <- duckhtsbench::duckhts_bench_artifact_path(id)
  if (offline) {
    if (!file.exists(path)) stop("missing cached input: ", path)
    duckhtsbench::duckhts_bench_validate_identity(id, path)
    duckhtsbench::duckhts_bench_write_provenance(id, path)
  } else {
    duckhtsbench::duckhts_bench_fetch(id)
  }
}

stage_tabix_bed <- function(input, output) {
  id <- "tabix_split_bed"
  if (file.exists(output)) {
    valid <- tryCatch({
      duckhtsbench::duckhts_bench_validate_identity(id, output)
      TRUE
    }, error = function(error) FALSE)
    if (valid) {
      duckhtsbench::duckhts_bench_write_provenance(id, output)
      return(invisible(output))
    }
    if (unlink(c(output, paste0(output, ".provenance.tsv")), force = TRUE) != 0L) {
      stop("could not remove invalid tabix split artifact: ", output)
    }
  }

  con <- gzfile(input, "rt")
  lines <- readLines(con, n = 120000L)
  close(con)
  lines <- lines[!startsWith(lines, "#")]
  if (length(lines) < 100000L) stop("GENCODE input has fewer than 100000 records")
  lines <- lines[seq_len(100000L)]
  fields <- strsplit(lines, "\t", fixed = TRUE)
  if (!all(lengths(fields) == 9L)) stop("GENCODE GFF3 rows must have nine fields")
  bed <- vapply(fields, function(row) paste(c(
    row[[1L]], as.integer(row[[4L]]) - 1L, row[[5L]], row[[3L]],
    row[[6L]], row[[7L]], row[[2L]], row[[8L]], row[[9L]],
    row[[1L]], row[[3L]], row[[4L]], row[[5L]], row[[7L]],
    row[[8L]], row[[9L]]
  ), collapse = "\t"), character(1))
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(basename(output), ".partial-"), tmpdir = dirname(output))
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  writeLines(bed, temporary, useBytes = TRUE)
  duckhtsbench::duckhts_bench_validate_identity(id, temporary)
  if (!file.rename(temporary, output)) stop("could not publish tabix split artifact: ", output)
  duckhtsbench::duckhts_bench_write_provenance(id, output)
  invisible(output)
}

stage_tabix_bed(
  duckhtsbench::duckhts_bench_artifact_path(ids[[1L]]),
  duckhtsbench::duckhts_bench_artifact_path("tabix_split_bed")
)
cat("Staged", length(ids) + 1L, "tabix split inputs\n")
