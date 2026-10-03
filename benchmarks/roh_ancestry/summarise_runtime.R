cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-ancestry")
arms <- c("pooled", "single_AFR", "single_AMR", "single_EUR", "single_EAS", "ancestry")
rows <- lapply(arms, function(arm) {
  do.call(rbind, lapply(1:3, function(repetition) {
    path <- file.path(cache, sprintf("time_%s_%d.txt", arm, repetition))
    if (!file.exists(path)) stop("missing runtime record: ", path, call. = FALSE)
    lines <- readLines(path, warn = FALSE)
    elapsed_line <- grep("Elapsed (wall clock) time", lines, value = TRUE, fixed = TRUE)
    rss_line <- grep("Maximum resident set size", lines, value = TRUE)
    if (length(elapsed_line) != 1L || length(rss_line) != 1L) {
      stop("runtime record is incomplete: ", path, call. = FALSE)
    }
    elapsed <- sub("^.*: ", "", elapsed_line)
    parts <- as.numeric(strsplit(elapsed, ":", fixed = TRUE)[[1L]])
    elapsed_seconds <- if (length(parts) == 3L) {
      parts[[1L]] * 3600 + parts[[2L]] * 60 + parts[[3L]]
    } else if (length(parts) == 2L) {
      parts[[1L]] * 60 + parts[[2L]]
    } else stop("unrecognized elapsed time: ", elapsed, call. = FALSE)
    data.frame(chromosome = 20L, arm = arm, repetition = repetition,
      elapsed_seconds = elapsed_seconds,
      peak_rss_gib = as.numeric(sub("^.*: ", "", rss_line)) / 1024^2,
      threads = 4L, batch_size = if (arm == "ancestry") 16L else 64L,
      stringsAsFactors = FALSE)
  }))
})
result <- do.call(rbind, rows)
utils::write.csv(result, "benchmarks/results/roh-ancestry/chr20/runtime_replicates.csv",
                 row.names = FALSE)
print(result)
