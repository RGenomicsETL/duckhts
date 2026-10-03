cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-ancestry")
repo_results <- "benchmarks/results/roh-ancestry/chr20"
population_children <- read.table(file.path(cache, "pedigree.txt"), header = TRUE,
                                  stringsAsFactors = FALSE)
population_children <- population_children[
  population_children$FatherID != "0" & population_children$MotherID != "0" &
    population_children$Population %in% c("ACB", "ASW", "CLM", "MXL", "PEL", "PUR",
                                          "YRI", "ESN", "CEU", "CHS"),
  c("SampleID", "Population"), drop = FALSE]
truth <- read.csv(file.path(cache, "chr20.truth_intervals.csv"), stringsAsFactors = FALSE)
children <- population_children$SampleID
arms <- c("pooled", "single_population", "ancestry_tuned")
thresholds <- c(1000000L, 2000000L, 5000000L)
mode_for_population <- c(ACB = "single_AFR", ASW = "single_AFR", ESN = "single_AFR",
                         YRI = "single_AFR", CLM = "single_AMR", MXL = "single_AMR",
                         PEL = "single_AMR", PUR = "single_AMR", CEU = "single_EUR",
                         CHS = "single_EAS")
canonical_segments <- function(x) {
  x <- x[order(x$sample, x$chrom, x$start, x$end), , drop = FALSE]
  x[c("sample", "chrom", "start", "end", "length", "n_markers", "quality")]
}
read_replicates <- function(arm) {
  paths <- file.path(cache, sprintf("chr20.arm_%s_%d.csv", arm, 1:3))
  if (!all(file.exists(paths))) stop("three completed arm repetitions are required")
  result <- lapply(paths, read.csv, stringsAsFactors = FALSE)
  reference <- canonical_segments(result[[1L]])
  if (!all(vapply(result[-1L], function(x) {
    isTRUE(all.equal(reference, canonical_segments(x), tolerance = 1e-12,
                     check.attributes = FALSE))
  }, logical(1L)))) stop("repeated arm output differs")
  result
}
replicates <- setNames(lapply(c("pooled", "ancestry", "single_AFR", "single_AMR",
                                "single_EUR", "single_EAS"), read_replicates),
                       c("pooled", "ancestry", "single_AFR", "single_AMR",
                         "single_EUR", "single_EAS"))
eligible_ids <- list(
  ancestry = population_children$SampleID,
  single_AFR = population_children$SampleID[population_children$Population %in%
    c("ACB", "ASW", "ESN", "YRI")],
  single_AMR = population_children$SampleID[population_children$Population %in%
    c("CLM", "MXL", "PEL", "PUR")],
  single_EUR = population_children$SampleID[population_children$Population == "CEU"],
  single_EAS = population_children$SampleID[population_children$Population == "CHS"]
)
for (arm in names(eligible_ids)) {
  for (rep in seq_len(3L)) {
    observed_ids <- replicates[[arm]][[rep]]$sample
    if (anyNA(observed_ids) || any(!observed_ids %in% eligible_ids[[arm]])) {
      stop("an arm output contains an ineligible child", call. = FALSE)
    }
  }
}
segments <- list(
  pooled = replicates$pooled[[1L]],
  ancestry_tuned = replicates$ancestry[[1L]]
)
single_parts <- lapply(names(mode_for_population), function(population) {
  source <- mode_for_population[[population]]
  data <- replicates[[source]][[1L]]
  ids <- population_children$SampleID[population_children$Population == population]
  data[data$sample %in% ids, , drop = FALSE]
})
names(single_parts) <- names(mode_for_population)
segments$single_population <- do.call(rbind, single_parts)
segments <- lapply(segments, function(x) {
  x <- x[!is.na(x$length), , drop = FALSE]
  x[order(x$sample, x$start), , drop = FALSE]
})
for (arm in names(segments)) {
  x <- segments[[arm]]
  for (id in unique(x$sample)) {
    one <- x[x$sample == id, , drop = FALSE]
    if (nrow(one) > 1L && any(one$start[-1L] <= one$end[-nrow(one)])) {
      stop("overlapping ROH calls prevent unique interval coverage")
    }
  }
}
covered_bp <- function(start, end, intervals) {
  if (!nrow(intervals)) return(0)
  sum(pmax(0, pmin(end, intervals$end) - pmax(start, intervals$start) + 1))
}
child_metrics <- do.call(rbind, lapply(arms, function(arm) {
  data <- segments[[arm]]
  do.call(rbind, lapply(children, function(id) {
    calls <- data[data$sample == id, , drop = FALSE]
    truth_calls <- truth[truth$sample_id == id, c("start", "end"), drop = FALSE]
    if (nrow(calls)) {
      support_bp <- vapply(seq_len(nrow(calls)), function(i) {
        covered_bp(calls$start[[i]], calls$end[[i]], truth_calls)
      }, numeric(1L))
      supported <- support_bp / calls$length >= 0.9
    } else {
      support_bp <- numeric()
      supported <- logical()
    }
    row <- data.frame(sample_id = id, arm = arm, total_roh_bp = sum(calls$length),
                      froh_chr20 = sum(calls$length) / 64444167,
                      supported_count_1mb = sum(calls$length >= thresholds[[1L]] & supported),
                      supported_length_1mb = sum(calls$length[calls$length >= thresholds[[1L]] & supported]),
                      truth_bp = sum(truth_calls$end - truth_calls$start + 1),
                      truth_covered_bp = if (nrow(truth_calls)) {
                        sum(vapply(seq_len(nrow(truth_calls)), function(i) {
                          covered_bp(truth_calls$start[[i]], truth_calls$end[[i]],
                                     calls[c("start", "end")])
                        }, numeric(1L)))
                      } else 0)
    for (index in seq_along(thresholds)) {
      name <- c("1mb", "2mb", "5mb")[[index]]
      selected <- calls$length >= thresholds[[index]]
      row[[paste0("unsupported_count_", name)]] <- sum(selected & !supported)
      row[[paste0("unsupported_length_", name)]] <- sum(calls$length[selected & !supported])
    }
    row
  }))
}))
child_metrics <- merge(child_metrics, population_children, by.x = "sample_id", by.y = "SampleID")
summary_rows <- do.call(rbind, lapply(arms, function(arm) {
  do.call(rbind, lapply(unique(population_children$Population), function(population) {
    x <- child_metrics[child_metrics$arm == arm & child_metrics$Population == population, , drop = FALSE]
    unsupported_count <- vapply(c("1mb", "2mb", "5mb"), function(n) {
      mean(x[[paste0("unsupported_count_", n)]])
    }, numeric(1L))
    unsupported_length <- vapply(c("1mb", "2mb", "5mb"), function(n) {
      mean(x[[paste0("unsupported_length_", n)]])
    }, numeric(1L))
    data.frame(Population = population, Arm = arm,
      `Unsupported ROH count / child (>=1 Mb; >=2 Mb; >=5 Mb)` = paste(formatC(unsupported_count, digits = 3, format = "f"), collapse = "; "),
      `Unsupported ROH length / child (>=1 Mb; >=2 Mb; >=5 Mb)` = paste(formatC(unsupported_length, digits = 0, format = "f"), collapse = "; "),
      `Supported count / child` = formatC(mean(x$supported_count_1mb), digits = 3, format = "f"),
      `Supported length / child` = formatC(mean(x$supported_length_1mb), digits = 0, format = "f"),
      FROH = formatC(mean(x$froh_chr20), digits = 6, format = "f"),
      `Truth-run length covered` = if (sum(x$truth_bp) > 0) {
        formatC(sum(x$truth_covered_bp) / sum(x$truth_bp) * 100, digits = 2, format = "f")
      } else "NA",
      check.names = FALSE, stringsAsFactors = FALSE)
  }))
}))
write.csv(child_metrics[order(child_metrics$Population, child_metrics$arm, child_metrics$sample_id), ],
          file.path(cache, "chr20.child_metrics.csv"), row.names = FALSE)
write.csv(summary_rows[order(summary_rows$Population, summary_rows$Arm), ],
          file.path(repo_results, "population_arm_summary.csv"), row.names = FALSE,
          quote = TRUE)
write.csv(child_metrics, file.path(repo_results, "child_metrics.csv"), row.names = FALSE)
