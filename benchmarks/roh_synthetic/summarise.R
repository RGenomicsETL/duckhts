# Score the three ROH arms on the synthetic chr20 children against the planted
# truth, as benchmark_roh_synthetic_truth.Rmd declares.
#
# Coordinates are 1-based and inclusive. Every quantity is a count of bases per
# child; populations are summed only after each child's values are kept.
options(scipen = 999)
artifact <- duckhtsbench::duckhts_bench_artifact_path
cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-synthetic")
results <- "benchmarks/results/roh-synthetic"
dir.create(results, recursive = TRUE, showWarnings = FALSE)

chromosome_bp <- 64444167
window_bp <- 100000
primary_min_call_bp <- 500000
min_call_bp <- c(0, primary_min_call_bp, 1000000)
margin <- 0.01
bootstrap_seed <- 20261005
bootstrap_draws <- 2000L
admixed <- c("ACB", "ASW", "CLM", "MXL", "PEL", "PUR")
controls <- c("YRI", "ESN", "CEU", "CHS")
single_arm <- c(ACB = "single_AFR", ASW = "single_AFR", ESN = "single_AFR",
                YRI = "single_AFR", CLM = "single_AMR", MXL = "single_AMR",
                PEL = "single_AMR", PUR = "single_AMR", CEU = "single_EUR",
                CHS = "single_EAS")
arms <- c("pooled", "single_population", "ancestry_tuned")

truth <- utils::read.csv(artifact("roh_synthetic_chr20_truth"), stringsAsFactors = FALSE)
children <- unique(truth$sample_id)
pedigree <- utils::read.table(artifact("roh_ancestry_pedigree"), header = TRUE,
                              stringsAsFactors = FALSE)
population <- stats::setNames(pedigree$Population[match(children, pedigree$SampleID)], children)
if (anyNA(population) || !all(population %in% names(single_arm))) {
  stop("every synthetic child needs a pedigree population", call. = FALSE)
}

# The source children's 100 kb windows: called genotypes and heterozygous calls.
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:", config = list(threads = 2L))
windows <- DBI::dbGetQuery(con, sprintf(
  "SELECT sample_id, window_id, sites, heterozygotes FROM read_parquet(%s) ORDER BY sample_id, window_id",
  as.character(DBI::dbQuoteString(con, artifact("roh_ancestry_chr20_truth_windows")))))
DBI::dbDisconnect(con, shutdown = TRUE)
window_count <- (chromosome_bp - 1) %/% window_bp + 1
if (nrow(windows) != length(children) * window_count || !all(children %in% windows$sample_id)) {
  stop("the source windows do not cover every child and window", call. = FALSE)
}
windows$start <- windows$window_id * window_bp + 1
windows$end <- pmin((windows$window_id + 1) * window_bp, chromosome_bp)
windows$class <- ifelse(windows$sites == 0, "no_record",
                        ifelse(windows$heterozygotes <= 1, "low_heterozygosity", "other"))
windows <- split(windows[c("start", "end", "class")], windows$sample_id)

read_arm <- function(arm) {
  path <- file.path(cache, sprintf("chr20.arm_%s.csv", arm))
  if (!file.exists(path)) stop("missing arm output: ", path, call. = FALSE)
  calls <- utils::read.csv(path, stringsAsFactors = FALSE)
  calls <- calls[!is.na(calls$length), c("sample", "start", "end", "length"), drop = FALSE]
  if (!all(calls$sample %in% children) || any(calls$end - calls$start + 1 != calls$length)) {
    stop("unexpected calls in arm ", arm, call. = FALSE)
  }
  calls
}
raw <- lapply(c(pooled = "pooled", ancestry = "ancestry", single_AFR = "single_AFR",
                single_AMR = "single_AMR", single_EUR = "single_EUR",
                single_EAS = "single_EAS"), read_arm)
single <- do.call(rbind, lapply(names(single_arm), function(p) {
  calls <- raw[[single_arm[[p]]]]
  calls[calls$sample %in% children[population == p], , drop = FALSE]
}))
calls_by_arm <- list(pooled = raw$pooled, single_population = single,
                     ancestry_tuned = raw$ancestry)

# Bases of each [from, to] covered by intervals that do not overlap each other.
covered_bp <- function(starts, ends, from, to) {
  vapply(seq_along(from), function(i) {
    sum(pmax(0, pmin(to[[i]], ends) - pmax(from[[i]], starts) + 1))
  }, numeric(1L))
}
# The parts of the intervals that lie outside every cut.
outside_cuts <- function(starts, ends, cut_starts, cut_ends) {
  for (i in seq_along(cut_starts)) {
    left_end <- pmin(ends, cut_starts[[i]] - 1)
    right_start <- pmax(starts, cut_ends[[i]] + 1)
    keep_left <- left_end >= starts
    keep_right <- right_start <= ends
    new_starts <- c(starts[keep_left], right_start[keep_right])
    new_ends <- c(left_end[keep_left], ends[keep_right])
    starts <- new_starts
    ends <- new_ends
  }
  list(starts = starts, ends = ends)
}

planted_lengths <- sort(unique(truth$length_bp))
child_rows <- list()
for (arm in arms) {
  arm_calls <- split(calls_by_arm[[arm]], calls_by_arm[[arm]]$sample)
  for (id in children) {
    calls <- arm_calls[[id]]
    if (is.null(calls)) calls <- calls_by_arm[[arm]][0L, , drop = FALSE]
    calls <- calls[order(calls$start), , drop = FALSE]
    if (nrow(calls) > 1L && any(calls$start[-1L] <= calls$end[-nrow(calls)])) {
      stop("overlapping ROH calls prevent unique base counts", call. = FALSE)
    }
    planted <- truth[truth$sample_id == id, , drop = FALSE]
    child_windows <- windows[[id]]
    for (minimum in min_call_bp) {
      selected <- calls[calls$length >= minimum, , drop = FALSE]
      segment_covered <- covered_bp(selected$start, selected$end, planted$start, planted$end)
      outside <- outside_cuts(selected$start, selected$end, planted$start, planted$end)
      window_outside <- covered_bp(outside$starts, outside$ends,
                                   child_windows$start, child_windows$end)
      by_class <- function(class) sum(window_outside[child_windows$class == class])
      row <- data.frame(
        sample_id = id, Population = population[[id]], arm = arm, min_call_bp = minimum,
        calls = nrow(selected), called_bp = sum(selected$length),
        planted_bp = sum(planted$length_bp), planted_covered_bp = sum(segment_covered),
        called_low_heterozygosity_bp = by_class("low_heterozygosity"),
        called_no_record_bp = by_class("no_record"),
        called_unsupported_bp = by_class("other"), stringsAsFactors = FALSE)
      for (length_bp in planted_lengths) {
        row[[sprintf("planted_covered_bp_%dmb", length_bp / 1e6)]] <-
          sum(segment_covered[planted$length_bp == length_bp])
      }
      if (row$called_bp != row$planted_covered_bp + row$called_low_heterozygosity_bp +
          row$called_no_record_bp + row$called_unsupported_bp) {
        stop("called bases do not sum over their categories", call. = FALSE)
      }
      child_rows[[length(child_rows) + 1L]] <- row
    }
  }
}
child_metrics <- do.call(rbind, child_rows)
child_metrics <- child_metrics[order(child_metrics$min_call_bp, child_metrics$Population,
                                     child_metrics$arm, child_metrics$sample_id), ]
utils::write.csv(child_metrics, file.path(results, "child_metrics.csv"), row.names = FALSE)

ratio <- function(numerator, denominator) {
  if (sum(denominator) == 0) NA_real_ else sum(numerator) / sum(denominator)
}
population_order <- c(admixed, controls)
segments_of_length <- table(truth$length_bp) / length(children)
summary_rows <- list()
for (minimum in min_call_bp) for (p in population_order) for (arm in arms) {
  x <- child_metrics[child_metrics$min_call_bp == minimum & child_metrics$Population == p &
                       child_metrics$arm == arm, , drop = FALSE]
  row <- data.frame(
    min_call_bp = minimum, Population = p, arm = arm, children = nrow(x),
    sensitivity = ratio(x$planted_covered_bp, x$planted_bp), stringsAsFactors = FALSE)
  for (length_bp in planted_lengths) {
    name <- sprintf("%dmb", length_bp / 1e6)
    row[[paste0("sensitivity_", name)]] <- ratio(
      x[[paste0("planted_covered_bp_", name)]],
      rep(length_bp * segments_of_length[[as.character(length_bp)]], nrow(x)))
  }
  row$false_discovery <- ratio(x$called_unsupported_bp, x$called_bp)
  row$outside_planted <- ratio(x$called_bp - x$planted_covered_bp, x$called_bp)
  row$called_bp_per_child <- mean(x$called_bp)
  row$unsupported_bp_per_child <- mean(x$called_unsupported_bp)
  row$low_heterozygosity_bp_per_child <- mean(x$called_low_heterozygosity_bp)
  row$no_record_bp_per_child <- mean(x$called_no_record_bp)
  summary_rows[[length(summary_rows) + 1L]] <- row
}
utils::write.csv(do.call(rbind, summary_rows), file.path(results, "population_arm_summary.csv"),
                 row.names = FALSE)

# Paired bootstrap over the children of one population: the same resampled
# children score both arms, so the interval is of the difference between arms.
set.seed(bootstrap_seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
         sample.kind = "Rejection")
interval <- function(values) {
  if (anyNA(values)) return(c(NA_real_, NA_real_))
  unname(stats::quantile(values, c(0.025, 0.975), type = 7))
}
criterion_rows <- list()
primary <- child_metrics[child_metrics$min_call_bp == primary_min_call_bp, , drop = FALSE]
for (p in population_order) {
  of_arm <- function(arm) {
    x <- primary[primary$Population == p & primary$arm == arm, , drop = FALSE]
    x[match(children[population == p], x$sample_id), , drop = FALSE]
  }
  ancestry <- of_arm("ancestry_tuned")
  draws <- matrix(sample.int(nrow(ancestry), nrow(ancestry) * bootstrap_draws, replace = TRUE),
                  ncol = bootstrap_draws)
  for (comparison in c("pooled", "single_population")) {
    comparator <- of_arm(comparison)
    difference <- function(rows, numerator, denominator) {
      ratio(ancestry[[numerator]][rows], ancestry[[denominator]][rows]) -
        ratio(comparator[[numerator]][rows], comparator[[denominator]][rows])
    }
    all_rows <- seq_len(nrow(ancestry))
    sensitivity_interval <- interval(apply(draws, 2L, difference,
                                           "planted_covered_bp", "planted_bp"))
    discovery_interval <- interval(apply(draws, 2L, difference,
                                         "called_unsupported_bp", "called_bp"))
    no_sensitivity_loss <- sensitivity_interval[[1L]] > -margin
    better <- no_sensitivity_loss & discovery_interval[[2L]] < 0
    not_worse <- no_sensitivity_loss & discovery_interval[[2L]] < margin
    required <- if (p %in% admixed) {
      "better"
    } else if (comparison == "single_population") {
      "not worse"
    } else {
      "none"
    }
    criterion_rows[[length(criterion_rows) + 1L]] <- data.frame(
      Population = p, group = if (p %in% admixed) "admixed" else "control",
      comparison = comparison, children = nrow(ancestry),
      ancestry_sensitivity = ratio(ancestry$planted_covered_bp, ancestry$planted_bp),
      comparator_sensitivity = ratio(comparator$planted_covered_bp, comparator$planted_bp),
      sensitivity_difference = difference(all_rows, "planted_covered_bp", "planted_bp"),
      sensitivity_difference_low = sensitivity_interval[[1L]],
      sensitivity_difference_high = sensitivity_interval[[2L]],
      ancestry_false_discovery = ratio(ancestry$called_unsupported_bp, ancestry$called_bp),
      comparator_false_discovery = ratio(comparator$called_unsupported_bp, comparator$called_bp),
      false_discovery_difference = difference(all_rows, "called_unsupported_bp", "called_bp"),
      false_discovery_difference_low = discovery_interval[[1L]],
      false_discovery_difference_high = discovery_interval[[2L]],
      better = better, not_worse = not_worse, required = required,
      criterion_pass = switch(required, better = better, `not worse` = not_worse, none = NA),
      stringsAsFactors = FALSE)
  }
}
criterion <- do.call(rbind, criterion_rows)
utils::write.csv(criterion, file.path(results, "criterion_by_population.csv"), row.names = FALSE)

truth$Population <- population[truth$sample_id]
truth_rows <- do.call(rbind, lapply(population_order, function(p) {
  x <- truth[truth$Population == p, , drop = FALSE]
  data.frame(Population = p, children = length(unique(x$sample_id)), segments = nrow(x),
             planted_bp = sum(x$length_bp), records = sum(x$records),
             source_heterozygous = sum(x$source_heterozygous))
}))
utils::write.csv(truth_rows, file.path(results, "planted_truth_by_population.csv"),
                 row.names = FALSE)
utils::write.csv(truth[c("sample_id", "Population", "segment", "chrom", "start", "end",
                         "length_bp", "records", "source_heterozygous")],
                 file.path(results, "planted_truth.csv"), row.names = FALSE)

required <- criterion[criterion$required != "none", , drop = FALSE]
verdict <- if (any(required$criterion_pass %in% FALSE)) {
  "not held"
} else if (anyNA(required$criterion_pass)) {
  "not evaluable in every population"
} else {
  "held"
}
cat("synthetic-truth criterion:", verdict, "\n")
print(criterion[c("Population", "comparison", "sensitivity_difference",
                  "false_discovery_difference", "required", "criterion_pass")])
