# Genome-wide accuracy of the ROH arms against the truth of truth_genome.R,
# over the autosomes of the run, and the acceptance criterion per population.
# The truth marks runs of low heterozygosity, not autozygosity, so truth
# length and FROH compare arms within a population, not across populations.
# froh_assembly is called ROH length over the summed FAI length of the run's
# autosomes, gaps and centromeres included: no callable-base mask is
# registered, so it is a whole-assembly fraction and smaller than a FROH over
# callable bases. Every input must have been written for the run's
# chromosome set, as its marker records.
source("benchmarks/roh_ancestry/genome_common.R")
results_dir <- genome_results
chromosomes <- genome_chromosomes()
pedigree <- read.table(duckhtsbench::duckhts_bench_artifact_path("roh_ancestry_pedigree"),
                       header = TRUE, stringsAsFactors = FALSE)
children <- readLines(genome_children_list)
population_children <- pedigree[
  pedigree$FatherID != "0" & pedigree$MotherID != "0" &
    pedigree$SampleID %in% children, c("SampleID", "Population"), drop = FALSE]
if (nrow(population_children) != 377L || anyDuplicated(population_children$SampleID)) {
  stop("pedigree does not uniquely identify all 377 children", call. = FALSE)
}
# truth_genome.R checked these intervals' chr20 part against the amended chr20
# truth when it wrote them; the marker check confirms they are that file.
truth_path <- file.path(genome_cache, "autosome.truth_intervals.csv")
genome_check_output(truth_path, chromosomes)
truth <- read.csv(truth_path, stringsAsFactors = FALSE)
if (anyDuplicated(truth[c("chromosome", "sample_id", "start", "end")])) {
  stop("truth intervals contain duplicate keys", call. = FALSE)
}
thresholds <- c(1000000L, 2000000L, 5000000L)
mode_for_population <- c(ACB = "single_AFR", ASW = "single_AFR", ESN = "single_AFR",
                         YRI = "single_AFR", CLM = "single_AMR", MXL = "single_AMR",
                         PEL = "single_AMR", PUR = "single_AMR", CEU = "single_EUR",
                         CHS = "single_EAS")
arm_names <- c("pooled", "single_AFR", "single_AMR", "single_EUR", "single_EAS",
               "ancestry")
canonical_segments <- function(x) {
  x <- x[order(x$sample, x$chrom, x$start, x$end), , drop = FALSE]
  x[c("sample", "chrom", "start", "end", "length", "n_markers", "quality")]
}
eligible_ids <- list(
  pooled = children,
  ancestry = children,
  single_AFR = population_children$SampleID[population_children$Population %in%
    c("ACB", "ASW", "ESN", "YRI")],
  single_AMR = population_children$SampleID[population_children$Population %in%
    c("CLM", "MXL", "PEL", "PUR")],
  single_EUR = population_children$SampleID[population_children$Population == "CEU"],
  single_EAS = population_children$SampleID[population_children$Population == "CHS"]
)
# Every arm decodes every autosome of the run, chr20 included.
read_repetition <- function(arm, repetition) {
  stem <- file.path(genome_cache, sprintf("autosome.arm_%s_%d", arm, repetition))
  genome_check_output(paste0(stem, ".csv"), chromosomes, paste0(stem, ".ok"))
  current <- read.csv(paste0(stem, ".csv"), stringsAsFactors = FALSE)
  observed <- current$sample
  if (anyNA(observed) || any(!observed %in% eligible_ids[[arm]])) {
    stop("arm output contains an ineligible child for ", arm, call. = FALSE)
  }
  if (any(!sub("^chr", "", current$chrom) %in% as.character(chromosomes))) {
    stop("arm output contains a chromosome outside the run for ", arm, call. = FALSE)
  }
  current
}
# The accuracy arms are deterministic, so one completed repetition suffices; any
# further completed repetitions must reproduce the first exactly.
completed_repetitions <- function(arm) {
  repetitions <- 1:3
  repetitions[file.exists(file.path(genome_cache,
    sprintf("autosome.arm_%s_%d.ok", arm, repetitions)))]
}
arm_segments <- setNames(vector("list", length(arm_names)), arm_names)
for (arm in arm_names) {
  repetitions <- completed_repetitions(arm)
  if (!length(repetitions) || repetitions[[1L]] != 1L) {
    stop("arm ", arm, " has no completed first repetition", call. = FALSE)
  }
  reference <- canonical_segments(read_repetition(arm, 1L))
  for (repetition in repetitions[-1L]) {
    current <- canonical_segments(read_repetition(arm, repetition))
    if (!isTRUE(all.equal(reference, current, tolerance = 1e-12,
                          check.attributes = FALSE))) {
      stop("repeated genome-wide arm output differs for ", arm, call. = FALSE)
    }
    rm(current)
    gc()
  }
  arm_segments[[arm]] <- reference
  rm(reference)
  gc()
}
runtime <- do.call(rbind, lapply(arm_names, function(arm) {
  data.frame(arm = arm, repetition = completed_repetitions(arm), stringsAsFactors = FALSE)
}))
runtime$time_file <- file.path(genome_cache,
  sprintf("autosome.arm_%s_%d.time.txt", runtime$arm, runtime$repetition))
if (any(!file.exists(runtime$time_file))) {
  stop("missing a completed autosome arm timing receipt", call. = FALSE)
}
time_values <- lapply(runtime$time_file, scan, what = double(), sep = "\t",
                       quiet = TRUE)
if (any(lengths(time_values) != 2L)) {
  stop("autosome arm timing receipt must contain elapsed seconds and peak RSS",
       call. = FALSE)
}
runtime$elapsed_seconds <- vapply(time_values, `[[`, numeric(1L), 1L)
runtime$peak_rss_kib <- vapply(time_values, `[[`, numeric(1L), 2L)
runtime$peak_rss_gib <- runtime$peak_rss_kib / 1024^2
runtime <- runtime[c("arm", "repetition", "elapsed_seconds", "peak_rss_kib",
                     "peak_rss_gib")]
write.csv(runtime, file.path(results_dir, "autosome_arm_runtime.csv"),
          row.names = FALSE)
segments <- list(
  pooled = arm_segments$pooled,
  ancestry_tuned = arm_segments$ancestry)
single_parts <- lapply(names(mode_for_population), function(population) {
  source <- mode_for_population[[population]]
  data <- arm_segments[[source]]
  ids <- population_children$SampleID[population_children$Population == population]
  data[data$sample %in% ids, , drop = FALSE]
})
names(single_parts) <- names(mode_for_population)
segments$single_population <- do.call(rbind, single_parts)
segments <- lapply(segments, function(x) {
  x <- x[!is.na(x$length), , drop = FALSE]
  x$chromosome <- as.integer(sub("^chr", "", as.character(x$chrom)))
  x[order(x$sample, x$chromosome, x$start), , drop = FALSE]
})
for (arm in names(segments)) {
  x <- segments[[arm]]
  chromosome_groups <- interaction(x$sample, x$chromosome, drop = TRUE,
                                   lex.order = TRUE)
  for (one in split(x, chromosome_groups)) {
    if (nrow(one) > 1L && any(one$start[-1L] <= one$end[-nrow(one)])) {
      stop("overlapping ROH calls prevent unique interval coverage", call. = FALSE)
    }
  }
}
assembly_bp <- sum(genome_lengths(chromosomes)$length_bp)
truth_by_child <- split(truth, truth$sample_id)
child_metrics <- do.call(rbind, lapply(names(segments), function(arm) {
  data <- segments[[arm]]
  calls_by_child <- split(data, data$sample)
  do.call(rbind, lapply(children, function(id) {
    calls <- calls_by_child[[id]]
    if (is.null(calls)) calls <- data[FALSE, , drop = FALSE]
    truth_calls <- truth_by_child[[id]]
    if (is.null(truth_calls)) truth_calls <- truth[FALSE, , drop = FALSE]
    call_indices_by_chromosome <- split(seq_len(nrow(calls)),
                                         as.character(calls$chromosome))
    truth_by_chromosome <- split(truth_calls, truth_calls$chromosome)
    support_bp <- numeric(nrow(calls))
    truth_covered_bp <- 0
    shared_chromosomes <- intersect(names(call_indices_by_chromosome),
                                    names(truth_by_chromosome))
    for (chromosome in shared_chromosomes) {
      call_indices <- call_indices_by_chromosome[[chromosome]]
      chromosome_calls <- calls[call_indices, , drop = FALSE]
      chromosome_truth <- truth_by_chromosome[[chromosome]]
      overlap_bp <- matrix(pmax(0,
        outer(chromosome_calls$end, chromosome_truth$end, pmin) -
          outer(chromosome_calls$start, chromosome_truth$start, pmax) + 1),
        nrow = nrow(chromosome_calls), ncol = nrow(chromosome_truth))
      support_bp[call_indices] <- rowSums(overlap_bp)
      truth_covered_bp <- truth_covered_bp + sum(overlap_bp)
    }
    supported <- support_bp / calls$length >= 0.9
    row <- data.frame(sample_id = id, arm = arm,
      total_roh_bp = sum(calls$length), froh_assembly = sum(calls$length) / assembly_bp,
      supported_count_1mb = sum(calls$length >= thresholds[[1L]] & supported),
      supported_length_1mb = sum(calls$length[calls$length >= thresholds[[1L]] & supported]),
      truth_bp = sum(truth_calls$end - truth_calls$start + 1),
      truth_covered_bp = truth_covered_bp)
    for (index in seq_along(thresholds)) {
      suffix <- c("1mb", "2mb", "5mb")[[index]]
      selected <- calls$length >= thresholds[[index]]
      row[[paste0("unsupported_count_", suffix)]] <- sum(selected & !supported)
      row[[paste0("unsupported_length_", suffix)]] <- sum(calls$length[selected & !supported])
    }
    row
  }))
}))
child_metrics <- merge(child_metrics, population_children,
                        by.x = "sample_id", by.y = "SampleID")
arm_order <- c("pooled", "single_population", "ancestry_tuned")
summary_rows <- do.call(rbind, lapply(arm_order, function(arm) {
  do.call(rbind, lapply(sort(unique(population_children$Population)), function(population) {
    x <- child_metrics[child_metrics$arm == arm &
      child_metrics$Population == population, , drop = FALSE]
    unsupported_count <- vapply(c("1mb", "2mb", "5mb"), function(suffix) {
      mean(x[[paste0("unsupported_count_", suffix)]])
    }, numeric(1L))
    unsupported_length <- vapply(c("1mb", "2mb", "5mb"), function(suffix) {
      mean(x[[paste0("unsupported_length_", suffix)]])
    }, numeric(1L))
    data.frame(Population = population, Arm = arm,
      `Unsupported ROH count / child (>=1 Mb; >=2 Mb; >=5 Mb)` =
        paste(formatC(unsupported_count, digits = 3, format = "f"), collapse = "; "),
      `Unsupported ROH length / child (>=1 Mb; >=2 Mb; >=5 Mb)` =
        paste(formatC(unsupported_length, digits = 0, format = "f"), collapse = "; "),
      `Supported count / child` = formatC(mean(x$supported_count_1mb), digits = 3, format = "f"),
      `Supported length / child` = formatC(mean(x$supported_length_1mb), digits = 0, format = "f"),
      `FROH (whole assembly)` = formatC(mean(x$froh_assembly), digits = 6, format = "f"),
      `Truth-run length covered` = if (sum(x$truth_bp) > 0) {
        formatC(sum(x$truth_covered_bp) / sum(x$truth_bp) * 100,
                digits = 2, format = "f")
      } else "NA",
      check.names = FALSE, stringsAsFactors = FALSE)
  }))
}))
write.csv(child_metrics[order(child_metrics$Population, child_metrics$arm,
  child_metrics$sample_id), ],
  file.path(results_dir, "autosome_child_metrics.csv"), row.names = FALSE)
write.csv(summary_rows[order(summary_rows$Population, summary_rows$Arm), ],
  file.path(results_dir, "autosome_population_arm_summary.csv"),
  row.names = FALSE, quote = TRUE)
admixed <- c("ACB", "ASW", "CLM", "MXL", "PEL", "PUR")
controls <- c("YRI", "ESN", "CEU", "CHS")
# As in criterion_chr20.R: a population without truth length has no truth
# coverage, so its coverage check is NA, not evaluable, and never counts as a
# pass. A failed count check still fails.
within_five_percent <- function(value, reference) {
  if (is.na(value) || is.na(reference)) return(NA)
  if (reference == 0) return(value == 0)
  abs(value - reference) <= 0.05 * abs(reference)
}
no_material_coverage_loss <- function(value, reference) {
  if (is.na(value) || is.na(reference)) return(NA)
  if (reference == 0) return(value == 0)
  value >= 0.95 * reference
}
populations <- sort(unique(child_metrics$Population))
criterion <- do.call(rbind, lapply(populations, function(population) {
  x <- child_metrics[child_metrics$Population == population, , drop = FALSE]
  average <- function(arm, field) mean(x[[field]][x$arm == arm])
  coverage <- function(arm) {
    selected <- x[x$arm == arm, , drop = FALSE]
    if (sum(selected$truth_bp) == 0) return(NA_real_)
    sum(selected$truth_covered_bp) / sum(selected$truth_bp)
  }
  ancestry_count <- average("ancestry_tuned", "unsupported_count_1mb")
  ancestry_coverage <- coverage("ancestry_tuned")
  if (population %in% admixed) {
    comparison_arms <- c("pooled", "single_population")
    count_pass <- vapply(comparison_arms, function(arm) {
      ancestry_count < average(arm, "unsupported_count_1mb")
    }, logical(1L))
    coverage_pass <- vapply(comparison_arms, function(arm) {
      no_material_coverage_loss(ancestry_coverage, coverage(arm))
    }, logical(1L))
    data.frame(Population = population, comparison = comparison_arms,
      ancestry_unsupported_count = ancestry_count,
      comparator_unsupported_count = vapply(comparison_arms, function(arm) {
        average(arm, "unsupported_count_1mb")
      }, numeric(1L)), ancestry_truth_coverage = ancestry_coverage,
      comparator_truth_coverage = vapply(comparison_arms, coverage, numeric(1L)),
      count_pass = count_pass, truth_coverage_pass = coverage_pass,
      criterion_pass = count_pass & coverage_pass)
  } else if (population %in% controls) {
    comparison <- "single_population"
    comparison_count <- average(comparison, "unsupported_count_1mb")
    comparison_coverage <- coverage(comparison)
    count_pass <- within_five_percent(ancestry_count, comparison_count)
    coverage_pass <- within_five_percent(ancestry_coverage, comparison_coverage)
    data.frame(Population = population, comparison = comparison,
      ancestry_unsupported_count = ancestry_count,
      comparator_unsupported_count = comparison_count,
      ancestry_truth_coverage = ancestry_coverage,
      comparator_truth_coverage = comparison_coverage,
      count_pass = count_pass, truth_coverage_pass = coverage_pass,
      criterion_pass = count_pass & coverage_pass)
  } else stop("unexpected population", call. = FALSE)
}))
write.csv(criterion, file.path(results_dir, "autosome_criterion_by_population.csv"),
          row.names = FALSE)
verdict <- if (any(criterion$criterion_pass %in% FALSE)) {
  "not held"
} else if (anyNA(criterion$criterion_pass)) {
  "not evaluable in every population"
} else {
  "held"
}
cat("children", nrow(child_metrics) / length(arm_order), "chromosomes", length(chromosomes),
    "truth runs", nrow(truth), "criterion", verdict, "\n")
print(criterion)

# The receipt of the genome-wide results: the chromosome set, the extension
# and the cache outputs the summaries read (each checked against its marker
# above or in summarise_genome_q.R), then the SHA-256 of every result table
# that truth_genome.R, verify_genome_site_sets.R, summarise_genome_q.R and this
# script wrote. The report reads a table only when its SHA-256 matches here.
arm_marker <- readLines(file.path(genome_cache, "autosome.arm_pooled_1.ok"))
q_path <- file.path(genome_cache, "autosome.q.csv")
genome_check_output(q_path, chromosomes)
cache_outputs <- c(truth_intervals = truth_path, q = q_path, stats::setNames(
  file.path(genome_cache, sprintf("autosome.arm_%s_1.csv", arm_names)),
  paste0("arm_", arm_names)))
result_tables <- c(
  "autosome_truth_window_distribution.csv", "autosome_truth_window_evidence.csv",
  "autosome_truth_by_child.csv", "autosome_truth_by_population.csv",
  "autosome_site_set_equivalence.csv", "autosome_q_by_child.csv",
  "autosome_q_summary_by_population.csv", "autosome_q_quality_by_population.csv",
  "autosome_q_site_counts_by_child.csv", "autosome_arm_runtime.csv",
  "autosome_child_metrics.csv", "autosome_population_arm_summary.csv",
  "autosome_criterion_by_population.csv")
receipt <- data.frame(
  field = c("chromosomes", "extension_sha256", "duckdb_r", "criterion",
            paste0("cache_", names(cache_outputs)), result_tables),
  value = c(paste(chromosomes, collapse = ","),
            sub("^extension=", "", arm_marker[startsWith(arm_marker, "extension=")]),
            sub("^duckdb_r=", "", arm_marker[startsWith(arm_marker, "duckdb_r=")]),
            verdict, vapply(cache_outputs, genome_sha256, character(1L)),
            vapply(file.path(results_dir, result_tables), genome_sha256, character(1L))),
  stringsAsFactors = FALSE)
utils::write.table(receipt, file.path(results_dir, "autosome_receipt.tsv"), sep = "\t",
                   row.names = FALSE, quote = FALSE)
