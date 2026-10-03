cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-ancestry")
results_dir <- "benchmarks/results/roh-ancestry"
q <- utils::read.csv(file.path(cache, "autosome.q.csv"), stringsAsFactors = FALSE)
pedigree <- read.table(file.path(cache, "pedigree.txt"), header = TRUE,
                       stringsAsFactors = FALSE)
children <- readLines(file.path(cache, "chr20.children.present.txt"))
metadata <- pedigree[pedigree$SampleID %in% children,
                     c("SampleID", "Population", "Superpopulation"), drop = FALSE]
if (nrow(metadata) != 377L || anyDuplicated(metadata$SampleID) ||
    nrow(q) != length(children) * 21L || anyDuplicated(q[c("sample_id", "group_id")]) ||
    length(unique(q$group_id)) != 21L || any(q$status != "ok") ||
    any(q$cor_pred < 0.4) || !setequal(unique(q$sample_id), children)) {
  stop("genome-wide q output failed child/group/status validation", call. = FALSE)
}
per_child_sum <- aggregate(proportion ~ sample_id, q, sum)
if (any(abs(per_child_sum$proportion - 1) > 1e-6)) {
  stop("genome-wide ancestry proportions do not sum to one", call. = FALSE)
}
diagnostics_consistent <- vapply(split(q, q$sample_id), function(sample_data) {
  length(unique(sample_data$used_variants)) == 1L &&
    length(unique(sample_data$cor_pred)) == 1L
}, logical(1L))
if (!all(diagnostics_consistent)) {
  stop("q diagnostics vary across ancestry groups within a child", call. = FALSE)
}
q <- merge(q, metadata, by.x = "sample_id", by.y = "SampleID", sort = FALSE)
q <- q[order(q$Population, q$group_id, q$sample_id), ]
utils::write.csv(q, file.path(results_dir, "autosome_q_by_child.csv"), row.names = FALSE)
q_means <- aggregate(proportion ~ Population + Superpopulation + group_id, q, mean)
q_means <- q_means[order(q_means$Population, q_means$group_id), ]
utils::write.csv(q_means,
  file.path(results_dir, "autosome_q_summary_by_population.csv"), row.names = FALSE)
quality <- do.call(rbind, lapply(sort(unique(q$Population)), function(population) {
  sample_data <- q[q$Population == population & q$group_id == unique(q$group_id)[[1L]], ]
  data.frame(Population = population, Superpopulation = unique(sample_data$Superpopulation),
    children = length(unique(sample_data$sample_id)),
    used_variants_min = min(sample_data$used_variants),
    used_variants_max = max(sample_data$used_variants),
    cor_pred_min = min(sample_data$cor_pred),
    cor_pred_median = median(sample_data$cor_pred),
    cor_pred_max = max(sample_data$cor_pred), stringsAsFactors = FALSE)
}))
utils::write.csv(quality,
  file.path(results_dir, "autosome_q_quality_by_population.csv"), row.names = FALSE)
site_counts <- unique(q[c("sample_id", "Population", "used_variants", "input_variants",
                           "dropped_variants", "reversed_variants", "flipped_variants",
                           "ambiguous_variants", "unmatched_variants", "invalid_variants",
                           "missing_variants", "cor_pred", "status")])
if (anyDuplicated(site_counts$sample_id)) stop("q site counts are not per-child", call. = FALSE)
utils::write.csv(site_counts,
  file.path(results_dir, "autosome_q_site_counts_by_child.csv"), row.names = FALSE)
cat("genome q", length(unique(q$sample_id)), "children;",
    length(unique(q$group_id)), "groups;", nrow(q), "proportions; cor_pred",
    range(site_counts$cor_pred), "; used variants", range(site_counts$used_variants), "\n")
