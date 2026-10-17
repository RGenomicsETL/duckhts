suppressPackageStartupMessages(library(tinytest))
setwd(tinytest::get_call_wd())
source("benchmarks/count_mixture_lowpass_genome/batch.R")
set.seed(332L)
selected <- list(sample = c("a", "b", "c"), other = c("d", "e"))
panels <- lapply(1:2, function(chrom) {
  n <- 10000L
  list(
    selected = selected, settings = genome_settings,
    sites = data.frame(chrom = as.character(chrom), pos = seq_len(n),
      f = runif(n, .01, .99), g = runif(n, .01, .99)),
    genotypes = setNames(lapply(unlist(selected), function(id) as.raw(sample(0:3, n, replace = TRUE))),
      unlist(selected))
  )
})
strata <- genome_strata()
for (relationship in genome_declaration$relationships) {
  for (depth in genome_declaration$depths) {
    for (error in genome_declaration$errors) {
      setting <- strata[strata$relationship == relationship & strata$depth == depth &
        strata$error == error & strata$share == .1, ]
      full <- genome_counts(panels, setting, 1L, covered_only = FALSE)
      sparse <- genome_counts(panels, setting, 1L, covered_only = TRUE)
      n <- full$counts$ref_count + full$counts$alt_count + full$counts$other_count
      expected <- full$counts[n > 0, ]
      rownames(expected) <- NULL
      expect_identical(sparse$counts, expected)
      expect_identical(sparse$qc, full$qc)
      expect_identical(sparse$identity, full$identity)
      expect_identical(sparse$qc[["sites"]], 20000)
      expect_identical(sparse$qc[["zero_depth_sites"]], as.numeric(sum(n == 0)))
      expect_identical(sparse$qc[["excluded_sites"]], as.numeric(sum(n > genome_settings$max_depth)))
    }
  }
}
