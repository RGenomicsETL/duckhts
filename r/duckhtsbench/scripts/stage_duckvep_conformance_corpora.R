#!/usr/bin/env Rscript
# Stage the registry-declared HPRC, Sniffles2, and dbVar DuckVEP corpora.

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("optparse is required to stage DuckVEP conformance corpora", call. = FALSE)
}
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("must be invoked as an R script", call. = FALSE)
package_dir <- normalizePath(
  file.path(dirname(sub("^--file=", "", script_arg[[1L]])), ".."), mustWork = TRUE
)
if (!nzchar(Sys.getenv("DUCKHTSBENCH_REGISTRY", unset = ""))) {
  Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(package_dir, "inst", "benchmark_registry.tsv"))
}
for (file in c("registry.R", "stage.R", "duckvep_corpora.R")) {
  source(file.path(package_dir, "R", file))
}

parser <- optparse::OptionParser(
  usage = "%prog [options]",
  option_list = list(
    optparse::make_option(
      "--corpus", default = "all",
      help = paste(
        "all, hprc-african4-chr22, sniffles2-chr22, or dbvar-chr22",
        "[default: %default]"
      )
    ),
    optparse::make_option(
      "--output-dir", dest = "output_dir", type = "character", default = NULL,
      help = "compatibility output root; default paths come from the artifact registry"
    ),
    optparse::make_option(
      "--bcftools", type = "character",
      default = Sys.getenv("BCFTOOLS_BIN", unset = Sys.which("bcftools"))
    ),
    optparse::make_option(
      "--tabix", type = "character",
      default = Sys.getenv("TABIX_BIN", unset = Sys.which("tabix"))
    ),
    optparse::make_option(
      "--curl", type = "character",
      default = Sys.getenv("CURL_BIN", unset = Sys.which("curl"))
    ),
    optparse::make_option(
      "--plan", action = "store_true", default = FALSE,
      help = "print the selected registry closure without staging"
    )
  )
)
options <- optparse::parse_args(parser)
definitions <- duckhts_bench_duckvep_corpus_definitions()
if (!options$corpus %in% c("all", names(definitions))) {
  stop("unknown DuckVEP corpus: ", options$corpus, call. = FALSE)
}
selected <- if (options$corpus == "all") definitions else definitions[options$corpus]
if (isTRUE(options$plan)) {
  ids <- unlist(lapply(selected, function(definition) {
    unlist(definition[c("source", "source_index", "manifest", "output", "output_index")],
      use.names = FALSE
    )
  }), use.names = FALSE)
  plan <- duckhts_bench_stage_plan("duckvep-conformance-corpora")
  print(plan[match(ids, plan$id), , drop = FALSE], row.names = FALSE)
  quit(status = 0L)
}

paths <- duckhts_bench_stage_duckvep_corpora(
  corpus = options$corpus,
  output_dir = options$output_dir,
  bcftools = options$bcftools,
  tabix = options$tabix,
  curl = options$curl
)
keys <- toupper(gsub("[^A-Za-z0-9]+", "_", names(paths)))
cat(paste0("DUCKVEP_CORPUS_", keys, "=", unname(paths)), sep = "\n")
cat("\n")
