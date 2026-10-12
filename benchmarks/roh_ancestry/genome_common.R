# Shared settings and cache identity of the genome-wide ROH ancestry drivers.
# Each driver runs from the repository root and sources this file.
genome_cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-ancestry")
genome_results <- "benchmarks/results/roh-ancestry"
genome_children_list <- file.path(genome_cache, "chr20.children.present.txt")

# The autosomes of a run: ROH_GENOME_CHROMOSOMES, autosome numbers separated
# by commas, or all 22 when it is unset. A run on fewer autosomes is a smoke
# run whose outputs cover only those autosomes. Every marker records the set,
# so a run on another set never reuses those outputs and the summaries refuse
# them.
genome_chromosomes <- function() {
  value <- Sys.getenv("ROH_GENOME_CHROMOSOMES", unset = "")
  if (!nzchar(value)) return(1:22)
  parts <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  chromosomes <- suppressWarnings(as.integer(parts))
  if (!length(parts) || anyNA(chromosomes) || any(as.character(chromosomes) != parts) ||
      any(chromosomes < 1L | chromosomes > 22L) || anyDuplicated(chromosomes)) {
    stop("ROH_GENOME_CHROMOSOMES must list distinct autosome numbers 1-22, ",
         "separated by commas", call. = FALSE)
  }
  sort(chromosomes)
}
genome_set_line <- function(chromosomes) paste0("chromosomes=", paste(chromosomes, collapse = ","))

# GRCh38 autosome lengths: lines 1-22 of the Ensembl 116 primary-assembly FASTA
# index (Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai, SHA-256
# 0998f61682f4041b11f0d156e1db6dae3e4c743e26643a3f45ea7faea70cb604), committed
# so the truth windows and the FROH denominator do not depend on an
# unregistered local file.
genome_fai <- "benchmarks/roh_ancestry/grch38_autosomes.fai"
genome_lengths <- function(chromosomes) {
  fai <- utils::read.delim(genome_fai, header = FALSE, stringsAsFactors = FALSE,
                           colClasses = c("character", "integer", "NULL", "NULL", "NULL"))
  lengths <- fai$V2[match(as.character(chromosomes), fai$V1)]
  if (anyNA(lengths)) stop("GRCh38 FAI lacks an autosome", call. = FALSE)
  data.frame(chromosome = chromosomes, length_bp = lengths)
}

# Registry paths of the per-chromosome inputs.
genome_artifact <- function(format, chromosome) {
  duckhtsbench::duckhts_bench_artifact_path(sprintf(format, chromosome))
}
genome_children_bcf <- function(chromosome) {
  genome_artifact("roh_ancestry_chr%d_children_bcf", chromosome)
}
genome_af_sites <- function(chromosome) genome_artifact("roh_ancestry_chr%d_af_sites", chromosome)
genome_reference_long <- function(chromosome) {
  genome_artifact("roh_ancestry_chr%d_reference_long", chromosome)
}

# Cache identity. An identity is one name=SHA-256 line per input file and one
# name=value line per setting. A marker holds the identity of an output and
# then output=<SHA-256 of the output>; a cached output is reused only when its
# marker holds exactly the current identity and the output's current SHA-256.
genome_sha256 <- function(path) {
  output <- system2("sha256sum", shQuote(path), stdout = TRUE)
  if (!is.null(attr(output, "status")) || length(output) != 1L) {
    stop("sha256sum failed: ", path, call. = FALSE)
  }
  sub(" .*$", "", output)
}
genome_identity <- function(files, values = character()) {
  missing <- files[!nzchar(files) | !file.exists(files)]
  if (length(missing)) {
    stop("missing inputs: ", paste(names(missing), missing, collapse = ", "), call. = FALSE)
  }
  c(paste0(names(files), "=", vapply(files, genome_sha256, character(1L))),
    paste0(names(values), "=", values))
}
genome_marker <- function(output) paste0(output, ".identity")
genome_current <- function(output, identity, marker = genome_marker(output)) {
  file.exists(output) && file.exists(marker) &&
    identical(readLines(marker), c(identity, paste0("output=", genome_sha256(output))))
}
genome_record <- function(output, identity, marker = genome_marker(output)) {
  writeLines(c(identity, paste0("output=", genome_sha256(output))), marker)
}
# A summary reads an output only when its marker records the run's
# chromosome set and the output's current SHA-256.
genome_check_output <- function(output, chromosomes, marker = genome_marker(output)) {
  lines <- if (file.exists(output) && file.exists(marker)) readLines(marker) else character()
  if (!identical(lines[startsWith(lines, "chromosomes=")], genome_set_line(chromosomes)) ||
      !identical(lines[startsWith(lines, "output=")],
                 paste0("output=", genome_sha256(output)))) {
    stop("output was not written for this chromosome set or changed since: ", output,
         call. = FALSE)
  }
}
