#' Stage the RefSeq E. coli K-12 MG1655 GenBank Benchmark Input
#'
#' Downloads the pinned NCBI RefSeq assembly archive
#' `GCF_000005845.2_ASM584v2_genomic.gbff.gz`, verifies it against the
#' registry identity (NCBI's published MD5 plus SHA-256 and byte size), and
#' derives the uncompressed `.gbff` that `benchmark_genbank_reader.Rmd` reads.
#' Network access occurs only in this explicit staging step; with
#' `fetch = FALSE` the compressed source must already be cached, which is how
#' the report renders.
#' @param fetch Whether to download a missing or invalid compressed source.
#' @return Named cache paths `gbff_gz` (the archive) and `gbff` (the input).
#' @export
duckhts_bench_stage_genbank <- function(fetch = TRUE) {
  ids <- c("genbank_ecoli_k12_gbff_gz", "genbank_ecoli_k12_gbff")
  plan <- duckhts_bench_stage_plan("genbank-reader")
  if (!identical(plan$id, ids)) stop("genbank-reader registry plan is incomplete", call. = FALSE)
  if (!identical(plan$transform, c("direct_download", "gunzip")) ||
      plan$locator[[2L]] != paste0("artifact:", ids[[1L]])) {
    stop("genbank-reader registry rows must be a direct download and its gunzip derivation", call. = FALSE)
  }
  paths <- stats::setNames(vapply(ids, duckhts_bench_artifact_path, character(1L)), c("gbff_gz", "gbff"))

  if (fetch) {
    duckhts_bench_fetch(ids[[1L]])
  } else {
    if (!file.exists(paths[["gbff_gz"]])) {
      stop("genbank source is not staged; run duckhts_bench_stage_genbank(): ", paths[["gbff_gz"]], call. = FALSE)
    }
    duckhts_bench_validate_identity(ids[[1L]], paths[["gbff_gz"]])
  }

  duckhts_bench_stage_gunzip(ids[[2L]], paths[["gbff_gz"]], paths[["gbff"]])
  invisible(paths)
}

#' Stage a Derived Artifact That Joins Gunzipped Record Files
#'
#' Decompresses each of `sources`, in order, into one `destination`, so a
#' set of gzipped record files whose records each end in `//` becomes one
#' record file. This is the registry's `gunzip_concatenate` transform, whose
#' locator lists the source artifacts in order; the growth dimension for
#' record-count scaling is then how many verified source parts a derived
#' artifact joins. A cached destination that matches its registered identity
#' is reused; otherwise it is rebuilt through a temporary file and validated
#' before publication, as [duckhts_bench_stage_gunzip()] does.
#' @param id Registry artifact identifier of the destination.
#' @param sources Staged, verified gzipped source paths, in the locator's order.
#' @param destination Cache path to publish.
#' @return The destination path, invisibly.
#' @export
duckhts_bench_stage_gunzip_concatenate <- function(id, sources, destination) {
  registry <- duckhts_bench_registry()
  row <- registry[registry$id == id, , drop = FALSE]
  if (nrow(row) != 1L) stop("unknown or non-unique benchmark artifact: ", id, call. = FALSE)
  if (row$transform != "gunzip_concatenate") {
    stop("registry transform for ", id, " must be gunzip_concatenate", call. = FALSE)
  }
  declared <- strsplit(row$locator, ";", fixed = TRUE)[[1L]]
  if (length(declared) != length(sources)) {
    stop("registry locator for ", id, " names ", length(declared), " sources; ",
         length(sources), " supplied", call. = FALSE)
  }
  if (file.exists(destination)) {
    valid <- tryCatch({
      duckhts_bench_validate_identity(id, destination)
      TRUE
    }, error = function(error) FALSE)
    if (valid) {
      duckhts_bench_write_provenance(id, destination)
      return(invisible(destination))
    }
    unlink(c(destination, paste0(destination, ".provenance.tsv")), force = TRUE)
  }
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(destination, ".partial-", Sys.getpid())
  unlink(temporary, force = TRUE)
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  output <- file(temporary, open = "wb")
  for (source in sources) {
    input <- gzfile(source, open = "rb")
    repeat {
      chunk <- readBin(input, what = "raw", n = 1048576L)
      if (!length(chunk)) break
      writeBin(chunk, output)
    }
    close(input)
  }
  close(output)
  duckhts_bench_validate_identity(id, temporary)
  if (!file.rename(temporary, destination)) {
    stop("could not publish the concatenated artifact: ", destination, call. = FALSE)
  }
  duckhts_bench_write_provenance(id, destination)
  invisible(destination)
}

#' Stage the RefSeq Plasmid Release Parts for Record-Count Scaling
#'
#' Downloads the first four `plasmid.N.genomic.gbff.gz` parts of the pinned
#' NCBI RefSeq release, each verified against NCBI's published MD5 for that
#' release, and derives three record files from them: part 1 alone, parts 1-2
#' joined, and parts 1-4 joined. The three grow the record count while every
#' other dimension (the records themselves) stays fixed, which is the scaling
#' input `benchmark_genbank_named_attributes.Rmd` reads. Network access occurs
#' only in this explicit staging step; with `fetch = FALSE` every part must
#' already be cached, which is how the report renders.
#' @param fetch Whether to download a missing or invalid part.
#' @return Named cache paths `parts` (the four archives) and `records_1`,
#'   `records_2`, `records_4` (the derived record files), invisibly.
#' @export
duckhts_bench_stage_genbank_plasmid <- function(fetch = TRUE) {
  part_ids <- paste0("genbank_plasmid_part", 1:4, "_gbff_gz")
  derived_ids <- c("genbank_plasmid_records_1", "genbank_plasmid_records_2", "genbank_plasmid_records_4")
  plan <- duckhts_bench_stage_plan("genbank-plasmid")
  if (!identical(plan$id, c(part_ids, derived_ids))) {
    stop("genbank-plasmid registry plan is incomplete", call. = FALSE)
  }
  if (!identical(plan$transform, c(rep("direct_download", 4L), rep("gunzip_concatenate", 3L)))) {
    stop("genbank-plasmid registry rows must be four direct downloads and three gunzip concatenations",
         call. = FALSE)
  }
  parts <- vapply(part_ids, duckhts_bench_artifact_path, character(1L))
  for (k in seq_along(part_ids)) {
    if (fetch) {
      duckhts_bench_fetch(part_ids[[k]])
    } else {
      if (!file.exists(parts[[k]])) {
        stop("genbank-plasmid part is not staged; run duckhts_bench_stage_genbank_plasmid(): ",
             parts[[k]], call. = FALSE)
      }
      duckhts_bench_validate_identity(part_ids[[k]], parts[[k]])
    }
  }
  derived <- vapply(derived_ids, duckhts_bench_artifact_path, character(1L))
  counts <- c(1L, 2L, 4L)
  for (k in seq_along(derived_ids)) {
    sources <- strsplit(plan$locator[[4L + k]], ";", fixed = TRUE)[[1L]]
    expected <- paste0("artifact:", part_ids[seq_len(counts[[k]])])
    if (!identical(sources, expected)) {
      stop("registry locator for ", derived_ids[[k]], " must join parts 1..", counts[[k]], call. = FALSE)
    }
    duckhts_bench_stage_gunzip_concatenate(derived_ids[[k]], unname(parts[seq_len(counts[[k]])]),
                                           derived[[k]])
  }
  invisible(list(parts = unname(parts), records_1 = derived[[1L]],
                 records_2 = derived[[2L]], records_4 = derived[[3L]]))
}

# Derive one registered gunzip artifact from its cached source. A cached output
# that still matches its registered identity is kept and re-receipted; anything
# else is rebuilt through a partial file so a failed derivation publishes nothing.
duckhts_bench_stage_gunzip <- function(id, source, destination) {
  if (file.exists(destination)) {
    valid <- tryCatch({
      duckhts_bench_validate_identity(id, destination)
      TRUE
    }, error = function(error) FALSE)
    if (valid) {
      duckhts_bench_write_provenance(id, destination)
      return(invisible(destination))
    }
    unlink(c(destination, paste0(destination, ".provenance.tsv")), force = TRUE)
  }
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(destination, ".partial-", Sys.getpid())
  unlink(temporary, force = TRUE)
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  duckhts_bench_gunzip(source, temporary)
  duckhts_bench_validate_identity(id, temporary)
  if (!file.rename(temporary, destination)) {
    stop("could not publish the uncompressed artifact: ", destination, call. = FALSE)
  }
  duckhts_bench_write_provenance(id, destination)
  invisible(destination)
}

duckhts_bench_gunzip <- function(source, destination) {
  input <- gzfile(source, open = "rb")
  on.exit(close(input), add = TRUE)
  output <- file(destination, open = "wb")
  on.exit(close(output), add = TRUE)
  repeat {
    chunk <- readBin(input, what = "raw", n = 1048576L)
    if (!length(chunk)) break
    writeBin(chunk, output)
  }
  invisible(destination)
}
