#!/usr/bin/env Rscript
# Network-free checks for the comparator and record-ID-only transformation.
source("test/duckvep/conformance/projection_differential.R")
source("test/duckvep/conformance/projection_fixtures.R")
expected <- data.frame(ID = c("a", "b"), transcript_id = "tx",
  cdna_start = c("4", NA_character_), reference_codons = c("Aaa", NA_character_))
stopifnot(nrow(duckvep_projection_compare(expected[2:1, ], expected)) == 0L)
for (field in names(expected)) {
  corrupt <- expected
  corrupt[[field]][[1L]] <- "wrong"
  stopifnot(nrow(duckvep_projection_compare(corrupt, expected)) > 0L)
}
for (field in c("cdna_start", "reference_codons")) {
  corrupt <- expected
  corrupt[[field]][[1L]] <- NA_character_
  stopifnot(nrow(duckvep_projection_compare(corrupt, expected)) == 1L)
  corrupt <- expected
  corrupt[[field]][[2L]] <- "value"
  stopifnot(nrow(duckvep_projection_compare(corrupt, expected)) == 1L)
}
stopifnot(nrow(duckvep_projection_compare(expected[-1L, ], expected)) == 1L,
  nrow(duckvep_projection_compare(rbind(expected, expected[1L, ]), expected)) > 0L,
  nrow(duckvep_projection_compare(expected, rbind(expected, expected[1L, ]))) > 0L)
local({
  path <- tempfile(fileext = ".json")
  on.exit(unlink(path))
  # Repeated biological alleles remain distinct physical records. Native
  # comparison must not inherit the corpus runner's distinct-allele sampling.
  records <- lapply(c("first", "second"), function(id) list(
    input = paste("chrDuck", 158, id, "C", "A", ".", ".", ".", sep = "\t"),
    transcript_consequences = list(list(transcript_id = "tx",
      consequence_terms = c("start_lost", "splice_region_variant")))))
  writeLines(vapply(records, function(x) jsonlite::toJSON(x, auto_unbox = TRUE), character(1L)), path)
  expected_so <- duckvep_projection_consequences_expected(path)
  stopifnot(nrow(expected_so) == 2L, identical(expected_so$ID, c("first", "second")),
    all(expected_so$consequence_terms == "splice_region_variant&start_lost"),
    nrow(duckvep_projection_compare(expected_so[2:1, ], expected_so)) == 0L)
  for (field in names(expected_so)) {
    corrupt <- expected_so
    corrupt[[field]][[1L]] <- "wrong"
    stopifnot(nrow(duckvep_projection_compare(corrupt, expected_so)) > 0L)
  }
  unresolved <- expected_so
  unresolved$status[[1L]] <- "unresolved"
  unresolved$reason[[1L]] <- "invalid_model_projection"
  stopifnot(nrow(duckvep_projection_compare(unresolved, expected_so)) == 2L,
    nrow(duckvep_projection_compare(expected_so[-1L, ], expected_so)) == 1L,
    nrow(duckvep_projection_compare(rbind(expected_so, expected_so[1L, ]), expected_so)) > 0L)
})
local({
  directory <- tempfile("projection-label-test-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  input <- file.path(directory, "original.vcf")
  output <- file.path(directory, "labelled.vcf")
  original <- c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    rep("chrDuck\t120\tduplicate\tA\tG\t.\t.\tCLASS=snv", 2L))
  writeLines(original, input)
  stopifnot(duckvep_projection_label_records(input, output) == 2L,
    identical(readLines(input), original))
  actual <- readLines(output)
  stopifnot(identical(actual[1:2], original[1:2]))
  fields <- strsplit(actual[-(1:2)], "\t", fixed = TRUE)
  stopifnot(!anyDuplicated(vapply(fields, `[`, character(1L), 3L)))
  fields <- lapply(fields, function(row) { row[[3L]] <- "duplicate"; row })
  stopifnot(identical(vapply(fields, paste, character(1L), collapse = "\t"), original[-(1:2)]))
})
message("projection comparator and record preservation: OK")
