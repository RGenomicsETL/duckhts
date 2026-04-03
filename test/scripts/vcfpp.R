#!/usr/bin/env Rscript

# Generate compact VCF fixtures for VCF-spec and DuckDB-mapping compliance.
#
# Sections:
#   1. Setup and argument parsing
#   2. Shared helpers
#   3. Mapping-compliance fixtures
#   4. Spec-compliance fixtures
#   5. Manifest and summary
#
# Notes:
#   - Uses vcfppR's writer API for reproducible fixture generation.
#   - Intentionally skips CSQ/ANN/BCSQ examples here because those go through
#     the dedicated VEP parser path and are tested separately.

args <- commandArgs(trailingOnly = TRUE)
full_args <- commandArgs(FALSE)
script_arg <- grep("^--file=", full_args, value = TRUE)
if (length(script_arg) != 1) {
  stop("Unable to determine script path from commandArgs()", call. = FALSE)
}
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

if (!requireNamespace("vcfppR", quietly = TRUE)) {
  stop(
    "vcfppR is required to generate the VCF compliance fixtures. Install vcfppR first.",
    call. = FALSE
  )
}

default_output_dirs <- c(
  file.path(repo_root, "test", "data"),
  file.path(repo_root, "r", "Rduckhts", "inst", "extdata")
)
output_dirs <- if (length(args) > 0) normalizePath(args, mustWork = FALSE) else default_output_dirs

for (output_dir in output_dirs) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# ---------------------------------------------------------------------------
# Section 2. Shared Helpers
# ---------------------------------------------------------------------------

tag_def <- function(id, number, type, description) {
  list(id = id, number = number, type = type, description = description)
}

write_fixture <- function(path,
                          contigs,
                          info_defs = list(),
                          format_defs = list(),
                          samples = character(),
                          records,
                          extra_lines = character()) {
  writer <- methods::new(vcfppR::vcfwriter, path, "VCFv4.2")
  on.exit(try(writer$close(), silent = TRUE), add = TRUE)

  writer$addLine("##source=duckhts/test/scripts/vcfpp.R")
  for (line in extra_lines) {
    writer$addLine(line)
  }
  for (contig in contigs) {
    writer$addContig(contig)
  }
  for (def in info_defs) {
    writer$addINFO(def$id, def$number, def$type, def$description)
  }
  for (def in format_defs) {
    writer$addFORMAT(def$id, def$number, def$type, def$description)
  }
  for (sample in samples) {
    writer$addSample(sample)
  }
  for (record in records) {
    writer$writeline(record)
  }

  writer$close()
  invisible(path)
}

render_fixture <- function(filename,
                           section,
                           purpose,
                           contigs,
                           info_defs = list(),
                           format_defs = list(),
                           samples = character(),
                           records,
                           extra_lines = character()) {
  written_paths <- character(length(output_dirs))
  for (i in seq_along(output_dirs)) {
    path <- file.path(output_dirs[[i]], filename)
    write_fixture(
      path = path,
      contigs = contigs,
      info_defs = info_defs,
      format_defs = format_defs,
      samples = samples,
      records = records,
      extra_lines = extra_lines
    )
    written_paths[[i]] <- path
  }

  data.frame(
    section = section,
    file = filename,
    purpose = purpose,
    outputs = paste(written_paths, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

manifest <- list()

# ---------------------------------------------------------------------------
# Section 3. Mapping-Compliance Fixtures
# ---------------------------------------------------------------------------

manifest[[length(manifest) + 1]] <- render_fixture(
  filename = "format_string_list.vcf",
  section = "mapping",
  purpose = "FORMAT string list example for Number=. fields such as DRAGEN LAA",
  contigs = c("chrM"),
  format_defs = list(
    tag_def("GT", "1", "String", "Genotype"),
    tag_def("LAA", ".", "String", "Local alt allele mapping")
  ),
  samples = c("SAMPLE1"),
  records = c(
    "chrM\t73\t.\tA\tC,G\t.\tPASS\t.\tGT:LAA\t1/2:1,2",
    "chrM\t150\t.\tT\tC\t.\tPASS\t.\tGT:LAA\t0/1:1",
    "chrM\t263\t.\tA\tG\t.\tPASS\t.\tGT:LAA\t0/0:."
  ),
  extra_lines = c("##duckhts_fixture_section=mapping_string_lists")
)

manifest[[length(manifest) + 1]] <- render_fixture(
  filename = "fixed_count_arrays.vcf",
  section = "mapping",
  purpose = "Exact fixed-cardinality INFO/FORMAT arrays for Number=2 and Number=4",
  contigs = c("chr1"),
  info_defs = list(
    tag_def("SB", "4", "Integer", "Strand counts")
  ),
  format_defs = list(
    tag_def("GT", "1", "String", "Genotype"),
    tag_def("HQ", "2", "Integer", "Haplotype qualities")
  ),
  samples = c("S1", "S2"),
  records = c(
    "chr1\t100\t.\tA\tC\t.\tPASS\tSB=10,20,30,40\tGT:HQ\t0/1:5,9\t1/1:7,11",
    "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT:HQ\t0/0:.,.\t./.:."
  ),
  extra_lines = c("##duckhts_fixture_section=mapping_fixed_counts")
)

manifest[[length(manifest) + 1]] <- render_fixture(
  filename = "mapping_number_families.vcf",
  section = "mapping",
  purpose = "Representative Number=1/A/R/G/. mappings across INFO and FORMAT fields",
  contigs = c("chr1"),
  info_defs = list(
    tag_def("IINT", "1", "Integer", "Scalar integer INFO"),
    tag_def("IFLT", "1", "Float", "Scalar float INFO"),
    tag_def("ISTR", "1", "String", "Scalar string INFO"),
    tag_def("FLAG", "0", "Flag", "Flag INFO"),
    tag_def("AINT", "A", "Integer", "Per-alt integer INFO"),
    tag_def("AFLT", "A", "Float", "Per-alt float INFO"),
    tag_def("RINT", "R", "Integer", "Per-allele integer INFO"),
    tag_def("GFLOAT", "G", "Float", "Per-genotype float INFO"),
    tag_def("VSTR", ".", "String", "Variable-length string INFO"),
    tag_def("VINT", ".", "Integer", "Variable-length integer INFO")
  ),
  format_defs = list(
    tag_def("GT", "1", "String", "Genotype"),
    tag_def("DP", "1", "Integer", "Read depth"),
    tag_def("AD", "R", "Integer", "Allelic depths"),
    tag_def("PL", "G", "Integer", "Phred-scaled genotype likelihoods"),
    tag_def("FT", "1", "String", "Sample filter"),
    tag_def("LAA", ".", "String", "Variable-length string FORMAT")
  ),
  samples = c("S1", "S2"),
  records = c(
    paste(
      "chr1\t100\trsA\tA\tC\t50\tPASS",
      "IINT=7;IFLT=1.5;ISTR=alpha;FLAG;AINT=2;AFLT=0.5;RINT=12,3;GFLOAT=0,10,20;VSTR=x,y;VINT=5,8",
      "GT:DP:AD:PL:FT:LAA",
      "0/1:12:9,3:0,10,20:PASS:1",
      "1/1:20:0,20:50,5,0:PASS:1",
      sep = "\t"
    ),
    paste(
      "chr1\t200\trsB\tG\tA,T\t99\tPASS",
      "IINT=11;IFLT=2.25;ISTR=beta;AINT=4,1;AFLT=0.40,0.10;RINT=20,7,2;GFLOAT=0,5,10,15,20,25;VSTR=red,blue,green;VINT=1,2,3",
      "GT:DP:AD:PL:FT:LAA",
      "1/2:18:9,5,4:0,10,20,30,40,50:LowQual:1,2",
      "0/2:15:10,0,5:0,15,25,35,45,55:PASS:2",
      sep = "\t"
    )
  ),
  extra_lines = c("##duckhts_fixture_section=mapping_number_families")
)

# ---------------------------------------------------------------------------
# Section 4. Spec-Compliance Fixtures
# ---------------------------------------------------------------------------

manifest[[length(manifest) + 1]] <- render_fixture(
  filename = "spec_standard_corrections.vcf",
  section = "spec",
  purpose = "Intentionally misdeclared standard tags that should be corrected by spec-aware bind logic",
  contigs = c("chr1"),
  info_defs = list(
    tag_def("SB", "1", "Integer", "Intentionally wrong: standard SB should be Number=4")
  ),
  format_defs = list(
    tag_def("GT", "1", "String", "Genotype"),
    tag_def("HQ", "1", "Integer", "Intentionally wrong: standard HQ should be Number=2"),
    tag_def("FT", ".", "String", "Intentionally wrong: standard FT should be Number=1")
  ),
  samples = c("S1", "S2"),
  records = c(
    "chr1\t300\t.\tC\tT\t.\tPASS\tSB=10,20,30,40\tGT:HQ:FT\t0/1:5,9:PASS\t1/1:7,11:LowQual",
    "chr1\t350\t.\tG\tA\t.\tPASS\t.\tGT:HQ:FT\t0/0:.,.:PASS\t./.:.:."
  ),
  extra_lines = c("##duckhts_fixture_section=spec_standard_corrections")
)

# ---------------------------------------------------------------------------
# Section 5. Manifest and Summary
# ---------------------------------------------------------------------------

manifest_df <- do.call(rbind, manifest)

for (output_dir in output_dirs) {
  manifest_path <- file.path(output_dir, "vcfpp_manifest.tsv")
  write.table(
    manifest_df[, c("section", "file", "purpose")],
    file = manifest_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

cat("Generated vcfppR fixtures:\n")
for (i in seq_len(nrow(manifest_df))) {
  cat(
    sprintf(
      "  [%s] %s - %s\n",
      manifest_df$section[[i]],
      manifest_df$file[[i]],
      manifest_df$purpose[[i]]
    )
  )
}
