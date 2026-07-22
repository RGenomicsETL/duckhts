library(tinytest)

test_htslib_contract <- function() {
  link <- if (.Platform$OS.type == "windows") "static" else "shared"
  config <- rduckhts_htslib_config(link)
  expect_identical(config$contract_version, 1L)
  expect_identical(config$htslib_version, "1.24")
  expect_identical(config$runtime_version, config$htslib_version)
  expect_identical(config$htslib_header_version, 102400L)
  expect_true(file.exists(file.path(config$include_dir, "htslib", "hts.h")))
  expect_true(file.exists(config$library_file))
  expect_true(nzchar(config$cppflags))
  expect_true(nzchar(config$ldflags))
  expect_true(is.list(config$features))
  expect_true(is.numeric(config$runtime_feature_bits))
  expect_true(nzchar(config$runtime_feature_string))

  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- DBI::dbConnect(drv)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  rduckhts_load(con)
  info <- rduckhts_htslib_info(con)
  expect_identical(as.character(info$version[[1L]]), "1.24")
  expect_true(info$feature_bits[[1L]] > 0)
  expect_true(grepl("libcurl=", info$feature_string[[1L]], fixed = TRUE))

  contract_path <- system.file("htslib_config.R", package = "Rduckhts")
  contract <- new.env(parent = baseenv())
  sys.source(contract_path, envir = contract)
  relocated <- tempfile("rduckhts_htslib_")
  dir.create(file.path(relocated, "lib"), recursive = TRUE)
  expect_true(file.copy(config$include_dir, relocated, recursive = TRUE))
  relocated_library <- file.path(relocated, "lib", basename(config$library_file))
  if (!file.symlink(config$library_file, relocated_library)) {
    expect_true(file.copy(config$library_file, relocated_library))
  }
  relocated_config <- contract$htslib_config(link, root = relocated)
  expect_identical(relocated_config$include_dir, file.path(relocated, "include"))
  expect_identical(relocated_config$library_file, relocated_library)

  relocated_header <- file.path(relocated, "include", "htslib", "hts.h")
  header <- readLines(relocated_header, warn = FALSE)
  header <- sub("#define HTS_VERSION 102400", "#define HTS_VERSION 102399", header,
                fixed = TRUE)
  writeLines(header, relocated_header)
  expect_error(
    contract$htslib_config(link, root = relocated),
    pattern = "header/receipt version mismatch"
  )

  rtinycc_artifact_api <- requireNamespace("Rtinycc", quietly = TRUE) &&
    all(c(
      "tcc_state", "tcc_include_paths", "tcc_lib_paths",
      "tcc_add_library", "tcc_compile_string", "tcc_output_file"
    ) %in% getNamespaceExports("Rtinycc"))
  if (
    .Platform$OS.type != "unix" ||
      !isTRUE(config$shared_available) ||
      !rtinycc_artifact_api
  ) {
    return(invisible(NULL))
  }

  consumer_code <- paste(c(
    "#include <stdio.h>",
    "#include <string.h>",
    "#include <htslib/hts.h>",
    "#include <htslib/sam.h>",
    "#include <htslib/vcf.h>",
    "static int alignment_ok(const char *path, const char *reference) {",
    "    samFile *file = sam_open(path, \"r\");",
    "    sam_hdr_t *header; bam1_t *record; int result;",
    "    if (!file) return 0;",
    "    if (reference && hts_set_fai_filename(file, reference) != 0) return 0;",
    "    header = sam_hdr_read(file); record = bam_init1();",
    "    result = header && record && sam_read1(file, header, record) >= 0;",
    "    bam_destroy1(record); sam_hdr_destroy(header); sam_close(file);",
    "    return result;",
    "}",
    "static int variant_ok(const char *path) {",
    "    htsFile *file = bcf_open(path, \"r\");",
    "    bcf_hdr_t *header; bcf1_t *record; int result;",
    "    if (!file) return 0;",
    "    header = bcf_hdr_read(file); record = bcf_init();",
    "    result = header && record && bcf_read(file, header, record) == 0;",
    "    bcf_destroy(record); bcf_hdr_destroy(header); bcf_close(file);",
    "    return result;",
    "}",
    "int main(int argc, char **argv) {",
    "    if (argc != 7 || strcmp(hts_version(), argv[1]) != 0) return 10;",
    "    if (!alignment_ok(argv[2], NULL)) return 11;",
    "    if (!alignment_ok(argv[3], argv[4])) return 12;",
    "    if (!variant_ok(argv[5]) || !variant_ok(argv[6])) return 13;",
    "    puts(hts_version()); return 0;",
    "}"
  ), collapse = "\n")

  executable <- tempfile("rduckhts_consumer_")
  state <- Rtinycc::tcc_state(
    output = "exe",
    include_path = c(Rtinycc::tcc_include_paths(), config$include_dir),
    lib_path = c(Rtinycc::tcc_lib_paths(), config$lib_dir)
  )
  expect_identical(Rtinycc::tcc_add_library(state, "hts"), 0L)
  expect_identical(Rtinycc::tcc_compile_string(state, consumer_code), 0L)
  expect_identical(Rtinycc::tcc_output_file(state, executable), 0L)

  fixture <- function(name) {
    system.file("extdata", name, package = "Rduckhts", mustWork = TRUE)
  }
  loader_path <- paste(
    c(config$lib_dir, Sys.getenv("LD_LIBRARY_PATH")),
    collapse = .Platform$path.sep
  )
  output <- system2(
    executable,
    shQuote(c(
      config$htslib_version,
      fixture("range.bam"),
      fixture("range.cram"),
      fixture("ce.fa"),
      fixture("vcf_file.bcf"),
      fixture("test_vep_tidy.vcf")
    )),
    stdout = TRUE,
    stderr = TRUE,
    env = paste0("LD_LIBRARY_PATH=", shQuote(loader_path))
  )
  expect_identical(attr(output, "status"), NULL)
  expect_identical(output[[length(output)]], config$htslib_version)
}

test_htslib_contract()
