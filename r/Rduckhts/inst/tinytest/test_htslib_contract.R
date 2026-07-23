library(tinytest)

test_htslib_contract <- function() {
  contract_path <- system.file("htslib_config.R", package = "Rduckhts")
  contract <- new.env(parent = baseenv())
  sys.source(contract_path, envir = contract)
  link <- contract$htslib_default_link()
  default_config <- rduckhts_htslib_config()
  expect_identical(default_config$link, link)
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

  expect_identical(contract$htslib_default_link(), link)
  expect_identical(contract$htslib_config()$link, link)
  expect_identical(contract$htslib_rpath(), contract$htslib_config()$rpath)
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

  rtinycc_memory_api <- requireNamespace("Rtinycc", quietly = TRUE) &&
    all(c(
      "tcc_state", "tcc_include_paths", "tcc_lib_paths",
      "tcc_add_library", "tcc_compile_string", "tcc_relocate",
      "tcc_call_symbol"
    ) %in% getNamespaceExports("Rtinycc"))
  if (
    .Platform$OS.type != "unix" ||
      !isTRUE(config$shared_available) ||
      !rtinycc_memory_api
  ) {
    return(invisible(NULL))
  }

  fixture <- function(name) {
    system.file("extdata", name, package = "Rduckhts", mustWork = TRUE)
  }
  c_string <- function(value) {
    encodeString(value, quote = "\"")
  }
  consumer_code <- paste(c(
    "#include <htslib/hts.h>",
    "#include <htslib/sam.h>",
    "#include <htslib/vcf.h>",
    paste0(
      "static const char expected_version[] = ",
      c_string(config$htslib_version),
      ";"
    ),
    paste0("static const char bam_path[] = ", c_string(fixture("range.bam")), ";"),
    paste0("static const char cram_path[] = ", c_string(fixture("range.cram")), ";"),
    paste0("static const char reference_path[] = ", c_string(fixture("ce.fa")), ";"),
    paste0("static const char bcf_path[] = ", c_string(fixture("vcf_file.bcf")), ";"),
    paste0(
      "static const char vcf_path[] = ",
      c_string(fixture("test_vep_tidy.vcf")),
      ";"
    ),
    "static int text_equal(const char *left, const char *right) {",
    "    while (*left && *left == *right) { left++; right++; }",
    "    return *left == *right;",
    "}",
    "static int alignment_ok(const char *path, const char *reference) {",
    "    samFile *file = sam_open(path, \"r\");",
    "    sam_hdr_t *header; bam1_t *record; int result;",
    "    if (!file) return 0;",
    "    if (reference && hts_set_fai_filename(file, reference) != 0) {",
    "        sam_close(file); return 0;",
    "    }",
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
    "int check_consumer(void) {",
    "    if (!text_equal(hts_version(), expected_version)) return 10;",
    "    if (!alignment_ok(bam_path, NULL)) return 11;",
    "    if (!alignment_ok(cram_path, reference_path)) return 12;",
    "    if (!variant_ok(bcf_path) || !variant_ok(vcf_path)) return 13;",
    "    return 0;",
    "}"
  ), collapse = "\n")

  system_include_paths <- character()
  if (identical(Sys.info()[["sysname"]], "Darwin")) {
    sdk_root <- Sys.getenv("SDKROOT")
    if (!nzchar(sdk_root)) {
      xcrun <- Sys.which("xcrun")
      if (nzchar(xcrun)) {
        sdk_root <- suppressWarnings(system2(
          xcrun,
          c("--sdk", "macosx", "--show-sdk-path"),
          stdout = TRUE,
          stderr = FALSE
        ))
      }
    }
    sdk_root <- sdk_root[nzchar(sdk_root)][1L]
    sdk_include <- if (length(sdk_root)) {
      file.path(sdk_root, "usr", "include")
    } else {
      ""
    }
    expect_true(dir.exists(sdk_include))
    if (!dir.exists(sdk_include)) return(invisible(NULL))
    system_include_paths <- sdk_include
  }

  state <- Rtinycc::tcc_state(
    output = "memory",
    include_path = c(
      Rtinycc::tcc_include_paths(),
      system_include_paths,
      config$include_dir
    ),
    lib_path = c(Rtinycc::tcc_lib_paths(), config$lib_dir)
  )
  link_status <- Rtinycc::tcc_add_library(state, "hts")
  expect_identical(link_status, 0L)
  if (link_status != 0L) return(invisible(NULL))

  compile_status <- Rtinycc::tcc_compile_string(state, consumer_code)
  expect_identical(compile_status, 0L)
  if (compile_status != 0L) return(invisible(NULL))

  relocate_status <- Rtinycc::tcc_relocate(state)
  expect_identical(relocate_status, 0L)
  if (relocate_status != 0L) return(invisible(NULL))

  result <- Rtinycc::tcc_call_symbol(
    state,
    "check_consumer",
    return = "int"
  )
  expect_identical(result, 0L)
}

test_htslib_contract()
