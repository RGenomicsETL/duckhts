args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
    res <- list()
    i <- 1L
    while (i <= length(args)) {
        key <- args[[i]]
        if (startsWith(key, "--") && i + 1L <= length(args)) {
            res[[substring(key, 3L)]] <- args[[i + 1L]]
            i <- i + 2L
        } else {
            i <- i + 1L
        }
    }
    res
}

opts <- parse_args(args)
required <- c(
    "template", "output", "static-default", "static-dep-libs",
    "htslib-version", "htslib-version-number", "rduckhts-version",
    "duckdb-platform", "has-zlib", "has-bz2", "has-lzma",
    "has-libdeflate", "has-curl", "has-openssl", "has-plugins",
    "has-s3", "has-gcs"
)
missing <- required[vapply(required, function(name) {
    is.null(opts[[name]])
}, logical(1))]
if (length(missing)) {
    stop(
        "Missing required arguments for write_htslib_config: ",
        paste(missing, collapse = ", "),
        call. = FALSE
    )
}

if (!file.exists(opts$template)) {
    stop("htslib config template not found: ", opts$template, call. = FALSE)
}

as_r_bool <- function(value, name) {
    if (!value %in% c("yes", "no")) {
        stop(name, " must be yes or no", call. = FALSE)
    }
    if (identical(value, "yes")) "TRUE" else "FALSE"
}

escape_r_string <- function(value) {
    value <- gsub("\\\\", "\\\\\\\\", value)
    gsub('"', '\\\\"', value, fixed = TRUE)
}

replacements <- c(
    HTSLIB_STATIC_DEFAULT = opts[["static-default"]],
    HTSLIB_STATIC_DEP_LIBS = escape_r_string(opts[["static-dep-libs"]]),
    HTSLIB_VERSION = opts[["htslib-version"]],
    HTSLIB_VERSION_NUMBER = opts[["htslib-version-number"]],
    RDUCKHTS_VERSION = opts[["rduckhts-version"]],
    DUCKDB_PLATFORM = opts[["duckdb-platform"]],
    HTSLIB_HAS_ZLIB = as_r_bool(opts[["has-zlib"]], "has-zlib"),
    HTSLIB_HAS_BZ2 = as_r_bool(opts[["has-bz2"]], "has-bz2"),
    HTSLIB_HAS_LZMA = as_r_bool(opts[["has-lzma"]], "has-lzma"),
    HTSLIB_HAS_LIBDEFLATE = as_r_bool(opts[["has-libdeflate"]], "has-libdeflate"),
    HTSLIB_HAS_CURL = as_r_bool(opts[["has-curl"]], "has-curl"),
    HTSLIB_HAS_OPENSSL = as_r_bool(opts[["has-openssl"]], "has-openssl"),
    HTSLIB_HAS_PLUGINS = as_r_bool(opts[["has-plugins"]], "has-plugins"),
    HTSLIB_HAS_S3 = as_r_bool(opts[["has-s3"]], "has-s3"),
    HTSLIB_HAS_GCS = as_r_bool(opts[["has-gcs"]], "has-gcs")
)

text <- paste(readLines(opts$template, warn = FALSE), collapse = "\n")
for (name in names(replacements)) {
    token <- paste0("@", name, "@")
    parts <- strsplit(text, token, fixed = TRUE)[[1L]]
    if (length(parts) == 1L) {
        stop("Template is missing token: ", token, call. = FALSE)
    }
    text <- paste(parts, collapse = replacements[[name]])
}

out_dir <- dirname(opts$output)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(text, opts$output, useBytes = TRUE)
