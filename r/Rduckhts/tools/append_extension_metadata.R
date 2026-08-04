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

script_args <- commandArgs()
script_arg <- grep("^--file=", script_args, value = TRUE)
if (!length(script_arg)) {
    stop("append_extension_metadata.R must be run with Rscript", call. = FALSE)
}
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "R", "extension_metadata.R"))

opts <- parse_args(args)
required <- c(
    "library-file", "out-file", "extension-name", "duckdb-platform",
    "duckdb-version", "extension-version"
)
missing <- required[vapply(required, function(name) {
    is.null(opts[[name]]) || !nzchar(opts[[name]])
}, logical(1))]
if (length(missing)) {
    stop(
        "Missing required arguments for append_extension_metadata: ",
        paste(missing, collapse = ", "),
        call. = FALSE
    )
}

abi_type <- opts[["abi-type"]]
if (is.null(abi_type) || !nzchar(abi_type)) {
    abi_type <- "C_STRUCT"
}

duckhts_append_extension_metadata(
    library_file = opts[["library-file"]],
    out_file = opts[["out-file"]],
    extension_name = opts[["extension-name"]],
    duckdb_platform = opts[["duckdb-platform"]],
    duckdb_version = opts[["duckdb-version"]],
    extension_version = opts[["extension-version"]],
    abi_type = abi_type
)
