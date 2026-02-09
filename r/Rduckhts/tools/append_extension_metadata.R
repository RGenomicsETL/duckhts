args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
    res <- list()
    i <- 1
    while (i <= length(args)) {
        key <- args[[i]]
        if (startsWith(key, "--") && i + 1 <= length(args)) {
            res[[substring(key, 3)]] <- args[[i + 1]]
            i <- i + 2
        } else {
            i <- i + 1
        }
    }
    res
}

start_signature <- function() {
    sig <- as.raw(c(0, 147, 4, 16))
    sig <- c(sig, charToRaw("duckdb_signature"))
    sig <- c(sig, as.raw(c(128, 4)))
    sig
}

padded_byte_string <- function(x) {
    bytes <- charToRaw(x)
    if (length(bytes) > 32) {
        bytes <- bytes[1:32]
    }
    c(bytes, as.raw(rep(0, 32 - length(bytes))))
}

opts <- parse_args(args)

library_file <- opts[["library-file"]]
out_file <- opts[["out-file"]]
extension_name <- opts[["extension-name"]]
duckdb_platform <- opts[["duckdb-platform"]]
duckdb_version <- opts[["duckdb-version"]]
extension_version <- opts[["extension-version"]]
abi_type <- opts[["abi-type"]]

if (is.null(library_file) || is.null(out_file) || is.null(extension_name) ||
    is.null(duckdb_platform) || is.null(duckdb_version) || is.null(extension_version)) {
    stop("Missing required arguments for append_extension_metadata")
}

if (is.null(abi_type) || abi_type == "") {
    abi_type <- "C_STRUCT"
}

out_tmp <- paste0(out_file, ".tmp")

file.copy(library_file, out_tmp, overwrite = TRUE)

con <- file(out_tmp, open = "ab")
writeBin(start_signature(), con, useBytes = TRUE)
writeBin(padded_byte_string(""), con, useBytes = TRUE)
writeBin(padded_byte_string(""), con, useBytes = TRUE)
writeBin(padded_byte_string(""), con, useBytes = TRUE)
writeBin(padded_byte_string(abi_type), con, useBytes = TRUE)
writeBin(padded_byte_string(extension_version), con, useBytes = TRUE)
writeBin(padded_byte_string(duckdb_version), con, useBytes = TRUE)
writeBin(padded_byte_string(duckdb_platform), con, useBytes = TRUE)
writeBin(padded_byte_string("4"), con, useBytes = TRUE)
writeBin(as.raw(rep(0, 256)), con, useBytes = TRUE)
close(con)

file.rename(out_tmp, out_file)
