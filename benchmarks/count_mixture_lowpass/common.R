suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(duckdb))
suppressPackageStartupMessages(library(duckhtsbench))
source("benchmarks/count_mixture_lowpass/declaration.R")

lowpass_connect <- function(extension, threads = 1L) {
  temporary <- tempfile("lowpass-duckdb-")
  dir.create(temporary)
  con <- dbConnect(duckdb(shared_home = FALSE, config = list(
    allow_unsigned_extensions = "true", threads = as.character(threads),
    temp_directory = temporary)))
  dbExecute(con, paste("LOAD", dbQuoteString(con, normalizePath(extension))))
  dbExecute(con, paste("SET memory_limit=", dbQuoteString(con, lowpass_declaration$memory_limit)))
  con
}

lowpass_fit <- function(con, error, table = "counts") {
  dbGetQuery(con, sprintf(paste(
    "SELECT * FROM duckhts_count_mixture_lowpass(%s, %.17g, %.17g, 0,",
    "block_bases := 1000000) ORDER BY candidate, relationship_model"),
    dbQuoteString(con, table), error, lowpass_declaration$artefact_weight))
}

# Inspect the actual moment rows passed to the native aggregate, without
# reproducing the production SQL expressions in the benchmark oracle.
lowpass_moment_macro <- function(con) {
  body <- dbGetQuery(con, paste(
    "SELECT macro_definition FROM duckdb_functions()",
    "WHERE function_name='duckhts_count_mixture_lowpass'"))$macro_definition
  stopifnot(length(body) == 1L, grepl(", packed AS", body, fixed = TRUE))
  prefix <- strsplit(body, ", packed AS", fixed = TRUE)[[1L]][[1L]]
  dbExecute(con, paste0("CREATE TEMP MACRO benchmark_moments(counts_table, ",
    "block_bases := 1000000, max_depth := 8, max_candidates := 16) AS TABLE ",
    prefix, " SELECT smp, cand, block, moments.* FROM block_rows ORDER BY smp, cand, block"))
}

lowpass_moments <- function(con, block_bases = lowpass_declaration$block_bases) {
  dbGetQuery(con, sprintf("SELECT * FROM benchmark_moments('counts', block_bases := %d)", block_bases))
}

lowpass_relative_error <- function(actual, expected) {
  stopifnot(length(actual) == length(expected))
  max(abs(actual - expected) / pmax(1, abs(expected)))
}

lowpass_oracle <- function(counts, native, block_bases = lowpass_declaration$block_bases) {
  errors <- numeric()
  counter_ok <- TRUE
  for (candidate in unique(counts$candidate)) {
    x <- counts[counts$candidate == candidate, ]
    physical <- data.frame(chrom = x$chrom, index = floor((x$pos - 1) / block_bases))
    blocks <- unique(physical)
    blocks <- blocks[order(blocks$chrom, blocks$index), ]
    for (b in seq_len(nrow(blocks))) {
      d <- x[physical$chrom == blocks$chrom[[b]] & physical$index == blocks$index[[b]], ]
      n <- d$ref_count + d$alt_count
      usable <- complete.cases(d) & n <= lowpass_declaration$max_depth
      n[!usable] <- 0
      a <- d$alt_count
      a[!usable] <- 0
      f <- d$sample_af
      g <- d$candidate_af
      f[is.na(f)] <- 0
      g[is.na(g)] <- 0
      pairs <- n * (n - 1) / 2
      discordant <- a * (n - a)
      singles_basis <- cbind(1, f, g)
      pairs_basis <- cbind(1, f, g, f^2, f*g, g^2)
      expected <- c(
        as.numeric(crossprod(singles_basis, n * singles_basis)),
        colSums(a * singles_basis),
        as.numeric(crossprod(pairs_basis, pairs * pairs_basis)),
        colSums(discordant * pairs_basis), sum(pairs * f * (1-f)))
      row <- native[native$cand == candidate & native$block == b - 1L, ]
      stopifnot(nrow(row) == 1L)
      actual <- c(row$s[[1L]], row$y[[1L]], row$t[[1L]], row$z[[1L]], row$denominator)
      errors <- c(errors, lowpass_relative_error(actual, expected))
      counters <- c(nrow(d), sum(usable), sum(usable & n == 0), sum(n == 1),
                    sum(n), sum(d$other_count), sum(n > 1),
                    sum(d$ref_count + d$alt_count > lowpass_declaration$max_depth),
                    sum((n-a)*(n-a-1)/2), sum(discordant), sum(a*(a-1)/2),
                    sum(!complete.cases(d)))
      counter_ok <- counter_ok && identical(as.numeric(row$counts[[1L]]), as.numeric(counters))
    }
  }
  list(moment_error = max(errors), counters_exact = counter_ok)
}

lowpass_thread_check <- function(one, four, one_moments, four_moments) {
  fields <- c("s", "y", "t", "z", "denominator")
  a <- unlist(one_moments[fields], use.names = FALSE)
  b <- unlist(four_moments[fields], use.names = FALSE)
  finite <- is.finite(one$share) & is.finite(four$share)
  list(status_equal = identical(one$status, four$status),
       share_missing_equal = identical(is.na(one$share), is.na(four$share)),
       share_error = if (any(finite)) lowpass_relative_error(one$share[finite], four$share[finite]) else 0,
       thread_moment_error = lowpass_relative_error(a, b))
}
