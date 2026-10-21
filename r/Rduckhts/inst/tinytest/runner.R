if (requireNamespace("tinytest", quietly = TRUE)) {
  # HTSlib writes directly to the process stderr stream, outside R sinks.
  diagnostics_file <- tempfile("rduckhts-tinytest-stderr-")

  expression <- "tinytest::test_package('Rduckhts')"
  status <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", "-e", shQuote(expression)),
    stdout = "",
    stderr = diagnostics_file
  ))
  diagnostics <- readLines(diagnostics_file, warn = FALSE)

  if (!identical(status, 0L)) {
    if (length(diagnostics) > 0L) writeLines(diagnostics, stderr())
    unlink(diagnostics_file)
    stop("tinytest subprocess failed with status ", status, call. = FALSE)
  }

  expected_patterns <- c(
    "^\\[W::accept_mismatch\\] ",
    "^\\[W::bam_hdr_read\\] EOF marker is absent\\. The input is probably truncated$",
    "^\\[W::bcf_handle_decode_diagnostic\\] ",
    "^\\[W::bcf_hdr_check_sanity\\] ",
    "^\\[W::fai_get_val\\] Reference .* not found in FASTA file, returning empty sequence$",
    paste0("^\\[W::hts_reglist_create\\] Region '.*' specifies an unknown ",
           "reference name\\. Continue anyway$"),
    "^\\[W::sam_hdr_create\\] Ignored @HD type record without VN tag$",
    "^\\[W::sam_read1_sam\\] Parse error at line 4$",
    "^\\[W::vcf_parse\\] Contig '.*' is not defined in the header",
    "^\\[W::vcf_parse\\] Missing headers may cause later processing to fail$",
    "^\\[W::vcf_parse\\] POS is 0, this is probably wrong$"
  )
  consume_warning_message <- function(message, pending_prefixes) {
    message <- trimws(message)
    if (!nzchar(message)) {
      return(list(pending = pending_prefixes, warning = FALSE,
                  ordinary = FALSE, unexpected = FALSE))
    }
    if (length(pending_prefixes) == 0L) {
      return(list(pending = pending_prefixes, warning = FALSE,
                  ordinary = TRUE, unexpected = FALSE))
    }
    candidates <- paste(pending_prefixes, message)
    matching_prefix <- vapply(candidates, function(candidate) {
      any(vapply(expected_patterns, grepl, logical(1L), x = candidate))
    }, logical(1L))
    if (!any(matching_prefix)) {
      return(list(pending = pending_prefixes, warning = FALSE,
                  ordinary = FALSE, unexpected = TRUE))
    }
    used_prefix <- which(matching_prefix)[[1L]]
    list(pending = pending_prefixes[-used_prefix], warning = TRUE,
         ordinary = FALSE, unexpected = FALSE)
  }

  warning_prefix <- "\\[W::[^]]+\\]"
  expected_htslib_warning <- rep(FALSE, length(diagnostics))
  unexpected_htslib_warning <- rep(FALSE, length(diagnostics))
  pending_warning_prefixes <- character()

  for (index in seq_along(diagnostics)) {
    line <- diagnostics[[index]]
    prefix_positions <- gregexpr(warning_prefix, line)[[1L]]
    prefix_lengths <- attr(prefix_positions, "match.length")
    line_warning <- FALSE
    line_ordinary <- FALSE
    cursor <- 1L

    if (prefix_positions[[1L]] >= 0L) {
      for (prefix_index in seq_along(prefix_positions)) {
        prefix_position <- prefix_positions[[prefix_index]]
        consumed <- consume_warning_message(
          substr(line, cursor, prefix_position - 1L),
          pending_warning_prefixes
        )
        pending_warning_prefixes <- consumed$pending
        line_warning <- line_warning || consumed$warning
        line_ordinary <- line_ordinary || consumed$ordinary
        unexpected_htslib_warning[[index]] <-
          unexpected_htslib_warning[[index]] || consumed$unexpected

        prefix_length <- prefix_lengths[[prefix_index]]
        prefix <- substr(line, prefix_position,
                         prefix_position + prefix_length - 1L)
        pending_warning_prefixes <- c(pending_warning_prefixes, prefix)
        line_warning <- TRUE
        cursor <- prefix_position + prefix_length
      }
    }

    consumed <- consume_warning_message(
      substr(line, cursor, nchar(line)), pending_warning_prefixes
    )
    pending_warning_prefixes <- consumed$pending
    line_warning <- line_warning || consumed$warning
    line_ordinary <- line_ordinary || consumed$ordinary
    unexpected_htslib_warning[[index]] <-
      unexpected_htslib_warning[[index]] || consumed$unexpected
    expected_htslib_warning[[index]] <- line_warning && !line_ordinary
  }

  if (any(unexpected_htslib_warning) || length(pending_warning_prefixes) > 0L) {
    writeLines(diagnostics, stderr())
    unlink(diagnostics_file)
    stop("tinytest emitted an unexpected HTSlib warning", call. = FALSE)
  }
  if (!all(expected_htslib_warning)) {
    writeLines(diagnostics[!expected_htslib_warning], stderr())
  }
  unlink(diagnostics_file)
} else {
  message("tinytest not available; skipping tests")
}
