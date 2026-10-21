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
  warning_prefix <- "\\[W::[^]]+\\]"
  leading_warning_prefixes <- paste0("^(", warning_prefix, "[[:space:]]*)+")
  expected_htslib_warning <- rep(FALSE, length(diagnostics))
  unexpected_htslib_warning <- rep(FALSE, length(diagnostics))
  pending_warning_prefixes <- character()

  for (index in seq_along(diagnostics)) {
    line <- diagnostics[[index]]
    prefix_text <- regmatches(
      line, regexpr(leading_warning_prefixes, line)
    )
    if (length(prefix_text) > 0L) {
      prefixes <- regmatches(
        prefix_text, gregexpr(warning_prefix, prefix_text)
      )[[1L]]
      pending_warning_prefixes <- c(pending_warning_prefixes, prefixes)
      message <- sub(leading_warning_prefixes, "", line)
      expected_htslib_warning[[index]] <- TRUE
    } else if (length(pending_warning_prefixes) > 0L) {
      message <- line
      expected_htslib_warning[[index]] <- TRUE
    } else {
      next
    }

    if (!nzchar(message)) next
    candidates <- paste(pending_warning_prefixes, message)
    matching_prefix <- vapply(candidates, function(candidate) {
      any(vapply(expected_patterns, grepl, logical(1L), x = candidate))
    }, logical(1L))
    if (!any(matching_prefix)) {
      unexpected_htslib_warning[[index]] <- TRUE
      next
    }
    pending_warning_prefixes <- pending_warning_prefixes[-which(matching_prefix)[[1L]]]
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
