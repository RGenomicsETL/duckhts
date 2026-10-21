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
    "^(\\[W::vcf_parse\\] )+Contig '.*' is not defined in the header",
    "^(\\[W::vcf_parse\\] )+Missing headers may cause later processing to fail$",
    "^(\\[W::vcf_parse\\] )+POS is 0, this is probably wrong$"
  )
  expected_htslib_warning <- rep(FALSE, length(diagnostics))
  for (pattern in expected_patterns) {
    expected_htslib_warning <- expected_htslib_warning |
      grepl(pattern, diagnostics)
  }

  if (!identical(status, 0L)) {
    if (length(diagnostics) > 0L) writeLines(diagnostics, stderr())
    unlink(diagnostics_file)
    stop("tinytest subprocess failed with status ", status, call. = FALSE)
  }

  unexpected_htslib_warning <- grepl("^\\[W::", diagnostics) &
    !expected_htslib_warning
  if (any(unexpected_htslib_warning)) {
    writeLines(diagnostics[!expected_htslib_warning], stderr())
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
