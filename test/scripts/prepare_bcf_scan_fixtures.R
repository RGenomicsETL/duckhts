# Maintainer-only, network-free derivation from vcfpp.R's committed source.
# Run from the repository root after vcfpp.R. Originally encoded with bcftools
# 1.23.1-70-g6dbd8fef / HTSlib 1.22.1 and bgzip 1.19; no external data is required.
stopifnot(file.exists("test/data/bcf_scan_contigs.vcf"),
          nzchar(Sys.which("bcftools")), nzchar(Sys.which("bgzip")))

prepare_bcf_scan_fixtures <- function() {
  directory <- tempfile("bcf-scan-fixtures-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  source <- "test/data/bcf_scan_contigs.vcf"
  lines <- readLines(source)
  records <- lines[!startsWith(lines, "#")]
  stopifnot(length(records) == 5L, sum(duplicated(records)) == 1L)
  run <- function(args) stopifnot(system2("bcftools", shQuote(args)) == 0L)
  outputs <- source
  for (header in c("full", "partial", "none", "empty")) {
    selected <- switch(header,
      full = lines,
      partial = lines[!startsWith(lines, "##contig=<ID=chr3,") &
                      lines != "##contig=<ID=chr3>"],
      none = lines[!startsWith(lines, "##contig=")],
      empty = lines[startsWith(lines, "#")])
    text <- file.path(directory, "input.vcf")
    writeLines(selected, text)
    path <- paste0("test/data/bcf_scan_contigs.", header, ".vcf.gz")
    stopifnot(system2("bgzip", c("-c", shQuote(text)), stdout = path) == 0L)
    for (format in c("tbi", "csi")) {
      index <- paste0(path, ".index.", format) # Deliberately not autodiscovered.
      run(c("index", "-f", if (format == "tbi") "-t" else "-c", "-o", index, path))
      outputs <- c(outputs, index)
    }
    decoded <- system2("bcftools", c("view", "-H", shQuote(path)), stdout = TRUE)
    stopifnot(identical(decoded, if (header == "empty") character() else records))
    outputs <- c(outputs, path)
  }
  bcf <- "test/data/bcf_scan_contigs.bcf"
  run(c("view", "--no-version", "-Ob", "-o", bcf, source))
  run(c("index", "-f", bcf))
  stopifnot(identical(system2("bcftools", c("view", "-H", bcf), stdout = TRUE), records))
  outputs <- c(outputs, bcf, paste0(bcf, ".csi"))
  stopifnot(all(file.copy(outputs, "r/Rduckhts/inst/extdata", overwrite = TRUE)))
}

prepare_bcf_scan_fixtures()
