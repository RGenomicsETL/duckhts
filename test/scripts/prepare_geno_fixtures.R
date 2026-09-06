# Maintainer-only, network-free derivation. Valid GT/PS source is written by
# vcfpp.R. Preserve deliberately invalid PS declarations/widths as separate
# witnesses, not by changing the valid source or weakening the reader oracle.
prepare_geno_fixtures <- function() {
  stopifnot(nzchar(Sys.which("bcftools")), nzchar(Sys.which("bgzip")))
  run <- function(args) stopifnot(system2("bcftools", shQuote(args)) == 0L)
  source <- "test/data/geno_calls.vcf"
  stopifnot(file.exists(source))
  lines <- readLines(source)
  bcf <- sub("vcf$", "bcf", source)
  vcf <- paste0(source, ".gz")
  run(c("view", "--no-version", "-Ob", "-o", bcf, source))
  run(c("index", "-f", bcf))
  stopifnot(system2("bgzip", c("-c", shQuote(source)), stdout = vcf) == 0L)
  run(c("index", "-f", "-t", vcf))
  outputs <- c(source, bcf, paste0(bcf, ".csi"), vcf, paste0(vcf, ".tbi"))
  expected <- system2("bcftools", c("view", "-H", source), stdout = TRUE)
  for (path in c(bcf, vcf)) {
    stopifnot(identical(system2("bcftools", c("view", "-H", path), stdout = TRUE), expected))
  }
  header <- lines[startsWith(lines, "#")]
  for (kind in c("ps_type", "ps_number", "ps_width", "gt_allele")) {
    selected <- header
    record <- "chrG\t10\tbad\tA\tC\t.\tPASS\t.\tGT:PS\t0|1:10\t1|1:20"
    if (kind == "ps_type") {
      selected <- sub("ID=PS,Number=1,Type=Integer", "ID=PS,Number=1,Type=String", selected, fixed = TRUE)
    } else if (kind == "ps_number") {
      selected <- sub("ID=PS,Number=1,Type=Integer", "ID=PS,Number=.,Type=Integer", selected, fixed = TRUE)
    } else if (kind == "ps_width") {
      record <- sub("0|1:10", "0|1:10,11", record, fixed = TRUE)
    } else record <- sub("0|1:10", "0|2:10", record, fixed = TRUE)
    text <- paste0("test/data/geno_", kind, ".vcf")
    writeLines(c(selected, record), text)
    binary <- sub("vcf$", "bcf", text)
    run(c("view", "--no-version", "-Ob", "-o", binary, text))
    outputs <- c(outputs, text, binary)
  }
  stopifnot(all(file.copy(outputs, "r/Rduckhts/inst/extdata", overwrite = TRUE)))
  # Explicit VCF 4.4 phase-prefix witnesses stay textual: the locally installed
  # bcftools encoder may implement an older VCF version than bundled HTSlib.
  stopifnot(file.copy("test/data/geno_vcf44.vcf", "r/Rduckhts/inst/extdata", overwrite = TRUE))
}
prepare_geno_fixtures()
