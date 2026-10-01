# Stream a pedigree-selected, biallelic-SNV BCF for one autosome.
stage_roh_children <- function(chromosome, source_vcf, pedigree, output,
                               bcftools = "/usr/local/bin/bcftools",
                               expected_source_sha256 = NULL) {
  if (length(chromosome) != 1L || !is.numeric(chromosome) ||
      !is.finite(chromosome) || chromosome != as.integer(chromosome) ||
      chromosome < 1L || chromosome > 22L) {
    stop("`chromosome` must be an autosome number from 1 through 22", call. = FALSE)
  }
  if (length(source_vcf) != 1L || !nzchar(source_vcf) ||
      length(pedigree) != 1L || !file.exists(pedigree) ||
      length(output) != 1L || !nzchar(output) ||
      length(bcftools) != 1L || !file.exists(bcftools)) {
    stop("source, pedigree, output, and bcftools must be supplied", call. = FALSE)
  }
  expected <- c(ACB = 20L, ASW = 13L, CLM = 35L, MXL = 32L, PEL = 35L,
                PUR = 35L, YRI = 56L, ESN = 43L, CEU = 57L, CHS = 51L)
  pedigree_data <- read.table(pedigree, header = TRUE, sep = "", quote = "",
                              comment.char = "", colClasses = "character",
                              stringsAsFactors = FALSE)
  children <- pedigree_data[
    pedigree_data$FatherID != "0" & pedigree_data$MotherID != "0" &
      pedigree_data$Population %in% names(expected), , drop = FALSE]
  observed <- table(factor(children$Population, levels = names(expected)))
  if (!identical(as.integer(observed), unname(expected))) {
    stop("pedigree child counts differ from the registered expectations", call. = FALSE)
  }
  source_hash <- if (file.exists(source_vcf)) {
    unname(digest::digest(file = source_vcf, algo = "sha256"))
  } else NA_character_
  if (!is.null(expected_source_sha256) &&
      !identical(source_hash, expected_source_sha256)) {
    stop("local VCF checksum differs from the expected source identity", call. = FALSE)
  }
  sample_file <- paste0(output, ".samples.txt")
  temporary <- paste0(output, ".partial-", Sys.getpid(), ".bcf")
  temporary_index <- paste0(temporary, ".csi")
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(c(temporary, temporary_index, sample_file), force = TRUE), add = TRUE)
  utils::write.table(children$SampleID, sample_file, row.names = FALSE,
                     col.names = FALSE, quote = FALSE)
  query <- function(arguments) {
    value <- system2(bcftools, shQuote(arguments), stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(value, "status"))) {
      stop("bcftools command failed: ", paste(arguments, collapse = " "),
           call. = FALSE)
    }
    value
  }
  header_samples <- query(c("query", "-l", source_vcf))
  absent <- setdiff(children$SampleID, header_samples)
  if (length(absent)) {
    stop("eligible children absent from VCF header: ", paste(absent, collapse = ", "),
         call. = FALSE)
  }
  command <- paste(
    shQuote(bcftools), "view --force-samples -S", shQuote(sample_file),
    "-m2 -M2 -v snps -Ou", shQuote(source_vcf), "|",
    shQuote(bcftools), "view -c1 -Ob -o", shQuote(temporary))
  status <- system2("/bin/bash", c("-o", "pipefail", "-c", shQuote(command)))
  if (status != 0L) stop("bcftools failed to build the children BCF", call. = FALSE)
  status <- system2(bcftools, shQuote(c("index", "-f", temporary)))
  if (status != 0L || !file.exists(temporary_index)) {
    stop("bcftools failed to index the children BCF", call. = FALSE)
  }
  output_samples <- query(c("query", "-l", temporary))
  if (!identical(output_samples, children$SampleID)) {
    stop("children BCF sample order differs from the pedigree selection", call. = FALSE)
  }
  count_text <- system2(bcftools, shQuote(c("index", "-n", temporary)),
                       stdout = TRUE, stderr = FALSE)
  if (!is.null(attr(count_text, "status")) || length(count_text) != 1L) {
    stop("could not read the children BCF record count", call. = FALSE)
  }
  record_count <- as.integer(count_text)
  if (is.na(record_count)) stop("children BCF record count is not an integer", call. = FALSE)
  output_sha256 <- unname(digest::digest(file = temporary, algo = "sha256"))
  if (!file.rename(temporary, output)) {
    stop("could not publish the children BCF", call. = FALSE)
  }
  if (!file.rename(temporary_index, paste0(output, ".csi"))) {
    unlink(output, force = TRUE)
    stop("could not publish the children BCF index", call. = FALSE)
  }
  list(chromosome = as.integer(chromosome), source = source_vcf,
       source_sha256 = source_hash, source_samples = length(header_samples),
       samples = length(output_samples),
       records = record_count, output = output, output_sha256 = output_sha256)
}
