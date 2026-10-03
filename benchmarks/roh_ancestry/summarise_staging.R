chromosomes <- 1:22
receipts <- lapply(chromosomes, function(chromosome) {
  path <- file.path("benchmarks/results/roh-ancestry", sprintf("chr%d", chromosome),
                    "input_receipt.tsv")
  if (!file.exists(path)) stop("missing receipt: ", path, call. = FALSE)
  receipt <- utils::read.delim(path, stringsAsFactors = FALSE, quote = "",
                               check.names = FALSE, na.strings = "")
  source <- receipt[grepl("source_vcf", receipt$artifact), , drop = FALSE]
  output <- receipt[grepl("children_bcf", receipt$artifact), , drop = FALSE]
  if (nrow(source) != 1L || nrow(output) != 1L) {
    stop("receipt lacks one source and one children BCF row: ", path, call. = FALSE)
  }
  data.frame(chromosome = chromosome, source_url = source$locator,
    source_last_modified = source$remote_last_modified, source_etag = source$remote_etag,
    source_samples = source$samples, child_samples = output$samples,
    child_records = output$records, output_sha256 = output$sha256,
    stringsAsFactors = FALSE)
})
result <- do.call(rbind, receipts)
result <- result[order(result$chromosome), ]
utils::write.csv(result,
  "benchmarks/results/roh-ancestry/autosome_staging_summary.csv", row.names = FALSE)
print(result)
