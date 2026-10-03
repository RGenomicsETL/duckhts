test_roh_ancestry_staging <- function() {
  bcftools <- Sys.which("bcftools")
  if (!nzchar(bcftools)) {
    stop("bcftools is required for the network-free staging test", call. = FALSE)
  }
  source("benchmarks/roh_ancestry_stage.R")
  directory <- tempfile("roh-ancestry-stage-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  populations <- c(ACB = 20L, ASW = 13L, CLM = 35L, MXL = 32L, PEL = 35L,
                   PUR = 35L, YRI = 56L, ESN = 43L, CEU = 57L, CHS = 51L)
  pedigree <- do.call(rbind, lapply(names(populations), function(population) {
    count <- populations[[population]]
    data.frame(FamilyID = rep(population, count),
      SampleID = sprintf("%s_%03d", population, seq_len(count)),
      FatherID = "father", MotherID = "mother", Sex = "0",
      Population = population, Superpopulation = "test",
      stringsAsFactors = FALSE)
  }))
  pedigree_path <- file.path(directory, "pedigree.txt")
  utils::write.table(pedigree, pedigree_path, sep = " ", row.names = FALSE,
                     quote = FALSE)
  vcf_path <- file.path(directory, "source.vcf")
  samples <- pedigree$SampleID
  header <- paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                    "INFO", "FORMAT", samples), collapse = "\t")
  record <- paste(c("chr20", "10", "rs1", "A", "C", ".", "PASS",
                    "AF=0.1;AF_AFR=0.2;AF_AMR=0.3;AF_EAS=0.4;AF_EUR=0.5",
                    "GT", rep("0|1", length(samples))), collapse = "\t")
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=chr20,length=100>",
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
    '##INFO=<ID=AF_AFR,Number=A,Type=Float,Description="AFR allele frequency">',
    '##INFO=<ID=AF_AMR,Number=A,Type=Float,Description="AMR allele frequency">',
    '##INFO=<ID=AF_EAS,Number=A,Type=Float,Description="EAS allele frequency">',
    '##INFO=<ID=AF_EUR,Number=A,Type=Float,Description="EUR allele frequency">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', header, record),
    vcf_path)
  output <- file.path(directory, "children.bcf")
  source_hash <- unname(digest::digest(file = vcf_path, algo = "sha256"))
  registry <- data.frame(
    id = c("roh_ancestry_chr20_source", "roh_ancestry_chr20_children_bcf"),
    locator = c(vcf_path, "derived:source+pedigree"),
    cache_relpath = c("source.vcf", "children.bcf"),
    stringsAsFactors = FALSE)
  result <- stage_roh_children_from_registry(
    20L, registry, pedigree_path, output, bcftools = bcftools,
    expected_source_sha256 = source_hash)
  stopifnot(result$records == 1L, result$source_samples == 377L,
            result$samples == 377L,
    identical(result$output_sha256, unname(digest::digest(file = output, algo = "sha256"))),
    identical(system2(bcftools, shQuote(c("query", "-l", output)),
                      stdout = TRUE), samples),
    identical(system2(bcftools, shQuote(c("query", "-f", "%INFO/AF\\t%INFO/AF_AFR\\n",
                                          output)), stdout = TRUE), "0.1\t0.2"),
    file.exists(paste0(output, ".csi")),
    !file.exists(paste0(output, ".samples.txt")))

  # The registered record and sample counts are checked before publishing: a
  # wrong count is an error that leaves nothing at the destination, and the right
  # counts pass.
  registered <- registry
  registered$supplier_identity <- c("", "records=2;samples=377;sha256=0000")
  rejected <- file.path(directory, "rejected.bcf")
  failed <- tryCatch({
    stage_roh_children_from_registry(20L, registered, pedigree_path, rejected,
                                     bcftools = bcftools,
                                     expected_source_sha256 = source_hash)
    FALSE
  }, error = function(error) grepl("differs from its registered identity", conditionMessage(error)))
  stopifnot(failed, !file.exists(rejected), !file.exists(paste0(rejected, ".csi")))
  registered$supplier_identity[[2L]] <- paste0("records=1;samples=377;sha256=", result$output_sha256)
  # The SHA-256 is the evaluated file's identity, not a reproducibility check.
  accepted <- file.path(directory, "accepted.bcf")
  verified <- stage_roh_children_from_registry(20L, registered, pedigree_path, accepted,
                                               bcftools = bcftools,
                                               expected_source_sha256 = source_hash)
  stopifnot(verified$records == 1L, verified$samples == 377L, file.exists(accepted))
}

test_roh_ancestry_staging()
