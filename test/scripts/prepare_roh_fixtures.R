# Regenerates the runs-of-homozygosity fixtures in test/data (roh_fixture.vcf.gz,
# roh_af.tsv.gz, roh_map_chr1.txt, roh_map_chr2.txt and their indexes).
# Run from the repository root with bgzip and tabix on PATH:
#   Rscript test/scripts/prepare_roh_fixtures.R
# The output is deterministic. The expected segments in test/sql/roh.test come
# from test/scripts/roh_bcftools_expected.sh, which runs `bcftools roh` on these files.

out_dir <- "test/data"
set.seed(318)

contigs <- data.frame(name = c("chr1", "chr2"), n = c(330L, 240L), length = c(4500000L, 3500000L))
samples <- c("S1", "S2", "S3", "S4")

# Planted homozygous stretches: sample -> contig -> site index ranges.
planted <- list(
  S1 = list(chr1 = list(120:230)),
  S2 = list(chr1 = list(40:75), chr2 = list(60:170)),
  S3 = list(),
  S4 = list(chr1 = list(200:330), chr2 = list(1:240))
)
# Sites with a record-level anomaly, by contig and site index.
missing_af_site <- list(chr1 = 150L, chr2 = 100L)
zero_af_site <- list(chr1 = 160L)
# Sites left out of roh_af.tsv.gz, so --AF-file sees fewer sites than --AF-tag.
af_file_omitted <- list(chr1 = c(100L, 101L), chr2 = 10L)
# Gaps above 100 kb, to exercise genetic-map interpolation across several nodes.
big_gaps <- list(chr1 = c(60L, 135L, 210L, 290L), chr2 = c(30L, 120L, 200L))

in_planted <- function(sample, contig, index) {
  ranges <- planted[[sample]][[contig]]
  any(vapply(ranges, function(r) index %in% r, logical(1)))
}

phred_likelihoods <- function(called) {
  near <- sample(18:45, 1)
  far <- near + sample(10:40, 1)
  switch(as.character(called),
    "0" = c(0L, near, far),
    "1" = c(near, 0L, far),
    "2" = c(far, near, 0L))
}

vcf_lines <- c(
  "##fileformat=VCFv4.2",
  paste0("##contig=<ID=", contigs$name, ",length=", contigs$length, ">"),
  "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate allele frequency\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">",
  paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples),
        collapse = "\t"))
af_lines <- character()
map_lines <- list()

gt_text <- function(dosage) c("0/0", "0/1", "1/1")[dosage + 1L]
pl_text <- function(pl) paste(pl, collapse = ",")

for (c_index in seq_len(nrow(contigs))) {
  contig <- contigs$name[c_index]
  n <- contigs$n[c_index]
  gaps <- sample(2000:12000, n, replace = TRUE)
  gaps[big_gaps[[contig]]] <- sample(150000:350000, length(big_gaps[[contig]]))
  positions <- 20000L + cumsum(as.integer(gaps))
  stopifnot(max(positions) < contigs$length[c_index] - 50000L)
  freq <- round(runif(n, 0.15, 0.5), 3)
  ref_alleles <- sample(c("A", "C", "G", "T"), n, replace = TRUE)
  alt_alleles <- vapply(ref_alleles, function(r) sample(setdiff(c("A", "C", "G", "T"), r), 1), "")

  for (i in seq_len(n)) {
    pos <- positions[i]
    ref <- ref_alleles[i]
    alt <- alt_alleles[i]
    info <- paste0("AF=", formatC(freq[i], format = "f", digits = 3))
    if (identical(missing_af_site[[contig]], i)) info <- "."
    if (identical(zero_af_site[[contig]], i)) info <- "AF=0.000"

    calls <- character(length(samples))
    for (s in seq_along(samples)) {
      sample_name <- samples[s]
      if (in_planted(sample_name, contig, i)) {
        allele <- rbinom(1, 1, freq[i])
        dosage <- 2L * allele
        if (runif(1) < 0.04) dosage <- 1L # genotyping error inside the run
      } else {
        dosage <- as.integer(sum(rbinom(2, 1, freq[i])))
      }
      pl <- phred_likelihoods(dosage)
      if (runif(1) < 0.05) pl <- c(0L, 3L, 9L) # poorly resolved site
      calls[s] <- paste0(gt_text(dosage), ":", pl_text(pl))
    }

    if (contig == "chr1" && i == 170L) {
      calls[1] <- "./.:0,50,100" # S1: GT missing, PL informative
      calls[2] <- "0/1:0,0,0" # S2: uninformative PL
      calls[3] <- "./.:." # S3: nothing
      calls[4] <- "1/1:300,200,0" # S4: PL above the 255 cap
    }
    if (contig == "chr1" && i == 190L) { # biallelic indel
      ref <- "AC"
      alt <- "A"
    }
    vcf_lines <- c(vcf_lines, paste(c(contig, pos, ".", ref, alt, ".", "PASS", info,
                                      "GT:PL", calls), collapse = "\t"))
    if (info != "." && !(i %in% af_file_omitted[[contig]])) {
      af_file_value <- round(freq[i] * 0.8 + 0.05, 3)
      if (identical(zero_af_site[[contig]], i)) af_file_value <- 0
      af_lines <- c(af_lines, paste(contig, pos, paste0(ref, ",", alt),
                                    formatC(af_file_value, format = "f", digits = 3), sep = "\t"))
    }

    if (contig == "chr1" && i == 175L) { # multiallelic record, skipped by bcftools roh
      vcf_lines <- c(vcf_lines, paste(c(contig, pos + 1L, ".", "A", "C,G", ".", "PASS", "AF=0.300,0.200",
                                        "GT:PL", rep("0/1:10,0,20,30,40,50", length(samples))),
                                      collapse = "\t"))
    }
    if (contig == "chr1" && i == 180L) { # record without ALT, skipped by default
      vcf_lines <- c(vcf_lines, paste(c(contig, pos + 1L, ".", "A", ".", ".", "PASS", ".",
                                        "GT", rep("0/0", length(samples))), collapse = "\t"))
    }
  }

  # Sparse IMPUTE2-style map (position, rate, cM). Nodes cover most of the contig; the tail is
  # left uncovered so the last sites fall beyond the map.
  node_positions <- seq(10000L, max(positions) - 40000L, by = 90000L)
  rates <- round(runif(length(node_positions), 0.2, 3.0), 4)
  cm <- round(cumsum(c(0, diff(node_positions) / 1e6 * rates[-1])), 6)
  map_lines[[contig]] <- c("position COMBINED_rate(cM/Mb) Genetic_Map(cM)",
                           paste(node_positions, rates, formatC(cm, format = "f", digits = 6)))
}

vcf_path <- file.path(out_dir, "roh_fixture.vcf")
writeLines(vcf_lines, vcf_path)
system2("bgzip", c("-f", vcf_path))
system2("tabix", c("-f", "-p", "vcf", paste0(vcf_path, ".gz")))

af_path <- file.path(out_dir, "roh_af.tsv")
writeLines(af_lines, af_path)
system2("bgzip", c("-f", af_path))
system2("tabix", c("-f", "-s1", "-b2", "-e2", paste0(af_path, ".gz")))

for (contig in names(map_lines)) {
  writeLines(map_lines[[contig]], file.path(out_dir, paste0("roh_map_", contig, ".txt")))
}

# Degenerate inputs: no records, and a single site.
small_header <- c(
  "##fileformat=VCFv4.2",
  "##contig=<ID=chr1,length=100000>",
  "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate allele frequency\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">",
  paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"),
        collapse = "\t"))
writeLines(small_header, file.path(out_dir, "roh_empty.vcf"))
writeLines(c(small_header,
             paste(c("chr1", "5000", ".", "A", "C", ".", "PASS", "AF=0.050", "GT:PL", "1/1:60,30,0"),
                   collapse = "\t")),
           file.path(out_dir, "roh_single_site.vcf"))

# GT-only VCF (no FORMAT/PL): one homozygous stretch between heterozygous flanks.
gt_only_header <- c(
  "##fileformat=VCFv4.2",
  "##contig=<ID=chr1,length=200000>",
  "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate allele frequency\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1", "S2"),
        collapse = "\t"))
gt_only_calls <- function(i) {
  s1 <- if (i >= 20 && i <= 70) c("0/0", "1/1")[1 + i %% 2] else "0/1"
  s2 <- if (i %% 7 == 0) "./." else "0/1"
  c(s1, s2)
}
writeLines(c(gt_only_header,
             vapply(1:90, function(i) {
               paste(c("chr1", 1000L * i, ".", "G", "T", ".", "PASS", "AF=0.300", "GT", gt_only_calls(i)),
                     collapse = "\t")
             }, "")),
           file.path(out_dir, "roh_gt_only.vcf"))
