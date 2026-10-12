# Regenerates test/data/roh_repeated.vcf.gz and roh_repeated_af.tsv.gz (with
# their indexes), the fixture for duckhts_roh_af_table(site_rule := ...) in
# test/sql/roh.test, and copies them to r/Rduckhts/inst/extdata.
# Run from the repository root: Rscript test/scripts/prepare_roh_repeated_fixture.R
# It uses bgzip and tabix from the RBCFTools package (HTSlib 1.24). The output
# is deterministic.
#
# Positions 500 to 10000 hold repeated records that bcftools roh --AF-file and
# the matching rule resolve differently. Every record there is homozygous for
# REF with strong PL, and a heterozygous site with strong PL lies between two
# such positions, so that with fast transitions each used position is a run of
# one marker and the runs name the sites used.

bgzip <- RBCFTools::bgzip_path()
tabix <- RBCFTools::tabix_path()
out_dir <- "test/data"
package_dir <- "r/Rduckhts/inst/extdata"
stopifnot(dir.exists(out_dir), dir.exists(package_dir))

hom <- "0/0:0,60,120"
het <- "0/1:60,0,60"
hom3 <- "0/0:0,60,120,60,120,120" # three alleles
record <- function(pos, id, ref, alt, sample) {
  paste("chr1", pos, id, ref, alt, "50", "PASS", ".", "GT:PL", sample, sep = "\t")
}
# Records at the test positions, in file order, and the frequency rows of the
# position (in file order; 4000 has two).
tests <- list(
  list(500, list(c("flank1", "C", "T")), list(c("C,T", "0.3"))),
  list(1000, list(c("A_indel", "A", "AT"), c("A_snp", "A", "G")), list(c("A,G", "0.3"))),
  list(2000, list(c("B_snp_T", "A", "T"), c("B_snp_G", "A", "G")), list(c("A,G", "0.3"))),
  list(3000, list(c("C_ins1", "A", "AT"), c("C_ins2", "A", "ATT")), list(c("A,G", "0.3"))),
  list(4000, list(c("D_snp", "A", "G")), list(c("A,T", "0.2"), c("A,G", "0.3"))),
  list(5000, list(c("E_multi", "A", "G,T"), c("E_snp", "A", "G")), list(c("A,G", "0.3"))),
  list(6000, list(c("F_snp", "A", "G"), c("F_ins", "A", "AT")), list(c("A,AT", "0.3"))),
  list(7000, list(c("G_mnp", "AC", "GT"), c("G_snp", "A", "G")), list(c("A,G", "0.3"))),
  list(8000, list(c("H_other", "ACG", "AT"), c("H_snp", "A", "G")), list(c("A,G", "0.3"))),
  list(9000, list(c("flank2", "C", "T")), list(c("C,T", "0.3"))),
  list(10000, list(c("I_symbolic", "A", "<*>,G"), c("I_snp", "A", "G")), list(c("A,G", "0.3"))),
  # Not in the frequency file: roh.test lists it with alleles A,<*>,G, which
  # bcftools would decode through the second ALT and DuckHTS skips.
  list(11000, list(c("J_symbolic", "A", "<*>,G")), list())
)

lines <- c(
  "##fileformat=VCFv4.2",
  "##contig=<ID=chr1,length=100000>",
  "##ALT=<ID=*,Description=\"Any other allele\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"
)
frequencies <- character()
for (k in seq_along(tests)) {
  test <- tests[[k]]
  pos <- test[[1]]
  if (k > 1) {
    separator <- pos - if (pos <= 1000) 250 else 500
    lines <- c(lines, record(separator, paste0("sep", separator), "C", "T", het))
    frequencies <- c(frequencies, paste("chr1", separator, "C,T", "0.5", sep = "\t"))
  }
  for (r in test[[2]]) {
    lines <- c(lines, record(pos, r[1], r[2], r[3], if (grepl(",", r[3])) hom3 else hom))
  }
  for (f in test[[3]]) frequencies <- c(frequencies, paste("chr1", pos, f[1], f[2], sep = "\t"))
}

text <- tempfile(fileext = ".vcf")
writeLines(lines, text)
vcf <- file.path(out_dir, "roh_repeated.vcf.gz")
stopifnot(system2(bgzip, c("-c", shQuote(text)), stdout = vcf) == 0L,
          system2(tabix, c("-f", "-p", "vcf", shQuote(vcf))) == 0L)
af_text <- tempfile(fileext = ".tsv")
writeLines(frequencies, af_text)
af <- file.path(out_dir, "roh_repeated_af.tsv.gz")
stopifnot(system2(bgzip, c("-c", shQuote(af_text)), stdout = af) == 0L,
          system2(tabix, c("-f", "-s1", "-b2", "-e2", shQuote(af))) == 0L)
unlink(c(text, af_text))
outputs <- c(vcf, paste0(vcf, ".tbi"), af, paste0(af, ".tbi"))
stopifnot(all(file.copy(outputs, package_dir, overwrite = TRUE)))
