# Maintainer-only, network-free reconstruction of the fixtures of
# duckhts_bam_mismatch_counts. The reference and the SAM records built below
# are the source authority. The samtools, bcftools and tabix binaries bundled
# by the RBCFTools package supply the BAM, CRAM, VCF and BCF encodings and
# their indexes; tools on PATH are not used. Run from the repository root.
# Generated with RBCFTools 1.24-1.1.1.9000: samtools 1.24, bcftools 1.24, HTSlib 1.24.
stopifnot(file.exists("src/bam_mismatch_counts.c"), requireNamespace("RBCFTools", quietly = TRUE),
          utils::packageVersion("RBCFTools") >= "1.24.1.1.0")
tools <- c(samtools = RBCFTools::samtools_path(), bcftools = RBCFTools::bcftools_path(),
           tabix = RBCFTools::tabix_path())
stopifnot(all(file.exists(tools)))

data_dir <- "test/data"
package_dir <- "r/Rduckhts/inst/extdata"
bases <- c("A", "C", "G", "T")

# A fixed sequence: no random numbers, so every R version writes the same bases.
contig_bases <- function(length, shift) {
  index <- seq_len(length)
  bases[((index * 7L + (index %/% 3L) * 5L + shift) %% 4L) + 1L]
}
reference <- list(ref1 = contig_bases(260L, 0L), ref2 = contig_bases(80L, 1L))

# The stored sequence of an alignment that starts at pos1: reference bases for
# M, the given bases for I and S, nothing for D and H. `substitutions` are
# zero-based offsets in the stored sequence whose base is replaced by the next
# base of A, C, G, T, so the read differs from the reference there.
aligned_sequence <- function(contig, pos1, operations, substitutions = integer(), unknown = integer()) {
  sequence <- character()
  at <- pos1
  for (operation in operations) {
    code <- operation[[1L]]
    length <- as.integer(operation[[2L]])
    if (code == "M") {
      sequence <- c(sequence, reference[[contig]][seq.int(at, length.out = length)])
      at <- at + length
    } else if (code %in% c("I", "S")) {
      sequence <- c(sequence, rep_len(c("T", "G"), length))
    } else if (code == "D") {
      at <- at + length
    }
  }
  for (offset in substitutions) {
    sequence[[offset + 1L]] <- bases[(match(sequence[[offset + 1L]], bases) %% 4L) + 1L]
  }
  sequence[unknown + 1L] <- "N"
  paste(sequence, collapse = "")
}

cigar_text <- function(operations) paste(vapply(operations, function(o) paste0(o[[2L]], o[[1L]]), ""), collapse = "")
sam_line <- function(name, flag, contig, pos1, mapq, operations, quality, mate = c("*", "0", "0"), ...) {
  sequence <- aligned_sequence(contig, pos1, operations, ...)
  qualities <- if (is.na(quality)) "*" else strrep(intToUtf8(quality + 33L), nchar(sequence))
  paste(name, flag, contig, pos1, mapq, cigar_text(operations), mate[[1L]], mate[[2L]], mate[[3L]],
        sequence, qualities, sep = "\t")
}
m <- function(length) list("M", length)

records <- c(
  # A proper pair: the first mate matches, the second (reverse) has one mismatch.
  sam_line("pair1", 99L, "ref1", 11L, 60L, list(m(20L)), 40L, mate = c("=", "41", "50")),
  sam_line("pair1", 147L, "ref1", 41L, 60L, list(m(20L)), 35L, mate = c("=", "11", "-50"), substitutions = 3L),
  # A soft clip and an insertion; one mismatch in the last three aligned bases.
  sam_line("clip_ins", 0L, "ref1", 71L, 60L, list(list("S", 5L), m(10L), list("I", 2L), m(8L)), 30L,
           substitutions = 23L),
  # A deletion on the reverse strand; one mismatch next to it and one far from it.
  sam_line("deletion", 16L, "ref1", 101L, 60L, list(m(10L), list("D", 3L), m(10L)), 20L,
           substitutions = c(8L, 17L)),
  # Hard clips shift the cycle: the stored bases are cycles 4 to 15.
  sam_line("hardclip", 0L, "ref1", 131L, 60L, list(list("H", 3L), m(12L), list("H", 2L)), 25L),
  # No base qualities.
  sam_line("noqual", 0L, "ref1", 151L, 60L, list(m(12L)), NA),
  # Over mask records: an SNV at 173, a deletion at 180-181 and a symbolic
  # deletion at 186-188. The read has an N, a mismatch at the masked SNV and a
  # mismatch at position 190, which only the default flank masks.
  sam_line("masked", 0L, "ref1", 171L, 60L, list(m(20L)), 10L, substitutions = c(2L, 19L), unknown = 0L),
  # Left out by the default filters: a duplicate and a low mapping quality.
  sam_line("duplicate", 1024L, "ref1", 201L, 60L, list(m(20L)), 40L),
  sam_line("lowmapq", 0L, "ref1", 221L, 5L, list(m(20L)), 40L),
  # A second contig, known to the mask header and without mask records.
  sam_line("other_contig", 0L, "ref2", 11L, 60L, list(m(20L)), 40L, substitutions = 0L),
  paste("unmapped", 4L, "*", 0L, 0L, "*", "*", 0L, 0L, "ACGT", "IIII", sep = "\t"))

prefix <- file.path(data_dir, "bam_mismatch")
fasta <- paste0(prefix, ".fa")
writeLines(unlist(lapply(names(reference), function(contig) {
  c(paste0(">", contig), paste(reference[[contig]], collapse = ""))
})), fasta)
run <- function(tool, arguments) {
  status <- system2(tools[[tool]], shQuote(arguments))
  if (status != 0L) stop(tool, " failed: ", paste(arguments, collapse = " "), call. = FALSE)
}
run("samtools", c("faidx", fasta))

sam <- tempfile(fileext = ".sam")
writeLines(c("@HD\tVN:1.6\tSO:coordinate",
             sprintf("@SQ\tSN:%s\tLN:%d", names(reference), lengths(reference)), records), sam)
bam <- paste0(prefix, ".bam")
run("samtools", c("view", "--no-PG", "-b", "-o", bam, sam))
run("samtools", c("index", bam))
cram <- paste0(prefix, ".cram")
run("samtools", c("view", "--no-PG", "-C", "-T", fasta, "-o", cram, sam))
run("samtools", c("index", cram))

# The mask. `contigs` are the header contig lines; `names` rename the contig
# of the records, for the fixture whose names differ from the alignment's.
write_mask <- function(output, contigs, record_contig) {
  ref_at <- function(pos1, length) paste(reference$ref1[seq.int(pos1, length.out = length)], collapse = "")
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    sprintf("##contig=<ID=%s,length=%d>", contigs, lengths(reference)),
    "##ALT=<ID=DEL,Description=\"Deletion\">",
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    paste(record_contig, 173L, ".", ref_at(173L, 1L), "T", ".", ".", ".", sep = "\t"),
    paste(record_contig, 180L, ".", ref_at(180L, 2L), ref_at(180L, 1L), ".", ".", ".", sep = "\t"),
    paste(record_contig, 186L, ".", ref_at(186L, 1L), "<DEL>", ".", ".", "END=188", sep = "\t")), vcf)
  run("bcftools", c("view", "--no-version", "-Oz", "-o", output, vcf))
  run("tabix", c("-f", "-p", "vcf", output))
}
mask <- paste0(prefix, ".mask.vcf.gz")
write_mask(mask, names(reference), "ref1")
run("bcftools", c("view", "--no-version", "-Ob", "-o", paste0(prefix, ".mask.bcf"), mask))
run("bcftools", c("index", "-f", paste0(prefix, ".mask.bcf")))
write_mask(paste0(prefix, ".other_names.vcf.gz"), paste0("chr", names(reference)), "chrref1")

# The R package tests read the BAM, the reference and the BCF mask.
package_files <- paste0(prefix, c(".fa", ".fa.fai", ".bam", ".bam.bai", ".mask.bcf", ".mask.bcf.csi"))
stopifnot(all(file.copy(package_files, package_dir, overwrite = TRUE)))
