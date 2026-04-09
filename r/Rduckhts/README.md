
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Rduckhts: DuckDB HTS File Reader Extension for R

[![CRAN
Status](https://www.r-pkg.org/badges/version/Rduckhts)](https://cran.r-project.org/package=Rduckhts)[![R-universe
version](https://RGenomicsETL.r-universe.dev/Rduckhts/badges/version)](https://RGenomicsETL.r-universe.dev/Rduckhts)

`Rduckhts` provides an R interface to a [DuckDB](https://duckdb.org/)
`HTS` (High Throughput Sequencing) file reader extension. This enables
reading common bioinformatics file formats such as `VCF`/`BCF`,
`SAM`/`BAM`/`CRAM`, `FASTA`, `FASTQ`, `GFF`, `GTF`, and tabix-indexed
files directly from `R` using `SQL` queries via
[`duckhts`](https://github.com/RGenomicsETL/duckhts).

## How it works

Following [RBCFTools](https://github.com/RGenomicsETL/RBCFTools), tables
are created and returned instead of data frames. `VCF`/`BCF`,
`SAM`/`BAM`/`CRAM`, `FASTA`, `FASTQ`, `GFF`, `GTF`, and `tabix` formats
can be queried. We support region queries for indexed files, and we
target Linux, macOS, and RTools.
[`htslib`](https://github.com/samtools/htslib) 1.23 is bundled so build
dependencies stay minimal. The extensnion is built by adapting the
generic extension infracstructure by using only makefiles unlike the
submitted communtity extension
[`duckhts`](https://github.com/RGenomicsETL/duckhts).

## Installation

The package can be installed from r-universe

``` r
# Install 'Rduckhts' in R:
install.packages('Rduckhts', repos = c('https://rgenomicsetl.r-universe.dev', 'https://cloud.r-project.org'))
# When on CRAN
install.packages("Rduckhts")
```

## System Requirements

Installation requires `htslib` dependencies such ad zlib and libbz2, and
optionally for full functionally liblzma, libcurl, and openssl. The
package requires GNU make. On Windows’s Rtools, `htslib` plugins are not
enabled.

## Quick Start

The extension is loaded with
`rduckhts_load(con, extension_path = NULL)`. The wrappers break down
into:

- readers: `rduckhts_bcf()`, `rduckhts_bam()`, `rduckhts_fasta()`,
  `rduckhts_fastq()`, `rduckhts_gff()`, `rduckhts_gtf()`,
  `rduckhts_tabix()`, `rduckhts_bed()`
- reference helpers: `rduckhts_fasta_index()`, `rduckhts_fasta_nuc()`
- compression/indexing: `rduckhts_bgzip()`, `rduckhts_bgunzip()`,
  `rduckhts_bam_index()`, `rduckhts_bcf_index()`,
  `rduckhts_tabix_index()`
- metadata helpers: `rduckhts_hts_header()`, `rduckhts_hts_index()`,
  `rduckhts_hts_index_spans()`, `rduckhts_hts_index_raw()`

Start with one reader, then materialize tables and compose the richer
helpers around them.

``` r
library(DBI)
library(duckdb)
library(Rduckhts)


fasta_path <- system.file("extdata", "ce.fa", package = "Rduckhts")
fastq_r1 <- system.file("extdata", "r1.fq", package = "Rduckhts")
fastq_r2 <- system.file("extdata", "r2.fq", package = "Rduckhts")
con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
rduckhts_load(con)
#> [1] TRUE

rduckhts_fasta(con, "sequences", fasta_path, overwrite = TRUE)
rduckhts_fastq(con, "reads", fastq_r1, mate_path = fastq_r2, overwrite = TRUE)

dbGetQuery(con, "SELECT COUNT(*) AS n FROM sequences")
#>   n
#> 1 7
dbGetQuery(con, "SELECT COUNT(*) AS n FROM reads")
#>    n
#> 1 10
```

## Function Catalog

Use `rduckhts_functions()` inside R to inspect the generated extension
catalog.

## Extension Function Catalog

This section is generated from `functions.yaml`.

### Readers

| Function      | Kind  | Returns | R helper               | Description                                                                                                                                                                                                                                                                                                                                                                                                     |
|---------------|-------|---------|------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `read_bcf`    | table | table   | `rduckhts_bcf`         | Read VCF and BCF variant data with typed INFO, FORMAT, typed CSQ/ANN/BCSQ subfields, optional tidy sample output, and optional bcftools-style CSQ type overrides.                                                                                                                                                                                                                                               |
| `read_bam`    | table | table   | `rduckhts_bam`         | Read SAM, BAM, and CRAM alignments with optional typed SAMtags and auxiliary tag maps. Use sequence_encoding := ‘nt16’ to return SEQ as UTINYINT\[\] and quality_representation := ‘phred’ to return QUAL as UTINYINT\[\] instead of VARCHAR.                                                                                                                                                                   |
| `read_fasta`  | table | table   | `rduckhts_fasta`       | Read FASTA records or indexed FASTA regions as sequence rows. Use sequence_encoding := ‘nt16’ to return SEQUENCE as UTINYINT\[\] (htslib nt16 4-bit codes) instead of VARCHAR.                                                                                                                                                                                                                                  |
| `read_bed`    | table | table   | `rduckhts_bed`         | Read BED3-BED12 interval files with canonical typed columns and optional tabix-backed region filtering.                                                                                                                                                                                                                                                                                                         |
| `fasta_nuc`   | table | table   | `rduckhts_fasta_nuc`   | Compute bedtools nuc-style nucleotide composition for supplied BED intervals or generated fixed-width bins over a FASTA reference.                                                                                                                                                                                                                                                                              |
| `read_fastq`  | table | table   | `rduckhts_fastq`       | Read single-end, paired-end, or interleaved FASTQ files with optional legacy quality decoding. By default, FASTQ qualities are interpreted as modern Phred+33 input. Use sequence_encoding := ‘nt16’ to return SEQUENCE as UTINYINT\[\] and quality_representation := ‘phred’ to return QUALITY as UTINYINT\[\] instead of VARCHAR. input_quality_encoding accepts ‘phred33’, ‘auto’, ‘phred64’, or ‘solexa64’. |
| `read_gff`    | table | table   | `rduckhts_gff`         | Read GFF annotations with optional parsed attribute maps and indexed region filtering.                                                                                                                                                                                                                                                                                                                          |
| `read_gtf`    | table | table   | `rduckhts_gtf`         | Read GTF annotations with optional parsed attribute maps and indexed region filtering.                                                                                                                                                                                                                                                                                                                          |
| `read_tabix`  | table | table   | `rduckhts_tabix`       | Read generic tabix-indexed text data with optional header handling and type inference.                                                                                                                                                                                                                                                                                                                          |
| `fasta_index` | table | table   | `rduckhts_fasta_index` | Build a FASTA index (.fai) and return a single row with columns success (BOOLEAN) and index_path (VARCHAR).                                                                                                                                                                                                                                                                                                     |

### Metadata

| Function                  | Kind        | Returns | R helper                           | Description                                                                                                                   |
|---------------------------|-------------|---------|------------------------------------|-------------------------------------------------------------------------------------------------------------------------------|
| `detect_quality_encoding` | table       | table   | `rduckhts_detect_quality_encoding` | Inspect a FASTQ file’s observed quality ASCII range and report compatible legacy encodings with a heuristic guessed encoding. |
| `read_hts_header`         | table       | table   | `rduckhts_hts_header`              | Inspect HTS headers in parsed, raw, or combined form across supported formats.                                                |
| `read_hts_index`          | table       | table   | `rduckhts_hts_index`               | Inspect high-level HTS index metadata such as sequence names and mapped counts.                                               |
| `read_hts_index_spans`    | table       | table   | `rduckhts_hts_index_spans`         | Expand index metadata into span and chunk rows suitable for low-level index inspection.                                       |
| `read_hts_index_raw`      | table_macro | table   | `rduckhts_hts_index_raw`           | Return the raw on-disk HTS index blob together with basic identifying metadata.                                               |

### Compression

| Function  | Kind  | Returns | R helper           | Description                                                                           |
|-----------|-------|---------|--------------------|---------------------------------------------------------------------------------------|
| `bgzip`   | table | table   | `rduckhts_bgzip`   | Compress a plain file to BGZF and return the created output path and byte counts.     |
| `bgunzip` | table | table   | `rduckhts_bgunzip` | Decompress a BGZF-compressed file and return the created output path and byte counts. |

### Indexing

| Function      | Kind  | Returns | R helper               | Description                                                                                        |
|---------------|-------|---------|------------------------|----------------------------------------------------------------------------------------------------|
| `bam_index`   | table | table   | `rduckhts_bam_index`   | Build a BAM or CRAM index and report the written index path and format.                            |
| `bcf_index`   | table | table   | `rduckhts_bcf_index`   | Build a TBI or CSI index for a VCF or BCF file and report the written index path and format.       |
| `tabix_index` | table | table   | `rduckhts_tabix_index` | Build a tabix index for a BGZF-compressed text file using a preset or explicit coordinate columns. |

### Variants

| Function             | Kind        | Returns | R helper            | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
|----------------------|-------------|---------|---------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `bcftools_liftover`  | scalar      | STRUCT  | `rduckhts_liftover` | Row-oriented liftover kernel intended to mirror bcftools +liftover semantics as closely as possible while returning one STRUCT per input row with fields: src_chrom, src_pos, src_ref, src_alt, dest_chrom, dest_pos, dest_end, dest_ref, dest_alt, mapped, reverse_complemented, swap, reject_reason, and note. Set no_left_align := true to skip post-liftover left-alignment of lifted indels (mirrors –no-left-align in bcftools +liftover).                                           |
| `duckdb_liftover`    | table_macro | table   | `rduckhts_liftover` | DuckDB-specific wrapper over bcftools_liftover that takes either a table name or a derived-table expression plus column-name strings for chrom/pos/ref/alt and returns the lifted table. The no_left_align parameter mirrors –no-left-align in bcftools +liftover.                                                                                                                                                                                                                         |
| `bcftools_score`     | table       | table   | `rduckhts_score`    | Compute polygenic scores from one genotype BCF/VCF and one summary-statistics file with bcftools +score-compatible GT/DS/HDS/AP/GP/AS dosage semantics, sample subsetting, and region/target/FILTER-string controls.                                                                                                                                                                                                                                                                       |
| `bcftools_munge_row` | scalar      | STRUCT  |                     | Normalize one summary-statistics row into GWAS-VCF-style fields (chrom/pos/ref/alt/effect metrics), resolving REF/ALT orientation against a FASTA reference and applying swap-aware sign/frequency/count transforms. The output flag `alleles_swapped` means REF/ALT orientation was swapped to match the FASTA reference.                                                                                                                                                                 |
| `duckdb_munge`       | table_macro | table   | `rduckhts_munge`    | DuckDB macro wrapper over bcftools_munge_row that maps source columns (via preset or explicit map) and returns normalized GWAS-VCF-style rows with lean outputs and explicit `alleles_swapped` semantics. Output columns: chrom, pos, id, ref, alt, alleles_swapped, filter, ns, ez, nc, es, se, lp, af, ac, ne (16 columns). For METAL meta-analysis output with SI/I2/CQ/ED columns, use duckdb_munge_metal.                                                                             |
| `duckdb_munge_metal` | table_macro | table   | `rduckhts_munge`    | Extended munge macro with METAL meta-analysis output columns. Same as duckdb_munge but additionally emits: si (imputation info, from INFO input), i2 (Cochran’s I² heterogeneity, from HET_I2), cq (Cochran’s Q -log10 p, from HET_LP or -log10(HET_P)), and ed (effect direction string, from DIRE; +/- flipped on allele swap). The R wrapper rduckhts_munge() auto-dispatches to this macro when metal keys (INFO, HET_I2, HET_P, HET_LP, DIRE) are present in the resolved column map. |

### Sequence UDFs

| Function          | Kind   | Returns      | R helper | Description                                                                                           |
|-------------------|--------|--------------|----------|-------------------------------------------------------------------------------------------------------|
| `seq_revcomp`     | scalar | VARCHAR      |          | Compute the reverse complement of a DNA sequence using A, C, G, T, and N bases.                       |
| `seq_canonical`   | scalar | VARCHAR      |          | Return the lexicographically smaller of a sequence and its reverse complement.                        |
| `seq_hash_2bit`   | scalar | UBIGINT      |          | Encode a short DNA sequence as a 2-bit unsigned integer hash.                                         |
| `seq_encode_4bit` | scalar | UTINYINT\[\] |          | Encode an IUPAC DNA sequence as a list of 4-bit base codes, preserving ambiguity symbols including N. |
| `seq_decode_4bit` | scalar | VARCHAR      |          | Decode a list of 4-bit IUPAC DNA base codes back into a sequence string.                              |
| `seq_gc_content`  | scalar | DOUBLE       |          | Compute GC fraction for a DNA sequence as a value between 0 and 1.                                    |
| `seq_kmers`       | table  | table        |          | Expand a sequence into positional k-mers with optional canonicalization.                              |

### SAM Flag UDFs

| Function                               | Kind   | Returns | R helper | Description                                                                                                                                                                        |
|----------------------------------------|--------|---------|----------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `sam_flag_bits`                        | scalar | STRUCT  |          | Decode a SAM flag into a struct of boolean bit fields using explicit SAM-oriented names such as `is_paired`, `is_proper_pair`, `is_next_segment_unmapped`, and `is_supplementary`. |
| `sam_flag_has`                         | scalar | BOOLEAN |          | Test whether any bits from the provided SAM flag mask are set in a flag value.                                                                                                     |
| `is_forward_aligned`                   | scalar | BOOLEAN |          | Test whether a mapped segment is aligned to the forward strand. Returns `NULL` for unmapped segments because SAM flag `0x10` does not define genomic strand when `0x4` is set.     |
| `is_paired`                            | scalar | BOOLEAN |          | Test whether the SAM flag indicates that the template has multiple segments in sequencing (`0x1`).                                                                                 |
| `is_proper_pair`                       | scalar | BOOLEAN |          | Test whether the SAM flag indicates that each segment is properly aligned according to the aligner (`0x2`).                                                                        |
| `is_unmapped`                          | scalar | BOOLEAN |          | Test whether the read itself is unmapped according to the SAM flag.                                                                                                                |
| `is_next_segment_unmapped`             | scalar | BOOLEAN |          | Test whether the next segment in the template is flagged as unmapped (`0x8`).                                                                                                      |
| `is_reverse_complemented`              | scalar | BOOLEAN |          | Test whether `SEQ` is stored reverse complemented (`0x10`); for mapped reads this corresponds to reverse-strand alignment.                                                         |
| `is_next_segment_reverse_complemented` | scalar | BOOLEAN |          | Test whether `SEQ` of the next segment in the template is stored reverse complemented (`0x20`).                                                                                    |
| `is_first_segment`                     | scalar | BOOLEAN |          | Test whether the read is marked as the first segment in the template.                                                                                                              |
| `is_last_segment`                      | scalar | BOOLEAN |          | Test whether the read is marked as the last segment in the template.                                                                                                               |
| `is_secondary`                         | scalar | BOOLEAN |          | Test whether the alignment is marked as secondary.                                                                                                                                 |
| `is_qc_fail`                           | scalar | BOOLEAN |          | Test whether the read failed vendor or pipeline quality checks.                                                                                                                    |
| `is_duplicate`                         | scalar | BOOLEAN |          | Test whether the alignment is flagged as a duplicate.                                                                                                                              |
| `is_supplementary`                     | scalar | BOOLEAN |          | Test whether the alignment is marked as supplementary.                                                                                                                             |

### CIGAR Utils

| Function                     | Kind   | Returns | R helper | Description                                                                                                         |
|------------------------------|--------|---------|----------|---------------------------------------------------------------------------------------------------------------------|
| `cigar_has_soft_clip`        | scalar | BOOLEAN |          | Test whether a CIGAR string contains any soft-clipped segment (`S`).                                                |
| `cigar_has_hard_clip`        | scalar | BOOLEAN |          | Test whether a CIGAR string contains any hard-clipped segment (`H`).                                                |
| `cigar_left_soft_clip`       | scalar | BIGINT  |          | Return the left-end soft-clipped length from a CIGAR string, or zero if the alignment does not start with `S`.      |
| `cigar_right_soft_clip`      | scalar | BIGINT  |          | Return the right-end soft-clipped length from a CIGAR string, or zero if the alignment does not end with `S`.       |
| `cigar_query_length`         | scalar | BIGINT  |          | Return the query-consuming length from a CIGAR string, counting `M`, `I`, `S`, `=`, and `X`.                        |
| `cigar_aligned_query_length` | scalar | BIGINT  |          | Return the aligned query length from a CIGAR string, counting `M`, `=`, and `X` but excluding clips and insertions. |
| `cigar_reference_length`     | scalar | BIGINT  |          | Return the reference-consuming length from a CIGAR string, counting `M`, `D`, `N`, `=`, and `X`.                    |
| `cigar_has_op`               | scalar | BOOLEAN |          | Test whether a CIGAR string contains at least one instance of the requested operator.                               |

## Common Workflows

### Region-aware variant and alignment queries

``` r
bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
bcf_index_path <- system.file("extdata", "vcf_file.bcf.csi", package = "Rduckhts")
bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
bam_index_path <- system.file("extdata", "range.bam.bai", package = "Rduckhts")

rduckhts_bcf(
  con, "variants_idx", bcf_path,
  region = "1:3000150-3000151",
  index_path = bcf_index_path,
  overwrite = TRUE
)
dbGetQuery(con, "SELECT CHROM, POS, REF, ALT FROM variants_idx")
#>   CHROM     POS REF ALT
#> 1     1 3000150   C   T
#> 2     1 3000151   C   T

rduckhts_bam(
  con, "bam_idx_reads", bam_path,
  region = "CHROMOSOME_I:1-1000",
  index_path = bam_index_path,
  overwrite = TRUE
)
dbGetQuery(con, "SELECT QNAME, FLAG, POS, MAPQ FROM bam_idx_reads")
#>                           QNAME FLAG POS MAPQ
#> 1 HS18_09653:4:1315:19857:61712  145 914   23
#> 2 HS18_09653:4:1308:11522:27107  161 934    0
```

### Interval + reference helpers

``` r
bed_path <- system.file("extdata", "targets.bed", package = "Rduckhts")
fai_path <- tempfile("duckhts_readme_", fileext = ".fai")
rduckhts_fasta_index(con, fasta_path, index_path = fai_path)
#>   success                                        index_path
#> 1    TRUE /tmp/RtmpDtZah6/duckhts_readme_112d8f4646b7fb.fai

rduckhts_bed(con, "targets", bed_path, overwrite = TRUE)
dbGetQuery(con, "SELECT chrom, start, \"end\", name, block_count FROM targets")
#>            chrom start end    name block_count
#> 1   CHROMOSOME_I     0  10 target1           2
#> 2   CHROMOSOME_I    10  20 target2           1
#> 3  CHROMOSOME_II     0   8 target3          NA
#> 4 CHROMOSOME_III     0   6 target4           1

rduckhts_fasta_nuc(con, fasta_path, bed_path = bed_path, index_path = fai_path)
#>            chrom start end pct_at pct_gc num_a num_c num_g num_t num_n
#> 1   CHROMOSOME_I     0  10  0.400  0.600     2     4     2     2     0
#> 2   CHROMOSOME_I    10  20  0.500  0.500     4     3     2     1     0
#> 3  CHROMOSOME_II     0   8  0.375  0.625     2     4     1     1     0
#> 4 CHROMOSOME_III     0   6  0.500  0.500     2     2     1     1     0
#>   num_other seq_len
#> 1         0      10
#> 2         0      10
#> 3         0       8
#> 4         0       6
rduckhts_fasta_nuc(con, fasta_path, bin_width = 10, region = "CHROMOSOME_I:1-20", index_path = fai_path)
#>          chrom start end pct_at pct_gc num_a num_c num_g num_t num_n num_other
#> 1 CHROMOSOME_I     0  10    0.4    0.6     2     4     2     2     0         0
#> 2 CHROMOSOME_I    10  20    0.5    0.5     4     3     2     1     0         0
#>   seq_len
#> 1      10
#> 2      10
unlink(fai_path)
```

### Liftover score-style rows

``` r
lift_src <- tempfile("duckhts_liftover_src_", fileext = ".fa")
lift_dst <- tempfile("duckhts_liftover_dst_", fileext = ".fa")
lift_chain <- tempfile("duckhts_liftover_", fileext = ".chain")

writeLines(c(
  ">chrF",
  "ACGTACGTAA",
  ">chrR",
  "AACCGGTTAA"
), lift_src)
writeLines(c(
  ">chrLiftF",
  "ACGTACGTAA",
  ">chrLiftR",
  "TTAACCGGTT"
), lift_dst)
writeLines(c(
  "chain 100 chrF 10 + 0 10 chrLiftF 10 + 0 10 1",
  "10",
  "",
  "chain 100 chrR 10 + 0 10 chrLiftR 10 - 0 10 2",
  "10"
), lift_chain)

rduckhts_fasta_index(con, lift_src, index_path = paste0(lift_src, ".fai"))
#>   success                                                 index_path
#> 1    TRUE /tmp/RtmpDtZah6/duckhts_liftover_src_112d8f67b539a5.fa.fai
rduckhts_fasta_index(con, lift_dst, index_path = paste0(lift_dst, ".fai"))
#>   success                                                 index_path
#> 1    TRUE /tmp/RtmpDtZah6/duckhts_liftover_dst_112d8f37b8c8a3.fa.fai

lifted <- rduckhts_liftover(
  con,
  query = paste(
    "SELECT * FROM (VALUES",
    "('chrF', 2, 'C', 'T'),",
    "('chrR', 2, 'A', 'G'),",
    "('chrF', 11, 'A', 'T')",
    ") AS t(chrom, pos, ref, alt)"
  ),
  chain_path = lift_chain,
  dst_fasta_ref = lift_dst,
  ref_col = "ref",
  alt_col = "alt",
  src_fasta_ref = lift_src
)

lifted[, c(
  "src_chrom", "src_pos", "dest_chrom", "dest_pos",
  "dest_ref", "dest_alt", "mapped", "reverse_complemented",
  "reject_reason", "note"
)]
#>   src_chrom src_pos dest_chrom dest_pos dest_ref dest_alt mapped
#> 1      chrF       2   chrLiftF        2        C        T   TRUE
#> 2      chrR       2   chrLiftR        9        T        C   TRUE
#> 3      chrF      11   chrLiftF       10        A    AA,AT   TRUE
#>   reverse_complemented reject_reason   note
#> 1                FALSE          <NA>   <NA>
#> 2                 TRUE          <NA>   <NA>
#> 3                FALSE          <NA> Padded

unlink(c(lift_src, paste0(lift_src, ".fai"), lift_dst, paste0(lift_dst, ".fai"), lift_chain))
```

### Munge score-style rows

``` r
munge_fasta <- tempfile("duckhts_munge_", fileext = ".fa")
writeLines(c(
  ">chrF",
  "ACGTACGTAA"
), munge_fasta)
rduckhts_fasta_index(con, munge_fasta, index_path = paste0(munge_fasta, ".fai"))
#>   success                                         index_path
#> 1    TRUE /tmp/RtmpDtZah6/duckhts_munge_112d8ffd04e2a.fa.fai

munge_out <- rduckhts_munge(
  con,
  query = paste(
    "SELECT * FROM (VALUES",
    "('rs1', 2, 'chrF', 'A', 'C', 0.01, 1.10, 0.20, 0.98, 0.10, 0.01, 1000),",
    "('rs2', 2, 'chrF', 'C', 'A', 0.02, 0.90, -0.20, 0.98, 0.90, 0.01, 1000)",
    ") AS t(SNP, BP, CHR, A1, A2, P, OR_VALUE, BETA, INFO, FRQ, SE, N)"
  ),
  fasta_ref = munge_fasta,
  column_map = c(
    SNP = "SNP", BP = "BP", CHR = "CHR", A1 = "A1", A2 = "A2",
    P = "P", OR = "OR_VALUE", BETA = "BETA", INFO = "INFO", FRQ = "FRQ", SE = "SE", N = "N"
  )
)

munge_out[, c("chrom", "pos", "id", "ref", "alt", "alleles_swapped", "filter", "af", "es", "ns")]
#>   chrom pos  id ref alt alleles_swapped filter  af  es   ns
#> 1  chrF   2 rs2   C   A            TRUE   <NA> 0.1 0.2 1000
#> 2  chrF   2 rs1   C   A           FALSE   <NA> 0.1 0.2 1000

unlink(c(munge_fasta, paste0(munge_fasta, ".fai")))
```

### Polygenic risk scoring

`rduckhts_score()` computes per-sample polygenic risk scores (PRS) from
a genotype VCF/BCF and a GWAS summary statistics file, wrapping
`bcftools_score`.

``` r
vcf_path    <- system.file("extdata", "score_input.vcf",        package = "Rduckhts")
dosage_path <- system.file("extdata", "score_dosage.vcf",       package = "Rduckhts")
sumf_path   <- system.file("extdata", "score_summary.tsv",      package = "Rduckhts")
gwas_path   <- system.file("extdata", "score_gwas_summary.vcf", package = "Rduckhts")

# Hard-call (GT) PRS with PLINK-format summary statistics
# S1: 0×0.5 + 1×(−0.2) + 2×1.0 = 1.8
# S2: 1×0.5 + 2×(−0.2) + 0×1.0 = 0.1
gt_prs <- rduckhts_score(con, vcf_path, sumf_path, use = "GT", columns = "PLINK")
gt_prs[, c("SAMPLE", "score_summary")]
#>   SAMPLE score_summary
#> 1     S1    1.80000000
#> 2     S2    0.09999999

# Dosage-based PRS (DS field) for imputed genotypes
# S1: 0.1×0.5 + 0.8×(−0.2) + 1.8×1.0 = 1.69
# S2: 1.0×0.5 + 1.9×(−0.2) + 0.2×1.0 = 0.32
ds_prs <- rduckhts_score(con, dosage_path, sumf_path, use = "DS", columns = "PLINK")
ds_prs[, c("SAMPLE", "score_summary")]
#>   SAMPLE score_summary
#> 1     S1          1.69
#> 2     S2          0.32

# GWAS-VCF multi-PRS: each FORMAT/ES sample column becomes a separate PRS track
gwas_prs <- rduckhts_score(con, vcf_path, gwas_path, use = "GT")
gwas_prs[, c("SAMPLE", "PRS_A", "PRS_B")]
#>   SAMPLE      PRS_A PRS_B
#> 1     S1 1.80000000   1.0
#> 2     S2 0.09999999   0.3
```

### Compression + tabix round-trips

``` r
tmp_bed <- tempfile("duckhts_targets_", fileext = ".bed")
tmp_bgz <- paste0(tmp_bed, ".gz")
tmp_tbi <- paste0(tmp_bgz, ".tbi")
writeLines(c("chr1\t0\t10\ta", "chr1\t10\t20\tb"), tmp_bed)

rduckhts_bgzip(con, tmp_bed, output_path = tmp_bgz, keep = TRUE, overwrite = TRUE)
#>   success                                           output_path bytes_in
#> 1    TRUE /tmp/RtmpDtZah6/duckhts_targets_112d8f377a2741.bed.gz       25
#>   bytes_out
#> 1        84
rduckhts_tabix_index(con, tmp_bgz, preset = "bed", index_path = tmp_tbi, threads = 1)
#>   success                                                index_path
#> 1    TRUE /tmp/RtmpDtZah6/duckhts_targets_112d8f377a2741.bed.gz.tbi
#>   index_format
#> 1          TBI
rduckhts_bed(con, "targets_idx", tmp_bgz, region = "chr1:1-20", index_path = tmp_tbi, overwrite = TRUE)
dbGetQuery(con, "SELECT * FROM targets_idx")
#>   chrom start end name score strand thick_start thick_end item_rgb block_count
#> 1  chr1     0  10    a  <NA>   <NA>          NA        NA     <NA>          NA
#> 2  chr1    10  20    b  <NA>   <NA>          NA        NA     <NA>          NA
#>   block_sizes block_starts extra
#> 1        <NA>         <NA>  <NA>
#> 2        <NA>         <NA>  <NA>

unlink(c(tmp_bed, tmp_bgz, tmp_tbi))
```

## Sequence UDFs

The extension also exposes sequence utility UDFs directly in DuckDB SQL,
including 4-bit IUPAC DNA encode/decode helpers. These can be applied to
`SEQUENCE` columns from FASTA and FASTQ scans.

``` r
dbGetQuery(
  con,
  "SELECT
     NAME,
     seq_hash_2bit(substr(SEQUENCE, 1, 12)) AS hash_2bit_prefix,
     seq_encode_4bit(substr(SEQUENCE, 1, 16)) AS codes,
     seq_decode_4bit(seq_encode_4bit(substr(SEQUENCE, 1, 16))) AS roundtrip
   FROM sequences
   LIMIT 2"
)
#>            NAME hash_2bit_prefix                                          codes
#> 1  CHROMOSOME_I          9898352 4, 2, 2, 8, 1, 1, 4, 2, 2, 8, 1, 1, 4, 2, 2, 8
#> 2 CHROMOSOME_II          6038978 2, 2, 8, 1, 1, 4, 2, 2, 8, 1, 1, 4, 2, 2, 8, 1
#>          roundtrip
#> 1 GCCTAAGCCTAAGCCT
#> 2 CCTAAGCCTAAGCCTA

dbGetQuery(
  con,
  "SELECT
     NAME,
     MATE,
     seq_encode_4bit(substr(SEQUENCE, 1, 12)) AS codes,
     seq_decode_4bit(seq_encode_4bit(substr(SEQUENCE, 1, 12))) AS roundtrip
   FROM reads
   LIMIT 2"
)
#>                              NAME MATE                              codes
#> 1 HS25_09827:2:1201:1505:59795#49    1 2, 2, 4, 8, 8, 1, 4, 1, 4, 2, 1, 8
#> 2 HS25_09827:2:1201:1505:59795#49    2 1, 1, 4, 4, 1, 1, 1, 4, 1, 1, 4, 4
#>      roundtrip
#> 1 CCGTTAGAGCAT
#> 2 AAGGAAAGAAGG
```

### FASTA region queries

`read_fasta` supports indexed region queries via
`rduckhts_fasta(..., region = ...)`.

``` r
fai_path <- tempfile("duckhts_readme_", fileext = ".fai")
fai_info <- rduckhts_fasta_index(con, fasta_path, index_path = fai_path)
fai_info
#>   success                                        index_path
#> 1    TRUE /tmp/RtmpDtZah6/duckhts_readme_112d8f559811ec.fai

rduckhts_fasta(
  con, "fasta_region", fasta_path,
  region = "CHROMOSOME_I:1-25",
  overwrite = TRUE
)
dbGetQuery(con, "SELECT NAME, length(SEQUENCE) AS n FROM fasta_region")
#>           NAME  n
#> 1 CHROMOSOME_I 25
unlink(fai_path)
```

## Examples

### Region Queries

Region queries can use implicit sidecar indexes or an explicit
`index_path` for custom index names/locations.

``` r
bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
bcf_index_path <- system.file("extdata", "vcf_file.bcf.csi", package = "Rduckhts")
rduckhts_bcf(con, "variants", bcf_path, overwrite = TRUE)
variants <- dbGetQuery(con, "SELECT * FROM variants LIMIT 5")
variants
#>   CHROM     POS    ID  REF  ALT  QUAL FILTER INFO_TEST   INFO_DP4 INFO_AC
#> 1     1 3000150  <NA>    C    T  59.2   PASS        NA       NULL       2
#> 2     1 3000151  <NA>    C    T  59.2   PASS        NA       NULL       2
#> 3     1 3062915  id3D GTTT    G  12.9    q10        NA 1, 2, 3, 4       2
#> 4     1 3062915 idSNP    G T, C  12.6   test         5 1, 2, 3, 4    1, 1
#> 5     1 3106154  <NA> CAAA    C 342.0   PASS        NA       NULL       2
#>   INFO_AN INFO_INDEL INFO_STR FORMAT_TT_A FORMAT_GT_A FORMAT_GQ_A FORMAT_DP_A
#> 1       4      FALSE     <NA>        NULL         0/1         245          NA
#> 2       4      FALSE     <NA>        NULL         0/1         245          32
#> 3       4       TRUE     test        NULL         0/1         409          35
#> 4       3      FALSE     <NA>        0, 1         0/1         409          35
#> 5       4      FALSE     <NA>        NULL         0/1         245          32
#>                  FORMAT_GL_A FORMAT_TT_B FORMAT_GT_B FORMAT_GQ_B FORMAT_DP_B
#> 1                       NULL        NULL         0/1         245          NA
#> 2                       NULL        NULL         0/1         245          32
#> 3               -20, -5, -20        NULL         0/1         409          35
#> 4 -20, -5, -20, -20, -5, -20        0, 1           2         409          35
#> 5                       NULL        NULL         0/1         245          32
#>    FORMAT_GL_B
#> 1         NULL
#> 2         NULL
#> 3 -20, -5, -20
#> 4 -20, -5, -20
#> 5         NULL

rduckhts_bcf(
  con, "variants_idx", bcf_path,
  region = "1:3000150-3000151",
  index_path = bcf_index_path,
  overwrite = TRUE
)
dbGetQuery(con, "SELECT count(*) AS n FROM variants_idx")
#>   n
#> 1 2

# Span-oriented index view from the same file
index_spans_preview <- rduckhts_hts_index_spans(con, bcf_path, index_path = bcf_index_path)
head(index_spans_preview[, c("seqname", "tid", "index_type", "chunk_beg_vo", "chunk_end_vo")], 5)
#>   seqname tid index_type chunk_beg_vo chunk_end_vo
#> 1       1   0        CSI         1586         1713
#> 2       1   0        CSI         1713         1973
#> 3       1   0        CSI         1973         2109
#> 4       1   0        CSI         2109         2242
#> 5       1   0        CSI         2242         2372
```

### Remote VCF on S3

S3 files can be query when `htslib` is built with plugins enable. This
is not the case on RTools

``` r
# Enable htslib plugins for remote access (S3/GCS/HTTP)
setup_hts_env()

# Example S3 URL (1000 Genomes cohort VCF)
s3_base <- "s3://1000genomes-dragen-v3.7.6/data/cohorts/"
s3_path <- "gvcf-genotyper-dragen-3.7.6/hg19/3202-samples-cohort/"
s3_vcf_file <- "3202_samples_cohort_gg_chr22.vcf.gz"
s3_vcf_uri <- paste0(s3_base, s3_path, s3_vcf_file)

rduckhts_bcf(con, "s3_variants", s3_vcf_uri, region = "chr22:16050000-16050500", overwrite = TRUE)
dbGetQuery(con, "SELECT CHROM, COUNT(*) AS n FROM s3_variants GROUP BY CHROM")
#>   CHROM  n
#> 1 chr22 11
```

### FASTQ files

Three modes for fastq files, single, paired and interleaved

``` r
r1 <- system.file("extdata", "r1.fq", package = "Rduckhts")
r2 <- system.file("extdata", "r2.fq", package = "Rduckhts")
interleaved <- system.file("extdata", "interleaved.fq", package = "Rduckhts")
rduckhts_fastq(con, "paired_reads", r1, mate_path = r2, overwrite = TRUE)
rduckhts_fastq(con, "interleaved_reads", interleaved, interleaved = TRUE, overwrite = TRUE)
pairs <- dbGetQuery(con, "SELECT * FROM paired_reads WHERE MATE = 1 LIMIT 5")
pairs
#>                              NAME DESCRIPTION
#> 1 HS25_09827:2:1201:1505:59795#49        <NA>
#> 2 HS25_09827:2:1201:1559:70726#49        <NA>
#> 3 HS25_09827:2:1201:1564:39627#49        <NA>
#> 4 HS25_09827:2:1201:1565:91731#49        <NA>
#> 5 HS25_09827:2:1201:1624:69925#49        <NA>
#>                                                                                               SEQUENCE
#> 1 CCGTTAGAGCATTTGTTGAAAATGCTTTCCTTGCTCCATGTGATGACTCTGGTGCCCTTGTCAAAAGCCAGCTGGGCCTATTCGTGTGGGTCTGTTTCTG
#> 2 TTGTTAAAATGACCATACCCAAAGTGATCTACAGACTCAATACAATTTCTATTGAAATACCAATCACACTCTTCACAGAACTAGAAAAACAGTTCTAAAA
#> 3 ACGCGGCAATCCAATGTGTGAGTTGAGAAGCGGTGAGGAGGGAATCCTAATTTTATGAGCAGGTCAGGACCGTGGGAGATACCTGACACCTGAGATGGTA
#> 4 GACATGCCATAACATTCATGTTTTATGTGTACAAGTCAATGAATTTTAGTATATTTACAGAGTTGTATGACTGTCTCCACAATCTAATTTTAGGTTTCCA
#> 5 GCCAGCCTCCTTCTCAATGGTCTTTTTAAACATTATATGAAAACCAGACATTTACATTTGATTTCTTTTTCAATACTATACAGTTCTAAGAGAAAAAACA
#>                                                                                                QUALITY
#> 1 CABCFGDEEFFEFHGHGGFFGDIGIJFIFHHGHEIFGHBCGHDIFBE9GIAICGGICFIBFGGHGDGGGHE?GIGDFGGHEGIEJG>;FG<GGHACEFGH
#> 2 CABEFGFFGFHGGGGJGGFFGKIHHJFIEHHHGIEGGEHJGHDHFGHIGICIJEFIFGIF8GGHKFHGGFEI6GGGFIGHGGIE>EFCFHGGGHEJEAJE
#> 3 BACCFGBFGFHGGJGHGGFEGHIGIJHFEH:HHEHGHHBGGH9IAGHGFHIFJFFAFGIFDIGHKEIG<C>F,CGD66?7EFI5EEG>EGGGGD5=HH6E
#> 4 CABFFGFFJFHEGEGJGGDG?FIGHHHBGHHHGIIGHGHGGHDGHFHIDFCIKEGIFHGGII9HFFGGGEEIGGEEHGGEEGDEHFH>FGGGGHAFAHGE
#> 5 CABEFGFGIFGGGJGHGGFH?FDHGHDHGHEHHJCGHHFHDHDHFGHIGHIFFHGHFGGGI9GHF@IGGH;FICGEFEIHGGIEEFC:DEGGGBDJHHFF
#>   MATE                         PAIR_ID
#> 1    1 HS25_09827:2:1201:1505:59795#49
#> 2    1 HS25_09827:2:1201:1559:70726#49
#> 3    1 HS25_09827:2:1201:1564:39627#49
#> 4    1 HS25_09827:2:1201:1565:91731#49
#> 5    1 HS25_09827:2:1201:1624:69925#49
```

### FASTQ quality decoding and per-position histograms

FASTQ quality handling has two separate knobs:

- `input_quality_encoding` controls how incoming FASTQ ASCII is decoded.
  The default is modern `phred33`; use `phred64`, `solexa64`, or `auto`
  only for legacy data.
- `quality_representation` controls how qualities are returned to
  DuckDB: canonical Phred+33 text (`"string"`) or numeric Phred arrays
  (`"phred"`).

The flow is:

1.  Decode FASTQ ASCII using `input_quality_encoding`.
2.  Normalize to numeric Phred qualities.
3.  Return either numeric arrays or canonical Phred+33 text.

`BAM`/`CRAM` reads skip the text-decoding step because qualities are
already stored numerically.

``` r
legacy_fastq <- system.file("extdata", "legacy_phred64.fq", package = "Rduckhts")

rduckhts_detect_quality_encoding(con, legacy_fastq)
#>   format observed_ascii_min observed_ascii_max records_sampled
#> 1  fastq                104                104               1
#>       compatible_encodings guessed_encoding is_ambiguous
#> 1 phred33,phred64,solexa64          phred64         TRUE

quality_hist <- dbGetQuery(
  con,
  sprintf(
    "WITH q AS (
       SELECT NAME, QUALITY
       FROM read_fastq('%s', quality_representation := 'phred')
     ),
     expanded AS (
       SELECT
         NAME,
         generate_subscripts(QUALITY, 1) AS pos,
         unnest(QUALITY) AS q
       FROM q
     )
     SELECT pos, q AS phred, count(*) AS n_reads
     FROM expanded
     GROUP BY pos, phred
     ORDER BY pos, phred
     LIMIT 12",
    fastq_r1
  )
)
quality_hist
#>    pos phred n_reads
#> 1    1    33       1
#> 2    1    34       4
#> 3    2    32       5
#> 4    3    33       4
#> 5    3    34       1
#> 6    4    34       2
#> 7    4    36       2
#> 8    4    37       1
#> 9    5    37       5
#> 10   6    38       5
#> 11   7    33       1
#> 12   7    35       1
```

### GFF files

These can be open with or with attributes maps

``` r
gff_path <- system.file("extdata", "gff_file.gff.gz", package = "Rduckhts")
rduckhts_gff(con, "genes", gff_path, attributes_map = TRUE, overwrite = TRUE)
gene_annotations <- dbGetQuery(con, "SELECT seqname, start, \"end\" FROM genes WHERE feature = 'gene' LIMIT 5")
gene_annotations
#>   seqname   start     end
#> 1       X 2934816 2964270
```

### BAM/CRAM

When built with htslib codec, `CRAM` can be opened in addition to `BAM`
files. `index_path` can also be passed for region scans with
non-standard index names.

``` r
cram_path <- system.file("extdata", "range.cram", package = "Rduckhts")
ref_path <- system.file("extdata", "ce.fa", package = "Rduckhts")
bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
bam_index_path <- system.file("extdata", "range.bam.bai", package = "Rduckhts")

rduckhts_bam(con, "cram_reads", cram_path, reference = ref_path, overwrite = TRUE)
cram_reads <- dbGetQuery(con, "SELECT QNAME, FLAG, POS, MAPQ FROM cram_reads LIMIT 5")
cram_reads
#>                           QNAME FLAG  POS MAPQ
#> 1 HS18_09653:4:1315:19857:61712  145  914   23
#> 2 HS18_09653:4:1308:11522:27107  161  934    0
#> 3 HS18_09653:4:2314:14991:85680   83 1020   10
#> 4 HS18_09653:4:2108:14085:93656  147 1122   60
#> 5  HS18_09653:4:1303:4347:38100   83 1137   37

rduckhts_bam(
  con, "bam_idx_reads", bam_path,
  region = "CHROMOSOME_I:1-1000",
  index_path = bam_index_path,
  overwrite = TRUE
)
dbGetQuery(con, "SELECT count(*) AS n FROM bam_idx_reads")
#>   n
#> 1 2
```

### SAMtags + auxiliary tags

Standard SAMtags can be exposed as typed columns, and any remaining tags
are available via `AUXILIARY_TAGS`:

``` r
aux_path <- system.file("extdata", "aux_tags.sam.gz", package = "Rduckhts")
rduckhts_bam(con, "aux_reads", aux_path, standard_tags = TRUE, auxiliary_tags = TRUE, overwrite = TRUE)
dbGetQuery(con, "SELECT RG, NM, map_extract(AUXILIARY_TAGS, 'XZ') AS XZ FROM aux_reads LIMIT 1")
#>   RG NM  XZ
#> 1 x1  2 foo
```

### Tabix headers + types

Use `header = TRUE` to use the first non-meta row as column names, and
`auto_detect = TRUE` / `column_types` to control column typing:

``` r
tabix_header <- system.file("extdata", "header_tabix.tsv.gz", package = "Rduckhts")
tabix_meta <- system.file("extdata", "meta_tabix.tsv.gz", package = "Rduckhts")

rduckhts_tabix(con, "header_tabix", tabix_header, header = TRUE, overwrite = TRUE)
dbGetQuery(con, "SELECT chrom, pos FROM header_tabix LIMIT 2")
#>   chrom pos
#> 1  chr1   1
#> 2  chr1   2

rduckhts_tabix(con, "typed_tabix", tabix_meta, auto_detect = TRUE, overwrite = TRUE)
dbGetQuery(con, "SELECT typeof(column1) AS column1_type FROM typed_tabix LIMIT 1")
#>   column1_type
#> 1       BIGINT

rduckhts_tabix(con, "typed_tabix_explicit", tabix_header,
               header = TRUE,
               column_types = c("VARCHAR", "BIGINT", "VARCHAR"),
               overwrite = TRUE)
dbGetQuery(con, "SELECT pos + 1 AS pos_plus_one FROM typed_tabix_explicit LIMIT 1")
#>   pos_plus_one
#> 1            2
```

### HTS header and index metadata

Use metadata helpers to inspect parsed headers, raw header lines, index
summaries, span-oriented index views, and raw index blobs.

``` r
header_meta <- rduckhts_hts_header(con, bcf_path)
head(header_meta[, c("record_type", "id", "number", "value_type")], 5)
#>   record_type   id number value_type
#> 1  fileformat <NA>   <NA>       <NA>
#> 2      FILTER PASS   <NA>       <NA>
#> 3        INFO TEST      1    Integer
#> 4      FORMAT   TT      A    Integer
#> 5        INFO  DP4      4    Integer

header_raw <- rduckhts_hts_header(con, bcf_path, mode = "raw")
head(header_raw[, c("idx", "raw")], 5)
#>   idx
#> 1   0
#> 2   1
#> 3   2
#> 4   3
#> 5   4
#>                                                                                                                                                          raw
#> 1                                                                                                                                       ##fileformat=VCFv4.1
#> 2                                                                                                        ##FILTER=<ID=PASS,Description="All filters passed">
#> 3                                                                                           ##INFO=<ID=TEST,Number=1,Type=Integer,Description="Testing Tag">
#> 4 ##FORMAT=<ID=TT,Number=A,Type=Integer,Description="Testing Tag, with commas and \\"escapes\\" and escaped escapes combined with \\\\\\"quotes\\\\\\\\\\"">
#> 5                       ##INFO=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">

index_meta <- rduckhts_hts_index(con, bcf_path, index_path = bcf_index_path)
head(index_meta[, c("seqname", "mapped", "unmapped", "index_type")], 5)
#>   seqname mapped unmapped index_type
#> 1       1     11        0        CSI
#> 2       2      1        0        CSI
#> 3       3      1        0        CSI
#> 4       4      2        0        CSI

index_spans <- rduckhts_hts_index_spans(con, bcf_path, index_path = bcf_index_path)
head(index_spans[, c("seqname", "tid", "index_type", "chunk_beg_vo", "chunk_end_vo")], 5)
#>   seqname tid index_type chunk_beg_vo chunk_end_vo
#> 1       1   0        CSI         1586         1713
#> 2       1   0        CSI         1713         1973
#> 3       1   0        CSI         1973         2109
#> 4       1   0        CSI         2109         2242
#> 5       1   0        CSI         2242         2372

index_raw <- rduckhts_hts_index_raw(con, bcf_path, index_path = bcf_index_path)
head(index_raw, 1)
#> [1] index_type                                                       
#> [2] '/usr/local/lib/R/site-library/Rduckhts/extdata/vcf_file.bcf.csi'
#> [3] raw                                                              
#> <0 rows> (or 0-length row.names)
```

### Remote GTEx tabix example

GTEx eQTL matrices on EBI are tabix-indexed

``` r
gtex_url <- "http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/Brain_Cerebellar_Hemisphere.tsv.gz" 
rduckhts_tabix(con, "gtex_eqtl", gtex_url, region = "1:11868-14409",
                 header = TRUE, auto_detect = TRUE, overwrite = TRUE)
dbGetQuery(con, "SELECT * FROM gtex_eqtl LIMIT 5")
#>          variant r2    pvalue molecular_trait_object_id molecular_trait_id
#> 1 chr1_13550_G_A NA 0.0204520           ENSG00000188290    ENSG00000188290
#> 2 chr1_13550_G_A NA 0.0303633           ENSG00000230699    ENSG00000230699
#> 3 chr1_13550_G_A NA 0.1057900           ENSG00000177757    ENSG00000177757
#> 4 chr1_13550_G_A NA 0.1617190           ENSG00000241860    ENSG00000241860
#> 5 chr1_13550_G_A NA 0.1919580           ENSG00000198744    ENSG00000198744
#>         maf         gene_id median_tpm      beta       se  an ac chromosome
#> 1 0.0114286 ENSG00000188290  6.3960000  0.633986 0.270285 350  4          1
#> 2 0.0114286 ENSG00000230699  0.0674459 -0.980082 0.447861 350  4          1
#> 3 0.0114286 ENSG00000177757  1.2659000  0.631359 0.387738 350  4          1
#> 4 0.0114286 ENSG00000241860  0.1081970 -0.791695 0.562674 350  4          1
#> 5 0.0114286 ENSG00000198744 21.6284000 -0.592354 0.451705 350  4          1
#>   position ref alt type        rsid
#> 1    13550   G   A  SNP rs554008981
#> 2    13550   G   A  SNP rs554008981
#> 3    13550   G   A  SNP rs554008981
#> 4    13550   G   A  SNP rs554008981
#> 5    13550   G   A  SNP rs554008981
```

``` r
dbDisconnect(con, shutdown = TRUE)
```

## References

- DuckDB: <https://duckdb.org/>
- DuckDB Extension API: <https://duckdb.org/docs/extensions/overview>
- DuckDB extension template (C):
  <https://github.com/duckdb/extension-template-c>
- htslib: <https://github.com/samtools/htslib>
- RBCFTools: <https://github.com/RGenomicsETL/RBCFTools>

## License

GPL-3.
