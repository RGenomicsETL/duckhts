
# Rduckhts: DuckDB HTS File Reader Extension for R

[![CRAN
Status](https://www.r-pkg.org/badges/version/Rduckhts)](https://cran.r-project.org/package=Rduckhts)
[![R-CMD-check](https://github.com/RGenomicsETL/duckhts/workflows/CRAN-check/badge.svg)](https://github.com/RGenomicsETL/duckhts/actions)

We provide an R interface to the DuckDB HTS (High Throughput Sequencing)
file reader extension. We enable reading common bioinformatics file
formats such as VCF/BCF, SAM/BAM/CRAM, FASTA, FASTQ, GFF, GTF, and
tabix-indexed files directly from R using SQL queries via DuckDB.

## Key Features

We follow the RBCFTools table pattern and create DuckDB tables instead
of returning data frames. We support VCF/BCF, SAM/BAM/CRAM, FASTA,
FASTQ, GFF, GTF, and tabix. We support region queries for indexed files,
and we target Linux, macOS, and Windows (MinGW). We bundle htslib 1.23
so runtime dependencies stay minimal.

## Installation

``` r
remotes::install_github("RGenomicsETL/duckhts", subdir = "r/Rduckhts")
```

## Quick Start (runnable example)

``` r
library(DBI)
library(duckdb)
library(Rduckhts)

ext_path <- system.file("extdata", "duckhts.duckdb_extension", package = "Rduckhts")
if (!file.exists(ext_path)) {
    ext_path <- file.path("inst", "extdata", "duckhts.duckdb_extension")
}
if (!file.exists(ext_path)) {
    stop("DuckHTS extension not found at: ", ext_path)
}

tmp_dir <- tempdir()
fasta_path <- file.path(tmp_dir, "example.fa")
fastq_path <- file.path(tmp_dir, "example.fq")

writeLines(c(">seq1", "ACGTACGT", ">seq2", "TTGGCCAA"), fasta_path)
writeLines(c("@r1", "ACGT", "+", "IIII"), fastq_path)

con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
rduckhts_load(con, extension_path = ext_path)
```

    ## [1] TRUE

``` r
rduckhts_fasta(con, "sequences", fasta_path, overwrite = TRUE)
rduckhts_fastq(con, "reads", fastq_path, overwrite = TRUE)

dbGetQuery(con, "SELECT COUNT(*) AS n FROM sequences")
```

    ##   n
    ## 1 2

``` r
dbGetQuery(con, "SELECT COUNT(*) AS n FROM reads")
```

    ##   n
    ## 1 1

``` r
dbDisconnect(con, shutdown = TRUE)
```

## Available Functions

We load the extension with `rduckhts_load(con, extension_path = NULL)`.
We create tables with `rduckhts_bcf`, `rduckhts_bam`, `rduckhts_fasta`,
`rduckhts_fastq`, `rduckhts_gff`, `rduckhts_gtf`, and `rduckhts_tabix`
using the parameters documented in their help pages.

## Advanced Examples

### Region Queries

``` r
rduckhts_bcf(con, "chr1_variants", "file.vcf.gz", region = "chr1:1000000-2000000")
chr1_variants <- dbGetQuery(con, "SELECT * FROM chr1_variants WHERE POS BETWEEN 1500000 AND 1600000")
```

### Paired FASTQ

``` r
rduckhts_fastq(con, "paired_reads", "r1.fq", mate_path = "r2.fq")
rduckhts_fastq(con, "interleaved_reads", "interleaved.fq", interleaved = TRUE)
pairs <- dbGetQuery(con, "SELECT * FROM paired_reads WHERE MATE = 1")
```

### Attributes Map for GFF/GTF

``` r
rduckhts_gff(con, "genes", "annotations.gff3.gz", attributes_map = TRUE)
gene_annotations <- dbGetQuery(con, "SELECT seqname, start, end, attributes['gene_id'], attributes['gene_name'] FROM genes WHERE feature = 'gene'")
```

### CRAM with Reference

``` r
rduckhts_bam(con, "cram_reads", "reads.cram", reference = "reference.fa")
cram_reads <- dbGetQuery(con, "SELECT QNAME, FLAG, POS, MAPQ FROM cram_reads")
```

## System Requirements

We require zlib and libbz2, and we optionally use liblzma, libcurl, and
openssl during build. We build with cmake, GNU make, and a C compiler.
On Ubuntu/Debian we install build-essential, cmake, and the listed dev
packages. On macOS we install cmake and htslib. On Windows we use Rtools
and MinGW.

## Notes

We follow the RBCFTools pattern and create DuckDB tables rather than
returning data frames. We check for existing tables and require
`overwrite = TRUE` to replace them. We use region queries for indexed
files, we bundle htslib 1.23, and we enable plugins where supported.

## References

We reference DuckDB, htslib, and RBCFTools as the primary upstream
projects that inform this work.

## License

GPL-3.
