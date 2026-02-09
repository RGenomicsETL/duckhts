
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Rduckhts: DuckDB HTS File Reader Extension for R

[![CRAN
Status](https://www.r-pkg.org/badges/version/Rduckhts)](https://cran.r-project.org/package=Rduckhts)
[![R-CMD-check](https://github.com/RGenomicsETL/duckhts/workflows/R-CMD-check/badge.svg)](https://github.com/RGenomicsETL/duckhts/actions)

We provide an R interface to the [DuckDB](https://duckdb.org/) HTS (High
Throughput Sequencing) file reader extension. We enable reading common
bioinformatics file formats such as VCF/BCF, SAM/BAM/CRAM, FASTA, FASTQ,
GFF, GTF, and tabix-indexed files directly from R using SQL queries via
[DuckDB](https://duckdb.org/).

## Key Features

We follow the [RBCFTools](https://github.com/RGenomicsETL/RBCFTools)
table pattern and create [DuckDB](https://duckdb.org/) tables instead of
returning data frames. We support VCF/BCF, SAM/BAM/CRAM, FASTA, FASTQ,
GFF, GTF, and tabix. We support region queries for indexed files, and we
target Linux, macOS, and Windows (MinGW). We bundle
[htslib](https://github.com/samtools/htslib) 1.23 so runtime
dependencies stay minimal.

## Installation

We install from GitHub with
`remotes::install_github("RGenomicsETL/duckhts", subdir = "r/Rduckhts")`.

## Quick Start (runnable example)

``` r
library(DBI)
library(duckdb)
library(Rduckhts)

setup_hts_env()

ext_path <- system.file("extdata", "duckhts.duckdb_extension", package = "Rduckhts")
fasta_path <- system.file("extdata", "ce.fa", package = "Rduckhts")
fastq_r1 <- system.file("extdata", "r1.fq", package = "Rduckhts")
fastq_r2 <- system.file("extdata", "r2.fq", package = "Rduckhts")
con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
rduckhts_load(con, extension_path = ext_path)
#> [1] TRUE

rduckhts_fasta(con, "sequences", fasta_path, overwrite = TRUE)
rduckhts_fastq(con, "reads", fastq_r1, mate_path = fastq_r2, overwrite = TRUE)

dbGetQuery(con, "SELECT COUNT(*) AS n FROM sequences")
#>   n
#> 1 7
dbGetQuery(con, "SELECT COUNT(*) AS n FROM reads")
#>    n
#> 1 10

message("Loaded example data and created tables.")
```

We load the extension with `rduckhts_load(con, extension_path = NULL)`.
We create tables with `rduckhts_bcf`, `rduckhts_bam`, `rduckhts_fasta`,
`rduckhts_fastq`, `rduckhts_gff`, `rduckhts_gtf`, and `rduckhts_tabix`
using the parameters documented in their help pages.

## Advanced Examples

### Region Queries

``` r
bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
rduckhts_bcf(con, "variants", bcf_path, overwrite = TRUE)
variants <- dbGetQuery(con, "SELECT * FROM variants LIMIT 5")
variants
#>   CHROM     POS    ID  REF  ALT  QUAL FILTER INFO_TEST INFO_DP4 INFO_AC INFO_AN
#> 1     1 3000150  <NA>    C    T  59.2   PASS        NA       NA       2       4
#> 2     1 3000151  <NA>    C    T  59.2   PASS        NA       NA       2       4
#> 3     1 3062915  id3D GTTT    G  12.9    q10        NA        1       2       4
#> 4     1 3062915 idSNP    G T, C  12.6   test         5        1    1, 1       3
#> 5     1 3106154  <NA> CAAA    C 342.0   PASS        NA       NA       2       4
#>   INFO_INDEL INFO_STR FORMAT_TT_A FORMAT_GT_A FORMAT_GQ_A FORMAT_DP_A
#> 1      FALSE     <NA>        NULL         0/1         245          NA
#> 2      FALSE     <NA>        NULL         0/1         245          32
#> 3       TRUE     test        NULL         0/1         409          35
#> 4      FALSE     <NA>        0, 1         0/1         409          35
#> 5      FALSE     <NA>        NULL         0/1         245          32
#>                  FORMAT_GL_A FORMAT_TT_B FORMAT_GT_B FORMAT_GQ_B FORMAT_DP_B
#> 1                       NULL        NULL         0/1         245          NA
#> 2                       NULL        NULL         0/1         245          32
#> 3               -20, -5, -20        NULL         0/1         409          35
#> 4 -20, -5, -20, -20, -5, -20        0, 1          2|         409          35
#> 5                       NULL        NULL         0/1         245          32
#>    FORMAT_GL_B
#> 1         NULL
#> 2         NULL
#> 3 -20, -5, -20
#> 4 -20, -5, -20
#> 5         NULL
```

### Remote VCF on S3 (requires libcurl)

``` r
# Enable htslib plugins for remote access (S3/GCS/HTTP)
setup_hts_env()

# Example S3 URL (1000 Genomes cohort VCF)
s3_base <- "s3://1000genomes-dragen-v3.7.6/data/cohorts/"
s3_path <- "gvcf-genotyper-dragen-3.7.6/hg19/3202-samples-cohort/"
s3_vcf_file <- "3202_samples_cohort_gg_chr22.vcf.gz"
s3_vcf_uri <- paste0(s3_base, s3_path, s3_vcf_file)

# Query remote VCF directly with DuckDB + DuckHTS (region-scoped)
rduckhts_bcf(con, "s3_variants", s3_vcf_uri, region = "chr22:10000000-10550000", overwrite = TRUE)
dbGetQuery(con, "SELECT CHROM, COUNT(*) AS n FROM s3_variants GROUP BY CHROM")
#> [1] CHROM n    
#> <0 rows> (or 0-length row.names)
```

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

``` r
gff_path <- system.file("extdata", "gff_file.gff.gz", package = "Rduckhts")
rduckhts_gff(con, "genes", gff_path, attributes_map = TRUE, overwrite = TRUE)
gene_annotations <- dbGetQuery(con, "SELECT seqname, start, \"end\" FROM genes WHERE feature = 'gene' LIMIT 5")
gene_annotations
#>   seqname   start     end
#> 1       X 2934816 2964270
```

``` r
cram_path <- system.file("extdata", "range.cram", package = "Rduckhts")
ref_path <- system.file("extdata", "ce.fa", package = "Rduckhts")
rduckhts_bam(con, "cram_reads", cram_path, reference = ref_path, overwrite = TRUE)
cram_reads <- dbGetQuery(con, "SELECT QNAME, FLAG, POS, MAPQ FROM cram_reads LIMIT 5")
cram_reads
#>                           QNAME FLAG  POS MAPQ
#> 1 HS18_09653:4:1315:19857:61712  145  914   23
#> 2 HS18_09653:4:1308:11522:27107  161  934    0
#> 3 HS18_09653:4:2314:14991:85680   83 1020   10
#> 4 HS18_09653:4:2108:14085:93656  147 1122   60
#> 5  HS18_09653:4:1303:4347:38100   83 1137   37
```

## System Requirements

We require zlib and libbz2, and we optionally use liblzma, libcurl, and
openssl during build. We build with cmake, GNU make, and a C compiler.
On Ubuntu/Debian we install build-essential, cmake, and the listed dev
packages. On macOS we install cmake and
[htslib](https://github.com/samtools/htslib). On Windows we use Rtools
and MinGW.

## Notes

We follow the [RBCFTools](https://github.com/RGenomicsETL/RBCFTools)
pattern and create [DuckDB](https://duckdb.org/) tables rather than
returning data frames. We check for existing tables and require
`overwrite = TRUE` to replace them. We use region queries for indexed
files, we bundle [htslib](https://github.com/samtools/htslib) 1.23, and
we enable plugins where supported.

## References

- DuckDB: <https://duckdb.org/>
- DuckDB Extension API: <https://duckdb.org/docs/extensions/overview>
- DuckDB extension template (C):
  <https://github.com/duckdb/extension-template-c>
- htslib: <https://github.com/samtools/htslib>
- RBCFTools: <https://github.com/RGenomicsETL/RBCFTools>

## License

GPL-3.

``` r
dbDisconnect(con, shutdown = TRUE)
```
