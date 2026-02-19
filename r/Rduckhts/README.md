
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

The extension is loaded with `rduckhts_load(con, extension_path =
NULL)`. We can create tables with `rduckhts_bcf`, `rduckhts_bam`,
`rduckhts_fasta`, `rduckhts_fastq`, `rduckhts_gff`, `rduckhts_gtf`, and
`rduckhts_tabix` using the parameters documented in their help pages

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

## Examples

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

### Remote VCF on S3

S3 files can be query when `htslib` is built with plugins enable. This
is not the case on RTools

``` r
# Not run on CRAN because it requires network access.
# Enable htslib plugins for remote access (S3/GCS/HTTP)
setup_hts_env()

# Example S3 URL (1000 Genomes cohort VCF)
s3_base <- "s3://1000genomes-dragen-v3.7.6/data/cohorts/"
s3_path <- "gvcf-genotyper-dragen-3.7.6/hg19/3202-samples-cohort/"
s3_vcf_file <- "3202_samples_cohort_gg_chr22.vcf.gz"
s3_vcf_uri <- paste0(s3_base, s3_path, s3_vcf_file)

# Query remote VCF directly with DuckDB + DuckHTS (region-scoped)
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
files

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

Use the metadata helpers to inspect parsed headers and index summaries.

``` r
bcf_index_path <- system.file("extdata", "vcf_file.bcf.csi", package = "Rduckhts")

header_meta <- rduckhts_hts_header(con, bcf_path)
head(header_meta[, c("record_type", "id", "number", "value_type")], 5)
#>   record_type   id number value_type
#> 1  fileformat <NA>   <NA>       <NA>
#> 2      FILTER PASS   <NA>       <NA>
#> 3        INFO TEST      1    Integer
#> 4      FORMAT   TT      A    Integer
#> 5        INFO  DP4      4    Integer

index_meta <- rduckhts_hts_index(con, bcf_path, index_path = bcf_index_path)
head(index_meta[, c("seqname", "mapped", "unmapped", "index_type")], 5)
#>   seqname mapped unmapped index_type
#> 1       1     11        0        CSI
#> 2       2      1        0        CSI
#> 3       3      1        0        CSI
#> 4       4      2        0        CSI

bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
bam_index_path <- system.file("extdata", "range.bam.bai", package = "Rduckhts")
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
