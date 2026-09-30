# Rduckhts

[![CRAN
Status](https://www.r-pkg.org/badges/version/Rduckhts)](https://cran.r-project.org/package=Rduckhts)
[![R-universe
version](https://RGenomicsETL.r-universe.dev/Rduckhts/badges/version)](https://RGenomicsETL.r-universe.dev/Rduckhts)

Query sequencing files from R with SQL. Rduckhts bundles
[DuckHTS](https://github.com/RGenomicsETL/duckhts), a DuckDB extension
built on htslib. It reads the following formats as DuckDB tables:

- VCF/BCF
- BAM/CRAM
- FASTA/FASTQ
- GFF/GTF
- GenBank
- BED
- BigWig
- tabix-indexed text

Indexed files answer region queries without a full scan, and remote
files are read over HTTPS and S3. Results stay in DuckDB until you
collect them. htslib is bundled, so there is nothing else to install.

## Installation

``` r
install.packages("Rduckhts")

# Development version
install.packages("Rduckhts", repos = c("https://rgenomicsetl.r-universe.dev",
                                       "https://cloud.r-project.org"))
```

Rduckhts needs the `duckdb` R package, version 1.5.0 or newer.

## Quick start

Every reader has the same shape,
`rduckhts_<format>(con, table_name, path, ...)`. It creates a DuckDB
table that you query with SQL:

``` r
library(DBI)
library(Rduckhts)

con <- rduckhts_connect()

bcf <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
rduckhts_bcf(con, "variants", bcf)

dbGetQuery(con, "
  SELECT CHROM, count(*) AS variants, min(POS) AS first, max(POS) AS last
  FROM variants GROUP BY CHROM ORDER BY CHROM
")
#>   CHROM variants   first    last
#> 1     1       11 3000150 3184885
#> 2     2        1 3199812 3199812
#> 3     3        1 3212016 3212016
#> 4     4        2 3258448 3258501
```

## Common tasks

### Read one region of an indexed file

Pass `region` to any indexed reader. Only the matching records are
decoded:

``` r
rduckhts_bcf(con, "variants_region", bcf, region = "1:3000150-3000151")
dbGetQuery(con, "SELECT CHROM, POS, REF, ALT FROM variants_region")
#>   CHROM     POS REF ALT
#> 1     1 3000150   C   T
#> 2     1 3000151   C   T

bam <- system.file("extdata", "range.bam", package = "Rduckhts")
rduckhts_bam(con, "reads", bam, region = "CHROMOSOME_I:1-1000")
dbGetQuery(con, "SELECT QNAME, FLAG, POS, MAPQ, CIGAR FROM reads")
#>                           QNAME FLAG POS MAPQ    CIGAR
#> 1 HS18_09653:4:1315:19857:61712  145 914   23 78M1D22M
#> 2 HS18_09653:4:1308:11522:27107  161 934    0 58M1D42M
```

### Fetch requested alleles, exactly

[`rduckhts_geno_sites()`](https://rgenomicsetl.github.io/duckhts/Rduckhts/reference/rduckhts_geno_sites.md)
takes a table of requested alleles and returns each one’s genotypes from
a single indexed scan. A requested ALT inside a multi-allelic record
comes back with its index in `ALT`. Requests without an exact match stay
in the result with a `match_status` and no fabricated genotype.

``` r
geno <- system.file("extdata", "geno_sites.bcf", package = "Rduckhts")
dbWriteTable(con, "requests", temporary = TRUE, data.frame(
  request_id = 1:4, build = "GRCh38", chrom = c("chr1", "chr1", "chr1", "chr3"),
  pos = c(100L, 200L, 300L, 10L), ref = c("A", "T", "G", "C"),
  alt = c("G", "C", "A", "T")
))

fetched <- rduckhts_geno_sites(con, geno, "requests", source_build = "GRCh38",
                               format_fields = "DS")
fetched <- fetched[order(fetched$request_id), ]
rownames(fetched) <- NULL
fetched[, c("request_id", "chrom", "pos", "alt", "alt_index", "match_status")]
#>   request_id chrom pos alt alt_index       match_status
#> 1          1  chr1 100   G         2            matched
#> 2          1  chr1 100   G         2            matched
#> 3          2  chr1 200   C         1            matched
#> 4          3  chr1 300   A        NA allele_not_at_site
#> 5          4  chr3  10   T        NA             absent
```

The test file holds two records at chr1:100, so request 1 matches twice,
once per record. Request 3’s REF matches but its ALT is not at that
site. Request 4’s position is not in the file.

### Keep annotation attributes as columns

GFF, GTF and GenBank readers can promote chosen attributes to typed
columns:

``` r
gff <- system.file("extdata", "gff_attrs.gff3", package = "Rduckhts")
rduckhts_gff(con, "features", gff, attributes = c("ID", "Dbxref"))
dbGetQuery(con, "SELECT seqname, feature, start, \"end\", ID, Dbxref FROM features")
#>   seqname feature start end ID               Dbxref
#> 1    chr1    gene     1  10  g GeneID:1,HGNC:HGNC:1
```

### Read remote files

Remote paths work wherever a local path does, once
[`setup_hts_env()`](https://rgenomicsetl.github.io/duckhts/Rduckhts/reference/setup_hts_env.md)
has pointed htslib at its bundled network plugins. Call it once per
session before the first remote read. This reads a 100 kb slice of the
UCSC GRCh38 phyloP track and one region of a 1000 Genomes cohort VCF on
S3, fetching only the index and the blocks that the region needs:

``` r
setup_hts_env()

phylop <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw"
rduckhts_bigwig(con, "conservation", phylop, region = "chr22:20000000-20099999")
dbGetQuery(con, "
  SELECT count(*) AS intervals, round(min(VALUE), 3) AS min, round(max(VALUE), 3) AS max
  FROM conservation
")
#>   intervals     min   max
#> 1     96783 -10.787 9.602

cohort <- paste0(
  "s3://1000genomes-dragen-v3.7.6/data/cohorts/gvcf-genotyper-dragen-3.7.6/",
  "hg19/3202-samples-cohort/3202_samples_cohort_gg_chr22.vcf.gz"
)
rduckhts_bcf(con, "cohort", cohort, region = "chr22:16050000-16050500")
dbGetQuery(con, "SELECT count(*) AS variants FROM cohort")
#>   variants
#> 1       11
```

S3 access needs a build with htslib’s S3 support; Windows (Rtools)
builds do not include it.

### Write Parquet

Converters stream a file into Parquet and keep its header in the Parquet
metadata:

``` r
parquet <- tempfile(fileext = ".parquet")
rduckhts_bcf_convert_parquet(con, bcf, parquet, overwrite = TRUE)
dbGetQuery(con, sprintf("SELECT count(*) AS rows FROM read_parquet('%s')", parquet))
#>   rows
#> 1   15
```

### Work with a database file

`rduckhts_connect(dbdir = "file.duckdb")` works with database files too,
even read-only ones, and never changes the file’s catalog.

- **Other connections:** for another DBI or pool connection, call
  `rduckhts_load(con)`, then `rduckhts_install_macros(con)`.
- **SQL functions:** the extension’s functions (`read_bcf()`,
  `read_bam()` and the rest) are available directly on the connection.

## Learn more

- **Function reference:**
  [reference](https://rgenomicsetl.github.io/duckhts/Rduckhts/reference/),
  or
  [`rduckhts_functions()`](https://rgenomicsetl.github.io/duckhts/Rduckhts/reference/rduckhts_functions.md)
  in R for the SQL catalog.
- **[Worked
  examples](https://rgenomicsetl.github.io/duckhts/Rduckhts/articles/workflows.html):**
  - variant normalization and liftover;
  - polygenic scores;
  - Somalier-style relatedness;
  - coverage and bin counts;
  - FASTQ quality control;
  - interval indexes;
  - compression and indexing;
  - multi-file reading.
- **[Browser and
  webR](https://rgenomicsetl.github.io/duckhts/Rduckhts/articles/browser.html):**
  remote access and request settings in wasm builds.
- **SQL, the DuckDB command line and Python:** the [DuckHTS
  README](https://github.com/RGenomicsETL/duckhts).
- **Consequence annotation:** the
  [DuckVEP](https://github.com/RGenomicsETL/DuckVEP) extension and its
  Rduckvep package.
- **Building another package against the bundled htslib:**
  [`rduckhts_htslib_config()`](https://rgenomicsetl.github.io/duckhts/Rduckhts/reference/rduckhts_htslib_config.md).

## License

GPL (\>= 2). The bundled DuckHTS extension is MIT-licensed, and other
bundled third-party code keeps its own licence; see `inst/COPYRIGHT`.

## Credits

The GenBank reader and FASTA converter were contributed by [Ryan
Ward](https://github.com/ryandward) of [Nurture
Bio](https://github.com/Nurture-Bio). Thanks to all
[contributors](https://github.com/RGenomicsETL/duckhts/graphs/contributors),
and to the htslib, bcftools, DuckDB and RBCFTools projects that DuckHTS
builds on.
