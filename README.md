
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DuckHTS

A [DuckDB](https://duckdb.org/) extension (see the [DuckDB Extension
API](https://duckdb.org/docs/extensions/overview)) for reading
high-throughput sequencing (HTS) file formats using
[htslib](https://github.com/samtools/htslib).

Query VCF, BCF, BAM, CRAM, FASTA, FASTQ, GTF, GFF, and tabix-indexed
files directly using SQL.

Note: MSVC builds (windows_amd64/windows_arm64) are not supported. Use
MinGW/RTools for Windows.

## Functions

| Function                                                                      | Description                 | Schema                                                                                                          |
|-------------------------------------------------------------------------------|-----------------------------|-----------------------------------------------------------------------------------------------------------------|
| `read_bcf(path, [region, tidy_format])`                                       | Read VCF/BCF files          | CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO\_*, FORMAT\_*                                                      |
| `read_bam(path, [region, reference, standard_tags, auxiliary_tags])`          | Read SAM/BAM/CRAM files     | QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL, READ_GROUP_ID, SAMPLE_ID (+ SAMtags / AUX) |
| `read_fasta(path)`                                                            | Read FASTA files            | NAME, DESCRIPTION, SEQUENCE                                                                                     |
| `read_fastq(path, [mate_path, interleaved])`                                  | Read FASTQ files            | NAME, DESCRIPTION, SEQUENCE, QUALITY (+ MATE, PAIR_ID when paired/interleaved)                                  |
| `read_gff(path, [region, attributes_map])`                                    | Read GFF3 files             | seqname, source, feature, start, end, score, strand, frame, attributes (+ attributes_map MAP when enabled)      |
| `read_gtf(path, [region, attributes_map])`                                    | Read GTF files              | seqname, source, feature, start, end, score, strand, frame, attributes (+ attributes_map MAP when enabled)      |
| `read_tabix(path, [region, header, header_names, auto_detect, column_types])` | Read any tabix-indexed file | column0, column1, … (auto-detected)                                                                             |

Notes:

- `read_fastq` with `mate_path` requires exact QNAME pairing; mismatches
  raise an error.
- `read_bcf` region queries on files missing `##contig` headers return
  an empty result (with a warning) instead of erroring.
- `read_tabix` column inference now respects tabix header/meta
  configuration (meta char + line_skip), so header lines are not
  mistaken as data.
- `read_tabix` supports `header := true` to use the first non-meta line
  as column names, and `header_names := [...]` to override names.
- `read_tabix` supports `auto_detect := true` to infer numeric types and
  `column_types := [...]` to set types explicitly.
- `read_bam` supports `standard_tags := true` for typed SAMtags columns
  and `auxiliary_tags := true` for a string map of all remaining tags.

## Examples

``` sql
LOAD 'duckhts';

-- Read a VCF file (tidy FORMAT columns)
SELECT CHROM, POS, REF, ALT, QUAL, FORMAT_GT
FROM read_bcf('test/data/formatcols.vcf.gz', tidy_format := true)
WHERE CHROM = '1' AND POS > 1000000;

-- Region query on an indexed VCF (requires .tbi or .csi index)
SELECT * FROM read_bcf('test/data/vcf_file.bcf', region := '1:3000150-3000151');

-- VEP annotations (CSQ/BCSQ/ANN when present)
SELECT CHROM, POS, VEP_Allele, VEP_Consequence
FROM read_bcf('test/data/test_vep.vcf')
LIMIT 1;

-- Convert a VCF slice to Parquet using DuckDB COPY
COPY (
  SELECT *
  FROM read_bcf('test/data/vcf_file.bcf', region := '1:3000000-3000500')
) TO 'variants.parquet' (FORMAT PARQUET);

-- Read a BAM file
SELECT QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, READ_GROUP_ID, SAMPLE_ID
FROM read_bam('test/data/range.bam')
WHERE FLAG & 4 = 0;  -- mapped reads only

-- Region query on an indexed BAM
SELECT count(*) FROM read_bam('test/data/range.bam', region := 'CHROMOSOME_I:1-1000');

-- CRAM with explicit reference
SELECT count(*) FROM read_bam('test/data/range.cram', reference := 'test/data/ce.fa');

-- Read FASTA sequences
SELECT NAME, length(SEQUENCE) as seq_length
FROM read_fasta('test/data/ce.fa');

-- Read FASTQ and compute average quality
SELECT NAME, length(SEQUENCE) as read_length
FROM read_fastq('test/data/r1.fq');

-- Paired FASTQ (mate_path)
SELECT NAME, MATE, PAIR_ID
FROM read_fastq('test/data/r1.fq', mate_path := 'test/data/r2.fq');

-- Interleaved FASTQ
SELECT NAME, MATE, PAIR_ID
FROM read_fastq('test/data/interleaved.fq', interleaved := true);

-- Read GFF3 annotations
SELECT seqname, feature, start, "end", attributes, attributes_map
FROM read_gff('test/data/gff_file.gff.gz', attributes_map := true)
WHERE feature = 'gene';

-- Read a tabix-indexed BED file
SELECT * FROM read_tabix('test/data/gff_file.gff.gz', region := 'X:2934816-2935190');
```

## Remote URLs and HTS_PATH

Remote URLs (S3/GCS/HTTP/S) are supported when htslib is built with
plugins enabled. htslib loads these plugins from the directory specified
by the `HTS_PATH` environment variable. Set `HTS_PATH` **before**
loading the extension so the plugins can be discovered.

Example (using the packaged htslib plugins directory):

``` bash
export HTS_PATH=/usr/local/lib/R/site-library/Rduckhts/duckhts_extension/htslib/libexec/htslib
duckdb -unsigned <<'SQL'
LOAD '/root/duckhts/build/release/duckhts.duckdb_extension';
SELECT CHROM, COUNT(*) AS n
FROM read_bcf('s3://1000genomes-dragen-v3.7.6/data/cohorts/gvcf-genotyper-dragen-3.7.6/hg19/3202-samples-cohort/3202_samples_cohort_gg_chr22.vcf.gz',
              region := 'chr22:16050000-16050500')
GROUP BY CHROM;
SQL
#> ┌─────────┬───────┐
#> │  CHROM  │   n   │
#> │ varchar │ int64 │
#> ├─────────┼───────┤
#> │ chr22   │    11 │
#> └─────────┴───────┘
```

Notes:

- On Windows (MinGW/RTools), plugins are typically disabled, so remote
  URLs will not work.
- If `HTS_PATH` is set after the extension is already loaded, restart
  the session and set it first.

If you don’t have htslib plugins installed locally, download the
prebuilt htslib plugin bundles from the GitHub Releases page and point
`HTS_PATH` at the extracted htslib/libexec/htslib directory:
<https://github.com/RGenomicsETL/duckhts/releases>

### S3 credentials and configuration

The htslib S3 plugin supports credentials embedded in the URL or
provided via environment variables or standard credentials files. For
AWS-style credentials, the most common variables are:

- `AWS_ACCESS_KEY_ID`
- `AWS_SECRET_ACCESS_KEY`
- `AWS_SESSION_TOKEN` (optional, for temporary credentials)
- `AWS_DEFAULT_REGION`
- `AWS_PROFILE` / `AWS_DEFAULT_PROFILE`
- `AWS_SHARED_CREDENTIALS_FILE` (override credentials file location)

You can also configure htslib-specific settings like
`HTS_S3_ADDRESS_STYLE`, `HTS_S3_HOST`, and `HTS_S3_S3CFG` for
non-default S3 endpoints or path-style access.

See the htslib S3 plugin documentation for full details, URL syntax, and
short‑lived credentials support:
<https://www.htslib.org/doc/htslib-s3-plugin.html>

## Building

### Environment setup

Run the one-time configure step to create the Python venv and detect
platform settings:

``` bash
make configure
```

Note: MSVC builds (windows_amd64/windows_arm64) are not supported. Use
MinGW/RTools for Windows.

### Prerequisites

- C compiler (GCC or Clang)
- CMake ≥ 3.5
- Make
- Python 3 + venv
- Git
- [htslib](https://github.com/samtools/htslib) build dependencies: zlib,
  libbz2, liblzma, libdeflate, libcurl, libcrypto (OpenSSL)

On Debian/Ubuntu:

``` bash
sudo apt install build-essential cmake python3 python3-venv git \
    zlib1g-dev libbz2-dev liblzma-dev libdeflate-dev libcurl4-openssl-dev libssl-dev
```

On macOS:

``` bash
brew install cmake htslib xz libdeflate
```

### Vendor [htslib](https://github.com/samtools/htslib)

``` bash
./scripts/vendor_htslib.sh
```

This downloads and verifies [htslib](https://github.com/samtools/htslib)
1.23 into `third_party/htslib/`.

### Build

``` bash
make configure    # one-time setup (Python venv, platform detection)
make release      # build optimised extension
```

The build runs [htslib](https://github.com/samtools/htslib)’s Makefile
(`make lib-static`) in-tree.

The extension binary is written to
`build/release/duckhts.duckdb_extension`.

### Debug build

``` bash
make debug
```

## Loading

``` sql
-- Unsigned extensions must be loaded with -unsigned flag:
-- duckdb -unsigned

LOAD '/path/to/duckhts.duckdb_extension';
```

## Testing

SQL tests live in `test/sql/` using [DuckDB](https://duckdb.org/)’s
SQLLogicTest format.

Before running tests for the first time, prepare the indexed test data:

``` bash
./test/scripts/prepare_test_data.sh   # requires samtools, bcftools, bgzip, tabix
```

This copies files from the vendored
[htslib](https://github.com/samtools/htslib) test suite into
`test/data/` and builds the required indexes (BAI, CSI, TBI, FAI) so
region queries work.

Then run:

``` bash
make test_release
```

## R demo

The R package lives under r/Rduckhts and provides helpers to load the
extension and create [DuckDB](https://duckdb.org/) tables from HTS
files. See its README for R-specific usage:
[r/Rduckhts/README.Rmd](r/Rduckhts/README.Rmd).

``` r
library(DBI)
library(duckdb)

drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
con <- dbConnect(drv, dbdir = ":memory:")
ext_path <- normalizePath("build/release/duckhts.duckdb_extension", mustWork = FALSE)
dbExecute(con, sprintf("LOAD '%s'", ext_path))
#> [1] 0

dbGetQuery(con, "
  SELECT *
  FROM read_bcf('test/data/formatcols.vcf.gz', tidy_format := true)
  LIMIT 5
")
#>   CHROM POS ID REF ALT QUAL FILTER SAMPLE_ID  FORMAT_S
#> 1     1 100  a   A   T   NA   PASS        S1         a
#> 2     1 100  a   A   T   NA   PASS        S²   bbbbbbb
#> 3     1 100  a   A   T   NA   PASS        S3 ccccccccc

parquet_path <- tempfile(fileext = ".parquet")
dbExecute(con, sprintf(
  "COPY (SELECT * FROM read_bcf('test/data/formatcols.vcf.gz', tidy_format := true)) TO '%s' (FORMAT PARQUET)",
  parquet_path
))
#> [1] 3
file.exists(parquet_path)
#> [1] TRUE

dbGetQuery(con, "
  SELECT NAME, SEQUENCE, QUALITY, MATE, PAIR_ID
  FROM read_fastq('test/data/r1.fq', mate_path := 'test/data/r2.fq')
  LIMIT 5
")
#>                              NAME
#> 1 HS25_09827:2:1201:1505:59795#49
#> 2 HS25_09827:2:1201:1505:59795#49
#> 3 HS25_09827:2:1201:1559:70726#49
#> 4 HS25_09827:2:1201:1559:70726#49
#> 5 HS25_09827:2:1201:1564:39627#49
#>                                                                                               SEQUENCE
#> 1 CCGTTAGAGCATTTGTTGAAAATGCTTTCCTTGCTCCATGTGATGACTCTGGTGCCCTTGTCAAAAGCCAGCTGGGCCTATTCGTGTGGGTCTGTTTCTG
#> 2 AAGGAAAGAAGGGAGGGAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAAGTAGGAAGAATTCATCTACCCAATT
#> 3 TTGTTAAAATGACCATACCCAAAGTGATCTACAGACTCAATACAATTTCTATTGAAATACCAATCACACTCTTCACAGAACTAGAAAAACAGTTCTAAAA
#> 4 TTTTCTTTTATTAATTTTATACTTACATTTAAGTCTTTATTCCATTTTGAGTCAATGTTTGTATATGATGAGAGATAGGGGTCTAGTTTCATACTTCTAC
#> 5 ACGCGGCAATCCAATGTGTGAGTTGAGAAGCGGTGAGGAGGGAATCCTAATTTTATGAGCAGGTCAGGACCGTGGGAGATACCTGACACCTGAGATGGTA
#>                                                                                                QUALITY
#> 1 CABCFGDEEFFEFHGHGGFFGDIGIJFIFHHGHEIFGHBCGHDIFBE9GIAICGGICFIBFGGHGDGGGHE?GIGDFGGHEGIEJG>;FG<GGHACEFGH
#> 2 <CBB>DCHFEFBHAGCGACF7CJI8HBIIEFGFEBG?DCGA?ACFGGI=BEDG?EFEHFFFEHFD?HG+DFH>FFHGFBFE4F@I3HF@>A5F?GFH<EG
#> 3 CABEFGFFGFHGGGGJGGFFGKIHHJFIEHHHGIEGGEHJGHDHFGHIGICIJEFIFGIF8GGHKFHGGFEI6GGGFIGHGGIE>EFCFHGGGHEJEAJE
#> 4 ;CBCEFDHDGFGHDGDIGEF@EJIIGEEIECGFHGFHGGGHHHHGGKIFFEHGEGHFIEFFHHGDHHGJEGF?FBHFFGCHHFFII>GCFCFFGGCEBF?
#> 5 BACCFGBFGFHGGJGHGGFEGHIGIJHFEH:HHEHGHHBGGH9IAGHGFHIFJFFAFGIFDIGHKEIG<C>F,CGD66?7EFI5EEG>EGGGGD5=HH6E
#>   MATE                         PAIR_ID
#> 1    1 HS25_09827:2:1201:1505:59795#49
#> 2    2 HS25_09827:2:1201:1505:59795#49
#> 3    1 HS25_09827:2:1201:1559:70726#49
#> 4    2 HS25_09827:2:1201:1559:70726#49
#> 5    1 HS25_09827:2:1201:1564:39627#49

dbGetQuery(con, "
  SELECT QNAME, RNAME, POS, READ_GROUP_ID, SAMPLE_ID
  FROM read_bam('test/data/rg.sam.gz')
  LIMIT 5
")
#>   QNAME RNAME POS READ_GROUP_ID SAMPLE_ID
#> 1    a1    xx   1            x1        x1
#> 2    b1    xx   1            x2        x2
#> 3    c1    xx   1          <NA>      <NA>
#> 4    a2    xx  11            x1        x1
#> 5    b2    xx  11            x2        x2
```

### SAMtags + auxiliary tags

Standard SAMtags can be surfaced as typed columns and non-standard tags
captured in a map for ad hoc access:

``` r
dbGetQuery(con, "
  SELECT RG, NM, map_extract(AUXILIARY_TAGS, 'XZ') AS XZ
  FROM read_bam('test/data/aux_tags.sam.gz', standard_tags := true, auxiliary_tags := true)
  LIMIT 1
")
#>   RG NM  XZ
#> 1 x1  2 foo

dbGetQuery(con, "
  SELECT seqname, feature, start, \"end\", attributes_map
  FROM read_gff('test/data/gff_file.gff.gz', attributes_map := true)
  WHERE feature = 'gene'
  LIMIT 5
")
#>   seqname feature   start     end
#> 1       X    gene 2934816 2964270
#>                                                              attributes_map
#> 1 ID, Name, biotype, OTTHUMG00000137358, OTTHUMG00000137358, protein_coding

dbGetQuery(con, "
  SELECT column0, column1
  FROM read_tabix('test/data/meta_tabix.tsv.gz')
  LIMIT 2
")
#>   column0 column1
#> 1    chr1       1
#> 2    chr1       2

dbGetQuery(con, "
  SELECT chrom, pos
  FROM read_tabix('test/data/header_tabix.tsv.gz', header := true)
  LIMIT 2
")
#>   chrom pos
#> 1  chr1   1
#> 2  chr1   2

dbGetQuery(con, "
  SELECT typeof(column1) AS column1_type
  FROM read_tabix('test/data/meta_tabix.tsv.gz', auto_detect := true)
  LIMIT 1
")
#>   column1_type
#> 1       BIGINT

dbGetQuery(con, "
  SELECT pos + 1 AS pos_plus_one
  FROM read_tabix('test/data/header_tabix.tsv.gz', header := true,
                  column_types := ['VARCHAR','BIGINT','VARCHAR'])
  LIMIT 1
")
#>   pos_plus_one
#> 1            2

dbDisconnect(con, shutdown = TRUE)
```

Rendering this document requires a built extension at
`build/release/duckhts.duckdb_extension`.

## Project Structure

    src/
      duckhts.c          # Extension entry point
      bcf_reader.c       # VCF/BCF reader (read_bcf)
      bam_reader.c       # SAM/BAM/CRAM reader (read_bam)
      seq_reader.c       # FASTA/FASTQ reader (read_fasta, read_fastq)
      tabix_reader.c     # Tabix/GTF/GFF reader (read_tabix, read_gtf, read_gff)
      vep_parser.c       # VEP/CSQ annotation parser
      include/
        vcf_types.h
        vep_parser.h
    third_party/
      htslib/            # Vendored htslib 1.23 (built automatically)
    test/
      sql/               # SQL logic tests
    duckdb_capi/
      duckdb.h           # DuckDB C API headers
      duckdb_extension.h
    r/
      Rduckhts/          # R package harness

## References

- DuckDB: <https://duckdb.org/>
- DuckDB Extension API: <https://duckdb.org/docs/extensions/overview>
- DuckDB extension template (C):
  <https://github.com/duckdb/extension-template-c>
- htslib: <https://github.com/samtools/htslib>
- RBCFTools: <https://github.com/RGenomicsETL/RBCFTools>

## License

MIT
