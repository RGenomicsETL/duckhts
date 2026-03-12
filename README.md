
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

## Extension Function Catalog

This section is generated from `functions.yaml`.

### Readers

| Function      | Kind  | Returns | R helper               | Description                                                                             |
|---------------|-------|---------|------------------------|-----------------------------------------------------------------------------------------|
| `read_bcf`    | table | table   | `rduckhts_bcf`         | Read VCF and BCF variant data with typed INFO, FORMAT, and optional tidy sample output. |
| `read_bam`    | table | table   | `rduckhts_bam`         | Read SAM, BAM, and CRAM alignments with optional typed SAMtags and auxiliary tag maps.  |
| `read_fasta`  | table | table   | `rduckhts_fasta`       | Read FASTA records or indexed FASTA regions as sequence rows.                           |
| `read_fastq`  | table | table   | `rduckhts_fastq`       | Read single-end, paired-end, or interleaved FASTQ files.                                |
| `read_gff`    | table | table   | `rduckhts_gff`         | Read GFF annotations with optional parsed attribute maps and indexed region filtering.  |
| `read_gtf`    | table | table   | `rduckhts_gtf`         | Read GTF annotations with optional parsed attribute maps and indexed region filtering.  |
| `read_tabix`  | table | table   | `rduckhts_tabix`       | Read generic tabix-indexed text data with optional header handling and type inference.  |
| `fasta_index` | table | table   | `rduckhts_fasta_index` | Build a FASTA index and return the index path used by the operation.                    |

### Metadata

| Function               | Kind        | Returns | R helper                   | Description                                                                             |
|------------------------|-------------|---------|----------------------------|-----------------------------------------------------------------------------------------|
| `read_hts_header`      | table       | table   | `rduckhts_hts_header`      | Inspect HTS headers in parsed, raw, or combined form across supported formats.          |
| `read_hts_index`       | table       | table   | `rduckhts_hts_index`       | Inspect high-level HTS index metadata such as sequence names and mapped counts.         |
| `read_hts_index_spans` | table_macro | table   | `rduckhts_hts_index_spans` | Expand index metadata into span and chunk rows suitable for low-level index inspection. |
| `read_hts_index_raw`   | table_macro | table   | `rduckhts_hts_index_raw`   | Return the raw on-disk HTS index blob together with basic identifying metadata.         |

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

| Function                       | Kind   | Returns | R helper | Description                                                             |
|--------------------------------|--------|---------|----------|-------------------------------------------------------------------------|
| `is_segmented`                 | scalar | BOOLEAN |          | Test whether the SAM flag marks a read as part of a segmented template. |
| `is_properly_aligned`          | scalar | BOOLEAN |          | Test whether the SAM flag indicates a properly aligned read pair.       |
| `is_properly_segmented`        | scalar | BOOLEAN |          | Alias for is_properly_aligned(flag) using segmented-read terminology.   |
| `is_unmapped`                  | scalar | BOOLEAN |          | Test whether the read itself is unmapped according to the SAM flag.     |
| `is_mate_unmapped`             | scalar | BOOLEAN |          | Test whether the mate read is flagged as unmapped.                      |
| `is_reverse_complemented`      | scalar | BOOLEAN |          | Test whether the read is aligned to the reverse strand.                 |
| `is_mate_reverse_complemented` | scalar | BOOLEAN |          | Test whether the mate read is aligned to the reverse strand.            |
| `is_first_segment`             | scalar | BOOLEAN |          | Test whether the read is marked as the first segment in the template.   |
| `is_last_segment`              | scalar | BOOLEAN |          | Test whether the read is marked as the last segment in the template.    |
| `is_secondary`                 | scalar | BOOLEAN |          | Test whether the alignment is marked as secondary.                      |
| `is_qc_fail`                   | scalar | BOOLEAN |          | Test whether the read failed vendor or pipeline quality checks.         |
| `is_duplicate`                 | scalar | BOOLEAN |          | Test whether the alignment is flagged as a duplicate.                   |
| `is_supplementary`             | scalar | BOOLEAN |          | Test whether the alignment is marked as supplementary.                  |

`read_fastq` with `mate_path` requires exact QNAME pairing. `read_bam`
supports typed `standard_tags` and `auxiliary_tags` maps. `read_tabix`
supports header-aware parsing (`header`, `header_names`) and optional
type inference (`auto_detect`, `column_types`). Region lists in
comma-separated form are supported by `read_bam`, `read_bcf`,
`read_fasta`, `read_gff`, `read_gtf`, and `read_tabix`. `read_bam`
multi-region queries are deduplicated by htslib, while
`read_bcf`/`read_fasta`/`read_gff`/`read_gtf`/`read_tabix` chain regions
and can return duplicates for overlaps.

## Examples

The examples below run directly against bundled local test files and
show the main reader APIs using the `R` `DBI` and `duckdb` packages.

``` r
library(DBI)
library(duckdb)

drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
con <- dbConnect(drv, dbdir = ":memory:")
ext_path <- normalizePath("build/release/duckhts.duckdb_extension", mustWork = FALSE)
dbExecute(con, sprintf("LOAD '%s'", ext_path))
#> [1] 0

dbGetQuery(con, "
  SELECT CHROM, POS, REF, ALT, SAMPLE_ID
  FROM read_bcf('test/data/formatcols.vcf.gz', tidy_format := true)
  LIMIT 3
")
#>   CHROM POS REF ALT SAMPLE_ID
#> 1     1 100   A   T        S1
#> 2     1 100   A   T        S²
#> 3     1 100   A   T        S3

dbGetQuery(con, "
  SELECT count(*) AS n
  FROM read_bam('test/data/range.bam', region := 'CHROMOSOME_I:1-1000')
")
#>   n
#> 1 2

dbGetQuery(con, "
  SELECT * FROM fasta_index('test/data/ce.fa')
")
#>   success index_path
#> 1    TRUE

dbGetQuery(con, "
  SELECT NAME, length(SEQUENCE) AS seq_length
  FROM read_fasta('test/data/ce.fa', region := 'CHROMOSOME_I:1-25')
")
#>           NAME seq_length
#> 1 CHROMOSOME_I         25

dbGetQuery(con, "
  SELECT NAME, MATE, PAIR_ID
  FROM read_fastq('test/data/interleaved.fq', interleaved := true)
  LIMIT 3
")
#>                              NAME MATE                         PAIR_ID
#> 1 HS25_09827:2:1201:1505:59795#49    1 HS25_09827:2:1201:1505:59795#49
#> 2 HS25_09827:2:1201:1505:59795#49    2 HS25_09827:2:1201:1505:59795#49
#> 3 HS25_09827:2:1201:1559:70726#49    1 HS25_09827:2:1201:1559:70726#49

dbDisconnect(con, shutdown = TRUE)
```

## Remote URLs and HTS_PATH

Remote URLs (S3/GCS/HTTP/S) can work in two htslib build modes:

1.  Dynamic plugin mode (`ENABLE_PLUGINS`): remote handlers are loaded
    from `HTS_PATH`.
2.  Static-handler mode (plugins disabled): handlers are compiled into
    `libhts` and `HTS_PATH` is not needed.

Use `HTS_PATH` only when you want dynamic plugin discovery (for example,
to point at an external htslib plugin directory).

Example (works in static-handler mode and plugin mode):

``` bash
# Not run by default because it requires network access and a built extension.
extension_path=build/release/duckhts.duckdb_extension
duckdb -unsigned <<SQL
LOAD '${extension_path}';
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

If you need dynamic plugin mode, set `HTS_PATH` before loading the
extension, for example:

``` bash
export HTS_PATH=$(Rscript --quiet -e 'cat(Rduckhts:::duckhts_htslib_plugins_dir(),sep="")')
```

If `HTS_PATH` is changed after loading, restart the session and reload
the extension.

If you don’t have htslib plugins installed locally, download the
prebuilt binaries from the r-universe-binaries GitHub release and point
`HTS_PATH` at the extracted htslib/libexec/htslib directory inside the
package bundle.
<https://github.com/RGenomicsETL/duckhts/releases/tag/r-universe-binaries>

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

The R package lives under `r/Rduckhts` and provides helpers to load the
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

dbGetQuery(con, "
  SELECT idx, raw
  FROM read_hts_header('test/data/formatcols.vcf.gz', mode := 'raw')
  LIMIT 3
")
#>   idx                                                 raw
#> 1   0                                ##fileformat=VCFv4.3
#> 2   1 ##FILTER=<ID=PASS,Description="All filters passed">
#> 3   2                                     ##contig=<ID=1>

dbGetQuery(con, "
  SELECT seqname, tid, index_type, chunk_beg_vo, chunk_end_vo
  FROM read_hts_index_spans('test/data/formatcols.vcf.gz')
  LIMIT 3
")
#>   seqname tid index_type chunk_beg_vo chunk_end_vo
#> 1       1   0        CSI           NA           NA

dbGetQuery(con, "
  SELECT index_type, octet_length(raw) AS raw_bytes
  FROM read_hts_index_raw('test/data/formatcols.vcf.gz')
")
#>   index_type raw_bytes
#> 1        CSI        30
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
