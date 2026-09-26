# Detect FASTQ Quality Encoding

Inspects a FASTQ file's observed quality ASCII range and reports
compatible legacy encodings with a heuristic guessed encoding.

## Usage

``` r
rduckhts_detect_quality_encoding(con, path, max_records = 10000)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- path:

  Path to the FASTQ file

- max_records:

  Maximum number of records to inspect

## Value

A data frame with the detected quality encoding summary
