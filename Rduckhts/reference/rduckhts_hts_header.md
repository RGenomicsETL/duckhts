# Read HTS Header Metadata

Reads file header records from HTS-supported formats using the DuckHTS
extension.

## Usage

``` r
rduckhts_hts_header(con, path, format = NULL, mode = NULL)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- path:

  Path to input HTS file

- format:

  Optional format hint (e.g., "auto", "vcf", "bcf", "bam", "cram",
  "tabix")

- mode:

  Header output mode: "parsed" (default), "raw", or "both"

## Value

A data frame with parsed header metadata.
