# Read HTS Index Metadata

Reads index metadata from HTS-supported index files via DuckHTS.

## Usage

``` r
rduckhts_hts_index(con, path, format = NULL, index_path = NULL)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- path:

  Path to input HTS file

- format:

  Optional format hint (e.g., "auto", "vcf", "bcf", "bam", "cram",
  "tabix")

- index_path:

  Optional explicit path to index file

## Value

A data frame with index metadata.
