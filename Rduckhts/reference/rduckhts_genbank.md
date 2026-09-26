# Create GenBank Feature Table

Creates a DuckDB table from a GenBank flat file using the DuckHTS
extension. Features are emitted in `read_gff`'s column shape, so a
GenBank record substitutes for a GFF without a schema change. Locations
built with `join()` or [`order()`](https://rdrr.io/r/base/order.html)
give one row per segment in biological order, `complement(...)` sets
strand `"-"`, and the GFF3 phase of each CDS segment is carried from
`/codon_start` across segments. `Parent` links a feature to the gene
sharing its `/locus_tag` wherever that gene appears in the record,
repeated qualifiers become one key with comma-joined values, and
valueless qualifiers read `true`. Records stream one at a time, and a
record without a terminating `//` or with a malformed location is an
error naming the feature and line.

## Usage

``` r
rduckhts_genbank(
  con,
  table_name = NULL,
  path,
  attributes_map = FALSE,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- table_name:

  Name for the created table, or `NULL` to create the `genbank_data`
  view

- path:

  Path to the GenBank flat file, optionally bgzipped

- attributes_map:

  Logical. If `TRUE`, add a parsed `MAP(VARCHAR, VARCHAR)` column
  alongside the raw attribute string

- overwrite:

  Logical. If TRUE, overwrites an existing table

## Value

Invisible TRUE on success
