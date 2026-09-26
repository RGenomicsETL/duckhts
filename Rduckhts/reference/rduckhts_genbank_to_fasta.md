# Write GenBank Sequence as FASTA

Writes the ORIGIN sequence of each record in a GenBank flat file to a
FASTA file using the DuckHTS extension. Records are written under the
same name `rduckhts_genbank` reports as `seqname`, so feature
coordinates land on the contig of that name; the DEFINITION follows the
name without its trailing period, as in NCBI's FASTA export. The FASTA
is written to a temporary file beside `output_path` and renamed into
place only after the input has been read to a clean end, so a failure
never leaves a partial output and an existing file is never lost.
Records without an ORIGIN block are skipped; zero written records is an
error.

## Usage

``` r
rduckhts_genbank_to_fasta(
  con,
  path,
  output_path = NULL,
  line_width = 70,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- path:

  Path to the GenBank flat file, optionally bgzipped

- output_path:

  Optional explicit output path for the FASTA file, one non-empty
  string; defaults to `path` with `.fa` appended

- line_width:

  Sequence characters per output line, a whole number between 1 and
  `.Machine$integer.max`

- overwrite:

  Logical. If `TRUE`, replace an existing output file; otherwise an
  existing output is an error

## Value

A data frame with columns `success`, `output_path` and `records_written`
