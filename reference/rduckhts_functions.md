# List DuckHTS Extension Functions

Returns the package-bundled function catalog generated from the
top-level `functions.yaml` manifest in the duckhts repository.

## Usage

``` r
rduckhts_functions(category = NULL, kind = NULL)
```

## Arguments

- category:

  Optional function category filter.

- kind:

  Optional function kind filter such as `"scalar"`, `"table"`, or
  `"table_macro"`.

## Value

A data frame describing the extension functions, including the DuckDB
function name, kind, category, signature, return type, optional R helper
wrapper, short description, and example SQL.

## Examples

``` r
catalog <- rduckhts_functions()
subset(catalog, category == "Sequence UDFs", select = c("name", "description"))
#>               name
#> 13     seq_revcomp
#> 14   seq_canonical
#> 15   seq_hash_2bit
#> 16 seq_encode_4bit
#> 17 seq_decode_4bit
#> 18  seq_gc_content
#> 19       seq_kmers
#>                                                                                              description
#> 13                       Compute the reverse complement of a DNA sequence using A, C, G, T, and N bases.
#> 14                        Return the lexicographically smaller of a sequence and its reverse complement.
#> 15                                         Encode a short DNA sequence as a 2-bit unsigned integer hash.
#> 16 Encode an IUPAC DNA sequence as a list of 4-bit base codes, preserving ambiguity symbols including N.
#> 17                              Decode a list of 4-bit IUPAC DNA base codes back into a sequence string.
#> 18                                    Compute GC fraction for a DNA sequence as a value between 0 and 1.
#> 19                              Expand a sequence into positional k-mers with optional canonicalization.
subset(rduckhts_functions(kind = "table"), select = c("name", "r_wrapper"))
#>               name            r_wrapper
#> 1         read_bcf         rduckhts_bcf
#> 2         read_bam         rduckhts_bam
#> 3       read_fasta       rduckhts_fasta
#> 4       read_fastq       rduckhts_fastq
#> 5         read_gff         rduckhts_gff
#> 6         read_gtf         rduckhts_gtf
#> 7       read_tabix       rduckhts_tabix
#> 8      fasta_index rduckhts_fasta_index
#> 9  read_hts_header  rduckhts_hts_header
#> 10  read_hts_index   rduckhts_hts_index
#> 11       seq_kmers                     
```
