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
#> 32     seq_revcomp
#> 33   seq_canonical
#> 34   seq_hash_2bit
#> 35 seq_encode_4bit
#> 36 seq_decode_4bit
#> 37  seq_gc_content
#> 38       seq_kmers
#>                                                                                              description
#> 32                       Compute the reverse complement of a DNA sequence using A, C, G, T, and N bases.
#> 33                        Return the lexicographically smaller of a sequence and its reverse complement.
#> 34                                         Encode a short DNA sequence as a 2-bit unsigned integer hash.
#> 35 Encode an IUPAC DNA sequence as a list of 4-bit base codes, preserving ambiguity symbols including N.
#> 36                              Decode a list of 4-bit IUPAC DNA base codes back into a sequence string.
#> 37                                    Compute GC fraction for a DNA sequence as a value between 0 and 1.
#> 38                              Expand a sequence into positional k-mers with optional canonicalization.
subset(rduckhts_functions(kind = "table"), select = c("name", "r_wrapper"))
#>                         name                        r_wrapper
#> 1                   read_bcf                     rduckhts_bcf
#> 2                   read_bam                     rduckhts_bam
#> 3                 read_fasta                   rduckhts_fasta
#> 4                   read_bed                     rduckhts_bed
#> 5                  fasta_nuc               rduckhts_fasta_nuc
#> 6                 read_fastq                   rduckhts_fastq
#> 7    detect_quality_encoding rduckhts_detect_quality_encoding
#> 8                   read_gff                     rduckhts_gff
#> 9                   read_gtf                     rduckhts_gtf
#> 10                read_tabix                   rduckhts_tabix
#> 11               fasta_index             rduckhts_fasta_index
#> 12                     bgzip                   rduckhts_bgzip
#> 13                   bgunzip                 rduckhts_bgunzip
#> 14                 bam_index               rduckhts_bam_index
#> 15                 bcf_index               rduckhts_bcf_index
#> 16               tabix_index             rduckhts_tabix_index
#> 17            bam_bin_counts          rduckhts_bam_bin_counts
#> 18  duckhts_bam_bed_coverage        rduckhts_bam_bed_coverage
#> 19          duckhts_mosdepth                rduckhts_mosdepth
#> 20 duckhts_samtools_idxstats       rduckhts_samtools_idxstats
#> 21           read_hts_header              rduckhts_hts_header
#> 22            read_hts_index               rduckhts_hts_index
#> 23      read_hts_index_spans         rduckhts_hts_index_spans
#> 24            bcftools_score                   rduckhts_score
#> 25                 seq_kmers                                 
```
