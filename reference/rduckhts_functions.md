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
#> 44     seq_revcomp
#> 45   seq_canonical
#> 46   seq_hash_2bit
#> 47 seq_encode_4bit
#> 48 seq_decode_4bit
#> 49  seq_gc_content
#> 50       seq_kmers
#>                                                                                              description
#> 44                       Compute the reverse complement of a DNA sequence using A, C, G, T, and N bases.
#> 45                        Return the lexicographically smaller of a sequence and its reverse complement.
#> 46                                         Encode a short DNA sequence as a 2-bit unsigned integer hash.
#> 47 Encode an IUPAC DNA sequence as a list of 4-bit base codes, preserving ambiguity symbols including N.
#> 48                              Decode a list of 4-bit IUPAC DNA base codes back into a sequence string.
#> 49                                    Compute GC fraction for a DNA sequence as a value between 0 and 1.
#> 50                              Expand a sequence into positional k-mers with optional canonicalization.
subset(rduckhts_functions(kind = "table"), select = c("name", "r_wrapper"))
#>                              name                        r_wrapper
#> 1                        read_bcf                     rduckhts_bcf
#> 2                        read_bam                     rduckhts_bam
#> 3                     read_pileup                  rduckhts_pileup
#> 4                      read_fasta                   rduckhts_fasta
#> 5                        read_bed                     rduckhts_bed
#> 6                       fasta_nuc               rduckhts_fasta_nuc
#> 7       duckhts_cgranges_overlaps                                 
#> 8  duckhts_cgranges_overlaps_bulk                                 
#> 9                      read_fastq                   rduckhts_fastq
#> 10        detect_quality_encoding rduckhts_detect_quality_encoding
#> 11                       read_gff                     rduckhts_gff
#> 12                       read_gtf                     rduckhts_gtf
#> 13                     read_tabix                   rduckhts_tabix
#> 14                    fasta_index             rduckhts_fasta_index
#> 15                          bgzip                   rduckhts_bgzip
#> 16                        bgunzip                 rduckhts_bgunzip
#> 17                      bam_index               rduckhts_bam_index
#> 18                      bcf_index               rduckhts_bcf_index
#> 19                    tabix_index             rduckhts_tabix_index
#> 20                 bam_bin_counts          rduckhts_bam_bin_counts
#> 21       duckhts_bam_bed_coverage        rduckhts_bam_bed_coverage
#> 22               duckhts_mosdepth                rduckhts_mosdepth
#> 23      duckhts_samtools_idxstats       rduckhts_samtools_idxstats
#> 24                read_hts_header              rduckhts_hts_header
#> 25                 read_hts_index               rduckhts_hts_index
#> 26           read_hts_index_spans         rduckhts_hts_index_spans
#> 27                 bcftools_score                   rduckhts_score
#> 28                      seq_kmers                                 
```
