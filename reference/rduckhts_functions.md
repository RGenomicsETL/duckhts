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
#> 88     seq_revcomp
#> 89   seq_canonical
#> 90   seq_hash_2bit
#> 91 seq_encode_4bit
#> 92 seq_decode_4bit
#> 93  seq_gc_content
#> 94       seq_kmers
#>                                                                                              description
#> 88                       Compute the reverse complement of a DNA sequence using A, C, G, T, and N bases.
#> 89                        Return the lexicographically smaller of a sequence and its reverse complement.
#> 90                                         Encode a short DNA sequence as a 2-bit unsigned integer hash.
#> 91 Encode an IUPAC DNA sequence as a list of 4-bit base codes, preserving ambiguity symbols including N.
#> 92                              Decode a list of 4-bit IUPAC DNA base codes back into a sequence string.
#> 93                                    Compute GC fraction for a DNA sequence as a value between 0 and 1.
#> 94                              Expand a sequence into positional k-mers with optional canonicalization.
subset(rduckhts_functions(kind = "table"), select = c("name", "r_wrapper"))
#>                              name                        r_wrapper
#> 1               duckhts_simd_info               rduckhts_simd_info
#> 2        duckhts_simd_kernel_info        rduckhts_simd_kernel_info
#> 3        duckhts_simd_set_backend        rduckhts_simd_set_backend
#> 4                        read_bcf                     rduckhts_bcf
#> 5                     read_bcf_v2                                 
#> 6               read_bcf_appender                                 
#> 7                        read_bam                     rduckhts_bam
#> 8                     read_pileup                  rduckhts_pileup
#> 9                      read_fasta                   rduckhts_fasta
#> 10                       read_bed                     rduckhts_bed
#> 11                      fasta_nuc               rduckhts_fasta_nuc
#> 12      duckhts_cgranges_overlaps                                 
#> 13 duckhts_cgranges_overlaps_bulk                                 
#> 14                     read_fastq                   rduckhts_fastq
#> 15        detect_quality_encoding rduckhts_detect_quality_encoding
#> 16                       read_gff                     rduckhts_gff
#> 17                       read_gtf                     rduckhts_gtf
#> 18                     read_tabix                   rduckhts_tabix
#> 19                    fasta_index             rduckhts_fasta_index
#> 20                          bgzip                   rduckhts_bgzip
#> 21                        bgunzip                 rduckhts_bgunzip
#> 22                      bam_index               rduckhts_bam_index
#> 23                      bcf_index               rduckhts_bcf_index
#> 24                    tabix_index             rduckhts_tabix_index
#> 25                 bam_bin_counts          rduckhts_bam_bin_counts
#> 26       duckhts_bam_bed_coverage        rduckhts_bam_bed_coverage
#> 27               duckhts_mosdepth                rduckhts_mosdepth
#> 28      duckhts_samtools_idxstats       rduckhts_samtools_idxstats
#> 29                read_hts_header              rduckhts_hts_header
#> 30                 read_hts_index               rduckhts_hts_index
#> 31           read_hts_index_spans         rduckhts_hts_index_spans
#> 32                 bcftools_score                   rduckhts_score
#> 33                      seq_kmers                                 
```
