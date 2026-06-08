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
#> 85     seq_revcomp
#> 86   seq_canonical
#> 87   seq_hash_2bit
#> 88 seq_encode_4bit
#> 89 seq_decode_4bit
#> 90  seq_gc_content
#> 91       seq_kmers
#>                                                                                              description
#> 85                       Compute the reverse complement of a DNA sequence using A, C, G, T, and N bases.
#> 86                        Return the lexicographically smaller of a sequence and its reverse complement.
#> 87                                         Encode a short DNA sequence as a 2-bit unsigned integer hash.
#> 88 Encode an IUPAC DNA sequence as a list of 4-bit base codes, preserving ambiguity symbols including N.
#> 89                              Decode a list of 4-bit IUPAC DNA base codes back into a sequence string.
#> 90                                    Compute GC fraction for a DNA sequence as a value between 0 and 1.
#> 91                              Expand a sequence into positional k-mers with optional canonicalization.
subset(rduckhts_functions(kind = "table"), select = c("name", "r_wrapper"))
#>                              name                        r_wrapper
#> 1               duckhts_simd_info               rduckhts_simd_info
#> 2        duckhts_simd_kernel_info        rduckhts_simd_kernel_info
#> 3        duckhts_simd_set_backend        rduckhts_simd_set_backend
#> 4                        read_bcf                     rduckhts_bcf
#> 5                        read_bam                     rduckhts_bam
#> 6                     read_pileup                  rduckhts_pileup
#> 7                      read_fasta                   rduckhts_fasta
#> 8                        read_bed                     rduckhts_bed
#> 9                       fasta_nuc               rduckhts_fasta_nuc
#> 10      duckhts_cgranges_overlaps                                 
#> 11 duckhts_cgranges_overlaps_bulk                                 
#> 12                     read_fastq                   rduckhts_fastq
#> 13        detect_quality_encoding rduckhts_detect_quality_encoding
#> 14                       read_gff                     rduckhts_gff
#> 15                       read_gtf                     rduckhts_gtf
#> 16                     read_tabix                   rduckhts_tabix
#> 17                    fasta_index             rduckhts_fasta_index
#> 18                          bgzip                   rduckhts_bgzip
#> 19                        bgunzip                 rduckhts_bgunzip
#> 20                      bam_index               rduckhts_bam_index
#> 21                      bcf_index               rduckhts_bcf_index
#> 22                    tabix_index             rduckhts_tabix_index
#> 23                 bam_bin_counts          rduckhts_bam_bin_counts
#> 24       duckhts_bam_bed_coverage        rduckhts_bam_bed_coverage
#> 25               duckhts_mosdepth                rduckhts_mosdepth
#> 26      duckhts_samtools_idxstats       rduckhts_samtools_idxstats
#> 27                read_hts_header              rduckhts_hts_header
#> 28                 read_hts_index               rduckhts_hts_index
#> 29           read_hts_index_spans         rduckhts_hts_index_spans
#> 30                 bcftools_score                   rduckhts_score
#> 31                      seq_kmers                                 
```
