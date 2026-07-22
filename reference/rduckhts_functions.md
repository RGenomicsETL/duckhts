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
#>                name
#> 102     seq_revcomp
#> 103   seq_canonical
#> 104   seq_hash_2bit
#> 105 seq_encode_4bit
#> 106 seq_decode_4bit
#> 107  seq_gc_content
#> 108       seq_kmers
#>                                                                                                                                                                                                                                                                                                                                                                                                             description
#> 102 Compute the reverse complement of a DNA sequence using A, C, G, T, and N bases. Overloaded: accepts either a VARCHAR text sequence (returns VARCHAR) or a UTINYINT[] of htslib nt16 codes as produced by read_bam(sequence_encoding := 'nt16') (returns UTINYINT[]); the nt16 overload is bit-identical to the text path after decoding, so BAM pipelines can reverse-complement without leaving the nt16 encoding.
#> 103                                          Return the lexicographically smaller of a sequence and its reverse complement. Overloaded: accepts either a VARCHAR text sequence (returns VARCHAR) or a UTINYINT[] of htslib nt16 codes as produced by read_bam(sequence_encoding := 'nt16') (returns UTINYINT[]); the nt16 overload compares by decoded base order and is bit-identical to the text path after decoding.
#> 104                                                                                                                                                                                  Encode a short DNA sequence as a 2-bit unsigned integer hash. Overloaded to also accept a UTINYINT[] of htslib nt16 codes (from read_bam(sequence_encoding := 'nt16')); non-ACGT codes yield NULL, bit-identical to the text path.
#> 105                                                                                                                                                                                                                                                                                                               Encode an IUPAC DNA sequence as a list of 4-bit base codes, preserving ambiguity symbols including N.
#> 106                                                                                                                                                                                                                                                                                                                                            Decode a list of 4-bit IUPAC DNA base codes back into a sequence string.
#> 107                                        Compute GC fraction for a DNA sequence as a value between 0 and 1. Overloaded: accepts either a VARCHAR text sequence or a UTINYINT[] of htslib nt16 codes as produced by read_bam(sequence_encoding := 'nt16'); the nt16 overload classifies codes directly and is bit-identical to the text path, so BAM pipelines can compute GC without decoding sequences back to text.
#> 108                                                                                                                                                                                                                                                                                                                                            Expand a sequence into positional k-mers with optional canonicalization.
subset(rduckhts_functions(kind = "table"), select = c("name", "r_wrapper"))
#>                              name                        r_wrapper
#> 1               duckhts_simd_info               rduckhts_simd_info
#> 2        duckhts_simd_kernel_info        rduckhts_simd_kernel_info
#> 3        duckhts_simd_set_backend        rduckhts_simd_set_backend
#> 4              duckvep_model_load                                 
#> 5                duckvep_annotate                                 
#> 6                duckvep_so_terms                                 
#> 7                        read_bcf                     rduckhts_bcf
#> 8                     read_bcf_v2                                 
#> 9               read_bcf_appender                                 
#> 10                       read_bam                     rduckhts_bam
#> 11                    read_pileup                  rduckhts_pileup
#> 12                     read_fasta                   rduckhts_fasta
#> 13                       read_bed                     rduckhts_bed
#> 14                      fasta_nuc               rduckhts_fasta_nuc
#> 15      duckhts_cgranges_overlaps                                 
#> 16 duckhts_cgranges_overlaps_bulk                                 
#> 17                     read_fastq                   rduckhts_fastq
#> 18                    read_bigwig                  rduckhts_bigwig
#> 19        detect_quality_encoding rduckhts_detect_quality_encoding
#> 20                       read_gff                     rduckhts_gff
#> 21                       read_gtf                     rduckhts_gtf
#> 22                     read_tabix                   rduckhts_tabix
#> 23                    fasta_index             rduckhts_fasta_index
#> 24                          bgzip                   rduckhts_bgzip
#> 25                        bgunzip                 rduckhts_bgunzip
#> 26                      bam_index               rduckhts_bam_index
#> 27                      bcf_index               rduckhts_bcf_index
#> 28                    tabix_index             rduckhts_tabix_index
#> 29                 bam_bin_counts          rduckhts_bam_bin_counts
#> 30       duckhts_bam_bed_coverage        rduckhts_bam_bed_coverage
#> 31               duckhts_mosdepth                rduckhts_mosdepth
#> 32      duckhts_samtools_idxstats       rduckhts_samtools_idxstats
#> 33                read_hts_header              rduckhts_hts_header
#> 34                 read_hts_index               rduckhts_hts_index
#> 35           read_hts_index_spans         rduckhts_hts_index_spans
#> 36                 bcftools_score                   rduckhts_score
#> 37                      seq_kmers                                 
```
