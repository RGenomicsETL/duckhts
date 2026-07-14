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
#> 95      seq_revcomp
#> 96    seq_canonical
#> 97    seq_hash_2bit
#> 98  seq_encode_4bit
#> 99  seq_decode_4bit
#> 100  seq_gc_content
#> 101       seq_kmers
#>                                                                                                                                                                                                                                                                                                                                                                                                             description
#> 95  Compute the reverse complement of a DNA sequence using A, C, G, T, and N bases. Overloaded: accepts either a VARCHAR text sequence (returns VARCHAR) or a UTINYINT[] of htslib nt16 codes as produced by read_bam(sequence_encoding := 'nt16') (returns UTINYINT[]); the nt16 overload is bit-identical to the text path after decoding, so BAM pipelines can reverse-complement without leaving the nt16 encoding.
#> 96                                           Return the lexicographically smaller of a sequence and its reverse complement. Overloaded: accepts either a VARCHAR text sequence (returns VARCHAR) or a UTINYINT[] of htslib nt16 codes as produced by read_bam(sequence_encoding := 'nt16') (returns UTINYINT[]); the nt16 overload compares by decoded base order and is bit-identical to the text path after decoding.
#> 97                                                                                                                                                                                   Encode a short DNA sequence as a 2-bit unsigned integer hash. Overloaded to also accept a UTINYINT[] of htslib nt16 codes (from read_bam(sequence_encoding := 'nt16')); non-ACGT codes yield NULL, bit-identical to the text path.
#> 98                                                                                                                                                                                                                                                                                                                Encode an IUPAC DNA sequence as a list of 4-bit base codes, preserving ambiguity symbols including N.
#> 99                                                                                                                                                                                                                                                                                                                                             Decode a list of 4-bit IUPAC DNA base codes back into a sequence string.
#> 100                                        Compute GC fraction for a DNA sequence as a value between 0 and 1. Overloaded: accepts either a VARCHAR text sequence or a UTINYINT[] of htslib nt16 codes as produced by read_bam(sequence_encoding := 'nt16'); the nt16 overload classifies codes directly and is bit-identical to the text path, so BAM pipelines can compute GC without decoding sequences back to text.
#> 101                                                                                                                                                                                                                                                                                                                                            Expand a sequence into positional k-mers with optional canonicalization.
subset(rduckhts_functions(kind = "table"), select = c("name", "r_wrapper"))
#>                              name                        r_wrapper
#> 1               duckhts_simd_info               rduckhts_simd_info
#> 2        duckhts_simd_kernel_info        rduckhts_simd_kernel_info
#> 3        duckhts_simd_set_backend        rduckhts_simd_set_backend
#> 4              duckvep_model_load                                 
#> 5                        read_bcf                     rduckhts_bcf
#> 6                     read_bcf_v2                                 
#> 7               read_bcf_appender                                 
#> 8                        read_bam                     rduckhts_bam
#> 9                     read_pileup                  rduckhts_pileup
#> 10                     read_fasta                   rduckhts_fasta
#> 11                       read_bed                     rduckhts_bed
#> 12                      fasta_nuc               rduckhts_fasta_nuc
#> 13      duckhts_cgranges_overlaps                                 
#> 14 duckhts_cgranges_overlaps_bulk                                 
#> 15                     read_fastq                   rduckhts_fastq
#> 16        detect_quality_encoding rduckhts_detect_quality_encoding
#> 17                       read_gff                     rduckhts_gff
#> 18                       read_gtf                     rduckhts_gtf
#> 19                     read_tabix                   rduckhts_tabix
#> 20                    fasta_index             rduckhts_fasta_index
#> 21                          bgzip                   rduckhts_bgzip
#> 22                        bgunzip                 rduckhts_bgunzip
#> 23                      bam_index               rduckhts_bam_index
#> 24                      bcf_index               rduckhts_bcf_index
#> 25                    tabix_index             rduckhts_tabix_index
#> 26                 bam_bin_counts          rduckhts_bam_bin_counts
#> 27       duckhts_bam_bed_coverage        rduckhts_bam_bed_coverage
#> 28               duckhts_mosdepth                rduckhts_mosdepth
#> 29      duckhts_samtools_idxstats       rduckhts_samtools_idxstats
#> 30                read_hts_header              rduckhts_hts_header
#> 31                 read_hts_index               rduckhts_hts_index
#> 32           read_hts_index_spans         rduckhts_hts_index_spans
#> 33                 bcftools_score                   rduckhts_score
#> 34                      seq_kmers                                 
```
