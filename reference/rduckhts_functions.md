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
#> 43     seq_revcomp
#> 44   seq_canonical
#> 45   seq_hash_2bit
#> 46 seq_encode_4bit
#> 47 seq_decode_4bit
#> 48  seq_gc_content
#> 49       seq_kmers
#>                                                                                              description
#> 43                       Compute the reverse complement of a DNA sequence using A, C, G, T, and N bases.
#> 44                        Return the lexicographically smaller of a sequence and its reverse complement.
#> 45                                         Encode a short DNA sequence as a 2-bit unsigned integer hash.
#> 46 Encode an IUPAC DNA sequence as a list of 4-bit base codes, preserving ambiguity symbols including N.
#> 47                              Decode a list of 4-bit IUPAC DNA base codes back into a sequence string.
#> 48                                    Compute GC fraction for a DNA sequence as a value between 0 and 1.
#> 49                              Expand a sequence into positional k-mers with optional canonicalization.
subset(rduckhts_functions(kind = "table"), select = c("name", "r_wrapper"))
#>                              name                        r_wrapper
#> 1                        read_bcf                     rduckhts_bcf
#> 2                        read_bam                     rduckhts_bam
#> 3                      read_fasta                   rduckhts_fasta
#> 4                        read_bed                     rduckhts_bed
#> 5                       fasta_nuc               rduckhts_fasta_nuc
#> 6       duckhts_cgranges_overlaps                                 
#> 7  duckhts_cgranges_overlaps_bulk                                 
#> 8                      read_fastq                   rduckhts_fastq
#> 9         detect_quality_encoding rduckhts_detect_quality_encoding
#> 10                       read_gff                     rduckhts_gff
#> 11                       read_gtf                     rduckhts_gtf
#> 12                     read_tabix                   rduckhts_tabix
#> 13                    fasta_index             rduckhts_fasta_index
#> 14                          bgzip                   rduckhts_bgzip
#> 15                        bgunzip                 rduckhts_bgunzip
#> 16                      bam_index               rduckhts_bam_index
#> 17                      bcf_index               rduckhts_bcf_index
#> 18                    tabix_index             rduckhts_tabix_index
#> 19                 bam_bin_counts          rduckhts_bam_bin_counts
#> 20       duckhts_bam_bed_coverage        rduckhts_bam_bed_coverage
#> 21               duckhts_mosdepth                rduckhts_mosdepth
#> 22      duckhts_samtools_idxstats       rduckhts_samtools_idxstats
#> 23                read_hts_header              rduckhts_hts_header
#> 24                 read_hts_index               rduckhts_hts_index
#> 25           read_hts_index_spans         rduckhts_hts_index_spans
#> 26                 bcftools_score                   rduckhts_score
#> 27                      seq_kmers                                 
```
