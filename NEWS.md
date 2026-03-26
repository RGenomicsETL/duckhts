# DuckHTS Extension News

## duckhts 0.1.3.9001 (2026-03-13)

- harden liftover error handling and diagnostics: `duckdb_liftover(...)` now raises explicit SQL errors for null/empty `chrom` and non-positive `pos`, liftover load failures surface the underlying chain/FASTA cause, macro registration no longer fails silently, and new SQL error-stress tests cover 1,000,000-row invalid-input paths
- add a `liftover(...)` table macro for score-style variant rows (`chrom`, `pos`, optional `ref`/`alt`) backed by a new `liftover_variant(...)` scalar kernel that uses UCSC chain files plus destination/source FASTA references to return lifted coordinates, lifted alleles, reverse-complement/swap indicators, and warning strings such as `IFFY` and `UNMAPPED`
- source the generated community extension descriptor `version` field from the repo-level `description.yml` instead of duplicating it in `functions.yaml`
- add `quality_representation := 'phred'` to `read_bam(...)` and `read_fastq(...)` so base qualities can be returned as `UTINYINT[]` raw Phred values instead of SAM/FASTQ text
- add `input_quality_encoding := 'phred33' | 'auto' | 'phred64' | 'solexa64'` to `read_fastq(...)`; default to modern `phred33`, with optional legacy decoding and canonical Phred output on read
- add `detect_quality_encoding(...)` to inspect FASTQ quality ASCII ranges and report compatible encodings plus a heuristic guessed encoding
- add `sequence_encoding := 'nt16'` parameter to `read_bam(...)`, `read_fasta(...)`, and `read_fastq(...)` to return sequence data as `UTINYINT[]` using htslib nt16 4-bit codes (`=ACMGRSVTWYHKDBN` → 0-15) instead of decoded `VARCHAR` strings; defaults to `'string'` for backward compatibility
- refactor `seq_encode_4bit(...)` / `seq_decode_4bit(...)` to use htslib's `seq_nt16_table[]` and `seq_nt16_str[]` directly instead of custom switch tables; `U` (RNA) now encodes as `T` (code 8) and code 0 (`=`) is accepted; **breaking**: unknown characters (e.g. `!`) now map to `N` (code 15) instead of returning NULL, matching htslib behavior — all encoding paths (UDF + reader `sequence_encoding := 'nt16'`) are now unified on the same shared code
- stop advertising unsupported `attributes_map := TRUE` on generic `read_tabix(...)`; parsed attribute maps remain available on `read_gff(...)` and `read_gtf(...)`
- centralize CSQ/ANN/BCSQ subfield typing for `read_bcf(...)` with builtin rules derived from `bcftools +split-vep`, conservative ANN defaults, and an `additional_csq_column_types := ...` override parameter using bcftools-style `PATTERN TYPE` entries
- add BGZF compression and decompression table functions: `bgzip(...)` and `bgunzip(...)`, both defaulting to preserving the source file unless `keep := FALSE` is requested
- add HTS index builders: `bam_index(...)`, `bcf_index(...)`, and `tabix_index(...)`
- add HTS metadata readers: `read_hts_header(...)`, `read_hts_index(...)`, `read_hts_index_spans(...)`, and `read_hts_index_raw(...)`
- add interval readers/helpers: `read_bed(...)` for BED3-BED12 input and `fasta_nuc(...)` for bedtools nuc-style FASTA interval composition over BED intervals or fixed-width bins
- add sequence helpers: `seq_encode_4bit(...)`, `seq_decode_4bit(...)`, `seq_gc_content(...)`, and `seq_kmers(...)`
- add SAM flag decoders: `sam_flag_bits(...)` for bulk struct output and `sam_flag_has(...)` for direct mask checks
- rename ambiguous SAM flag helpers to explicit/spec-aligned names: `is_paired(...)`, `is_proper_pair(...)`, `is_next_segment_unmapped(...)`, and `is_next_segment_reverse_complemented(...)`
- add `is_forward_aligned(...)` as a strand helper with SAM-spec semantics for mapped segments
- add `CIGAR utils` scalar helpers for soft clips, operator detection, and query/reference length calculations
- extend `read_bam(...)` with `standard_tags := TRUE` typed SAM tag columns and `auxiliary_tags := TRUE` for the remaining tags as a map
- improve `read_tabix(...)`, `read_gff(...)`, and `read_gtf(...)` with header-based column names, basic type inference, explicit column type overrides, and tabix metadata-aware parsing
- harden region-query behavior for files with incomplete contig metadata by returning empty results with a warning instead of failing
- add SQL coverage for BED reading, FASTA composition metrics, the new metadata readers, sequence/SAM-flag helpers, typed tabix parsing, and BGZF/index round-trips

## duckhts 0.1.1.9000 (2026-02-10)

- validate paired FASTQ mate files by exact QNAME match and equal record counts
- detect odd-length interleaved FASTQ input and raise an error
- return empty results (with warning) for region queries when contig headers are missing
- add non-conforming VCF fixture without ##contig for region query tests
- generic tabix reader now respects tabix header/meta configuration when inferring columns
- read_tabix supports header-based column names via header := true and header_names
- read_bam supports standard_tags (typed SAMtags columns) and auxiliary_tags (map of remaining tags)
- standard_tags + auxiliary_tags demo added to README R examples
- read_tabix supports auto_detect for basic numeric inference and column_types overrides
- added SQL/R demos and tests for tabix type inference and column overrides
