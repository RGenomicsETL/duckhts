# DuckHTS Extension News

## duckhts 0.1.3.9001 (2026-03-13)

- add BGZF compression and decompression table functions: `bgzip(...)` and `bgunzip(...)`, both defaulting to preserving the source file unless `keep := FALSE` is requested
- add HTS index builders: `bam_index(...)`, `bcf_index(...)`, and `tabix_index(...)`
- add HTS metadata readers: `read_hts_header(...)`, `read_hts_index(...)`, `read_hts_index_spans(...)`, and `read_hts_index_raw(...)`
- add sequence helpers: `seq_encode_4bit(...)`, `seq_decode_4bit(...)`, `seq_gc_content(...)`, and `seq_kmers(...)`
- extend `read_bam(...)` with `standard_tags := TRUE` typed SAM tag columns and `auxiliary_tags := TRUE` for the remaining tags as a map
- improve `read_tabix(...)`, `read_gff(...)`, and `read_gtf(...)` with header-based column names, basic type inference, explicit column type overrides, and tabix metadata-aware parsing
- harden region-query behavior for files with incomplete contig metadata by returning empty results with a warning instead of failing
- add SQL coverage for the new metadata readers, sequence/SAM-flag helpers, typed tabix parsing, and BGZF/index round-trips

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
