# DuckHTS Extension News

## duckhts 0.1.1.9000 (2026-02-10)

- validate paired FASTQ mate files by exact QNAME match and equal record counts
- detect odd-length interleaved FASTQ input and raise an error
- return empty results (with warning) for region queries when contig headers are missing
- add non-conforming VCF fixture without ##contig for region query tests
- generic tabix reader now respects tabix header/meta configuration when inferring columns
- read_tabix supports header-based column names via header := true and header_names
- read_bam supports standard_tags (typed SAMtags columns) and auxiliary_tags (map of remaining tags)
- standard_tags + auxiliary_tags demo added to README R examples
