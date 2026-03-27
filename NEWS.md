# DuckHTS Extension News

## duckhts 0.1.3.9002 (2026-03-27)

- fix `bcftools_score(...)` CNT column naming to match upstream: when `q_score_thr` and `counts` are both active, count columns are now named `<prs>_CNT_p<thr>` (upstream pattern) instead of `<prs>_p<thr>_CNT` (**breaking** column name change for q_score_thr + counts queries)
- fix `bcftools_score(...)` threshold boundary precision to match upstream: parse `q_score_thr` values with `strtof`→double promotion before `-log10`, reproducing the exact float→double comparison asymmetry in upstream `bcftools +score`; markers at exact P-value boundaries may now be excluded where they were previously included
- add `scripts/score_conformance.sh` end-to-end conformance test that runs both DuckHTS and upstream `bcftools +score`, UNPIVOTs to long format, and FULL OUTER JOINs with numeric tolerance to catch any output divergence
- fix `bcftools_score(...)` wrong-result bugs: remove incorrect METAL `Zscore → P` and SSF `standard_error → P` column mappings that produced wrong PRS scores, add haploid GT support for chrX male samples, and fix memory leak when summary markers are skipped due to missing effect sizes
- fix `duckdb_liftover(...)` wrong-result bugs: skip reverse complement for symbolic alleles (`<DEL>`, `*`, multi-alt) that cannot be complemented, and detect insertions (1-bp ref, multi-bp alt) as indels so they route through the indel liftover path instead of the SNP path
- fix `duckdb_munge(...)` wrong-result bugs: emit rows with `filter := 'MissingContig'` instead of aborting the query when FASTA ref fetch fails for unknown chromosomes/scaffolds, and propagate NAN correctly when swapping AC with missing NS (matches upstream arithmetic)
- harden `bcftools_score(...)` upstream parity: port `bcf_hdr_name2id_flexible()` for flexible chromosome name resolution (chr prefix strip/prepend, 23→X, 24→Y, 26/MT→chrM aliases), port numerically stable `-log10(p)` parsing for very small p-values via mantissa/exponent splitting, and handle `NA`/`.` as missing values in summary stats fields (BETA, OR, P, LP)
- add comprehensive `bcftools_score(...)` test coverage: DS/HDS/AP/GP dosage modes with real FORMAT values, OR-to-beta conversion, PLINK2 preset with LOG10_P, custom columns_file mapping, allele mismatch (zero-match), missing genotype (`./. `) handling, NA in summary stats, chr-prefix flexible matching, auto-detection priority, and small p-value threshold precision
- tighten `bcftools_score(...)` TSV matching parity: enforce marker-ID column requirements for `use_variant_id := true`, keep SNP-ID fallback behavior explicit when CHR/BP are unavailable, and add SQL coverage for rsID-vs-CHR/BP matching paths
- slim munge output schema at the extension level by removing legacy `si`, `i2`, `cq`, and `ed` fields, and clarify orientation semantics by renaming `swapped` to `alleles_swapped` in `bcftools_munge_row(...)` / `duckdb_munge(...)`
- add `bcftools_munge_row(...)` and `duckdb_munge(...)` entries to `functions.yaml` so munge APIs are documented in the generated function catalog and community extension descriptor
- add `benchmark_munge.Rmd` and `make bench-munge` to benchmark DuckHTS munge against `bcftools +munge` with normalized output-group parity checks
- harden liftover error handling and diagnostics: `duckdb_liftover(...)` and `bcftools_liftover(...)` now raise explicit SQL errors for invalid required inputs, liftover load failures surface the underlying chain/FASTA cause, macro registration no longer fails silently, and new SQL error-stress tests cover 1,000,000-row invalid-input paths
- align liftover result semantics more closely with `bcftools +liftover` by splitting the old `warning` field into `reject_reason` for rejected rows and `note` for emitted rows with extra annotations, using upstream-style reject names such as `MissingContig`, `UnmappedAnchors`, `MismatchAnchors`, and `MissingFasta`
- add README liftover examples and expand impossible-liftover coverage so unmapped rows and invalid-input error paths are exercised in both SQL and package-level tests
- add a `liftover(...)` table macro for score-style variant rows (`chrom`, `pos`, optional `ref`/`alt`) backed by a new `liftover_variant(...)` scalar kernel that uses UCSC chain files plus destination/source FASTA references to return lifted coordinates, lifted alleles, reverse-complement/swap indicators, and liftover status annotations
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
