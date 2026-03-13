## Extension Function Catalog

This section is generated from `functions.yaml`.

### Readers

| Function | Kind | Returns | R helper | Description |
| --- | --- | --- | --- | --- |
| `read_bcf` | table | table | `rduckhts_bcf` | Read VCF and BCF variant data with typed INFO, FORMAT, and optional tidy sample output. |
| `read_bam` | table | table | `rduckhts_bam` | Read SAM, BAM, and CRAM alignments with optional typed SAMtags and auxiliary tag maps. |
| `read_fasta` | table | table | `rduckhts_fasta` | Read FASTA records or indexed FASTA regions as sequence rows. |
| `read_bed` | table | table | `rduckhts_bed` | Read BED3-BED12 interval files with canonical typed columns and optional tabix-backed region filtering. |
| `fasta_nuc` | table | table | `rduckhts_fasta_nuc` | Compute bedtools nuc-style nucleotide composition for supplied BED intervals or generated fixed-width bins over a FASTA reference. |
| `read_fastq` | table | table | `rduckhts_fastq` | Read single-end, paired-end, or interleaved FASTQ files. |
| `read_gff` | table | table | `rduckhts_gff` | Read GFF annotations with optional parsed attribute maps and indexed region filtering. |
| `read_gtf` | table | table | `rduckhts_gtf` | Read GTF annotations with optional parsed attribute maps and indexed region filtering. |
| `read_tabix` | table | table | `rduckhts_tabix` | Read generic tabix-indexed text data with optional header handling and type inference. |
| `fasta_index` | table | table | `rduckhts_fasta_index` | Build a FASTA index and return the index path used by the operation. |

### Compression

| Function | Kind | Returns | R helper | Description |
| --- | --- | --- | --- | --- |
| `bgzip` | table | table | `rduckhts_bgzip` | Compress a plain file to BGZF and return the created output path and byte counts. |
| `bgunzip` | table | table | `rduckhts_bgunzip` | Decompress a BGZF-compressed file and return the created output path and byte counts. |

### Indexing

| Function | Kind | Returns | R helper | Description |
| --- | --- | --- | --- | --- |
| `bam_index` | table | table | `rduckhts_bam_index` | Build a BAM or CRAM index and report the written index path and format. |
| `bcf_index` | table | table | `rduckhts_bcf_index` | Build a TBI or CSI index for a VCF or BCF file and report the written index path and format. |
| `tabix_index` | table | table | `rduckhts_tabix_index` | Build a tabix index for a BGZF-compressed text file using a preset or explicit coordinate columns. |

### Metadata

| Function | Kind | Returns | R helper | Description |
| --- | --- | --- | --- | --- |
| `read_hts_header` | table | table | `rduckhts_hts_header` | Inspect HTS headers in parsed, raw, or combined form across supported formats. |
| `read_hts_index` | table | table | `rduckhts_hts_index` | Inspect high-level HTS index metadata such as sequence names and mapped counts. |
| `read_hts_index_spans` | table_macro | table | `rduckhts_hts_index_spans` | Expand index metadata into span and chunk rows suitable for low-level index inspection. |
| `read_hts_index_raw` | table_macro | table | `rduckhts_hts_index_raw` | Return the raw on-disk HTS index blob together with basic identifying metadata. |

### Sequence UDFs

| Function | Kind | Returns | R helper | Description |
| --- | --- | --- | --- | --- |
| `seq_revcomp` | scalar | VARCHAR |  | Compute the reverse complement of a DNA sequence using A, C, G, T, and N bases. |
| `seq_canonical` | scalar | VARCHAR |  | Return the lexicographically smaller of a sequence and its reverse complement. |
| `seq_hash_2bit` | scalar | UBIGINT |  | Encode a short DNA sequence as a 2-bit unsigned integer hash. |
| `seq_encode_4bit` | scalar | UTINYINT[] |  | Encode an IUPAC DNA sequence as a list of 4-bit base codes, preserving ambiguity symbols including N. |
| `seq_decode_4bit` | scalar | VARCHAR |  | Decode a list of 4-bit IUPAC DNA base codes back into a sequence string. |
| `seq_gc_content` | scalar | DOUBLE |  | Compute GC fraction for a DNA sequence as a value between 0 and 1. |
| `seq_kmers` | table | table |  | Expand a sequence into positional k-mers with optional canonicalization. |

### SAM Flag UDFs

| Function | Kind | Returns | R helper | Description |
| --- | --- | --- | --- | --- |
| `sam_flag_bits` | scalar | STRUCT |  | Decode a SAM flag into a struct of boolean bit fields using explicit SAM-oriented names such as `is_paired`, `is_proper_pair`, `is_next_segment_unmapped`, and `is_supplementary`. |
| `sam_flag_has` | scalar | BOOLEAN |  | Test whether any bits from the provided SAM flag mask are set in a flag value. |
| `is_forward_aligned` | scalar | BOOLEAN |  | Test whether a mapped segment is aligned to the forward strand. Returns `NULL` for unmapped segments because SAM flag `0x10` does not define genomic strand when `0x4` is set. |
| `is_paired` | scalar | BOOLEAN |  | Test whether the SAM flag indicates that the template has multiple segments in sequencing (`0x1`). |
| `is_proper_pair` | scalar | BOOLEAN |  | Test whether the SAM flag indicates that each segment is properly aligned according to the aligner (`0x2`). |
| `is_unmapped` | scalar | BOOLEAN |  | Test whether the read itself is unmapped according to the SAM flag. |
| `is_next_segment_unmapped` | scalar | BOOLEAN |  | Test whether the next segment in the template is flagged as unmapped (`0x8`). |
| `is_reverse_complemented` | scalar | BOOLEAN |  | Test whether `SEQ` is stored reverse complemented (`0x10`); for mapped reads this corresponds to reverse-strand alignment. |
| `is_next_segment_reverse_complemented` | scalar | BOOLEAN |  | Test whether `SEQ` of the next segment in the template is stored reverse complemented (`0x20`). |
| `is_first_segment` | scalar | BOOLEAN |  | Test whether the read is marked as the first segment in the template. |
| `is_last_segment` | scalar | BOOLEAN |  | Test whether the read is marked as the last segment in the template. |
| `is_secondary` | scalar | BOOLEAN |  | Test whether the alignment is marked as secondary. |
| `is_qc_fail` | scalar | BOOLEAN |  | Test whether the read failed vendor or pipeline quality checks. |
| `is_duplicate` | scalar | BOOLEAN |  | Test whether the alignment is flagged as a duplicate. |
| `is_supplementary` | scalar | BOOLEAN |  | Test whether the alignment is marked as supplementary. |

### CIGAR Utils

| Function | Kind | Returns | R helper | Description |
| --- | --- | --- | --- | --- |
| `cigar_has_soft_clip` | scalar | BOOLEAN |  | Test whether a CIGAR string contains any soft-clipped segment (`S`). |
| `cigar_has_hard_clip` | scalar | BOOLEAN |  | Test whether a CIGAR string contains any hard-clipped segment (`H`). |
| `cigar_left_soft_clip` | scalar | BIGINT |  | Return the left-end soft-clipped length from a CIGAR string, or zero if the alignment does not start with `S`. |
| `cigar_right_soft_clip` | scalar | BIGINT |  | Return the right-end soft-clipped length from a CIGAR string, or zero if the alignment does not end with `S`. |
| `cigar_query_length` | scalar | BIGINT |  | Return the query-consuming length from a CIGAR string, counting `M`, `I`, `S`, `=`, and `X`. |
| `cigar_aligned_query_length` | scalar | BIGINT |  | Return the aligned query length from a CIGAR string, counting `M`, `=`, and `X` but excluding clips and insertions. |
| `cigar_reference_length` | scalar | BIGINT |  | Return the reference-consuming length from a CIGAR string, counting `M`, `D`, `N`, `=`, and `X`. |
| `cigar_has_op` | scalar | BOOLEAN |  | Test whether a CIGAR string contains at least one instance of the requested operator. |

