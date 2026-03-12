## Extension Function Catalog

This section is generated from `functions.yaml`.

### Readers

| Function | Kind | Returns | R helper | Description |
| --- | --- | --- | --- | --- |
| `read_bcf` | table | table | `rduckhts_bcf` | Read VCF and BCF variant data with typed INFO, FORMAT, and optional tidy sample output. |
| `read_bam` | table | table | `rduckhts_bam` | Read SAM, BAM, and CRAM alignments with optional typed SAMtags and auxiliary tag maps. |
| `read_fasta` | table | table | `rduckhts_fasta` | Read FASTA records or indexed FASTA regions as sequence rows. |
| `read_fastq` | table | table | `rduckhts_fastq` | Read single-end, paired-end, or interleaved FASTQ files. |
| `read_gff` | table | table | `rduckhts_gff` | Read GFF annotations with optional parsed attribute maps and indexed region filtering. |
| `read_gtf` | table | table | `rduckhts_gtf` | Read GTF annotations with optional parsed attribute maps and indexed region filtering. |
| `read_tabix` | table | table | `rduckhts_tabix` | Read generic tabix-indexed text data with optional header handling and type inference. |
| `fasta_index` | table | table | `rduckhts_fasta_index` | Build a FASTA index and return the index path used by the operation. |

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
| `is_segmented` | scalar | BOOLEAN |  | Test whether the SAM flag marks a read as part of a segmented template. |
| `is_properly_aligned` | scalar | BOOLEAN |  | Test whether the SAM flag indicates a properly aligned read pair. |
| `is_properly_segmented` | scalar | BOOLEAN |  | Alias for is_properly_aligned(flag) using segmented-read terminology. |
| `is_unmapped` | scalar | BOOLEAN |  | Test whether the read itself is unmapped according to the SAM flag. |
| `is_mate_unmapped` | scalar | BOOLEAN |  | Test whether the mate read is flagged as unmapped. |
| `is_reverse_complemented` | scalar | BOOLEAN |  | Test whether the read is aligned to the reverse strand. |
| `is_mate_reverse_complemented` | scalar | BOOLEAN |  | Test whether the mate read is aligned to the reverse strand. |
| `is_first_segment` | scalar | BOOLEAN |  | Test whether the read is marked as the first segment in the template. |
| `is_last_segment` | scalar | BOOLEAN |  | Test whether the read is marked as the last segment in the template. |
| `is_secondary` | scalar | BOOLEAN |  | Test whether the alignment is marked as secondary. |
| `is_qc_fail` | scalar | BOOLEAN |  | Test whether the read failed vendor or pipeline quality checks. |
| `is_duplicate` | scalar | BOOLEAN |  | Test whether the alignment is flagged as a duplicate. |
| `is_supplementary` | scalar | BOOLEAN |  | Test whether the alignment is marked as supplementary. |

