# Future DuckHTS SIMD Kernel Proposals

Date: 2026-05-29

This note captures code-grounded SIMD opportunities for later review. The intent is to treat SIMD as a general data-parallel kernel layer for packed sequence data, byte classification, quality-score math, integer reductions, delimiter parsing, and allele string comparison.

## Candidate hotspot map

### 1. BAM 4-bit packed sequence counts

Primary target:

- `src/bam_bin_counts.c::count_read_gc_bases`

Current scalar shape:

```c
static inline uint64_t count_read_gc_bases(const bam1_t *rec) {
    uint64_t gc = 0;
    const uint8_t *seq = bam_get_seq(rec);
    int len = rec->core.l_qseq;
    for (int i = 0; i < len; i++) {
        int code = bam_seqi(seq, i);
        if (code == 2 || code == 4) gc++;
    }
    return gc;
}
```

Why this matters:

- BAM stores two nt16 bases per byte, so this avoids text decoding.
- `C=2`, `G=4`, `A=1`, `T=8`, `N=15` are easy to classify in vector lanes.
- This can feed `bam_bin_counts(..., stats := 'gc')` and future QC summaries.

Potential internal result type:

```c
typedef struct {
    uint64_t a;
    uint64_t c;
    uint64_t g;
    uint64_t t;
    uint64_t n;
    uint64_t other;
    uint64_t gc;
    uint64_t called;
} duckhts_simd_bam_nt16_counts_t;
```

Potential kernel:

```c
void duckhts_simd_bam_nt16_counts(
    const uint8_t *packed_seq,
    int32_t n_bases,
    duckhts_simd_bam_nt16_counts_t *out
);
```

Potential manifest entry:

```c
DUCKHTS_SIMD_KERNEL(BAM_NT16_COUNTS,
                    bam_nt16_counts,
                    BAM_NT16_COUNTS,
                    "bam_nt16_counts")
```

Integration sketch:

```c
const duckhts_simd_dispatch_table_t *simd = duckhts_simd_dispatch_snapshot();

duckhts_simd_bam_nt16_counts_t s;
duckhts_simd_bam_nt16_counts_with_table(
    simd,
    bam_get_seq(rec),
    rec->core.l_qseq,
    &s
);

gc_bases = s.gc;
read_bases = s.called;
```

---

### 2. Unify `fasta_nuc(...)` with SIMD dispatch

Primary target:

- `src/interval_udf.c::fasta_nuc_scan`
- private nucleotide counter currently reached through `count_nucleotides(...)`

Opportunity:

- Replace duplicated local scalar/AVX2 nucleotide counting with the shared DuckHTS SIMD dispatch layer.
- Extend the existing `seq_base_counts` idea into a richer base histogram used by both `seq_gc_content(...)` and `fasta_nuc(...)`.

Potential result type:

```c
typedef struct {
    uint64_t a;
    uint64_t c;
    uint64_t g;
    uint64_t t;
    uint64_t n;
    uint64_t other;
    uint64_t gc;
    uint64_t called;
    int invalid;
} duckhts_simd_base_hist_t;
```

Potential manifest entry:

```c
DUCKHTS_SIMD_KERNEL(SEQ_BASE_HIST,
                    base_hist,
                    BASE_HIST,
                    "seq_base_hist")
```

Beneficiaries:

- `seq_gc_content(...)`
- `fasta_nuc(...)`
- future reference-window summaries

---

### 3. Sequence and quality one-pass stats

Primary targets:

- `src/bam_reader.c` sequence and quality extraction around `bam_get_seq(...)` and `bam_get_qual(...)`
- `src/quality_encoding.c::update_ascii_range`
- `src/bam_bed_coverage.c::accumulate_record_baseq_stats`
- future FASTQ/QC summaries

Current scalar quality range code:

```c
static void update_ascii_range(duckhts_quality_detect_result *out,
                               const char *text, size_t len) {
    size_t i;
    for (i = 0; i < len; i++) {
        unsigned char c = (unsigned char)text[i];
        if ((int)c < out->observed_ascii_min) out->observed_ascii_min = (int)c;
        if ((int)c > out->observed_ascii_max) out->observed_ascii_max = (int)c;
    }
}
```

Potential quality stats type:

```c
typedef struct {
    uint64_t sum;
    uint64_t count;
    uint32_t min;
    uint32_t max;
    uint64_t q20;
    uint64_t q30;
    uint64_t below_threshold;
} duckhts_simd_qual_stats_t;
```

Potential kernels:

```c
void duckhts_simd_qual_phred_stats(
    const uint8_t *qual,
    size_t len,
    uint8_t threshold,
    duckhts_simd_qual_stats_t *out
);

void duckhts_simd_qual_ascii_stats(
    const char *qual,
    size_t len,
    int ascii_offset,
    uint8_t threshold,
    duckhts_simd_qual_stats_t *out
);
```

Potential SQL-facing later feature:

```sql
seq_read_stats(SEQ, QUAL)
```

Potential output fields:

- `read_len`
- `num_a`, `num_c`, `num_g`, `num_t`, `num_n`, `num_other`
- `gc_bases`, `called_bases`, `gc_fraction`
- `qual_sum`, `qual_min`, `qual_max`
- `q20_bases`, `q30_bases`, `lowq_bases`, `mean_qual`

---

## `duckdb_mosdepth_simd` design sketch

Primary target:

- `src/mosdepth_table.c`

Important existing stages:

1. Build a difference array with `add_coverage_delta(...)`.
2. Prefix-sum the difference array into per-base depth.
3. Scan depth for cumulative depth, min/max, distribution, thresholds, quantization, and RLE output.

Keep scalar initially:

- htslib iteration
- CIGAR interpretation
- overlap handling
- pair correction

These paths are branchy and compatibility-sensitive. The first SIMD layer should operate after coverage has been materialized as a contiguous `int32_t *coverage` array.

### Proposed internal depth summary type

```c
typedef struct {
    uint64_t sum_depth;
    uint64_t length;
    uint32_t min_depth;
    uint32_t max_depth;
    uint64_t covered_bases;
} duckhts_simd_depth_summary_t;
```

### Proposed depth kernels

```c
void duckhts_simd_depth_prefix_sum_i32(
    int32_t *coverage_delta,
    int64_t len
);

void duckhts_simd_depth_summary_u32(
    const int32_t *coverage,
    int64_t len,
    uint32_t min_depth,
    duckhts_simd_depth_summary_t *out
);

void duckhts_simd_depth_threshold_counts_u32(
    const int32_t *coverage,
    int64_t len,
    const int64_t *thresholds,
    size_t n_thresholds,
    uint64_t *counts
);

void duckhts_simd_depth_rle_boundaries_i32(
    const int32_t *coverage,
    int64_t len,
    /* output change offsets or callback */
);
```

Potential manifest entries:

```c
DUCKHTS_SIMD_KERNEL(DEPTH_I32_PREFIX_SUM,
                    depth_i32_prefix_sum,
                    DEPTH_I32_PREFIX_SUM,
                    "depth_i32_prefix_sum")

DUCKHTS_SIMD_KERNEL(DEPTH_U32_SUMMARY,
                    depth_u32_summary,
                    DEPTH_U32_SUMMARY,
                    "depth_u32_summary")

DUCKHTS_SIMD_KERNEL(DEPTH_U32_THRESHOLDS,
                    depth_u32_thresholds,
                    DEPTH_U32_THRESHOLDS,
                    "depth_u32_thresholds")

DUCKHTS_SIMD_KERNEL(DEPTH_I32_RLE,
                    depth_i32_rle,
                    DEPTH_I32_RLE,
                    "depth_i32_rle")
```

### Mosdepth integration sketch

```c
const duckhts_simd_dispatch_table_t *simd =
    duckhts_simd_dispatch_snapshot();

/* 1. Convert delta array to depth array. */
duckhts_simd_depth_prefix_sum_i32_with_table(simd, coverage, chrom_len);

/* 2. Compute global summary. */
duckhts_simd_depth_summary_t summary;
duckhts_simd_depth_summary_u32_with_table(
    simd,
    coverage,
    chrom_len,
    min_depth,
    &summary
);

/* 3. Optional thresholds. */
duckhts_simd_depth_threshold_counts_u32_with_table(
    simd,
    coverage,
    chrom_len,
    thresholds->data,
    thresholds->count,
    threshold_counts
);

/* 4. Optional per-base/quantized RLE. */
duckhts_simd_depth_rle_i32_with_table(...);
```

### SIMD-friendly mosdepth operations

#### Summary math

High-confidence first target:

- sum depths
- min depth
- max depth
- count bases with depth >= threshold

This is pure integer reduction over contiguous arrays.

#### Threshold counts

For thresholds like `1,5,10,20,30`, compare each loaded depth vector against each threshold and accumulate per-threshold counts.

#### Quantized output

Current scalar shape calls `quantize_bucket(...)` per base and emits a run when the bucket changes. SIMD can classify blocks and quickly identify all-same-bucket spans, while keeping actual BED emission scalar.

#### Per-base RLE

Current scalar shape compares `coverage[i]` to the previous depth. SIMD can compare adjacent vectors to find change masks, then only scalar-handle the marked positions.

#### Prefix sum

A vectorized `int32_t` prefix sum is possible but more delicate than summary reductions. Do this after summary and threshold kernels are validated.

---

## `duckhts_bam_bed_coverage(...)` SIMD targets

Primary target:

- `src/bam_bed_coverage.c::flush_completed_tiles`

Current scalar reduction shape:

```c
for (j = 0; j < tile->len; j++) {
    uint32_t dpre = tile->depth_pre[j];
    uint32_t dpost = tile->depth_post[j];
    uint32_t dfwd = tile->depth_fwd_post ? tile->depth_fwd_post[j] : 0;
    uint32_t drev = tile->depth_rev_post ? tile->depth_rev_post[j] : 0;
    if (dpre >= (uint32_t)bind->min_depth) {
        region->covbases_pre++;
        region->summed_cov_pre += dpre;
    }
    ...
}
```

This maps directly to the proposed depth summary kernel:

```c
depth_u32_summary(depth_pre, len, min_depth, &summary_pre);
depth_u32_summary(depth_post, len, min_depth, &summary_post);
```

Outputs map to:

- `covbases_pre`
- `summed_cov_pre`
- `covbases_post`
- `summed_cov_post`
- forward/reverse post-filter summaries

Secondary targets:

- `src/bam_bed_coverage.c::accumulate_record_baseq_stats`
- `src/bam_bed_coverage.c::accumulate_record_on_tile`

For quality-filtered increments:

- if `min_baseq == 0`, use a fast range-increment path
- if quality exists and `min_baseq > 0`, vector-scan quality bytes
  - all-pass block: range increment
  - all-fail block: skip
  - mixed block: scalar fallback

---

## Parser SIMD targets

### VEP/CSQ parser

Primary target:

- `src/vep_parser.c`

Current delimiter-heavy operations include:

- `strchr(token, '|')`
- comma splitting
- field counting
- whitespace trimming

Potential kernels:

```c
void duckhts_simd_count_byte(
    const char *text,
    size_t len,
    char needle,
    uint64_t *count
);

void duckhts_simd_find_delims(
    const char *text,
    size_t len,
    const char *delims,
    size_t n_delims,
    uint32_t *positions,
    size_t max_positions,
    uint32_t *n_positions
);
```

Beneficiaries:

- CSQ/ANN/BCSQ splitting
- transcript count preallocation
- field validation
- reducing repeated `strchr(...)` scans

### Tabix, BED, GFF/GTF, summary-stat parsing

Likely additional targets:

- `src/tabix_reader.c`
- `src/interval_udf.c` BED parsing helpers
- score/summary-stat parsing code paths

A generic delimiter scanner can handle tabs, commas, semicolons, pipes, and newlines with one reusable byte-classification kernel.

---

## Allele normalization SIMD targets

Primary targets:

- `src/liftover_udf.c::scalar_is_left_aligned`
- `src/liftover_udf.c::scalar_trim_right`
- related scalar prefix/suffix comparisons
- future `bcftools norm`-style and DuckVEP allele normalization paths

Potential kernels:

```c
size_t duckhts_simd_common_prefix2(
    const char *a,
    size_t alen,
    const char *b,
    size_t blen
);

size_t duckhts_simd_common_suffix2(
    const char *a,
    size_t alen,
    const char *b,
    size_t blen
);

int duckhts_simd_dna_validate(
    const char *s,
    size_t len
);
```

Beneficiaries:

- left/right trimming
- left alignment
- allele validation
- VariantKey-style input validation
- future VEP allele normalization

---

## Suggested implementation order

### Priority 1: BAM 4-bit packed base counts

Deliverable:

- `bam_nt16_counts`

Use in:

- `bam_bin_counts(..., stats := 'gc')`

Why first:

- clear current hotspot
- avoids decoding to text
- straightforward scalar oracle
- easy real-BAM benchmark

### Priority 2: unified base histogram

Deliverable:

- `seq_base_hist`

Use in:

- `seq_gc_content(...)`
- `fasta_nuc(...)`

Why:

- removes duplicate local nucleotide-counting logic
- reuses existing dispatch infrastructure
- easy cross-backend validation

### Priority 3: quality stats kernels

Deliverables:

- `qual_phred_stats`
- `qual_ascii_stats` or `qual_ascii_range`

Use in:

- `detect_quality_encoding(...)`
- future `seq_read_stats(...)`
- BAM/FASTQ QC paths
- `duckhts_bam_bed_coverage(...)`

### Priority 4: mosdepth depth-array summaries

Deliverables:

- `depth_u32_summary`
- `depth_u32_threshold_counts`
- later `depth_i32_rle_changes`
- later `depth_i32_prefix_sum`

Use in:

- `duckhts_mosdepth(...)`
- `duckhts_bam_bed_coverage(...)`

### Priority 5: delimiter scanner

Deliverables:

- `count_byte_u8`
- `find_delims_u8`

Use in:

- `vep_parser.c`
- `tabix_reader.c`
- BED/GFF/GTF/summary-stat parsing

---

## Design rules

1. Keep the dispatch table per-kernel and immutable, following the existing SIMD dispatch model.
2. Always provide scalar fallback.
3. Keep optional backends compile- and runtime-gated.
4. Snapshot dispatch once per DuckDB vector/chunk or per worker phase, not once per row/base.
5. Keep R wrappers thin; backend inventories and kernel availability remain extension-owned.
6. Benchmark each kernel with both synthetic microbenchmarks and at least one real file workload.
7. Add function catalog, README, NEWS, SQL tests, and R tests only when a public SQL/R surface changes.

## Kernel families

Avoid one large generic SIMD API. Keep typed families:

1. ASCII/base kernels
   - `seq_base_counts`
   - `seq_base_hist`
   - delimiter scan
   - prefix/suffix compare

2. BAM packed 4-bit kernels
   - `bam_nt16_counts`

3. Quality byte kernels
   - min/max/sum/Q20/Q30/below-threshold

4. Depth integer kernels
   - prefix/sum/min/max/threshold/RLE
