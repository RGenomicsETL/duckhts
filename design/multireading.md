# Multi-File Reading in duckhts

## Part 1: Macro-Based Multi-File API (Works Today)

The `SET VARIABLE` + `query()` + `string_agg` + `UNION ALL BY NAME` pattern
provides comprehensive multi-file reading with no extension changes. All
examples below are verified against DuckDB v1.5.1 with the duckhts extension.

### Core Pattern

```sql
-- 1. Build the UNION ALL query string from a glob
SET VARIABLE q = (
    SELECT string_agg(
        'SELECT ''' || file || ''' AS filename, * FROM read_bam(''' || file || ''')',
        ' UNION ALL BY NAME '
    ) FROM glob('samples/*.bam') g(file)
);

-- 2. Execute
SELECT * FROM query(getvariable('q'));
```

DuckDB pushes filters and projections into each UNION ALL branch. Each
`read_bam` call gets its own contig-parallel scan internally.

### Per-File Parameters

The generated SQL string is arbitrary text -- each UNION ALL arm can have
completely different parameters. This is not a limitation of the pattern;
it is the pattern's main strength.

**Different region per file** (via a lookup table):

```sql
CREATE TABLE file_regions AS VALUES
    ('sample1.bam', 'chr1:1000-2000'),
    ('sample2.bam', 'chr17:7500000-7600000'),
    ('sample3.bam', 'chr13:32300000-32400000');

SET VARIABLE q = (
    SELECT string_agg(
        'SELECT ''' || column0 || ''' AS filename, * FROM read_bam('''
        || column0 || ''', region := ''' || column1 || ''')',
        ' UNION ALL BY NAME '
    ) FROM file_regions
);
SELECT * FROM query(getvariable('q'));
```

**Different mate_path per FASTQ pair:**

```sql
SET VARIABLE q = (
    SELECT string_agg(
        'SELECT ''' || file || ''' AS filename, * FROM read_fastq('''
        || file || ''', mate_path := '''
        || replace(file, '_R1.fq.gz', '_R2.fq.gz') || ''')',
        ' UNION ALL BY NAME '
    ) FROM glob('samples/*_R1.fq.gz') g(file)
);
SELECT * FROM query(getvariable('q'));
```

**Different tidy_format per BCF file:**

```sql
SET VARIABLE q = (
    SELECT string_agg(
        'SELECT ''' || file || ''' AS filename, * FROM read_bcf('''
        || file || ''', tidy_format := true)',
        ' UNION ALL BY NAME '
    ) FROM glob('cohort/*.bcf') g(file)
);
SELECT * FROM query(getvariable('q'));
```

Any named parameter for any reader can be embedded this way: `region`,
`index_path`, `mate_path`, `tidy_format`, `header_names`, `standard_tags`,
`auxiliary_tags`, `sequence_encoding`, `quality_representation`, `interleaved`,
`auto_detect`, `column_types`, `attributes_map`, `additional_csq_column_types`.

### Schema Reconciliation by Format

`UNION ALL BY NAME` matches columns by name, filling NULLs for columns absent
in a given arm. This handles heterogeneous schemas across files.

| Format | Schema stability | `UNION ALL` | `UNION ALL BY NAME` |
|---|---|---|---|
| BAM (default) | Fixed 13 cols | Works | Works |
| BAM (`standard_tags`) | Fixed ~70 cols (SAM spec superset) | Works | Works |
| BAM (`auxiliary_tags`) | Fixed -- all tags in one MAP column | Works | Works |
| BED | Fixed 13-col superset, NULLs for missing | Works | Works |
| BCF/VCF (wide, default) | Varies per file (samples, INFO, FORMAT) | **Fails** (column count mismatch) | Works (sparse NULLs) |
| BCF/VCF (`tidy_format`) | Mostly stable (samples become rows) | **Fails** (type mismatch on positional) | **Works** -- the correct multi-file VCF shape |
| `read_tabix` (no header) | `column0..N` positional | Fails if different counts | Works but semantically wrong (`column3` means different things) |
| `read_tabix` (`header=true`) | Real column names from `#header` | Fails | Works correctly by name matching |
| FASTA/FASTQ | Fixed schema | Works | Works |

**Key takeaway:** always use `UNION ALL BY NAME` for multi-file reads.
Plain `UNION ALL` only works when all files have identical schemas.

### Benchmarks

All measured on 20-core machine, extension loaded via `rduckhts_load()`.

| Files | Data | Time | Rows | Throughput |
|---|---|---|---|---|
| 1 | 316 MB exome BAM | ~3s | 3.25M | -- |
| 100 (distinct via glob) | 100 BAMs | 13.4s | -- | -- |
| 1000 (same BAM, simulating scale) | 316 MB x 1000 | 2m 27s | 3.25B | 22M rows/s |

Query plan overhead is linear in N (N scan nodes, N bind calls). The
practical ceiling is around 1000 files before plan compilation becomes
the bottleneck rather than I/O.

### What Does Not Work in SQL (Verified)

These are DuckDB-level constraints on C API table functions, not duckhts bugs:

| Approach | Error |
|---|---|
| `LATERAL read_bam(col)` | "does not support lateral join column parameters" |
| `read_bam(unnest([...]))` | "UNNEST not supported here" |
| `read_bam(['a.bam', 'b.bam'])` | "No function matches ... VARCHAR[]" |
| `read_bam((SELECT ...))` | "cannot contain subqueries" |
| `read_bam(glob(...))` | "glob is a table function" |

The `SET VARIABLE` + `query()` pattern works around all of these.

### Caveats

- **Two statements required.** `SET VARIABLE` then `SELECT ... query()`. Cannot
  be wrapped in a single TABLE macro (macros cannot contain `SET`).
- **Empty glob** sets the variable to NULL, causing a parse error in `query()`.
  Guard with `COALESCE` or check file count first.
- **Plan overhead** is O(N) in file count. Measured acceptable up to ~1000 files.

---

## Part 2: C-Level Multi-File Reader -- Research Questions

A native C-level multi-file reader (`read_bam('*.bam')` accepting globs or
lists directly) would eliminate the two-statement ceremony, reduce plan
overhead from O(N) to O(1), and enable DuckDB-managed 2D parallelism across
files and contigs. The macro-based API proves the feature is useful; the
question is whether the C implementation is tractable given htslib's resource
model.

### The Core Problem: htslib File Handle Proliferation

Each open HTS file requires:

- **`htsFile*`** -- wraps a BGZF file descriptor with its own decompression
  buffer (~64 KB default), optional thread pool allocation, and format-specific
  state.
- **`hts_idx_t*`** -- the loaded index. BAI indexes are typically 5-50 MB in
  memory. CSI indexes can be larger. These are either `mmap`'d or read into
  heap. Each thread in the current contig-claiming model loads its own copy.
- **`hts_itr_t*`** -- iterator state for region queries. Small, but one per
  active (file, contig) work unit.

For N files with T threads, the naive model opens `N * T` file handles and
`N * T` index copies. At 1000 files and 16 threads, that's 16,000 `htsFile*`
handles and 16,000 index loads.

### Specific Concerns

**File descriptor limits.** Default `ulimit -n` is often 1024. Each `htsFile*`
consumes one fd (or more for remote URLs). A pool of 16,000 handles blows
past the limit. Raising `ulimit` is a sysadmin concern, not an extension
concern.

**Index memory.** A 30x WGS BAM index is ~10 MB. 1000 files x 16 threads x
10 MB = 160 GB of index memory. Even with shared (not per-thread) index
loads: 1000 x 10 MB = 10 GB. The current per-thread index model
(`sam_index_load3` in each thread's `local_init`) does not share.

**htslib threading.** The current readers use `hts_set_threads(fp, 2)` per
file handle for BGZF decompression. With N open handles, that's `2N`
decompression threads competing with DuckDB's thread pool. htslib's thread
pool (`hts_tpool`) can be shared across handles, but the C API for that
requires careful lifecycle management.

**Handle caching / recycling.** Threads should not hold open handles to files
they are not currently scanning. A work-unit model where threads claim
(file, contig) pairs needs an LRU cache or explicit open/close cycle when
switching files. Opening and closing htslib handles is not free (involves
header parsing, BGZF state setup, and potentially index loading).

**Bind-time probing.** The current model opens the file at bind time to
detect schema. For N files, bind becomes O(N) in file opens. Options:
first-file-only probing (lazy, detects mismatches at scan time) or full
probe (expensive for remote files).

### Current Reader Parallelism Model (for Reference)

BAM and BCF use a 1D atomic counter over contigs within a single file:

```
global->current_contig = 0  (atomic)
Thread claims: next = __sync_fetch_and_add(&global->current_contig, 1)
If next >= n_contigs: done
Otherwise: sam_itr_queryi(idx, next, 0, HTS_POS_MAX)
```

Capped at `min(n_contigs, 16)` threads. Each thread opens its own `htsFile*`
and loads its own index copy.

Extending this to 2D (files x contigs) is conceptually straightforward --
the work-unit queue becomes `[(file_0, contig_0), (file_0, contig_1), ...,
(file_1, contig_0), ...]` -- but the resource management questions above
must be answered first.

### Per-Reader Index Summary

| Reader | Index Type | Parallel? | Index Required? |
|---|---|---|---|
| `read_bam` | .bai/.csi/.crai | Yes (contig-claiming) | No (degrades to seq), Yes for `region` |
| `read_bcf` | .tbi/.csi | Yes (contig-claiming) | No (degrades to seq), Yes for `region` |
| `read_tabix` | .tbi | No | No (degrades), Yes for `region` |
| `read_gff/gtf` | .tbi | No | No (degrades), Yes for `region` |
| `read_fasta` | .fai | No | Only for `region` |
| `read_fastq` | none | No | No |
| `read_bed` | .tbi | No | Only for `region` |

### Named Parameters per Reader

These are the parameters that a C-level multi-file reader would need to
handle. In the macro API, they are per-file (arbitrary SQL per arm). In a
C-level reader, they would either be uniform across all files or require
a new parameter model (e.g., parallel lists).

- **read_bam**: `region`, `index_path`, `reference`, `sequence_encoding`, `quality_representation`, `standard_tags`, `auxiliary_tags`
- **read_bcf**: `region`, `index_path`, `tidy_format`, `additional_csq_column_types`
- **read_tabix**: `region`, `index_path`, `header`, `header_names`, `auto_detect`, `column_types`
- **read_gff/gtf**: same as tabix + `attributes_map`
- **read_fasta**: `region`, `index_path`, `sequence_encoding`
- **read_fastq**: `mate_path`, `sequence_encoding`, `quality_representation`, `input_quality_encoding`, `interleaved`
- **read_bed**: `region`, `index_path`

### What a C-Level Implementation Would Buy

- **Single-statement ergonomics**: `SELECT * FROM read_bam('*.bam')` instead
  of `SET VARIABLE` + `query()`.
- **O(1) plan overhead**: one scan node instead of N.
- **Managed parallelism**: DuckDB's scheduler assigns threads to (file, contig)
  work units, with proper handle pooling and backpressure.
- **`filename` column**: trivially added as a virtual column per row.
- **LIST(VARCHAR) input**: `read_bam(['a.bam', 'b.bam'])`.

### Open Design Questions

1. **Index sharing.** Can we load each file's index once and share across
   threads? `hts_idx_t*` is read-only after load -- this should be safe but
   needs verification against htslib's internal state.

2. **Handle pooling.** What's the right pool size? LRU eviction of
   `htsFile*` handles, or pre-open a bounded set and block when exhausted?

3. **Thread budget.** How to divide DuckDB's thread pool between files?
   With 16 threads and 100 files, do we scan 16 files in parallel (1 thread
   each) or 4 files in parallel (4 threads each for contig parallelism)?

4. **Schema merging for VCF.** Wide-format VCF files from different samples
   have different column sets. A C-level reader needs to implement
   `UNION ALL BY NAME` semantics internally (scan all headers at bind,
   compute union schema, NULL-fill missing columns per file).

5. **Remote file globbing.** POSIX `glob(3)` only works for local files.
   S3 glob requires `ListObjectsV2`. Plain HTTP has no directory listing.
   DuckDB's `FileSystem::GlobFiles()` handles all of these but is not
   exposed through the C API.

6. **Parameter uniformity.** In the macro API, each UNION ALL arm can have
   different parameters. A C-level `read_bam('*.bam', region := 'chr1')`
   applies `region` uniformly to all files. Per-file parameters would need
   a new mechanism (parallel lists, lookup table, naming templates). This
   is a real expressiveness regression vs. the macro API.
