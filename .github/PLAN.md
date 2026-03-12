# DuckHTS Implementation Plan

## ✅ Phase 0 — Repository Setup (COMPLETED)
- [x] Add vendoring layout under `third_party/`.
- [x] Add initial vendoring scripts under `scripts/`.
- [x] Add license capture workflow for upstream sources.

## ✅ Phase 1 — htslib Integration (COMPLETED)
- [x] Add minimal build of htslib into extension build.
- [x] Provide a small C wrapper layer for reading records.
- [x] Add basic `read_bcf` table function.
- [x] Implement header parsing and structured metadata capture.

## ✅ Phase 2 — SAM/BAM/CRAM Readers (COMPLETED)
- [x] Add `read_bam` table function (covers SAM/BAM/CRAM).
- [x] Implement indexed region filtering.
- [x] Add tag handling strategy for optional fields.
- [x] Preserve read group and reference metadata for round-trip conversion.

## ✅ Phase 3 — Additional Readers (COMPLETED)
- [x] Add `read_fasta` and `read_fastq` table functions.
- [x] Add `read_gff` and `read_gtf` table functions.
- [x] Add `read_tabix` generic tabix reader.
- [x] Implement VEP annotation parser.

## ✅ Phase 4 — Conformance & CI (COMPLETED)
- [x] Add conformance datasets (small, indexed).
- [x] Add tests in `test/sql/` for core read paths.
- [x] Add CI job to run conformance tests offline.

## ✅ Phase 5 — R Package Harness (COMPLETED)
- [x] Add an R package subdirectory for testing.
- [x] Provide `bootstrap.R` to vendor and build the extension from R.
- [x] Add smoke tests that load and query the extension from R.

## ✅ Phase 6 — R Package & CRAN Release (COMPLETED)
- [x] Build a self-contained R package around the extension.
- [x] Adapt the package to use CMake plus `configure` / `configure.win`.
- [x] Add cleanup scripts and simplify package bootstrapping.
- [x] Remove `vcpkg` dependency for CRAN compatibility.
- [x] Add comprehensive R package tests.
- [x] Complete CRAN submission / resubmission work.
- [x] Publish `Rduckhts` on CRAN.
- [x] Advance package development beyond the initial CRAN baseline; current repository version is `0.1.3-0.0.2.9001`.

## 📋 Phase 7 — Post-CRAN Hardening, Testing, Release Maintenance & Distribution
- [ ] Add more comprehensive edge-case tests for current readers.
- [ ] Add performance benchmarks for hot paths and reader regressions.
- [ ] Update `README.Rmd` with current CRAN and source installation instructions.
- [ ] Add DuckDB `COPY` / export examples for downstream analytics workflows.
- [ ] Add vignettes for common HTS workflows.
- [ ] Test across Linux, macOS, and Windows MinGW/RTools.
- [ ] Run `R CMD check --as-cran` on release candidates across supported platforms.
- [ ] Fix release-blocking NOTEs, WARNINGs, and ERRORs as they arise.
- [ ] Keep CRAN policy compliance and toolchain compatibility current.
- [ ] Keep DuckDB community-extension metadata aligned with packaged releases.
- [ ] Prepare release materials for CRAN and community-extensions as needed.

## 🔨 Phase 8 — BGZF Compression & Index Building

DuckDB's `COPY TO 'file.gz'` produces standard gzip — **not BGZF**. Tabix and all htslib region queries require BGZF block compression. This phase adds BGZF compression/decompression and index building as composable table functions, completing the read→process→write→index→re-read pipeline entirely within SQL.

### 8.1 BGZF Compression / Decompression

| Function | htslib API | Equivalent CLI | Purpose |
|----------|-----------|----------------|---------|
| `bgzip(path, ...)` | `bgzf_open` + `bgzf_write` + `bgzf_mt` | `bgzip` | Compress a plain file to BGZF (`.gz`) |
| `bgunzip(path, ...)` | `bgzf_open` + `bgzf_read` | `bgzip -d` | Decompress a BGZF file back to plain text |

Both return `TABLE(success BOOLEAN, output_path VARCHAR, bytes_in BIGINT, bytes_out BIGINT)`.

**`bgzip` parameters:**
- `path` (positional) — input plain-text file
- `output_path := NULL` — output file; default appends `.gz` to input path
- `threads := 4` — BGZF compression threads via `bgzf_mt(fp, threads, 128)`
- `level := -1` — compression level (0–9, -1 = zlib default ~6)
- `keep := FALSE` — if TRUE, keep the original uncompressed file; if FALSE, remove it after compression (like `bgzip` CLI default)
- `overwrite := FALSE` — if TRUE, allow overwriting an existing output file; if FALSE, error when output exists

**`bgunzip` parameters:**
- `path` (positional) — input `.gz` BGZF file
- `output_path := NULL` — output file; default strips `.gz` suffix
- `threads := 4` — decompression threads via `bgzf_mt`
- `keep := FALSE` — if TRUE, keep the compressed file; if FALSE, remove it
- `overwrite := FALSE` — if TRUE, allow overwriting an existing output file

**Implementation core (bgzip):**
```c
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN
#include <htslib/bgzf.h>
#include <stdio.h>
#include <string.h>

/* Buffer size for file I/O — 64 KiB matches BGZF block ceiling */
#define BGZIP_BUF_SIZE (64 * 1024)

typedef struct {
    char *output_path;
    int64_t bytes_in;
    int64_t bytes_out;
    int emitted;
} bgzip_bind_t;

static void bgzip_bind(duckdb_bind_info info) {
    duckdb_value pv = duckdb_bind_get_parameter(info, 0);
    char *path = duckdb_get_varchar(pv);
    duckdb_destroy_value(&pv);

    char *output_path = NULL;
    int threads = 4;
    int level = -1;
    int keep = 0;
    int overwrite = 0;

    /* Parse named parameters ... */
    duckdb_value v;
    v = duckdb_bind_get_named_parameter(info, "output_path");
    if (v && !duckdb_is_null_value(v)) output_path = duckdb_get_varchar(v);
    if (v) duckdb_destroy_value(&v);

    v = duckdb_bind_get_named_parameter(info, "threads");
    if (v && !duckdb_is_null_value(v)) threads = (int)duckdb_get_int64(v);
    if (v) duckdb_destroy_value(&v);

    v = duckdb_bind_get_named_parameter(info, "level");
    if (v && !duckdb_is_null_value(v)) level = (int)duckdb_get_int64(v);
    if (v) duckdb_destroy_value(&v);

    v = duckdb_bind_get_named_parameter(info, "keep");
    if (v && !duckdb_is_null_value(v)) keep = duckdb_get_bool(v);
    if (v) duckdb_destroy_value(&v);

    v = duckdb_bind_get_named_parameter(info, "overwrite");
    if (v && !duckdb_is_null_value(v)) overwrite = duckdb_get_bool(v);
    if (v) duckdb_destroy_value(&v);

    /* Default output: input_path + ".gz" */
    if (!output_path) {
        output_path = duckdb_malloc(strlen(path) + 4);
        snprintf(output_path, strlen(path) + 4, "%s.gz", path);
    }

    /* Overwrite guard */
    {
        struct stat st;
        if (!overwrite && stat(output_path, &st) == 0) {
            char err[512];
            snprintf(err, sizeof(err),
                "bgzip: output '%s' already exists (use overwrite := TRUE to replace)",
                output_path);
            duckdb_bind_set_error(info, err);
            duckdb_free(path); duckdb_free(output_path);
            return;
        }
    }

    /* Open input (plain) and output (BGZF) */
    FILE *in = fopen(path, "rb");
    if (!in) {
        char err[512];
        snprintf(err, sizeof(err), "bgzip: cannot open input %s", path);
        duckdb_bind_set_error(info, err);
        duckdb_free(path); duckdb_free(output_path);
        return;
    }

    char mode[8];
    snprintf(mode, sizeof(mode), "w%d", level >= 0 ? level : 6);
    BGZF *out = bgzf_open(output_path, mode);
    if (!out) {
        fclose(in);
        duckdb_bind_set_error(info, "bgzip: cannot open output for BGZF writing");
        duckdb_free(path); duckdb_free(output_path);
        return;
    }

    /* Enable multi-threaded compression */
    if (threads > 1) bgzf_mt(out, threads, 128);

    /* Stream compress */
    char buf[BGZIP_BUF_SIZE];
    int64_t total_in = 0, total_out = 0;
    size_t n;
    while ((n = fread(buf, 1, BGZIP_BUF_SIZE, in)) > 0) {
        if (bgzf_write(out, buf, n) < 0) {
            duckdb_bind_set_error(info, "bgzip: write error");
            fclose(in); bgzf_close(out);
            duckdb_free(path); duckdb_free(output_path);
            return;
        }
        total_in += n;
    }
    fclose(in);
    if (bgzf_close(out) < 0) {
        duckdb_bind_set_error(info, "bgzip: error finalising BGZF output");
        duckdb_free(path); duckdb_free(output_path);
        return;
    }

    /* Remove original if keep=FALSE */
    if (!keep) remove(path);

    /* ... set bind data with output_path, total_in, total_out ... */
    bgzip_bind_t *bind = duckdb_malloc(sizeof(bgzip_bind_t));
    bind->output_path = output_path;
    bind->bytes_in = total_in;
    bind->bytes_out = 0;  /* could stat output file */
    bind->emitted = 0;

    /* Output schema */
    duckdb_logical_type bt = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type vt = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bi = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_bind_add_result_column(info, "success",     bt);
    duckdb_bind_add_result_column(info, "output_path", vt);
    duckdb_bind_add_result_column(info, "bytes_in",    bi);
    duckdb_bind_add_result_column(info, "bytes_out",   bi);
    duckdb_destroy_logical_type(&bt);
    duckdb_destroy_logical_type(&vt);
    duckdb_destroy_logical_type(&bi);
    duckdb_bind_set_bind_data(info, bind, destroy_bgzip_bind);
    duckdb_free(path);
}
```

`bgunzip` follows the same pattern with `bgzf_open(path, "r")` + `bgzf_read` + `fwrite` to plain output.

### 8.2 Index Building

| Function | htslib API | Equivalent CLI | Output format |
|----------|-----------|----------------|---------------|
| `bam_index(path, ...)` | `sam_index_build3` | `samtools index` | see below |
| `bcf_index(path, ...)` | `bcf_index_build3` | `bcftools index` | see below |
| `tabix_index(path, ...)` | `tbx_index_build3` | `tabix` | see below |
| `fasta_index(path, ...)` | `fai_build3` | `samtools faidx` | ✅ already exists |

All return `TABLE(success BOOLEAN, index_path VARCHAR, index_format VARCHAR)` — the `index_format` column explicitly reports what was written.

### 8.3 Output Format Matrix

| Function | `min_shift` | Input | Index file produced | `index_format` value |
|----------|------------|-------|---------------------|----------------------|
| `bam_index` | `0` (default) | `.bam` | `<path>.bai` | `BAI` |
| `bam_index` | `14` | `.bam` | `<path>.csi` | `CSI` |
| `bam_index` | any | `.cram` | `<path>.crai` | `CRAI` |
| `bcf_index` | `0` (default for `.vcf.gz`) | `.vcf.gz` | `<path>.tbi` | `TBI` |
| `bcf_index` | `14` (default for `.bcf`) | `.bcf` | `<path>.csi` | `CSI` |
| `tabix_index` | `0` (default) | any `.gz` | `<path>.tbi` | `TBI` |
| `tabix_index` | `>0` | any `.gz` | `<path>.csi` | `CSI` |
| `fasta_index` | n/a | `.fa`/`.fasta` | `<path>.fai` | `FAI` |
| `fasta_index` | n/a | `.fa.gz` (bgzipped) | `<path>.fai` + `<path>.gzi` | `FAI+GZI` |

**Key rules:**
- `min_shift = 0` → classic format (BAI for BAM, TBI for VCF/tabix)
- `min_shift > 0` (typically 14) → CSI format (supports large chromosomes > 512 Mbp)
- CRAM always produces `.crai` regardless of `min_shift`
- `bcf_index` auto-detects: `.bcf` defaults to `min_shift=14` (CSI), `.vcf.gz` defaults to `min_shift=0` (TBI)

### 8.4 Index Builder Parameters

**`bam_index(path, ...)`:**
- `index_path := NULL` — explicit output path; NULL → auto-append `.bai`/`.csi`/`.crai`
- `min_shift := 0` — `0`=BAI (default), `14`=CSI
- `threads := 4` — htslib decompression threads for index building

**`bcf_index(path, ...)`:**
- `index_path := NULL` — explicit output path; NULL → auto-append `.tbi`/`.csi`
- `min_shift := NULL` — NULL=auto-detect from extension (`.bcf`→14, `.vcf.gz`→0), or explicit
- `threads := 4` — htslib threads

**`tabix_index(path, ...)`:**
- `index_path := NULL` — explicit output path; NULL → auto-append `.tbi`/`.csi`
- `preset := 'vcf'` — `'vcf'`/`'bed'`/`'gff'`/`'sam'` → maps to `tbx_conf_vcf` etc.
- `min_shift := 0` — `0`=TBI, `>0`=CSI
- `threads := 4`
- `seq_col := NULL` — custom: 1-based sequence name column (overrides preset)
- `start_col := NULL` — custom: 1-based start position column
- `end_col := NULL` — custom: 1-based end position column
- `comment_char := NULL` — header/comment character (default from preset)
- `skip_lines := NULL` — fixed number of header lines to skip

### 8.5 Implementation Sketch (`src/hts_index_builder.c`)

```c
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <string.h>

/* Shared bind data for all index builders */
typedef struct {
    char *index_path;
    char *index_format;  /* "BAI", "CSI", "CRAI", "TBI", "FAI", "FAI+GZI" */
    int emitted;
} index_build_bind_t;

static void destroy_index_build_bind(void *data) {
    index_build_bind_t *b = (index_build_bind_t *)data;
    if (!b) return;
    if (b->index_path) duckdb_free(b->index_path);
    if (b->index_format) duckdb_free(b->index_format);
    duckdb_free(b);
}

/* Shared output schema: success, index_path, index_format */
static void add_index_result_columns(duckdb_bind_info info) {
    duckdb_logical_type bt = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type vt = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_bind_add_result_column(info, "success", bt);
    duckdb_bind_add_result_column(info, "index_path", vt);
    duckdb_bind_add_result_column(info, "index_format", vt);
    duckdb_destroy_logical_type(&bt);
    duckdb_destroy_logical_type(&vt);
}

/* Shared init */
static void index_build_init(duckdb_init_info info) {
    index_build_bind_t *bind = (index_build_bind_t *)duckdb_init_get_bind_data(info);
    bind->emitted = 0;
}

/* Shared scan: emit single row */
static void index_build_scan(duckdb_function_info info, duckdb_data_chunk output) {
    index_build_bind_t *bind = (index_build_bind_t *)duckdb_function_get_bind_data(info);
    if (bind->emitted) { duckdb_data_chunk_set_size(output, 0); return; }
    bool *success = (bool *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 0));
    success[0] = true;
    duckdb_vector_assign_string_element(duckdb_data_chunk_get_vector(output, 1), 0,
                                         bind->index_path ? bind->index_path : "");
    duckdb_vector_assign_string_element(duckdb_data_chunk_get_vector(output, 2), 0,
                                         bind->index_format ? bind->index_format : "");
    bind->emitted = 1;
    duckdb_data_chunk_set_size(output, 1);
}

/* Determine output format string from context */
static const char *bam_index_format_str(const char *path, int min_shift) {
    if (strstr(path, ".cram")) return "CRAI";
    return (min_shift > 0) ? "CSI" : "BAI";
}

static const char *bcf_index_format_str(int min_shift) {
    return (min_shift > 0) ? "CSI" : "TBI";
}

static const char *tabix_index_format_str(int min_shift) {
    return (min_shift > 0) ? "CSI" : "TBI";
}

/* ---- bam_index bind ---- */
static void bam_index_bind(duckdb_bind_info info) {
    duckdb_value pv = duckdb_bind_get_parameter(info, 0);
    char *path = duckdb_get_varchar(pv);
    duckdb_destroy_value(&pv);
    if (!path || strlen(path) == 0) {
        duckdb_bind_set_error(info, "bam_index requires a file path");
        if (path) duckdb_free(path);
        return;
    }

    char *index_path = NULL;
    int min_shift = 0;
    int threads = 4;

    duckdb_value v;
    v = duckdb_bind_get_named_parameter(info, "index_path");
    if (v && !duckdb_is_null_value(v)) index_path = duckdb_get_varchar(v);
    if (v) duckdb_destroy_value(&v);
    v = duckdb_bind_get_named_parameter(info, "min_shift");
    if (v && !duckdb_is_null_value(v)) min_shift = (int)duckdb_get_int64(v);
    if (v) duckdb_destroy_value(&v);
    v = duckdb_bind_get_named_parameter(info, "threads");
    if (v && !duckdb_is_null_value(v)) threads = (int)duckdb_get_int64(v);
    if (v) duckdb_destroy_value(&v);

    int ret = sam_index_build3(path, index_path, min_shift, threads);
    if (ret != 0) {
        char err[512];
        snprintf(err, sizeof(err), "bam_index: failed to index %s (error %d)", path, ret);
        duckdb_bind_set_error(info, err);
        duckdb_free(path); if (index_path) duckdb_free(index_path);
        return;
    }

    const char *fmt = bam_index_format_str(path, min_shift);
    add_index_result_columns(info);
    index_build_bind_t *bind = duckdb_malloc(sizeof(index_build_bind_t));
    bind->index_path = index_path ? index_path : strdup_duckdb("");
    bind->index_format = strdup_duckdb(fmt);
    bind->emitted = 0;
    duckdb_bind_set_bind_data(info, bind, destroy_index_build_bind);
    duckdb_free(path);
}

/* ---- bcf_index bind ---- */
static void bcf_index_bind(duckdb_bind_info info) {
    duckdb_value pv = duckdb_bind_get_parameter(info, 0);
    char *path = duckdb_get_varchar(pv);
    duckdb_destroy_value(&pv);
    if (!path || strlen(path) == 0) {
        duckdb_bind_set_error(info, "bcf_index requires a file path");
        if (path) duckdb_free(path);
        return;
    }

    char *index_path = NULL;
    int threads = 4;
    /* Auto-detect: .bcf → CSI, .vcf.gz → TBI */
    int min_shift = (strstr(path, ".bcf") != NULL) ? 14 : 0;
    int min_shift_explicit = 0;

    duckdb_value v;
    v = duckdb_bind_get_named_parameter(info, "index_path");
    if (v && !duckdb_is_null_value(v)) index_path = duckdb_get_varchar(v);
    if (v) duckdb_destroy_value(&v);
    v = duckdb_bind_get_named_parameter(info, "min_shift");
    if (v && !duckdb_is_null_value(v)) {
        min_shift = (int)duckdb_get_int64(v);
        min_shift_explicit = 1;
    }
    if (v) duckdb_destroy_value(&v);
    v = duckdb_bind_get_named_parameter(info, "threads");
    if (v && !duckdb_is_null_value(v)) threads = (int)duckdb_get_int64(v);
    if (v) duckdb_destroy_value(&v);

    int ret = bcf_index_build3(path, index_path, min_shift, threads);
    if (ret != 0) {
        char err[512];
        snprintf(err, sizeof(err), "bcf_index: failed to index %s (error %d)", path, ret);
        duckdb_bind_set_error(info, err);
        duckdb_free(path); if (index_path) duckdb_free(index_path);
        return;
    }

    const char *fmt = bcf_index_format_str(min_shift);
    add_index_result_columns(info);
    index_build_bind_t *bind = duckdb_malloc(sizeof(index_build_bind_t));
    bind->index_path = index_path ? index_path : strdup_duckdb("");
    bind->index_format = strdup_duckdb(fmt);
    bind->emitted = 0;
    duckdb_bind_set_bind_data(info, bind, destroy_index_build_bind);
    duckdb_free(path);
}

/* ---- tabix_index bind ---- */
static void tabix_index_bind(duckdb_bind_info info) {
    duckdb_value pv = duckdb_bind_get_parameter(info, 0);
    char *path = duckdb_get_varchar(pv);
    duckdb_destroy_value(&pv);
    if (!path || strlen(path) == 0) {
        duckdb_bind_set_error(info, "tabix_index requires a file path");
        if (path) duckdb_free(path);
        return;
    }

    char *index_path = NULL;
    int min_shift = 0;
    int threads = 4;
    tbx_conf_t conf = tbx_conf_vcf;  /* default */

    duckdb_value v;
    v = duckdb_bind_get_named_parameter(info, "index_path");
    if (v && !duckdb_is_null_value(v)) index_path = duckdb_get_varchar(v);
    if (v) duckdb_destroy_value(&v);

    v = duckdb_bind_get_named_parameter(info, "preset");
    if (v && !duckdb_is_null_value(v)) {
        char *p = duckdb_get_varchar(v);
        if      (strcmp(p, "vcf") == 0) conf = tbx_conf_vcf;
        else if (strcmp(p, "bed") == 0) conf = tbx_conf_bed;
        else if (strcmp(p, "gff") == 0) conf = tbx_conf_gff;
        else if (strcmp(p, "sam") == 0) conf = tbx_conf_sam;
        duckdb_free(p);
    }
    if (v) duckdb_destroy_value(&v);

    /* Custom column overrides (1-based, like tabix CLI) */
    v = duckdb_bind_get_named_parameter(info, "seq_col");
    if (v && !duckdb_is_null_value(v)) conf.sc = (int32_t)duckdb_get_int64(v);
    if (v) duckdb_destroy_value(&v);
    v = duckdb_bind_get_named_parameter(info, "start_col");
    if (v && !duckdb_is_null_value(v)) conf.bc = (int32_t)duckdb_get_int64(v);
    if (v) duckdb_destroy_value(&v);
    v = duckdb_bind_get_named_parameter(info, "end_col");
    if (v && !duckdb_is_null_value(v)) conf.ec = (int32_t)duckdb_get_int64(v);
    if (v) duckdb_destroy_value(&v);
    v = duckdb_bind_get_named_parameter(info, "comment_char");
    if (v && !duckdb_is_null_value(v)) {
        char *cc = duckdb_get_varchar(v);
        if (cc && strlen(cc) > 0) conf.meta_char = cc[0];
        if (cc) duckdb_free(cc);
    }
    if (v) duckdb_destroy_value(&v);
    v = duckdb_bind_get_named_parameter(info, "skip_lines");
    if (v && !duckdb_is_null_value(v)) conf.line_skip = (int32_t)duckdb_get_int64(v);
    if (v) duckdb_destroy_value(&v);
    v = duckdb_bind_get_named_parameter(info, "min_shift");
    if (v && !duckdb_is_null_value(v)) min_shift = (int)duckdb_get_int64(v);
    if (v) duckdb_destroy_value(&v);
    v = duckdb_bind_get_named_parameter(info, "threads");
    if (v && !duckdb_is_null_value(v)) threads = (int)duckdb_get_int64(v);
    if (v) duckdb_destroy_value(&v);

    int ret = tbx_index_build3(path, index_path, min_shift, threads, &conf);
    if (ret != 0) {
        char err[512];
        snprintf(err, sizeof(err), "tabix_index: failed to index %s (error %d)", path, ret);
        duckdb_bind_set_error(info, err);
        duckdb_free(path); if (index_path) duckdb_free(index_path);
        return;
    }

    const char *fmt = tabix_index_format_str(min_shift);
    add_index_result_columns(info);
    index_build_bind_t *bind = duckdb_malloc(sizeof(index_build_bind_t));
    bind->index_path = index_path ? index_path : strdup_duckdb("");
    bind->index_format = strdup_duckdb(fmt);
    bind->emitted = 0;
    duckdb_bind_set_bind_data(info, bind, destroy_index_build_bind);
    duckdb_free(path);
}

/* ---- Registration ---- */
void register_bam_index_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "bam_index");
    duckdb_logical_type vt = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type it = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_table_function_add_parameter(tf, vt);
    duckdb_table_function_add_named_parameter(tf, "index_path", vt);
    duckdb_table_function_add_named_parameter(tf, "min_shift", it);
    duckdb_table_function_add_named_parameter(tf, "threads", it);
    duckdb_table_function_set_bind(tf, bam_index_bind);
    duckdb_table_function_set_init(tf, index_build_init);
    duckdb_table_function_set_function(tf, index_build_scan);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
    duckdb_destroy_logical_type(&vt);
    duckdb_destroy_logical_type(&it);
}

void register_bcf_index_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "bcf_index");
    duckdb_logical_type vt = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type it = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_table_function_add_parameter(tf, vt);
    duckdb_table_function_add_named_parameter(tf, "index_path", vt);
    duckdb_table_function_add_named_parameter(tf, "min_shift", it);
    duckdb_table_function_add_named_parameter(tf, "threads", it);
    duckdb_table_function_set_bind(tf, bcf_index_bind);
    duckdb_table_function_set_init(tf, index_build_init);
    duckdb_table_function_set_function(tf, index_build_scan);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
    duckdb_destroy_logical_type(&vt);
    duckdb_destroy_logical_type(&it);
}

void register_tabix_index_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "tabix_index");
    duckdb_logical_type vt = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type it = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_table_function_add_parameter(tf, vt);
    duckdb_table_function_add_named_parameter(tf, "index_path", vt);
    duckdb_table_function_add_named_parameter(tf, "preset", vt);
    duckdb_table_function_add_named_parameter(tf, "min_shift", it);
    duckdb_table_function_add_named_parameter(tf, "threads", it);
    duckdb_table_function_add_named_parameter(tf, "seq_col", it);
    duckdb_table_function_add_named_parameter(tf, "start_col", it);
    duckdb_table_function_add_named_parameter(tf, "end_col", it);
    duckdb_table_function_add_named_parameter(tf, "comment_char", vt);
    duckdb_table_function_add_named_parameter(tf, "skip_lines", it);
    duckdb_table_function_set_bind(tf, tabix_index_bind);
    duckdb_table_function_set_init(tf, index_build_init);
    duckdb_table_function_set_function(tf, index_build_scan);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
    duckdb_destroy_logical_type(&vt);
    duckdb_destroy_logical_type(&it);
}
```

### 8.6 Registration in `duckhts.c`

```c
/* hts_index_builder.c */
extern void register_bam_index_function(duckdb_connection connection);
extern void register_bcf_index_function(duckdb_connection connection);
extern void register_tabix_index_function(duckdb_connection connection);

/* bgzip.c */
extern void register_bgzip_function(duckdb_connection connection);
extern void register_bgunzip_function(duckdb_connection connection);

/* In DUCKDB_EXTENSION_ENTRYPOINT: */
register_bam_index_function(connection);
register_bcf_index_function(connection);
register_tabix_index_function(connection);
register_bgzip_function(connection);
register_bgunzip_function(connection);
```

### 8.7 Expected SQL Usage

```sql
-- Index a BAM file (produces .bai by default)
SELECT * FROM bam_index('sample.bam');
-- ┌─────────┬────────────────┬──────────────┐
-- │ success │ index_path     │ index_format │
-- │ true    │                │ BAI          │
-- └─────────┴────────────────┴──────────────┘

-- Index a BAM with CSI format (for large genomes)
SELECT * FROM bam_index('sample.bam', min_shift := 14);
-- index_format = 'CSI'

-- Index a CRAM (always produces .crai)
SELECT * FROM bam_index('sample.cram');
-- index_format = 'CRAI'

-- Index a VCF.gz (auto-detects TBI)
SELECT * FROM bcf_index('variants.vcf.gz');
-- index_format = 'TBI'

-- Index a BCF (auto-detects CSI)
SELECT * FROM bcf_index('variants.bcf');
-- index_format = 'CSI'

-- Force CSI for a VCF.gz
SELECT * FROM bcf_index('variants.vcf.gz', min_shift := 14);
-- index_format = 'CSI'

-- Index a BED.gz with tabix
SELECT * FROM tabix_index('regions.bed.gz', preset := 'bed');
-- index_format = 'TBI'

-- Index a custom TSV.gz with explicit columns
SELECT * FROM tabix_index('custom.tsv.gz',
    seq_col := 1, start_col := 2, end_col := 3,
    comment_char := '#');
-- index_format = 'TBI'

-- BGZF compress a plain file
SELECT * FROM bgzip('regions.bed');
-- ┌─────────┬──────────────────┬──────────┬───────────┐
-- │ success │ output_path      │ bytes_in │ bytes_out │
-- │ true    │ regions.bed.gz   │ 1048576  │ 204800    │
-- └─────────┴──────────────────┴──────────┴───────────┘

-- BGZF decompress
SELECT * FROM bgunzip('regions.bed.gz');

-- Full pipeline: export → bgzip → tabix → re-read
COPY (
    SELECT chrom, start, "end", count_total, count_fwd, count_rev
    FROM bam_bin_counts('sample.bam', bin_width := 50000)
    ORDER BY chrom, start
) TO 'bins.bed' (FORMAT CSV, DELIMITER '\t', HEADER FALSE);

SELECT * FROM bgzip('bins.bed');
SELECT * FROM tabix_index('bins.bed.gz', preset := 'bed');

-- Now queryable with region filtering:
SELECT * FROM read_tabix('bins.bed.gz', region := 'chr21');
```

### 8.8 `functions.yaml` Entries

```yaml
    {
      "name": "bam_index",
      "kind": "table",
      "category": "Indexing",
      "signature": "bam_index(path, index_path := NULL, min_shift := 0, threads := 4)",
      "returns": "table",
      "r_wrapper": "rduckhts_bam_index",
      "description": "Build a BAM/CRAM index. Output format: min_shift=0 → .bai (BAI), min_shift=14 → .csi (CSI), CRAM → .crai (CRAI).",
      "examples": ["SELECT * FROM bam_index('sample.bam');"]
    },
    {
      "name": "bcf_index",
      "kind": "table",
      "category": "Indexing",
      "signature": "bcf_index(path, index_path := NULL, min_shift := NULL, threads := 4)",
      "returns": "table",
      "r_wrapper": "rduckhts_bcf_index",
      "description": "Build a VCF/BCF index. Auto-detects format: .bcf → .csi (CSI), .vcf.gz → .tbi (TBI). Override with min_shift.",
      "examples": ["SELECT * FROM bcf_index('variants.vcf.gz');"]
    },
    {
      "name": "tabix_index",
      "kind": "table",
      "category": "Indexing",
      "signature": "tabix_index(path, preset := 'vcf', index_path := NULL, min_shift := 0, threads := 4, seq_col := NULL, start_col := NULL, end_col := NULL, comment_char := NULL, skip_lines := NULL)",
      "returns": "table",
      "r_wrapper": "rduckhts_tabix_index",
      "description": "Build a tabix index (.tbi) for a BGZF-compressed file. Presets: vcf, bed, gff, sam. Custom column layout via seq_col/start_col/end_col.",
      "examples": ["SELECT * FROM tabix_index('regions.bed.gz', preset := 'bed');"]
    },
    {
      "name": "bgzip",
      "kind": "table",
      "category": "Compression",
      "signature": "bgzip(path, output_path := NULL, threads := 4, level := -1, keep := FALSE)",
      "returns": "table",
      "r_wrapper": "rduckhts_bgzip",
      "description": "BGZF-compress a file using htslib. Required before tabix indexing. Output: .gz appended to input path by default.",
      "examples": ["SELECT * FROM bgzip('regions.bed');"]
    },
    {
      "name": "bgunzip",
      "kind": "table",
      "category": "Compression",
      "signature": "bgunzip(path, output_path := NULL, threads := 4, keep := FALSE)",
      "returns": "table",
      "r_wrapper": "rduckhts_bgunzip",
      "description": "Decompress a BGZF-compressed file. Output: .gz suffix stripped from input path by default.",
      "examples": ["SELECT * FROM bgunzip('regions.bed.gz');"]
    },
```

### 8.9 Conformance Tests (`test/sql/index_build.test`)

```
require duckhts

# --- bam_index: build BAI, verify format reported ---
query BTT
SELECT success, index_format, index_format FROM bam_index('__WORKING_DIRECTORY__/test/data/range.bam');
----
true	BAI	BAI

# --- bcf_index: build TBI for VCF.gz ---
query BTT
SELECT success, index_format, index_format FROM bcf_index('__WORKING_DIRECTORY__/test/data/formatcols.vcf.gz');
----
true	TBI	TBI

# --- tabix_index with bed preset ---
query BT
SELECT success, index_format FROM tabix_index('__WORKING_DIRECTORY__/test/data/some_file.bed.gz',
    preset := 'bed');
----
true	TBI

# --- verify indexed BAM is now queryable by region ---
query I
SELECT COUNT(*) > 0 FROM read_bam('__WORKING_DIRECTORY__/test/data/range.bam',
    region := '11');
----
true

# --- bgzip round-trip ---
# (requires a small test fixture in plain text)
```

### 8.10 New Source Files

| File | Purpose |
|------|---------|
| `src/hts_index_builder.c` | `bam_index`, `bcf_index`, `tabix_index` — shared init/scan, format-specific bind |
| `src/bgzip.c` | `bgzip`, `bgunzip` — BGZF compression/decompression via htslib |

## ⏸️ Phase 9 — Write Path (FUTURE)
- [ ] Implement write path for VCF/BCF if it becomes a clear project need.
- [ ] Reassess write-path scope after more reader/analytics feedback from downstream users.

## 🧬 Phase 10 — Coverage & Interval Primitives for Interoperable CNV Workflows (PROPOSED)

### 10.1 Problem Statement
- [ ] Build elegant, high-performance, and conformant coverage/range primitives for multiple downstream workflows rather than a single tool clone.
- [ ] Support WisecondorX-style fixed-bin read counting.
- [ ] Support NIPTeR / NIPTmer-style strand-aware bin counting and future GC-aware normalization.
- [ ] Support general CNV and fetal-CNV style workflows with stable, reusable building blocks.

### 10.2 Design Goals
- [ ] Keep primitives orthogonal:
  - fixed-bin read counting is not per-base pileup,
  - pileup overlap handling is not CNV bin semantics,
  - interval preprocessing is not BAM scan policy.
- [ ] Keep semantics explicit:
  - every counting API must declare its counting unit, overlap behavior, and filtering contract.
- [ ] Prefer composable building blocks over tool-specific one-offs.
- [ ] Prefer a small canonical API over compatibility macros and one-off wrappers.
- [ ] Use BED-compatible outputs as the primary interoperability contract.
- [ ] Reuse a canonical interval-annotation schema across FASTA/BED and BAM bin outputs where possible.
- [ ] Preserve SQL-first workflows where the existing readers/UDFs already suffice.
- [ ] Make conformance a first-class deliverable before signatures are frozen.
- [ ] Require an R wrapper for every public function added in Phase 10.

### 10.3 Current Constraints and Facts
- [ ] Keep public APIs DuckDB 1.4+ compatible.
- [ ] `third_party/cgranges` is vendored but not yet wired into `CMakeLists.txt` or runtime code.
- [ ] There are no implemented BAM coverage/binning/pileup table functions yet in `src/`.
- [ ] Existing SQL tests cover readers/UDFs but not coverage/binning semantics.
- [ ] Existing readers and scalar UDFs already provide a SQL-first baseline for filtered read-level and GC-related summaries.
- [ ] Upstream semantic references are available locally under `.sync/` and should be treated as the primary source for behavior review:
  - `.sync/WisecondorX`
  - `.sync/NIPTeR`
  - `.sync/Rsamtools`

### 10.4 Reference Semantics to Preserve
- [ ] WisecondorX compatibility:
  - count read starts per fixed-width bin,
  - for paired reads, skip improper pairs,
  - require `MAPQ >= 1`,
  - by default suppress duplicates heuristically by repeated `(read.pos, mate.pos)` or repeated `read.pos`; `--normdup` disables that suppression,
  - keep export compatible with indexed tabular outputs instead of NPZ-only exchange.
- [ ] NIPTeR / NIPTmer compatibility:
  - NIPTeR currently bins BAM `pos` values into fixed 50 kb bins with optional strand separation,
  - current NIPTeR code relies on `scanBam()` and does not expose an explicit proper-pair / duplicate / MAPQ filter in the binning API,
  - preserve forward/reverse counts separately when requested,
  - keep room for later GC-aware normalization columns.
- [ ] WisecondorX masking behavior:
  - reference-time mask removes bins with zero or very low aggregate coverage across the reference set,
  - reference preparation additionally removes anomalous bins by PCA-distance outlier filtering,
  - prediction-time `--maskrepeats` derives a distance cutoff iteratively from reference distances,
  - prediction-time `--blacklist` applies an extra BED-driven mask to final per-bin outputs.
- [ ] Rsamtools pileup compatibility:
  - treat pileup as a per-base path,
  - mirror `PileupParam` controls that matter downstream: `min_mapq`, `min_base_quality`, `distinguish_strands`, `include_deletions`, `include_insertions`,
  - note that Rsamtools pileup automatically excludes unmapped, secondary, QC-fail, and duplicate reads.
- [ ] samtools `bedcov` compatibility:
  - add a separate BED-guided summed-depth path,
  - keep default flag behavior aligned with `bedcov` (`UNMAP`, `SECONDARY`, `QCFAIL`, `DUP` excluded by default),
  - support `-g` / `-G` style flag inclusion/exclusion and `-j` style deletion / ref-skip handling.
- [ ] samtools `coverage` compatibility:
  - add a region-level summary path returning `numreads`, `covbases`, `coverage`, `meandepth`, `meanbaseq`, `meanmapq`,
  - support read-length, mapq, baseq, required-flag, excluded-flag, min-depth, depth-cap, and region restriction options.
- [ ] samtools `depth` compatibility:
  - add a per-position depth path with `-a` / `-aa` style zero-depth controls,
  - support BED-driven region restriction, overlap suppression (`-s`), deletion inclusion (`-J`), flag filters, read-length, mapq, and baseq options.

### 10.4A Consistency Notes to Verify Before Freezing APIs
- [ ] **Threading pattern**: the plan currently assumes contig-level partitioning plus per-handle `hts_set_threads(fp, 2)` broadly. Re-check this against htslib 1.23 and actual workloads before treating it as mandatory for every function.
- [ ] **Multi-region BAM iteration**: confirm the exact htslib 1.23 API and dedup semantics for multi-region iteration before standardizing on `sam_itr_regarray()` language across docs and code.
- [ ] **Flag defaults**: keep default excluded-flag sets explicit per compatibility target. Do not casually reuse `0x704` where supplementary (`0x800`) behavior also matters.
- [ ] **Pair / duplicate compatibility**: `require_flags` / `exclude_flags` are the preferred generic controls, but WisecondorX duplicate suppression is heuristic rather than pure flag-based. Mark any claimed compatibility as partial until tested against `.sync/WisecondorX`.
- [ ] **Coordinate conventions**: BED/bin outputs can remain 0-based half-open, but base-level outputs (`bam_depth`, `bam_pileup`) must have their coordinate contract checked against `samtools` / Rsamtools expectations before freezing docs or wrappers.
- [ ] **`bam_bedcov` naming**: verify whether the output metric should represent summed depth, covered bases above threshold, read count, or multiple columns; avoid implying `bedcov` semantics with the wrong column names.
- [ ] **`bam_depth` scope**: current test notes mention possible multi-file behavior while the proposed signature is single-path. Decide explicitly whether multi-file depth is in scope for v1.
- [ ] **Interval UDF shape**: examples like `interval_merge(SELECT ...)` are conceptual only. Confirm the actual DuckDB table-function interface before documenting them as callable SQL forms.

### 10.5 Revised Architecture

#### 10.5.1 Layer 0 — SQL-First Baseline
- [ ] Use `read_bam(...)`, existing SAM flag UDFs, and sequence helpers as the initial correctness baseline.
- [ ] Treat SQL-based aggregation as the parity oracle for later native kernels.

**WisecondorX-like bin counting in pure SQL (available today):**
```sql
-- read_start binning with proper-pair + MAPQ gate, 5 kb bins
SELECT RNAME AS chrom,
       (POS / 5000) * 5000 AS start,
       (POS / 5000) * 5000 + 5000 AS "end",
       COUNT(*) AS count_total,
       SUM(CASE WHEN NOT is_reverse_complemented(FLAG) THEN 1 ELSE 0 END) AS count_fwd,
       SUM(CASE WHEN is_reverse_complemented(FLAG) THEN 1 ELSE 0 END) AS count_rev
FROM read_bam('sample.bam')
WHERE is_properly_aligned(FLAG)
  AND NOT is_duplicate(FLAG)
  AND NOT is_secondary(FLAG)
  AND NOT is_supplementary(FLAG)
  AND MAPQ >= 1
GROUP BY RNAME, (POS / 5000)
ORDER BY RNAME, start;
```

**Per-read GC with filters (available today):**
```sql
SELECT RNAME, POS, seq_gc_content(SEQ) AS gc,
       LENGTH(SEQ) AS read_len
FROM read_bam('sample.bam')
WHERE is_properly_aligned(FLAG)
  AND NOT is_unmapped(FLAG)
  AND MAPQ >= 20;
```

#### 10.5.2 Layer 1 — Interval / Bin Primitives
- [ ] Add interval preprocessing / algebra helpers:
  - `read_bed(...)`
  - `fasta_nuc(...)`
  - `interval_merge(...)`
  - `interval_overlap(...)`
  - `interval_nearest(...)`
- [ ] Keep fixed-bin generation internal to `bam_bin_counts`; add a public bin generator only if a concrete downstream need remains after the counting/export APIs are stable.
- [ ] Use arithmetic binning for regular bins.
- [ ] Use `cgranges` for irregular interval overlap / contain / nearest operations and interval preprocessing.
- [ ] Do not hide interval merge/combine policy inside BAM counting functions.

**Build integration — add cgranges to CMakeLists.txt:**
```cmake
# In CMakeLists.txt, add to EXTENSION_SOURCES:
set(EXTENSION_SOURCES
        src/duckhts.c
        src/bcf_reader.c
        src/bam_reader.c
        src/seq_reader.c
        src/kmer_udf.c
        src/tabix_reader.c
        src/vep_parser.c
        src/hts_meta_reader.c
        src/interval_udf.c          # NEW: interval algebra UDFs
        src/bam_bin_counts.c         # NEW: bin-count kernel
        third_party/cgranges/cgranges.c
)
```

**cgranges internal adapter pattern (`src/interval_udf.c` sketch):**
```c
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include "../third_party/cgranges/cgranges.h"
#include <string.h>
#include <stdlib.h>

/* read_bed: parse BED / BED.GZ into canonical columns, preserving extra fields */
/* interval_merge / interval_overlap / interval_nearest: build a cgranges index
 * over canonical BED-like (chrom, start, end) inputs and emit merged/joined rows */
```

**Expected SQL usage for interval primitives:**
```sql
-- Read a BED file as canonical intervals
SELECT * FROM read_bed('targets.bed.gz', region := 'chr21');

-- Annotate BED intervals with nucleotide content from a FASTA
SELECT *
FROM fasta_nuc('reference.fa', 'targets.bed.gz');

-- Merge overlapping BED intervals before using as targets
SELECT * FROM interval_merge(
    SELECT chrom, start, "end" FROM read_tabix('targets.bed.gz')
);

-- Overlap-join reads against custom target regions
SELECT t.chrom, t.start, t.end, COUNT(*) AS n_reads
FROM read_bam('sample.bam') r
JOIN interval_overlap(
    (SELECT chrom, start, "end" FROM read_tabix('targets.bed.gz')),
    (SELECT RNAME AS chrom, POS AS start, POS + LENGTH(SEQ) AS "end" FROM read_bam('sample.bam'))
) o ON ...
GROUP BY t.chrom, t.start, t.end;
```

#### 10.5.3 Layer 2 — Native Counting Kernels
- [ ] Implement a dedicated CNV-oriented fixed-bin kernel:
  - `bam_bin_counts(path, bin_width, [...])`
- [ ] Implement a separate BED-guided summed-depth kernel:
  - `bam_bedcov(path, bed_path, [...])`
- [ ] Implement a region-summary kernel mirroring `samtools coverage`:
  - `bam_coverage(path, [...])`
- [ ] Implement separate base-depth / pileup kernels:
  - `bam_depth(path, [...])`
  - `bam_pileup(path, [...])`
- [ ] Do not add compatibility wrappers until these kernels and their R wrappers are stable.

**Kernel sketch summary:**
- `bam_bin_counts` follows the contig-parallel pattern already proven in `bam_reader.c`.
- The hot path computes `bin_id = pos / bin_width` and updates total / forward / reverse counts in one pass.
- Filtering is expressed with `require_flags`, `exclude_flags`, and `min_mapq`.
- Sparse output is the default; `emit_zero_bins := TRUE` is an explicit export-oriented mode.
- Region-restricted scans use the existing htslib index iterator path; full indexed scans use contig partitioning.

**Expected SQL usage for counting kernels:**
```sql
-- WisecondorX-like fixed-bin read-start counts
SELECT *
FROM bam_bin_counts(
    'sample.bam',
    bin_width := 5000,
    min_mapq := 1,
    require_flags := 0x2,   -- proper pair
    exclude_flags := 0x704  -- unmapped/secondary/qcfail/duplicate; supplementary handled separately if needed
);

-- NIPTeR-like strand split
SELECT chrom, start, "end", count_total, count_fwd, count_rev
FROM bam_bin_counts(
    'sample.bam',
    bin_width := 50000,
    strand_mode := 'split'
);

-- samtools bedcov-like summed coverage over targets
SELECT *
FROM bam_bedcov(
    'sample.bam',
    'targets.bed.gz',
    exclude_flags := 0x704
);

-- samtools coverage-like region summary
SELECT *
FROM bam_coverage(
    'sample.bam',
    min_mapq := 20,
    min_baseq := 10,
    region := 'chr21'
);
```

#### 10.5.4 Layer 3 — Export and Interoperability
- [ ] Keep table outputs canonical and BED-friendly.
- [ ] Add a writer/procedure layer for sorted `.bed.gz` / `.tsv.gz` plus `.tbi` / `.csi`.
- [ ] Ensure exported files round-trip through DuckHTS readers and remain consumable by external tools.

**Export round-trip pattern (COPY → bgzip → tabix_index → read_tabix):**
```sql
-- Step 1: Export bin counts to sorted plain BED (NOT .gz — DuckDB gzip ≠ BGZF)
COPY (
    SELECT chrom, start, "end", count_total, count_fwd, count_rev
    FROM bam_bin_counts('sample.bam', bin_width := 50000,
         strand_mode := 'split')
    ORDER BY chrom, start
) TO 'sample_bins.bed' (FORMAT CSV, DELIMITER '\t', HEADER FALSE);

-- Step 2: BGZF-compress (produces sample_bins.bed.gz, removes plain file)
SELECT * FROM bgzip('sample_bins.bed');

-- Step 3: Build tabix index (produces sample_bins.bed.gz.tbi)
SELECT * FROM tabix_index('sample_bins.bed.gz', preset := 'bed');

-- Step 4: Read back with region filtering — full round-trip
SELECT * FROM read_tabix('sample_bins.bed.gz',
    header_names := ['chrom','start','end','count_total','count_fwd','count_rev'],
    column_types := {'start': 'BIGINT', 'end': 'BIGINT',
                     'count_total': 'BIGINT', 'count_fwd': 'BIGINT',
                     'count_rev': 'BIGINT'},
    region := 'chr21');
```

### 10.6 Semantic Contract to Freeze Before Implementation

#### 10.6.1 Count Model
- [ ] Keep separate functions for distinct counting units:
  - `bam_bin_counts` => read-start counts in fixed bins
  - `bam_bedcov` => summed per-base depth across supplied BED intervals
  - `bam_coverage` => per-region summary statistics
  - `bam_depth` / `bam_pileup` => base-level outputs
- [ ] Avoid a wide `count_model` enum unless a new mode cannot be expressed by a separate function with a clearer contract.

#### 10.6.2 Filtering Contract
- [ ] Use `require_flags` and `exclude_flags` as the primary filtering controls for counting APIs.
- [ ] Keep `min_mapq` and `min_baseq` explicit where applicable.
- [ ] Prefer symbolic SAM flag names in docs and R wrappers; keep raw integer masks available as the exact low-level form for SQL power users and compatibility work.
- [ ] Add a higher-level duplicate or pair parameter only if a required upstream behavior cannot be represented cleanly with flags plus overlap settings.

#### 10.6.3 Strand Policy
- [ ] Support `strand_mode := 'combined' | 'split'`.
- [ ] Make `count_total`, `count_fwd`, and `count_rev` part of the stable bin-count contract.

#### 10.6.4 Overlap Policy
- [ ] Keep overlap handling explicit and local to the functions where it matters.
- [ ] `bam_bedcov` must document whether deletions and ref-skips contribute, matching `samtools bedcov -j`.
- [ ] `bam_depth` must document `-s` overlap suppression and `-J` deletion inclusion explicitly.
- [ ] `bam_pileup` remains separate for richer base-level Rsamtools-like output.

#### 10.6.5 Coordinate and Storage Contract
- [ ] Keep exported BED-like files in canonical `chrom`, `start`, `end` 0-based half-open form.
- [ ] Allow SQL outputs to expose convenience columns if useful, without changing the canonical on-disk contract.
- [ ] Require deterministic output sort order by `chrom`, `start`, `end`.
- [ ] Keep interval annotation columns reusable across functions:
  - `seq_len`
  - `pct_at`
  - `pct_gc`
  - `num_a`
  - `num_c`
  - `num_g`
  - `num_t`
  - `num_n`
  - `num_other`

#### 10.6.6 Region Semantics
- [ ] Keep region preprocessing separate from BAM scan policy.
- [ ] Preprocess irregular target sets with interval primitives rather than burying merge policy in BAM APIs.
- [ ] Document exact behavior for region strings, region files, merged targets, and overlapping intervals.

### 10.7 Proposed API Backlog (Prioritized)

#### 10.7.1 Interval / Bin Layer
- [ ] `read_bed(path, region := NULL, index_path := NULL, header := NULL, column_names := NULL, column_types := NULL)`
- [ ] `fasta_nuc(fasta_path, bed_path, stranded := FALSE, include_seq := FALSE, pattern := NULL, ignore_case := FALSE, full_header := FALSE, reference_index := NULL)`
- [ ] Stable output columns after the input BED columns:
  - `pct_at`
  - `pct_gc`
  - `num_a`
  - `num_c`
  - `num_g`
  - `num_t`
  - `num_n`
  - `num_other`
  - `seq_len`
- [ ] Optional output columns:
  - `seq` when `include_seq := TRUE`
  - `pattern_count` when `pattern` is provided
- [ ] `interval_merge(intervals, [...])`
- [ ] `interval_overlap(left, right, [mode := 'any'|'contain'|'within'])`
- [ ] `interval_nearest(left, right, [k := 1, max_distance := NULL])`

#### 10.7.2 Native Bin-Count Kernel
- [ ] `bam_bin_counts(path, bin_width := 5000, strand_mode := 'combined', min_mapq := 0, require_flags := NULL, exclude_flags := NULL, region := NULL, index_path := NULL, reference := NULL, emit_zero_bins := FALSE)`
- [ ] Stable output columns:
  - `chrom`
  - `start`
  - `end`
  - `bin_id`
  - `count_total`
  - `count_fwd`
  - `count_rev`
- [ ] Optional annotation columns, computed only when requested with a reference FASTA:
  - `pct_at`
  - `pct_gc`
  - `num_a`
  - `num_c`
  - `num_g`
  - `num_t`
  - `num_n`
  - `num_other`
  - `seq_len`
- [ ] Mask-related optional columns:
  - `is_masked_ref`
  - `is_masked_blacklist`
- [ ] Default behavior should favor sparse output; zero-filled bins are opt-in.

#### 10.7.3 BED Coverage and Pileup Kernels
- [ ] `bam_bedcov(path, bed_path, min_mapq := 0, require_flags := NULL, exclude_flags := DEFAULT_BEDCOV_EXCLUDES, count_deletions := TRUE, count_ref_skips := TRUE, region := NULL, index_path := NULL, reference := NULL)`
- [ ] Stable output columns:
  - `chrom`
  - `start`
  - `end`
  - `name` (nullable if BED has no name field)
  - `score` (nullable passthrough)
  - `strand` (nullable passthrough)
  - `bases_covered`
  - `read_count` (optional when requested)
- [ ] `bam_coverage(path, [bam_list := NULL, region := NULL, min_read_len := 0, min_mapq := 0, min_baseq := 0, require_flags := NULL, exclude_flags := DEFAULT_COVERAGE_EXCLUDES, max_depth := 1000000, min_depth := 1, index_path := NULL, reference := NULL])`
- [ ] Stable output columns:
  - `chrom`
  - `start`
  - `end`
  - `numreads`
  - `covbases`
  - `coverage`
  - `meandepth`
  - `meanbaseq`
  - `meanmapq`
- [ ] `bam_depth(path, [all_positions := FALSE, all_reference_positions := FALSE, bed_path := NULL, region := NULL, min_read_len := 0, min_mapq := 0, min_baseq := 0, include_deletions := FALSE, suppress_overlaps := FALSE, require_flags := NULL, exclude_flags := DEFAULT_DEPTH_EXCLUDES, index_path := NULL, reference := NULL])`
- [ ] Stable output columns:
  - `chrom`
  - `pos`
  - `depth`
- [ ] `bam_pileup(path, [min_mapq, min_baseq, require_flags, exclude_flags, include_deletions, include_insertions, distinguish_strands, region, index_path, reference])`
- [ ] Keep overlap handling in `bam_depth` / `bam_pileup`, not as a substitute for CNV bin semantics.

#### 10.7.4 Export / Interoperability Layer
- [ ] Add `bam_bin_export(...)` as a procedure/writer layer over canonical table outputs.
- [ ] Support:
  - sorted `.bed.gz` and `.tsv.gz`,
  - `.tbi` / `.csi` index generation,
  - deterministic schema suitable for DuckHTS and external tools.
- [ ] Canonical on-disk required fields:
  - `chrom`
  - `start`
  - `end`
  - `count_total`
- [ ] Recommended shared downstream fields:
  - `count_fwd`
  - `count_rev`
  - `bin_id`
- [ ] Reserve optional fields:
  - `gc_content`
  - `gc_count`
  - `at_count`
  - `n_count`
  - `mappability`
  - `quality_*`

#### 10.7.5 GC Utilities
- [ ] Document the SQL-first GC path already available:
  - `read_bam(...)` + `seq_gc_content(SEQ)` + flag/MAPQ filters
- [ ] Add native GC helpers only if benchmarks justify them:
  - `bam_gc_stats(...)`
  - or reference-aware GC columns folded into `bam_bin_counts(...)`

### 10.8 Performance Principles
- [ ] **Exploit DuckDB parallelism**:
  - Treat `duckdb_init_set_max_threads` plus contig-level partitioning as the primary full-scan fast path, not an unconditional rule for every function.
  - Each DuckDB thread opens its own `samFile*` + `hts_idx_t*` + `sam_hdr_t*` (no file handle sharing).
  - Contigs are claimed via `__sync_fetch_and_add` on an atomic counter when contig partitioning is actually used.
  - Re-check per-thread `hts_set_threads(fp, 2)` against oversubscription and small-region workloads before standardizing it.
  - For `bam_bin_counts`, per-contig bin arrays are independent — zero cross-thread synchronization during accumulation.
- [ ] Keep a fast path for the common CNV case:
  - fixed bins + read-start counting should avoid interval-tree overhead.
- [ ] Compute totals and strand splits in a single scan when requested together.
- [ ] Keep hot paths streaming and bounded-memory.
- [ ] Use index-backed region restriction where it helps, but preprocess overlapping/unsorted region sets explicitly.
- [ ] Benchmark native kernels against the existing SQL baseline.
- [ ] Do not bake downstream-specific quirks into generic defaults.

### 10.9 Conformance Matrix to Freeze Early
- [ ] Pin behavior across these axes before implementation:
  - counting unit,
  - pair handling,
  - duplicate handling,
  - improper / supplementary / secondary handling,
  - strandedness,
  - coordinate conventions,
  - region semantics,
  - export contract,
  - empty / edge cases.
- [ ] Add conformance fixtures for:
  - overlapping mates,
  - non-overlapping mates,
  - duplicate starts,
  - discordant pairs,
  - cross-contig mates,
  - empty regions and zero-count bins.

**Conformance reference table:**

| Aspect | WisecondorX | NIPTeR | Rsamtools pileup | samtools mpileup | DuckHTS `bam_bin_counts` | DuckHTS `bam_pileup` |
|--------|-------------|--------|------------------|------------------|--------------------------|----------------------|
| **Count unit** | read start | read start | pileup entry (base) | pileup entry (base) | read start | pileup entry |
| **Pair handling** | proper-pair gate | none (all reads) | none | none | via `require_flags` | via `require_flags` |
| **MAPQ filter** | ≥ 1 | implicit | configurable | configurable | `min_mapq` param | `min_mapq` param |
| **BASEQ filter** | none | none | configurable | configurable | n/a (bin-level) | `min_baseq` param |
| **Duplicate handling** | heuristic (pos+mate) | none | flag-based | flag-based | via `exclude_flags` or explicit compatibility mode if needed | `exclude_flags` |
| **Improper pairs** | excluded | included | included | included | via `require_flags` / `exclude_flags` | via `exclude_flags` |
| **Strand tracking** | no | yes (F/R) | optional | implicit | `strand_mode` | `distinguish_strands` |
| **Overlap (mates)** | n/a (read start) | n/a (read start) | double-counted | htslib-managed | n/a (read start) | verify exact `bam_pileup` contract before naming params |
| **Coordinates** | 0-based | 1-based R | 1-based R | 1-based | 0-based half-open | verify before freezing; likely base-position rather than BED-style |
| **Output format** | NPZ | R matrix | R data.frame | TSV | DuckDB table / .bed.gz | DuckDB table |

**Conformance test patterns (`test/sql/bam_bin_counts.test`):**
```
require duckhts

# --- WisecondorX-like parity: read-start, proper pair, MAPQ≥1 ---
query ITII
SELECT chrom, start, count_total, count_fwd
FROM bam_bin_counts('__WORKING_DIRECTORY__/test/data/range.bam',
     bin_width := 1000, require_flags := 0x2, min_mapq := 1)
WHERE chrom = 'CHROMOSOME_I' AND start = 0;
----
CHROMOSOME_I	0	2	1

# --- parity: SQL baseline must match native kernel ---
query I
SELECT COUNT(*) FROM (
    SELECT n.chrom, n.start, n.count_total AS native, s.ct AS sql_ct
    FROM bam_bin_counts('__WORKING_DIRECTORY__/test/data/range.bam',
         bin_width := 1000) n
    JOIN (
        SELECT RNAME AS chrom, (POS / 1000) * 1000 AS start, COUNT(*) AS ct
        FROM read_bam('__WORKING_DIRECTORY__/test/data/range.bam')
        WHERE NOT is_secondary(FLAG) AND NOT is_supplementary(FLAG)
          AND NOT is_unmapped(FLAG) AND NOT is_duplicate(FLAG)
        GROUP BY RNAME, (POS / 1000)
    ) s ON n.chrom = s.chrom AND n.start = s.start
    WHERE n.count_total != s.ct
);
----
0

# --- strand split: count_total = count_fwd + count_rev ---
query I
SELECT COUNT(*) FROM (
    SELECT * FROM bam_bin_counts('__WORKING_DIRECTORY__/test/data/range.bam',
         bin_width := 1000, strand_mode := 'split')
    WHERE count_total != count_fwd + count_rev
);
----
0

# --- empty region returns zero rows, not an error ---
query I
SELECT COUNT(*) FROM bam_bin_counts('__WORKING_DIRECTORY__/test/data/range.bam',
     bin_width := 1000, region := 'nonexistent_contig');
----
0
```

### 10.9A Extension-Level Test Plan

Keep extension tests split by feature family under `test/sql/` rather than continuing to grow `duckhts.test`.

**Recommended SQL test files:**
- `test/sql/index_build.test`
  - `bgzip`, `bgunzip`, `bam_index`, `bcf_index`, `tabix_index`
  - verify output rows, reported index format, overwrite behavior, and round-trip usability through existing readers
- `test/sql/read_bed.test`
  - plain BED, BED4/BED6 passthrough columns, BED with extra columns, header handling if supported, `.bed.gz` indexed region reads
- `test/sql/fasta_nuc.test`
  - known interval composition fixtures with exact counts for A/C/G/T/N/other, `%AT`, `%GC`, `seq_len`, optional sequence emission
- `test/sql/interval_udf.test`
  - `interval_merge`, `interval_overlap`, `interval_nearest` with sorted/unsorted inputs, touching intervals, empty inputs, and contig partitioning
- `test/sql/bam_bin_counts.test`
  - SQL baseline parity, strand invariants, region semantics, flag-filter semantics, proper-pair behavior, and optional annotation columns
- `test/sql/bam_bedcov.test`
  - parity against stable `samtools bedcov` expectations for deletions/ref-skips, counts, and threshold columns
- `test/sql/bam_coverage.test`
  - region summary columns and stable fixture expectations for `numreads`, `covbases`, `coverage`, `meandepth`, `meanbaseq`, `meanmapq`
- `test/sql/bam_depth.test`
  - `-a`/`-aa`-style semantics, `include_deletions`, overlap suppression, BED-guided restriction, and multi-file behavior if supported
- `test/sql/bam_pileup.test`
  - Rsamtools-style base-level behavior, minimum base/map quality, strand distinction, insertions/deletions toggles

**Extension test matrix by concern:**
- Schema contract:
  - exact column names/order for BED-compatible outputs
  - output types stable under optional arguments
- Region/index behavior:
  - indexed region queries
  - missing contig returns zero rows, not errors, once that contract is frozen
  - explicit `index_path` overrides
- Semantic parity:
  - SQL-vs-native parity for `bam_bin_counts`
  - fixture-based parity for `fasta_nuc`
  - stable output expectations mirroring documented `samtools` / `bedtools` behavior
- Flag and overlap behavior:
  - duplicates, secondary, supplementary, unmapped, QC-fail
  - proper/improper pair handling through `require_flags` / `exclude_flags`
  - overlapping mates, discordant pairs, cross-contig mates
- I/O round-trip:
  - `COPY`/writer output -> `bgzip` -> `tabix_index` -> `read_tabix`
  - index builders produce files consumable by existing readers

**Fixture plan:**
- Keep deterministic small fixtures under `test/data/`.
- Extend `test/scripts/prepare_test_data.sh` when new indexed files or derived fixtures are needed.
- Add tiny synthetic fixtures specifically for:
  - BED composition checks for `fasta_nuc`
  - duplicate/proper-pair/overlap edge cases for BAM depth/counting
  - masked/blacklisted interval joins
- Any fixture needed by package tests should also be copied into `r/Rduckhts/inst/extdata/` through the existing bootstrap/package flow.

### 10.9B R Package-Level Test Plan

The R package tests should validate the packaged product, not just the extension source tree.

**Recommended `tinytest` layout:**
- `r/Rduckhts/inst/tinytest/test_wrappers_basic.R`
  - exported wrapper existence
  - wrapper formals/signatures
  - argument validation and error messages
  - `rduckhts_functions()` catalog availability
- `r/Rduckhts/inst/tinytest/test_load_and_catalog.R`
  - `rduckhts_load()` against a fresh DuckDB connection
  - generated `functions.tsv` present and consistent with installed wrappers
- `r/Rduckhts/inst/tinytest/test_readers_integration.R`
  - current read-path wrappers (`bcf`, `bam`, `fasta`, `fastq`, `gff`, `gtf`, `tabix`, metadata helpers)
- `r/Rduckhts/inst/tinytest/test_indexes_bgzip.R`
  - `bgzip`, `bgunzip`, `bam_index`, `bcf_index`, `tabix_index`
  - temporary-file workflow from R into DuckDB and back to packaged readers
- `r/Rduckhts/inst/tinytest/test_phase10_coverage.R`
  - `read_bed`, `fasta_nuc`, `bam_bin_counts`, `bam_bedcov`, `bam_coverage`, `bam_depth`, `bam_pileup`
  - wrapper argument mapping to SQL named parameters
  - returned tables and data frames with expected columns
- `r/Rduckhts/inst/tinytest/test_installed_extdata.R`
  - required package fixtures exist after tarball install
  - indexed companion files (`.bai`, `.csi`, `.tbi`, `.fai`) are present where expected

**R package assertions should cover:**
- packaged extension loads without relying on the source checkout
- wrapper names and defaults stay aligned with `functions.yaml`
- installed `inst/extdata` contains every file used in examples/tests
- temporary output files created by SQL functions are visible and reusable from R
- DBI workflows succeed using only installed package assets

**Release-gate commands to plan around:**
- `cd r/Rduckhts && R CMD build .`
- install the produced tarball in a clean library
- `R -e "tinytest::test_package('Rduckhts')"`
- `R CMD check --as-cran <tarball>`

**README / documentation verification:**
- If Phase 10 adds user-facing functions, add at least one installed-package example path that exercises each major family.
- Prefer examples that use bundled `inst/extdata`.
- If `README.Rmd` rendering is blocked by the harness, record that explicitly and ask for user help before treating README verification as complete.

### 10.10 Implementation Sequence
- [ ] Implement the independent utility functions first:
  - `bgzip`
  - `bgunzip`
  - `bam_index`
  - `bcf_index`
  - `tabix_index`
- [ ] For each newly implemented public function above:
  - add/update the `functions.yaml` entry
  - run `python3 scripts/render_function_catalog.py`
  - wire the function into `src/duckhts.c`
  - wire any new `.c` file into top-level `CMakeLists.txt`
  - update `r/Rduckhts/R/bootstrap.R` so the file is copied into the package bundle
  - update `r/Rduckhts/configure` and `r/Rduckhts/configure.win` if they enumerate or compile sources explicitly
  - run the R bootstrap refresh and package tests
- [ ] Freeze the minimal public API and canonical BED-compatible output schemas for the remaining Phase 10 functions.
- [ ] Define consistent region semantics across existing readers before adding more region-aware APIs.
- [ ] Integrate `cgranges` and implement `read_bed` plus interval primitives.
- [ ] Implement `fasta_nuc` and align its annotation columns with optional `bam_bin_counts` bin annotations.
- [ ] Implement `bam_bin_counts` as the primary CNV-oriented fixed-bin kernel.
- [ ] Implement `bam_bedcov` for samtools-compatible BED-guided summed coverage.
- [ ] Implement `bam_coverage` and `bam_depth` to mirror `samtools coverage` and `samtools depth`.
- [ ] Implement separate `bam_pileup` for richer Rsamtools-like base-level output.
- [ ] Add export/index writer logic and round-trip tests.
- [ ] Document SQL-first GC recipes and defer native GC until benchmarks justify it.
- [ ] Add R wrappers for all public functions.
- [ ] Add the corresponding extension-level `.test` files and R `tinytest` files as each function family lands; do not defer all testing to the end.
- [ ] Update `functions.yaml`, regenerate catalogs/descriptors, bootstrap the R package copy, and add conformance + performance tests.

**New source files to add:**

| File | Purpose |
|------|---------|
| `src/hts_index_builder.c` | `bam_index`, `bcf_index`, `tabix_index` (Phase 8) |
| `src/bgzip.c` | `bgzip`, `bgunzip` BGZF compression/decompression (Phase 8) |
| `src/interval_udf.c` | `read_bed`, `fasta_nuc`, `interval_merge`, `interval_overlap`, `interval_nearest` |
| `src/bam_bin_counts.c` | `bam_bin_counts` fixed-bin read-start counting |
| `src/bam_bedcov.c` | `bam_bedcov` BED-guided summed depth |
| `src/bam_coverage.c` | `bam_coverage` samtools-style region summary |
| `src/bam_depth.c` | `bam_depth` samtools-style per-position depth |
| `src/bam_pileup.c` | `bam_pileup` Rsamtools-style base-level pileup |

**Registration and metadata notes:**
- Register only concrete public functions in `src/duckhts.c`; do not add compatibility macros unless they prove necessary after the base APIs stabilize.
- Every public function added here must get:
  - a `functions.yaml` entry,
  - an R wrapper in `Rduckhts`,
  - registration in the extension entrypoint,
  - build integration in top-level `CMakeLists.txt`,
  - package-side source integration in `r/Rduckhts/R/bootstrap.R`,
  - package build integration in `r/Rduckhts/configure` and `r/Rduckhts/configure.win`,
  - SQL tests and R tests,
  - README/function catalog documentation.

**Immediate implementation priority:**
1. `bgzip`
2. `bgunzip`
3. `bam_index`
4. `bcf_index`
5. `tabix_index`

These are largely independent of the BAM statistics API design and should be implemented first.

**Package sync checklist for new C sources:**
- [ ] Add the file to `CMakeLists.txt`
- [ ] Add the file to any source lists copied by `r/Rduckhts/R/bootstrap.R`
- [ ] Ensure `r/Rduckhts/configure` compiles/links it on Unix
- [ ] Ensure `r/Rduckhts/configure.win` compiles/links it on Windows
- [ ] Re-run `Rscript bootstrap.R ~/duckhts/` from `r/Rduckhts/`
- [ ] Rebuild/install the package tarball before running R tests

## 📝 Notes
- **Scope**: keep the extension focused on readers and analytics building blocks; downstream application modeling stays outside the extension.
- **Target**: maintain DuckDB 1.4+ C/C++ extension API compatibility.
- **Documentation**: primary user-facing documentation lives in `README.Rmd`; function metadata lives in `functions.yaml`.
- **Build**: keep the package self-contained and CRAN-compatible; no `vcpkg` in the R package path.
- **Interop goal**: one family of canonical tabular outputs should be reusable across WisecondorX-like, NIPTeR-like, and broader CNV workflows.

## 🔍 Review Feedback (2026-02-10)
- bcf_reader: region lookup currently errors for both “contig not found” and “no overlapping records” (TODO in `src/bcf_reader.c`); consider distinguishing to avoid false failures on empty regions.
- seq_reader: paired FASTQ path assumes reads are in lockstep but does not validate QNAME pairing; mismatched mates will silently pair (recommend add name check + test).
- seq_reader: interleaved mode toggles mate 1/2 regardless of QNAME suffix or odd record count; consider handling trailing unpaired read and/or validating suffixes.
- tabix_reader (generic): column count inferred from first non-# line only; files with variable columns or non-# meta-char may mis-bind schema (consider using tabix conf/meta-char when indexed, and add tests for varying columns).
- behavior consistency: bcf_reader errors on “region not found” while tabix_reader returns empty; decide on a consistent contract and document it in README/tests.
