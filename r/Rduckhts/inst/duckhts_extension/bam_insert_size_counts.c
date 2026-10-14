/* Insert-size semantics independently implemented from Fulcrum Genomics Riker
 * (MIT), v0.4.1 / 40d5e326. No read-name pairing or sequence/quality/tag decoding.
 * Each indexed worker owns its htslib handle and sparse histogram.
 *
 * Upstream notice: Copyright (c) 2026 Fulcrum Genomics LLC (MIT).
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN
#include "hts_io_tuning.h"
#include <htslib/sam.h>
#include <htslib/khash.h>
#include <limits.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

KHASH_MAP_INIT_INT64(isize, uint64_t)
#define ISIZE_MAX_BINS 1000000u
#define ISIZE_WINDOW 1000000
#define ISIZE_MAX_WORKERS 8

typedef struct {
    char *path;
    char *reference;
    char *index_path;
    bool duplicates;
    bool sequential;
    int threads;
    int decompression_threads;
    sam_hdr_t *header;
    uint64_t *task_offsets;
    uint64_t tasks;
} isize_bind_t;

typedef struct {
    const isize_bind_t *bind;
    atomic_uint_fast64_t next_task;
} isize_shared_t;

typedef struct {
    isize_shared_t *shared;
    khash_t(isize) *hist;
    const char *error;
} isize_worker_t;

typedef struct {
    khash_t(isize) *hist;
    khiter_t cursor;
} isize_result_t;

static void isize_bind_destroy(void *ptr) {
    isize_bind_t *bind = ptr;
    if (!bind) return;
    free(bind->path);
    free(bind->reference);
    free(bind->index_path);
    free(bind->task_offsets);
    sam_hdr_destroy(bind->header);
    free(bind);
}

static void isize_result_destroy(void *ptr) {
    isize_result_t *result = ptr;
    if (!result) return;
    kh_destroy(isize, result->hist);
    free(result);
}

static char *isize_string(duckdb_bind_info info, const char *name) {
    duckdb_value value = name ? duckdb_bind_get_named_parameter(info, name)
                             : duckdb_bind_get_parameter(info, 0);
    char *copy = NULL;
    if (value && !duckdb_is_null_value(value)) {
        char *text = duckdb_get_varchar(value);
        if (text) {
            copy = strdup(text);
            duckdb_free(text);
        }
        if (!copy) duckdb_bind_set_error(info, "insert-size counts: out of memory");
    }
    if (value) duckdb_destroy_value(&value);
    return copy;
}

static samFile *isize_open(const isize_bind_t *bind) {
    samFile *fp = sam_open(bind->path, "r");
    if (!fp) return NULL;
    duckhts_apply_remote_hts_tuning(fp, bind->path,
        bind->sequential ? DUCKHTS_HTS_IO_PROFILE_STREAMING : DUCKHTS_HTS_IO_PROFILE_INDEXED_REGION);
    /* Riker sam/pair_orientation.rs:43-51 needs CIGAR only for the reverse
     * read's end. CRAM can omit SEQ, QUAL, AUX, QNAME and MAPQ entirely. */
    if ((bind->reference && hts_set_opt(fp, CRAM_OPT_REFERENCE, bind->reference) < 0) ||
        (hts_get_format(fp)->format == cram &&
         hts_set_opt(fp, CRAM_OPT_REQUIRED_FIELDS,
                     SAM_FLAG | SAM_RNAME | SAM_POS | SAM_CIGAR | SAM_RNEXT | SAM_PNEXT | SAM_TLEN) < 0) ||
        (bind->decompression_threads && hts_set_threads(fp, bind->decompression_threads) < 0)) {
        sam_close(fp);
        return NULL;
    }
    return fp;
}

static void isize_bind(duckdb_bind_info info) {
    isize_bind_t *bind = calloc(1, sizeof(*bind));
    if (!bind) {
        duckdb_bind_set_error(info, "insert-size counts: out of memory");
        return;
    }
    duckdb_bind_set_bind_data(info, bind, isize_bind_destroy);
    bind->threads = 1;
    bind->path = isize_string(info, NULL);
    bind->reference = isize_string(info, "reference");
    bind->index_path = isize_string(info, "index_path");
    char *reference_index = isize_string(info, "reference_index_path");
    if (reference_index) {
        if (!bind->reference) {
            duckdb_bind_set_error(info, "reference_index_path requires reference");
            free(reference_index);
            return;
        }
        char *locator = duckhts_reference_locator(bind->reference, reference_index);
        free(reference_index);
        if (!locator) {
            duckdb_bind_set_error(info, "insert-size counts: out of memory");
            return;
        }
        free(bind->reference);
        bind->reference = locator;
    }
    if (!bind->path || !bind->path[0]) {
        duckdb_bind_set_error(info, "insert-size counts: path must not be NULL or empty");
        return;
    }
    duckdb_value value = duckdb_bind_get_named_parameter(info, "include_duplicates");
    if (value && !duckdb_is_null_value(value)) bind->duplicates = duckdb_get_bool(value);
    if (value) duckdb_destroy_value(&value);
    const char *options[] = {"threads", "decompression_threads"};
    for (int i = 0; i < 2; i++) {
        value = duckdb_bind_get_named_parameter(info, options[i]);
        if (value && !duckdb_is_null_value(value)) {
            int64_t n = duckdb_get_int64(value);
            if (n < (i == 0 ? 1 : 0) || n > (i == 0 ? ISIZE_MAX_WORKERS : INT_MAX)) {
                duckdb_bind_set_error(info, "threads must be 1..8; decompression_threads must be 0..INT_MAX");
                duckdb_destroy_value(&value);
                return;
            }
            if (i == 0) bind->threads = (int)n;
            else bind->decompression_threads = (int)n;
        }
        if (value) duckdb_destroy_value(&value);
    }
#if defined(__EMSCRIPTEN__)
    if (bind->threads != 1 || bind->decompression_threads != 0) {
        duckdb_bind_set_error(info, "insert-size counts: wasm requires threads=1 and decompression_threads=0");
        return;
    }
#endif
    char *mode = isize_string(info, "scan_mode");
    bool require_index = false;
    if (mode) {
        bind->sequential = strcmp(mode, "sequential") == 0;
        require_index = strcmp(mode, "indexed") == 0;
        if (!bind->sequential && !require_index && strcmp(mode, "auto") != 0) {
            duckdb_bind_set_error(info, "scan_mode must be auto, indexed or sequential");
            free(mode);
            return;
        }
        free(mode);
    }
    samFile *fp = isize_open(bind);
    if (!fp) {
        duckdb_bind_set_error(info, "insert-size counts: cannot open/configure alignment file");
        return;
    }
    bind->header = sam_hdr_read(fp);
    if (!bind->header) {
        sam_close(fp);
        duckdb_bind_set_error(info, "insert-size counts: cannot read alignment header");
        return;
    }
    if (!bind->sequential) {
        hts_idx_t *index = sam_index_load3(fp, bind->path, bind->index_path,
            duckhts_index_save_remote_flag(bind->path, bind->index_path) | HTS_IDX_SILENT_FAIL);
        if (!index && (require_index || bind->index_path)) {
            sam_close(fp);
            duckdb_bind_set_error(info, "insert-size counts: cannot load required index");
            return;
        }
        bind->sequential = index == NULL;
        hts_idx_destroy(index);
    }
    sam_close(fp);
    if (bind->sequential) bind->threads = 1;
    int nref = sam_hdr_nref(bind->header);
    bind->task_offsets = calloc((size_t)nref + 1, sizeof(uint64_t));
    if (!bind->task_offsets) {
        duckdb_bind_set_error(info, "insert-size counts: out of memory");
        return;
    }
    for (int tid = 0; tid < nref; tid++) {
        hts_pos_t length = sam_hdr_tid2len(bind->header, tid);
        bind->tasks += (uint64_t)(length / ISIZE_WINDOW + (length % ISIZE_WINDOW != 0));
        bind->task_offsets[tid + 1] = bind->tasks;
    }
    const char *names[] = {"pair_orientation", "insert_size", "read_pairs"};
    const duckdb_type types[] = {DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_BIGINT, DUCKDB_TYPE_UBIGINT};
    for (int i = 0; i < 3; i++) {
        duckdb_logical_type type = duckdb_create_logical_type(types[i]);
        duckdb_bind_add_result_column(info, names[i], type);
        duckdb_destroy_logical_type(&type);
    }
}

static bool isize_add(isize_worker_t *worker, uint64_t key, uint64_t count) {
    khiter_t slot = kh_get(isize, worker->hist, key);
    if (slot == kh_end(worker->hist)) {
        if (kh_size(worker->hist) >= ISIZE_MAX_BINS) {
            worker->error = "insert-size counts: one-million-bin histogram limit exceeded";
            return false;
        }
        int inserted;
        slot = kh_put(isize, worker->hist, key, &inserted);
        if (inserted < 0) {
            worker->error = "insert-size counts: out of memory growing histogram";
            return false;
        }
        kh_value(worker->hist, slot) = 0;
    }
    if (UINT64_MAX - kh_value(worker->hist, slot) < count) {
        worker->error = "insert-size counts: read-pair count overflow";
        return false;
    }
    kh_value(worker->hist, slot) += count;
    return true;
}

static bool isize_count(isize_worker_t *worker, const bam1_t *record) {
    const bam1_core_t *core = &record->core;
    uint16_t flag = core->flag;
    /* riker v0.4.1 commands/isize.rs:197-206: paired, both mapped,
     * primary, not first-segment or QC-fail, nonzero TLEN; duplicates optional.
     * MAPQ and proper-pair status are not predicates. Picard is context only. */
    if (!(flag & BAM_FPAIRED) ||
        (flag & (BAM_FUNMAP | BAM_FMUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL)) ||
        (flag & BAM_FREAD1) || core->isize == 0 ||
        (!worker->shared->bind->duplicates && (flag & BAM_FDUP))) return true;
    /* riker v0.4.1 sam/pair_orientation.rs:30-40: same reference required;
     * equal strands are TANDEM. Lines 43-59 compare mate start with inclusive
     * alignment end for reverse reads, or POS with POS + TLEN for forward
     * reads. Keep both positions zero-based here; equality is RF. */
    if (core->tid < 0 || core->mtid != core->tid || core->pos < 0 || core->mpos < 0) return true;
    unsigned orientation = 2;
    bool reverse = (flag & BAM_FREVERSE) != 0;
    bool mate_reverse = (flag & BAM_FMREVERSE) != 0;
    if (reverse != mate_reverse) {
        hts_pos_t positive = reverse ? core->mpos : core->pos;
        hts_pos_t negative = reverse ? bam_endpos(record) - 1 : core->pos + core->isize;
        /* riker v0.4.1 sam/pair_orientation.rs:56-59. */
        orientation = positive < negative ? 0 : 1;
    }
    uint64_t size = core->isize < 0 ? (uint64_t)(-(core->isize + 1)) + 1 : (uint64_t)core->isize;
    if (size > UINT32_MAX) {
        worker->error = "insert-size counts: absolute TLEN exceeds UINT32_MAX";
        return false;
    }
    /* riker v0.4.1 commands/isize.rs:211 uses unsigned absolute TLEN;
     * widen before negation so INT32_MIN is preserved. */
    return isize_add(worker, ((uint64_t)orientation << 32) | size, 1);
}

static void *isize_worker_main(void *ptr) {
    isize_worker_t *worker = ptr;
    const isize_bind_t *bind = worker->shared->bind;
    samFile *fp = isize_open(bind);
    sam_hdr_t *header = fp ? sam_hdr_read(fp) : NULL;
    bam1_t *record = bam_init1();
    hts_idx_t *index = NULL;
    hts_itr_t *iterator = NULL;
    if (!fp || !header || !record) {
        worker->error = "insert-size counts: cannot initialize worker";
        goto cleanup;
    }
    int status = -1;
    if (bind->sequential) {
        while ((status = sam_read1(fp, header, record)) >= 0) {
            if (!isize_count(worker, record)) goto cleanup;
        }
    } else {
        index = sam_index_load3(fp, bind->path, bind->index_path,
            duckhts_index_save_remote_flag(bind->path, bind->index_path));
        if (!index) {
            worker->error = "insert-size counts: cannot load worker index";
            goto cleanup;
        }
        uint64_t task;
        while ((task = atomic_fetch_add(&worker->shared->next_task, 1)) < bind->tasks) {
            int tid = 0;
            while (bind->task_offsets[tid + 1] <= task) tid++;
            hts_pos_t beg = (hts_pos_t)(task - bind->task_offsets[tid]) * ISIZE_WINDOW;
            hts_pos_t end = beg + ISIZE_WINDOW;
            hts_pos_t length = sam_hdr_tid2len(header, tid);
            if (end > length) end = length;
            iterator = sam_itr_queryi(index, tid, beg, end);
            if (!iterator) {
                worker->error = "insert-size counts: cannot create indexed iterator";
                goto cleanup;
            }
            while ((status = sam_itr_next(fp, iterator, record)) >= 0) {
                /* Overlapping alignments belong only to the window containing
                 * their start. This preserves exactly-once eligible-record counts. */
                if (record->core.pos < beg || record->core.pos >= end) continue;
                if (!isize_count(worker, record)) goto cleanup;
            }
            hts_itr_destroy(iterator);
            iterator = NULL;
            if (status < -1) break;
        }
    }
    if (status < -1) worker->error = "insert-size counts: alignment decode failed";
cleanup:
    hts_itr_destroy(iterator);
    hts_idx_destroy(index);
    bam_destroy1(record);
    sam_hdr_destroy(header);
    if (fp) sam_close(fp);
    return NULL;
}

static void isize_init(duckdb_init_info info) {
    const isize_bind_t *bind = duckdb_init_get_bind_data(info);
    isize_result_t *result = calloc(1, sizeof(*result));
    if (!result) {
        duckdb_init_set_error(info, "insert-size counts: out of memory");
        return;
    }
    duckdb_init_set_init_data(info, result, isize_result_destroy);
    duckdb_init_set_max_threads(info, 1);
    isize_shared_t shared = {.bind = bind};
    atomic_init(&shared.next_task, 0);
    isize_worker_t workers[ISIZE_MAX_WORKERS] = {0};
    pthread_t threads[ISIZE_MAX_WORKERS];
    int started = 0;
    const char *error = NULL;
    for (int i = 0; i < bind->threads; i++) {
        workers[i].shared = &shared;
        workers[i].hist = kh_init(isize);
        if (!workers[i].hist) {
            error = "insert-size counts: out of memory";
            break;
        }
        if (bind->threads == 1) {
            isize_worker_main(&workers[i]);
        } else {
#if !defined(__EMSCRIPTEN__)
            if (pthread_create(&threads[i], NULL, isize_worker_main, &workers[i]) != 0) {
                error = "insert-size counts: cannot create worker thread";
                break;
            }
            started++;
#endif
        }
    }
#if !defined(__EMSCRIPTEN__)
    for (int i = 0; i < started; i++) pthread_join(threads[i], NULL);
#endif
    for (int i = 0; i < bind->threads; i++) {
        if (workers[i].error) error = workers[i].error;
    }
    if (!error) {
        for (int i = 1; i < bind->threads && !error; i++) {
            for (khiter_t slot = kh_begin(workers[i].hist); slot != kh_end(workers[i].hist); slot++) {
                if (kh_exist(workers[i].hist, slot) &&
                    !isize_add(&workers[0], kh_key(workers[i].hist, slot), kh_value(workers[i].hist, slot))) {
                    error = workers[0].error;
                    break;
                }
            }
        }
    }
    result->hist = workers[0].hist;
    for (int i = 1; i < bind->threads; i++) kh_destroy(isize, workers[i].hist);
    if (error) duckdb_init_set_error(info, error);
}

static void isize_scan(duckdb_function_info info, duckdb_data_chunk output) {
    isize_result_t *result = duckdb_function_get_init_data(info);
    duckdb_vector orientation = duckdb_data_chunk_get_vector(output, 0);
    int64_t *sizes = duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 1));
    uint64_t *counts = duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 2));
    static const char *const names[] = {"FR", "RF", "TANDEM"};
    idx_t n = 0;
    while (result->cursor != kh_end(result->hist) && n < duckdb_vector_size()) {
        khiter_t slot = result->cursor++;
        if (!kh_exist(result->hist, slot)) continue;
        uint64_t key = kh_key(result->hist, slot);
        duckdb_vector_assign_string_element(orientation, n, names[key >> 32]);
        sizes[n] = (int64_t)(key & UINT32_MAX);
        counts[n] = kh_value(result->hist, slot);
        n++;
    }
    duckdb_data_chunk_set_size(output, n);
}

void register_duckhts_bam_insert_size_counts_function(duckdb_connection connection) {
    duckdb_table_function function = duckdb_create_table_function();
    duckdb_table_function_set_name(function, "duckhts_bam_insert_size_counts");
    duckdb_logical_type text = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type integer = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type boolean = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_table_function_add_parameter(function, text);
    duckdb_table_function_add_named_parameter(function, "reference", text);
    duckdb_table_function_add_named_parameter(function, "reference_index_path", text);
    duckdb_table_function_add_named_parameter(function, "index_path", text);
    duckdb_table_function_add_named_parameter(function, "scan_mode", text);
    duckdb_table_function_add_named_parameter(function, "include_duplicates", boolean);
    duckdb_table_function_add_named_parameter(function, "threads", integer);
    duckdb_table_function_add_named_parameter(function, "decompression_threads", integer);
    duckdb_table_function_set_bind(function, isize_bind);
    duckdb_table_function_set_init(function, isize_init);
    duckdb_table_function_set_function(function, isize_scan);
    duckdb_register_table_function(connection, function);
    duckdb_destroy_logical_type(&text);
    duckdb_destroy_logical_type(&integer);
    duckdb_destroy_logical_type(&boolean);
    duckdb_destroy_table_function(&function);
}
