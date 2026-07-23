/**
 * DuckHTS samtools idxstats-compatible table function.
 *
 * Compatibility target:
 *   samtools idxstats from samtools 1.23.1
 * Primary upstream reference:
 *   .sync/samtools/bam_index.c
 */

#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "include/hts_io_tuning.h"

#ifndef SAM_RNAME
#define SAM_RNAME SAM_POS
#endif

typedef struct {
    char *path;
    char *output_path;
    char *index_path;
    int threads;
    int overwrite;
    int used_index_fast_path;
    int emitted;
} samtools_idxstats_bind_t;

static char *append_suffix_idxstats(const char *path, const char *suffix) {
    size_t path_len;
    size_t suffix_len;
    char *out;
    if (!path || !suffix) return NULL;
    path_len = strlen(path);
    suffix_len = strlen(suffix);
    out = (char *)duckdb_malloc(path_len + suffix_len + 1);
    if (!out) return NULL;
    memcpy(out, path, path_len);
    memcpy(out + path_len, suffix, suffix_len + 1);
    return out;
}

static int file_exists_idxstats(const char *path) {
    FILE *fp;
    if (!path) return 0;
    fp = fopen(path, "rb");
    if (!fp) return 0;
    fclose(fp);
    return 1;
}

static int prepare_output_path_idxstats(const char *path, int overwrite, char *err, size_t errlen) {
    if (!path || !path[0]) {
        snprintf(err, errlen, "duckhts_samtools_idxstats: output must be non-empty");
        return -1;
    }
    if (!file_exists_idxstats(path)) return 0;
    if (!overwrite) {
        snprintf(err, errlen,
                 "duckhts_samtools_idxstats: output already exists (use overwrite := TRUE to replace)");
        return -1;
    }
    if (remove(path) != 0) {
        snprintf(err, errlen,
                 "duckhts_samtools_idxstats: failed to replace existing output '%s'",
                 path);
        return -1;
    }
    return 0;
}

static int slow_idxstats_to_file(samFile *fp, sam_hdr_t *header, FILE *out, int threads,
                                 char *err, size_t errlen) {
    int ret;
    int last_tid = -2;
    bam1_t *b = NULL;
    uint64_t (*count0)[2] = NULL;
    uint64_t (*counts)[2] = NULL;

    if (!fp || !header || !out) {
        snprintf(err, errlen, "duckhts_samtools_idxstats: internal error");
        return -1;
    }

    if (threads > 0 && hts_set_threads(fp, threads) != 0) {
        snprintf(err, errlen, "duckhts_samtools_idxstats: failed to configure decompression threads");
        return -1;
    }

    if (hts_set_opt(fp, CRAM_OPT_REQUIRED_FIELDS, SAM_RNAME | SAM_FLAG) != 0) {
        /* Keep going for BAM/other formats where this may be ignored. */
    }

    b = bam_init1();
    if (!b) {
        snprintf(err, errlen, "duckhts_samtools_idxstats: failed to allocate BAM record");
        return -1;
    }

    count0 = (uint64_t(*)[2])calloc((size_t)sam_hdr_nref(header) + 1, sizeof(*count0));
    if (!count0) {
        bam_destroy1(b);
        snprintf(err, errlen, "duckhts_samtools_idxstats: out of memory");
        return -1;
    }
    counts = count0 + 1;

    while ((ret = sam_read1(fp, header, b)) >= 0) {
        if (b->core.tid >= sam_hdr_nref(header) || b->core.tid < -1) {
            free(count0);
            bam_destroy1(b);
            snprintf(err, errlen, "duckhts_samtools_idxstats: encountered record with invalid reference id");
            return -1;
        }

        if (b->core.tid != last_tid) {
            if (last_tid >= -1) {
                if (counts[b->core.tid][0] + counts[b->core.tid][1]) {
                    free(count0);
                    bam_destroy1(b);
                    snprintf(err, errlen, "duckhts_samtools_idxstats: file is not position sorted");
                    return -1;
                }
            }
            last_tid = b->core.tid;
        }

        counts[b->core.tid][(b->core.flag & BAM_FUNMAP) ? 1 : 0]++;
    }

    if (ret != -1) {
        free(count0);
        bam_destroy1(b);
        snprintf(err, errlen, "duckhts_samtools_idxstats: failed while scanning input");
        return -1;
    }

    for (int i = 0; i < sam_hdr_nref(header); i++) {
        if (fprintf(out, "%s\t%" PRId64 "\t%" PRIu64 "\t%" PRIu64 "\n",
                    sam_hdr_tid2name(header, i),
                    (int64_t)sam_hdr_tid2len(header, i),
                    counts[i][0], counts[i][1]) < 0) {
            free(count0);
            bam_destroy1(b);
            snprintf(err, errlen, "duckhts_samtools_idxstats: failed to write output");
            return -1;
        }
    }
    if (fprintf(out, "*\t0\t%" PRIu64 "\t%" PRIu64 "\n", counts[-1][0], counts[-1][1]) < 0) {
        free(count0);
        bam_destroy1(b);
        snprintf(err, errlen, "duckhts_samtools_idxstats: failed to write output");
        return -1;
    }

    free(count0);
    bam_destroy1(b);
    return 0;
}

static int fast_idxstats_to_file(samFile *fp, sam_hdr_t *header, const char *path,
                                 const char *index_path, FILE *out, int *used_fast,
                                 char *err, size_t errlen) {
    hts_idx_t *idx = NULL;
    uint64_t mapped = 0, unmapped = 0;
    uint64_t n_no_coor = 0;

    if (!fp || !header || !path || !out || !used_fast) {
        snprintf(err, errlen, "duckhts_samtools_idxstats: internal error");
        return -1;
    }

    *used_fast = 0;
    idx = sam_index_load2(fp, path, index_path);
    if (!idx) return 1;

    for (int i = 0; i < sam_hdr_nref(header); ++i) {
        mapped = 0;
        unmapped = 0;
        hts_idx_get_stat(idx, i, &mapped, &unmapped);
        if (fprintf(out, "%s\t%" PRId64 "\t%" PRIu64 "\t%" PRIu64 "\n",
                    sam_hdr_tid2name(header, i),
                    (int64_t)sam_hdr_tid2len(header, i),
                    mapped, unmapped) < 0) {
            hts_idx_destroy(idx);
            snprintf(err, errlen, "duckhts_samtools_idxstats: failed to write output");
            return -1;
        }
    }

    n_no_coor = hts_idx_get_n_no_coor(idx);
    if (fprintf(out, "*\t0\t0\t%" PRIu64 "\n", n_no_coor) < 0) {
        hts_idx_destroy(idx);
        snprintf(err, errlen, "duckhts_samtools_idxstats: failed to write output");
        return -1;
    }

    hts_idx_destroy(idx);
    *used_fast = 1;
    return 0;
}

static int run_samtools_idxstats(samtools_idxstats_bind_t *bind, char *err, size_t errlen) {
    samFile *fp = NULL;
    sam_hdr_t *header = NULL;
    FILE *out = NULL;
    int rc = -1;

    if (prepare_output_path_idxstats(bind->output_path, bind->overwrite, err, errlen) != 0) {
        return -1;
    }

    fp = sam_open(bind->path, "r");
    if (!fp) {
        snprintf(err, errlen, "duckhts_samtools_idxstats: failed to open alignment file '%s'", bind->path);
        return -1;
    }
    duckhts_apply_remote_hts_tuning(fp, bind->path,
                                    DUCKHTS_HTS_IO_PROFILE_STREAMING);

    header = sam_hdr_read(fp);
    if (!header) {
        snprintf(err, errlen, "duckhts_samtools_idxstats: failed to read alignment header");
        goto cleanup;
    }

    out = fopen(bind->output_path, "wb");
    if (!out) {
        snprintf(err, errlen, "duckhts_samtools_idxstats: failed to open '%s': %s",
                 bind->output_path, strerror(errno));
        goto cleanup;
    }

    if (hts_get_format(fp)->format == bam) {
        int fast_rc = fast_idxstats_to_file(fp, header, bind->path, bind->index_path, out,
                                            &bind->used_index_fast_path, err, errlen);
        if (fast_rc == 0) {
            rc = 0;
            goto cleanup;
        }
        if (fast_rc < 0) {
            goto cleanup;
        }
        rewind(out);
#ifdef _WIN32
        _chsize(_fileno(out), 0);
#else
        if (ftruncate(fileno(out), 0) != 0) {
            snprintf(err, errlen, "duckhts_samtools_idxstats: failed to reset output after fast-path fallback");
            goto cleanup;
        }
#endif
        if (fseek(out, 0, SEEK_SET) != 0) {
            snprintf(err, errlen, "duckhts_samtools_idxstats: failed to reset output after fast-path fallback");
            goto cleanup;
        }
    }

    bind->used_index_fast_path = 0;
    if (slow_idxstats_to_file(fp, header, out, bind->threads, err, errlen) != 0) {
        goto cleanup;
    }

    rc = 0;

cleanup:
    if (out) {
        if (fclose(out) != 0 && rc == 0) {
            snprintf(err, errlen, "duckhts_samtools_idxstats: failed to close output");
            rc = -1;
        }
    }
    if (header) sam_hdr_destroy(header);
    if (fp) sam_close(fp);
    if (rc != 0 && bind->output_path) remove(bind->output_path);
    return rc;
}

static void destroy_samtools_idxstats_bind(void *data) {
    samtools_idxstats_bind_t *bind = (samtools_idxstats_bind_t *)data;
    if (!bind) return;
    if (bind->path) duckdb_free(bind->path);
    if (bind->output_path) duckdb_free(bind->output_path);
    if (bind->index_path) duckdb_free(bind->index_path);
    duckdb_free(bind);
}

static void samtools_idxstats_bind(duckdb_bind_info info) {
    duckdb_value path_val = duckdb_bind_get_parameter(info, 0);
    char *path = duckdb_get_varchar(path_val);
    char *output_path = NULL;
    char *index_path = NULL;
    samtools_idxstats_bind_t *bind = NULL;
    duckdb_value val;
    char err[512];

    memset(err, 0, sizeof(err));
    duckdb_destroy_value(&path_val);

    if (!path || path[0] == '\0') {
        duckdb_bind_set_error(info, "duckhts_samtools_idxstats requires a non-empty alignment path");
        if (path) duckdb_free(path);
        return;
    }

    bind = (samtools_idxstats_bind_t *)duckdb_malloc(sizeof(samtools_idxstats_bind_t));
    if (!bind) {
        duckdb_bind_set_error(info, "duckhts_samtools_idxstats: out of memory");
        duckdb_free(path);
        return;
    }
    memset(bind, 0, sizeof(*bind));
    bind->path = path;
    bind->threads = 0;
    bind->overwrite = 0;

    val = duckdb_bind_get_named_parameter(info, "output");
    if (val && !duckdb_is_null_value(val)) output_path = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "index_path");
    if (val && !duckdb_is_null_value(val)) index_path = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "threads");
    if (val && !duckdb_is_null_value(val)) bind->threads = (int)duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "overwrite");
    if (val && !duckdb_is_null_value(val)) bind->overwrite = duckdb_get_bool(val) ? 1 : 0;
    if (val) duckdb_destroy_value(&val);

    if (!output_path) output_path = append_suffix_idxstats(bind->path, ".idxstats.txt");
    if (!output_path) {
        destroy_samtools_idxstats_bind(bind);
        if (index_path) duckdb_free(index_path);
        duckdb_bind_set_error(info, "duckhts_samtools_idxstats: out of memory");
        return;
    }

    bind->output_path = output_path;
    bind->index_path = index_path;

    if (bind->threads < 0) {
        destroy_samtools_idxstats_bind(bind);
        duckdb_bind_set_error(info, "duckhts_samtools_idxstats: threads must be >= 0");
        return;
    }
    if (bind->index_path && bind->index_path[0] == '\0') {
        destroy_samtools_idxstats_bind(bind);
        duckdb_bind_set_error(info, "duckhts_samtools_idxstats: index_path must be non-empty when provided");
        return;
    }

    if (run_samtools_idxstats(bind, err, sizeof(err)) != 0) {
        destroy_samtools_idxstats_bind(bind);
        duckdb_bind_set_error(info, err);
        return;
    }

    {
        duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
        duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
        duckdb_bind_add_result_column(info, "success", bool_type);
        duckdb_bind_add_result_column(info, "path", varchar_type);
        duckdb_bind_add_result_column(info, "output_path", varchar_type);
        duckdb_bind_add_result_column(info, "used_index_fast_path", bool_type);
        duckdb_bind_add_result_column(info, "error_message", varchar_type);
        duckdb_destroy_logical_type(&bool_type);
        duckdb_destroy_logical_type(&varchar_type);
    }

    duckdb_bind_set_bind_data(info, bind, destroy_samtools_idxstats_bind);
}

static void samtools_idxstats_init(duckdb_init_info info) {
    samtools_idxstats_bind_t *bind = (samtools_idxstats_bind_t *)duckdb_init_get_bind_data(info);
    bind->emitted = 0;
}

static void samtools_idxstats_scan(duckdb_function_info info, duckdb_data_chunk output) {
    samtools_idxstats_bind_t *bind = (samtools_idxstats_bind_t *)duckdb_function_get_bind_data(info);
    duckdb_vector success_vec;
    duckdb_vector path_vec;
    duckdb_vector output_path_vec;
    duckdb_vector used_fast_vec;
    duckdb_vector error_vec;
    bool *success_data;
    bool *used_fast_data;

    if (bind->emitted) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    success_vec = duckdb_data_chunk_get_vector(output, 0);
    path_vec = duckdb_data_chunk_get_vector(output, 1);
    output_path_vec = duckdb_data_chunk_get_vector(output, 2);
    used_fast_vec = duckdb_data_chunk_get_vector(output, 3);
    error_vec = duckdb_data_chunk_get_vector(output, 4);

    success_data = (bool *)duckdb_vector_get_data(success_vec);
    used_fast_data = (bool *)duckdb_vector_get_data(used_fast_vec);

    success_data[0] = true;
    used_fast_data[0] = bind->used_index_fast_path ? true : false;
    duckdb_vector_assign_string_element(path_vec, 0, bind->path ? bind->path : "");
    duckdb_vector_assign_string_element(output_path_vec, 0, bind->output_path ? bind->output_path : "");
    duckdb_vector_assign_string_element(error_vec, 0, "");

    bind->emitted = 1;
    duckdb_data_chunk_set_size(output, 1);
}

void register_duckhts_samtools_idxstats_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);

    duckdb_table_function_set_name(tf, "duckhts_samtools_idxstats");
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "output", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "index_path", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "threads", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "overwrite", bool_type);
    duckdb_table_function_set_bind(tf, samtools_idxstats_bind);
    duckdb_table_function_set_init(tf, samtools_idxstats_init);
    duckdb_table_function_set_function(tf, samtools_idxstats_scan);
    duckdb_register_table_function(connection, tf);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_table_function(&tf);
}
