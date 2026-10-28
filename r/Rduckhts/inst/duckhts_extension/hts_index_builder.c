/**
 * DuckHTS Index Builders
 */

#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <htslib/sam.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>

#include "include/hts_io_tuning.h"

typedef enum {
    INDEX_BUILD_BAM,
    INDEX_BUILD_BCF,
    INDEX_BUILD_TABIX
} index_build_kind_t;

/* Bind records the request; init builds the index. DuckDB 2.0 binds a prepared
 * statement again at every EXECUTE, so work done at bind would run at PREPARE,
 * EXPLAIN and once more per execution
 * (https://github.com/RGenomicsETL/duckhts/issues/358). */
typedef struct {
    index_build_kind_t kind;
    char *path;
    char *requested_index_path;
    char *index_path;
    char *index_format;
    int min_shift;
    int threads;
    tbx_conf_t tbx_conf;
} index_build_bind_t;

typedef struct {
    int emitted;
} index_build_state_t;

static char *dup_string(const char *s) {
    if (!s) return NULL;
    size_t len = strlen(s) + 1;
    char *copy = (char *)duckdb_malloc(len);
    if (!copy) return NULL;
    memcpy(copy, s, len);
    return copy;
}

static int ends_with(const char *s, const char *suffix) {
    size_t s_len, suffix_len;
    if (!s || !suffix) return 0;
    s_len = strlen(s);
    suffix_len = strlen(suffix);
    if (suffix_len > s_len) return 0;
    return strcmp(s + s_len - suffix_len, suffix) == 0;
}

static char *append_suffix(const char *path, const char *suffix) {
    size_t path_len = strlen(path);
    size_t suffix_len = strlen(suffix);
    char *out = (char *)duckdb_malloc(path_len + suffix_len + 1);
    if (!out) return NULL;
    memcpy(out, path, path_len);
    memcpy(out + path_len, suffix, suffix_len + 1);
    return out;
}

static char *default_bam_index_path(const char *path, int min_shift) {
    if (ends_with(path, ".cram")) return append_suffix(path, ".crai");
    return append_suffix(path, min_shift > 0 ? ".csi" : ".bai");
}

static char *default_bcf_index_path(const char *path, int min_shift) {
    return append_suffix(path, min_shift > 0 ? ".csi" : ".tbi");
}

static char *default_tabix_index_path(const char *path, int min_shift) {
    return append_suffix(path, min_shift > 0 ? ".csi" : ".tbi");
}

static const char *index_build_name(index_build_kind_t kind) {
    switch (kind) {
    case INDEX_BUILD_BAM: return "bam_index";
    case INDEX_BUILD_BCF: return "bcf_index";
    case INDEX_BUILD_TABIX: return "tabix_index";
    }
    return "index";
}

static void destroy_index_build_bind(void *data) {
    index_build_bind_t *bind = (index_build_bind_t *)data;
    if (!bind) return;
    if (bind->path) duckdb_free(bind->path);
    if (bind->requested_index_path) duckdb_free(bind->requested_index_path);
    if (bind->index_path) duckdb_free(bind->index_path);
    if (bind->index_format) duckdb_free(bind->index_format);
    duckdb_free(bind);
}

static void add_result_columns(duckdb_bind_info info) {
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_bind_add_result_column(info, "success", bool_type);
    duckdb_bind_add_result_column(info, "index_path", varchar_type);
    duckdb_bind_add_result_column(info, "index_format", varchar_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_logical_type(&varchar_type);
}

/* Takes ownership of path and requested_index_path; returns NULL on allocation failure. */
static index_build_bind_t *new_index_build_bind(index_build_kind_t kind, char *path,
                                                char *requested_index_path, int min_shift,
                                                int threads) {
    index_build_bind_t *bind = (index_build_bind_t *)duckdb_malloc(sizeof(index_build_bind_t));
    const char *format;
    if (!bind) {
        duckdb_free(path);
        if (requested_index_path) duckdb_free(requested_index_path);
        return NULL;
    }
    memset(bind, 0, sizeof(*bind));
    bind->kind = kind;
    bind->path = path;
    bind->requested_index_path = requested_index_path;
    bind->min_shift = min_shift;
    bind->threads = threads;
    bind->tbx_conf = tbx_conf_vcf;

    if (requested_index_path) {
        bind->index_path = dup_string(requested_index_path);
    } else if (kind == INDEX_BUILD_BAM) {
        bind->index_path = default_bam_index_path(path, min_shift);
    } else if (kind == INDEX_BUILD_BCF) {
        bind->index_path = default_bcf_index_path(path, min_shift);
    } else {
        bind->index_path = default_tabix_index_path(path, min_shift);
    }
    if (kind == INDEX_BUILD_BAM && ends_with(path, ".cram")) {
        format = "CRAI";
    } else if (kind == INDEX_BUILD_BAM) {
        format = min_shift > 0 ? "CSI" : "BAI";
    } else {
        format = min_shift > 0 ? "CSI" : "TBI";
    }
    bind->index_format = dup_string(format);
    if (!bind->index_path || !bind->index_format) {
        destroy_index_build_bind(bind);
        return NULL;
    }
    return bind;
}

static int run_index_build(const index_build_bind_t *bind, char *err, size_t errlen) {
    const char *name = index_build_name(bind->kind);
    int ret;

    if (bind->kind == INDEX_BUILD_BAM) {
        htsFile *fp = hts_open(bind->path, "r");
        if (!fp) {
            snprintf(err, errlen, "%s: failed to open %s", name, bind->path);
            return -1;
        }
        duckhts_apply_remote_hts_tuning(fp, bind->path, DUCKHTS_HTS_IO_PROFILE_STREAMING);
        ret = sam_index_build4(fp, bind->path, bind->requested_index_path, bind->min_shift,
                               bind->threads);
        hts_close(fp);
    } else if (bind->kind == INDEX_BUILD_BCF) {
        ret = bcf_index_build3(bind->path, bind->requested_index_path, bind->min_shift,
                               bind->threads);
    } else {
        tbx_conf_t conf = bind->tbx_conf;
        ret = tbx_index_build3(bind->path, bind->requested_index_path, bind->min_shift,
                               bind->threads, &conf);
    }
    if (ret != 0) {
        snprintf(err, errlen, "%s: failed to build index for %s (error %d)", name, bind->path, ret);
        return -1;
    }
    return 0;
}

static void init_index_build(duckdb_init_info info) {
    const index_build_bind_t *bind = (const index_build_bind_t *)duckdb_init_get_bind_data(info);
    index_build_state_t *state = (index_build_state_t *)duckdb_malloc(sizeof(index_build_state_t));
    char err[512];

    if (!state) {
        duckdb_init_set_error(info, "Out of memory");
        return;
    }
    state->emitted = 0;
    duckdb_init_set_init_data(info, state, duckdb_free);
    duckdb_init_set_max_threads(info, 1);
    if (run_index_build(bind, err, sizeof(err)) != 0) {
        duckdb_init_set_error(info, err);
    }
}

static void scan_index_build(duckdb_function_info info, duckdb_data_chunk output) {
    const index_build_bind_t *bind = (const index_build_bind_t *)duckdb_function_get_bind_data(info);
    index_build_state_t *state = (index_build_state_t *)duckdb_function_get_init_data(info);
    duckdb_vector success_vec;
    duckdb_vector path_vec;
    duckdb_vector fmt_vec;
    bool *success_data;

    if (state->emitted) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    success_vec = duckdb_data_chunk_get_vector(output, 0);
    path_vec = duckdb_data_chunk_get_vector(output, 1);
    fmt_vec = duckdb_data_chunk_get_vector(output, 2);

    success_data = (bool *)duckdb_vector_get_data(success_vec);
    success_data[0] = true;
    duckdb_vector_assign_string_element(path_vec, 0, bind->index_path);
    duckdb_vector_assign_string_element(fmt_vec, 0, bind->index_format);
    state->emitted = 1;
    duckdb_data_chunk_set_size(output, 1);
}

static void bind_bam_index(duckdb_bind_info info) {
    duckdb_value path_val = duckdb_bind_get_parameter(info, 0);
    char *path = duckdb_get_varchar(path_val);
    char *index_path = NULL;
    index_build_bind_t *bind;
    duckdb_value val;
    int min_shift = 0;
    int threads = 4;

    duckdb_destroy_value(&path_val);
    if (!path || path[0] == '\0') {
        duckdb_bind_set_error(info, "bam_index requires a file path");
        if (path) duckdb_free(path);
        return;
    }

    val = duckdb_bind_get_named_parameter(info, "index_path");
    if (val && !duckdb_is_null_value(val)) index_path = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "min_shift");
    if (val && !duckdb_is_null_value(val)) min_shift = (int)duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "threads");
    if (val && !duckdb_is_null_value(val)) threads = (int)duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    bind = new_index_build_bind(INDEX_BUILD_BAM, path, index_path, min_shift, threads);
    if (!bind) {
        duckdb_bind_set_error(info, "bam_index: out of memory");
        return;
    }
    add_result_columns(info);
    duckdb_bind_set_bind_data(info, bind, destroy_index_build_bind);
}

static void bind_bcf_index(duckdb_bind_info info) {
    duckdb_value path_val = duckdb_bind_get_parameter(info, 0);
    char *path = duckdb_get_varchar(path_val);
    char *index_path = NULL;
    index_build_bind_t *bind;
    duckdb_value val;
    int min_shift = 0;
    int threads = 4;

    duckdb_destroy_value(&path_val);
    if (!path || path[0] == '\0') {
        duckdb_bind_set_error(info, "bcf_index requires a file path");
        if (path) duckdb_free(path);
        return;
    }

    min_shift = ends_with(path, ".bcf") ? 14 : 0;

    val = duckdb_bind_get_named_parameter(info, "index_path");
    if (val && !duckdb_is_null_value(val)) index_path = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "min_shift");
    if (val && !duckdb_is_null_value(val)) min_shift = (int)duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "threads");
    if (val && !duckdb_is_null_value(val)) threads = (int)duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    bind = new_index_build_bind(INDEX_BUILD_BCF, path, index_path, min_shift, threads);
    if (!bind) {
        duckdb_bind_set_error(info, "bcf_index: out of memory");
        return;
    }
    add_result_columns(info);
    duckdb_bind_set_bind_data(info, bind, destroy_index_build_bind);
}

static int parse_tbx_preset(const char *preset, tbx_conf_t *conf) {
    if (!preset || strcmp(preset, "vcf") == 0) {
        *conf = tbx_conf_vcf;
        return 0;
    }
    if (strcmp(preset, "bed") == 0) {
        *conf = tbx_conf_bed;
        return 0;
    }
    if (strcmp(preset, "gff") == 0) {
        *conf = tbx_conf_gff;
        return 0;
    }
    if (strcmp(preset, "sam") == 0) {
        *conf = tbx_conf_sam;
        return 0;
    }
    return -1;
}

static void bind_tabix_index(duckdb_bind_info info) {
    duckdb_value path_val = duckdb_bind_get_parameter(info, 0);
    char *path = duckdb_get_varchar(path_val);
    char *index_path = NULL;
    char *preset = NULL;
    index_build_bind_t *bind;
    duckdb_value val;
    tbx_conf_t conf = tbx_conf_vcf;
    int min_shift = 0;
    int threads = 4;

    duckdb_destroy_value(&path_val);
    if (!path || path[0] == '\0') {
        duckdb_bind_set_error(info, "tabix_index requires a file path");
        if (path) duckdb_free(path);
        return;
    }

    val = duckdb_bind_get_named_parameter(info, "index_path");
    if (val && !duckdb_is_null_value(val)) index_path = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "preset");
    if (val && !duckdb_is_null_value(val)) preset = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    if (parse_tbx_preset(preset, &conf) != 0) {
        duckdb_bind_set_error(info, "tabix_index: preset must be one of vcf, bed, gff, sam");
        duckdb_free(path);
        if (index_path) duckdb_free(index_path);
        if (preset) duckdb_free(preset);
        return;
    }

    val = duckdb_bind_get_named_parameter(info, "seq_col");
    if (val && !duckdb_is_null_value(val)) conf.sc = (int32_t)duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "start_col");
    if (val && !duckdb_is_null_value(val)) conf.bc = (int32_t)duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "end_col");
    if (val && !duckdb_is_null_value(val)) conf.ec = (int32_t)duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "comment_char");
    if (val && !duckdb_is_null_value(val)) {
        char *meta = duckdb_get_varchar(val);
        if (meta && meta[0] != '\0') conf.meta_char = meta[0];
        if (meta) duckdb_free(meta);
    }
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "skip_lines");
    if (val && !duckdb_is_null_value(val)) conf.line_skip = (int32_t)duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "min_shift");
    if (val && !duckdb_is_null_value(val)) min_shift = (int)duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "threads");
    if (val && !duckdb_is_null_value(val)) threads = (int)duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    bind = new_index_build_bind(INDEX_BUILD_TABIX, path, index_path, min_shift, threads);
    if (preset) duckdb_free(preset);
    if (!bind) {
        duckdb_bind_set_error(info, "tabix_index: out of memory");
        return;
    }
    bind->tbx_conf = conf;
    add_result_columns(info);
    duckdb_bind_set_bind_data(info, bind, destroy_index_build_bind);
}

void register_bam_index_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type int_type = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);

    duckdb_table_function_set_name(tf, "bam_index");
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "index_path", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "min_shift", int_type);
    duckdb_table_function_add_named_parameter(tf, "threads", int_type);
    duckdb_table_function_set_bind(tf, bind_bam_index);
    duckdb_table_function_set_init(tf, init_index_build);
    duckdb_table_function_set_function(tf, scan_index_build);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&int_type);
}

void register_bcf_index_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type int_type = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);

    duckdb_table_function_set_name(tf, "bcf_index");
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "index_path", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "min_shift", int_type);
    duckdb_table_function_add_named_parameter(tf, "threads", int_type);
    duckdb_table_function_set_bind(tf, bind_bcf_index);
    duckdb_table_function_set_init(tf, init_index_build);
    duckdb_table_function_set_function(tf, scan_index_build);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&int_type);
}

void register_tabix_index_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type int_type = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);

    duckdb_table_function_set_name(tf, "tabix_index");
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "preset", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "index_path", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "min_shift", int_type);
    duckdb_table_function_add_named_parameter(tf, "threads", int_type);
    duckdb_table_function_add_named_parameter(tf, "seq_col", int_type);
    duckdb_table_function_add_named_parameter(tf, "start_col", int_type);
    duckdb_table_function_add_named_parameter(tf, "end_col", int_type);
    duckdb_table_function_add_named_parameter(tf, "comment_char", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "skip_lines", int_type);
    duckdb_table_function_set_bind(tf, bind_tabix_index);
    duckdb_table_function_set_init(tf, init_index_build);
    duckdb_table_function_set_function(tf, scan_index_build);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&int_type);
}
