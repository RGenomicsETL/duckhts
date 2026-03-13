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

typedef struct {
    char *index_path;
    char *index_format;
    int emitted;
} index_build_bind_t;

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

static void destroy_index_build_bind(void *data) {
    index_build_bind_t *bind = (index_build_bind_t *)data;
    if (!bind) return;
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

static void init_index_build(duckdb_init_info info) {
    index_build_bind_t *bind = (index_build_bind_t *)duckdb_init_get_bind_data(info);
    bind->emitted = 0;
}

static void scan_index_build(duckdb_function_info info, duckdb_data_chunk output) {
    index_build_bind_t *bind = (index_build_bind_t *)duckdb_function_get_bind_data(info);
    duckdb_vector success_vec;
    duckdb_vector path_vec;
    duckdb_vector fmt_vec;
    bool *success_data;

    if (bind->emitted) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    success_vec = duckdb_data_chunk_get_vector(output, 0);
    path_vec = duckdb_data_chunk_get_vector(output, 1);
    fmt_vec = duckdb_data_chunk_get_vector(output, 2);

    success_data = (bool *)duckdb_vector_get_data(success_vec);
    success_data[0] = true;
    duckdb_vector_assign_string_element(path_vec, 0, bind->index_path ? bind->index_path : "");
    duckdb_vector_assign_string_element(fmt_vec, 0, bind->index_format ? bind->index_format : "");
    bind->emitted = 1;
    duckdb_data_chunk_set_size(output, 1);
}

static void bind_bam_index(duckdb_bind_info info) {
    duckdb_value path_val = duckdb_bind_get_parameter(info, 0);
    char *path = duckdb_get_varchar(path_val);
    char *index_path = NULL;
    const char *index_format;
    duckdb_value val;
    int min_shift = 0;
    int threads = 4;
    int ret;
    char err[512];

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

    ret = sam_index_build3(path, index_path, min_shift, threads);
    if (ret != 0) {
        snprintf(err, sizeof(err), "bam_index: failed to build index for %s (error %d)", path, ret);
        duckdb_bind_set_error(info, err);
        duckdb_free(path);
        if (index_path) duckdb_free(index_path);
        return;
    }

    if (!index_path) index_path = default_bam_index_path(path, min_shift);
    index_format = ends_with(path, ".cram") ? "CRAI" : (min_shift > 0 ? "CSI" : "BAI");

    add_result_columns(info);
    index_build_bind_t *bind = (index_build_bind_t *)duckdb_malloc(sizeof(index_build_bind_t));
    bind->index_path = index_path ? index_path : dup_string("");
    bind->index_format = dup_string(index_format);
    bind->emitted = 0;
    duckdb_bind_set_bind_data(info, bind, destroy_index_build_bind);
    duckdb_free(path);
}

static void bind_bcf_index(duckdb_bind_info info) {
    duckdb_value path_val = duckdb_bind_get_parameter(info, 0);
    char *path = duckdb_get_varchar(path_val);
    char *index_path = NULL;
    const char *index_format;
    duckdb_value val;
    int min_shift = 0;
    int threads = 4;
    int ret;
    char err[512];

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

    ret = bcf_index_build3(path, index_path, min_shift, threads);
    if (ret != 0) {
        snprintf(err, sizeof(err), "bcf_index: failed to build index for %s (error %d)", path, ret);
        duckdb_bind_set_error(info, err);
        duckdb_free(path);
        if (index_path) duckdb_free(index_path);
        return;
    }

    if (!index_path) index_path = default_bcf_index_path(path, min_shift);
    index_format = min_shift > 0 ? "CSI" : "TBI";

    add_result_columns(info);
    index_build_bind_t *bind = (index_build_bind_t *)duckdb_malloc(sizeof(index_build_bind_t));
    bind->index_path = index_path ? index_path : dup_string("");
    bind->index_format = dup_string(index_format);
    bind->emitted = 0;
    duckdb_bind_set_bind_data(info, bind, destroy_index_build_bind);
    duckdb_free(path);
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
    const char *index_format;
    duckdb_value val;
    tbx_conf_t conf = tbx_conf_vcf;
    int min_shift = 0;
    int threads = 4;
    int ret;
    char err[512];

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

    ret = tbx_index_build3(path, index_path, min_shift, threads, &conf);
    if (ret != 0) {
        snprintf(err, sizeof(err), "tabix_index: failed to build index for %s (error %d)", path, ret);
        duckdb_bind_set_error(info, err);
        duckdb_free(path);
        if (index_path) duckdb_free(index_path);
        if (preset) duckdb_free(preset);
        return;
    }

    if (!index_path) index_path = default_tabix_index_path(path, min_shift);
    index_format = min_shift > 0 ? "CSI" : "TBI";

    add_result_columns(info);
    index_build_bind_t *bind = (index_build_bind_t *)duckdb_malloc(sizeof(index_build_bind_t));
    bind->index_path = index_path ? index_path : dup_string("");
    bind->index_format = dup_string(index_format);
    bind->emitted = 0;
    duckdb_bind_set_bind_data(info, bind, destroy_index_build_bind);

    duckdb_free(path);
    if (preset) duckdb_free(preset);
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
