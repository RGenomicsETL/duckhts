#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <stdio.h>
#include <string.h>

#include "include/quality_encoding.h"

typedef struct {
    char *path;
    duckhts_file_access_t *access;
    int64_t max_records;
} quality_detect_bind_t;

typedef struct {
    duckhts_quality_detect_result detect;
    char *compatible_encodings;
    int emitted;
} quality_detect_init_t;

static char *dup_string(const char *s) {
    size_t len;
    char *copy;
    if (!s) return NULL;
    len = strlen(s) + 1;
    copy = (char *)duckdb_malloc(len);
    if (!copy) return NULL;
    memcpy(copy, s, len);
    return copy;
}

static void destroy_quality_detect_bind(void *data) {
    quality_detect_bind_t *bind = (quality_detect_bind_t *)data;
    if (!bind) return;
    if (bind->path) duckdb_free(bind->path);
    duckhts_file_access_destroy(bind->access);
    duckdb_free(bind);
}

static void destroy_quality_detect_init(void *data) {
    quality_detect_init_t *init = (quality_detect_init_t *)data;
    if (!init) return;
    if (init->compatible_encodings) duckdb_free(init->compatible_encodings);
    duckdb_free(init);
}

static void build_compatible_string(const duckhts_quality_detect_result *detect,
                                    char *buf, size_t buf_len) {
    int first = 1;
    buf[0] = '\0';
    if (detect->compatible_phred33) {
        snprintf(buf + strlen(buf), buf_len - strlen(buf), "%sphred33", first ? "" : ",");
        first = 0;
    }
    if (detect->compatible_phred64) {
        snprintf(buf + strlen(buf), buf_len - strlen(buf), "%sphred64", first ? "" : ",");
        first = 0;
    }
    if (detect->compatible_solexa64) {
        snprintf(buf + strlen(buf), buf_len - strlen(buf), "%ssolexa64", first ? "" : ",");
    }
}

static void detect_quality_encoding_bind(duckdb_bind_info info) {
    duckdb_value path_val = duckdb_bind_get_parameter(info, 0);
    duckdb_value max_records_val;
    char *path = duckdb_get_varchar(path_val);
    duckdb_logical_type varchar_type;
    duckdb_logical_type bigint_type;
    duckdb_logical_type bool_type;
    quality_detect_bind_t *bind;
    int64_t max_records = 10000;

    duckdb_destroy_value(&path_val);
    if (!path || path[0] == '\0') {
        duckdb_bind_set_error(info, "detect_quality_encoding requires a file path");
        if (path) duckdb_free(path);
        return;
    }

    max_records_val = duckdb_bind_get_named_parameter(info, "max_records");
    if (max_records_val && !duckdb_is_null_value(max_records_val)) {
        max_records = duckdb_get_int64(max_records_val);
    }
    if (max_records_val) duckdb_destroy_value(&max_records_val);

    bind = (quality_detect_bind_t *)duckdb_malloc(sizeof(quality_detect_bind_t));
    memset(bind, 0, sizeof(*bind));
    bind->path = path;
    bind->max_records = max_records;
    if (!duckhts_file_access_bind(info, "detect_quality_encoding", &bind->access)) {
        destroy_quality_detect_bind(bind);
        return;
    }

    varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_bind_add_result_column(info, "format", varchar_type);
    duckdb_bind_add_result_column(info, "observed_ascii_min", bigint_type);
    duckdb_bind_add_result_column(info, "observed_ascii_max", bigint_type);
    duckdb_bind_add_result_column(info, "records_sampled", bigint_type);
    duckdb_bind_add_result_column(info, "compatible_encodings", varchar_type);
    duckdb_bind_add_result_column(info, "guessed_encoding", varchar_type);
    duckdb_bind_add_result_column(info, "is_ambiguous", bool_type);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&bool_type);

    duckdb_bind_set_bind_data(info, bind, destroy_quality_detect_bind);
}

static void detect_quality_encoding_init(duckdb_init_info info) {
    quality_detect_bind_t *bind = (quality_detect_bind_t *)duckdb_init_get_bind_data(info);
    quality_detect_init_t *init = (quality_detect_init_t *)duckdb_malloc(sizeof(quality_detect_init_t));
    char compat[128];
    char err[512];
    memset(init, 0, sizeof(*init));
    if (duckhts_detect_fastq_quality_encoding(bind->access, bind->path, bind->max_records, &init->detect, err, sizeof(err)) != 0) {
        duckdb_init_set_error(info, err);
        duckdb_free(init);
        return;
    }
    build_compatible_string(&init->detect, compat, sizeof(compat));
    init->compatible_encodings = dup_string(compat);
    init->emitted = 0;
    duckdb_init_set_init_data(info, init, destroy_quality_detect_init);
}

static void detect_quality_encoding_scan(duckdb_function_info info, duckdb_data_chunk output) {
    quality_detect_init_t *bind = (quality_detect_init_t *)duckdb_function_get_init_data(info);
    int64_t *min_data;
    int64_t *max_data;
    int64_t *sampled_data;
    bool *ambig_data;

    if (bind->emitted) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    duckdb_vector_assign_string_element(duckdb_data_chunk_get_vector(output, 0), 0, "fastq");
    min_data = (int64_t *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 1));
    max_data = (int64_t *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 2));
    sampled_data = (int64_t *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 3));
    min_data[0] = bind->detect.observed_ascii_min;
    max_data[0] = bind->detect.observed_ascii_max;
    sampled_data[0] = bind->detect.records_sampled;
    duckdb_vector_assign_string_element(duckdb_data_chunk_get_vector(output, 4), 0,
                                        bind->compatible_encodings ? bind->compatible_encodings : "");
    duckdb_vector_assign_string_element(duckdb_data_chunk_get_vector(output, 5), 0,
                                        duckhts_quality_encoding_name(bind->detect.guessed_encoding));
    ambig_data = (bool *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 6));
    ambig_data[0] = bind->detect.is_ambiguous ? true : false;
    bind->emitted = 1;
    duckdb_data_chunk_set_size(output, 1);
}

void register_detect_quality_encoding_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);

    duckdb_table_function_set_name(tf, "detect_quality_encoding");
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "max_records", bigint_type);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);

    duckdb_table_function_set_bind(tf, detect_quality_encoding_bind);
    duckdb_table_function_set_init(tf, detect_quality_encoding_init);
    duckdb_table_function_set_function(tf, detect_quality_encoding_scan);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}
