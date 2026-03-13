/**
 * DuckHTS BGZF Compression / Decompression
 */

#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <htslib/bgzf.h>

#define BGZIP_IO_BUF_SIZE (64 * 1024)

typedef struct {
    char *path;
    int64_t bytes_in;
    int64_t bytes_out;
    int emitted;
} bgzip_bind_t;

static int ends_with(const char *s, const char *suffix) {
    size_t s_len, suffix_len;
    if (!s || !suffix) return 0;
    s_len = strlen(s);
    suffix_len = strlen(suffix);
    if (suffix_len > s_len) return 0;
    return strcmp(s + s_len - suffix_len, suffix) == 0;
}

static char *default_bgzip_output_path(const char *path) {
    size_t len = strlen(path);
    char *out = (char *)duckdb_malloc(len + 4);
    if (!out) return NULL;
    snprintf(out, len + 4, "%s.gz", path);
    return out;
}

static char *default_bgunzip_output_path(const char *path) {
    size_t len = strlen(path);
    char *out;
    if (ends_with(path, ".gz") && len > 3) {
        out = (char *)duckdb_malloc(len - 2);
        if (!out) return NULL;
        memcpy(out, path, len - 3);
        out[len - 3] = '\0';
        return out;
    }
    out = (char *)duckdb_malloc(len + 5);
    if (!out) return NULL;
    snprintf(out, len + 5, "%s.out", path);
    return out;
}

static int stat_file_size(const char *path, int64_t *size_out) {
    struct stat st;
    if (stat(path, &st) != 0) return -1;
    *size_out = (int64_t)st.st_size;
    return 0;
}

static void destroy_bgzip_bind(void *data) {
    bgzip_bind_t *bind = (bgzip_bind_t *)data;
    if (!bind) return;
    if (bind->path) duckdb_free(bind->path);
    duckdb_free(bind);
}

static void add_bgzip_result_columns(duckdb_bind_info info) {
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_bind_add_result_column(info, "success", bool_type);
    duckdb_bind_add_result_column(info, "output_path", varchar_type);
    duckdb_bind_add_result_column(info, "bytes_in", bigint_type);
    duckdb_bind_add_result_column(info, "bytes_out", bigint_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
}

static void bgzip_bind_common(duckdb_bind_info info, int decompress) {
    duckdb_value path_val = duckdb_bind_get_parameter(info, 0);
    char *input_path = duckdb_get_varchar(path_val);
    char *output_path = NULL;
    duckdb_value val;
    int threads = 4;
    int level = -1;
    int keep = 1;
    int overwrite = 0;
    int64_t bytes_in = 0;
    int64_t bytes_out = 0;
    FILE *plain_fp = NULL;
    BGZF *bgzf_fp = NULL;
    char io_buf[BGZIP_IO_BUF_SIZE];

    duckdb_destroy_value(&path_val);
    if (!input_path || input_path[0] == '\0') {
        duckdb_bind_set_error(info, decompress ? "bgunzip requires a file path" : "bgzip requires a file path");
        if (input_path) duckdb_free(input_path);
        return;
    }

    val = duckdb_bind_get_named_parameter(info, "output_path");
    if (val && !duckdb_is_null_value(val)) output_path = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "threads");
    if (val && !duckdb_is_null_value(val)) threads = (int)duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    if (!decompress) {
        val = duckdb_bind_get_named_parameter(info, "level");
        if (val && !duckdb_is_null_value(val)) level = (int)duckdb_get_int64(val);
        if (val) duckdb_destroy_value(&val);
    }

    val = duckdb_bind_get_named_parameter(info, "keep");
    if (val && !duckdb_is_null_value(val)) keep = duckdb_get_bool(val) ? 1 : 0;
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "overwrite");
    if (val && !duckdb_is_null_value(val)) overwrite = duckdb_get_bool(val) ? 1 : 0;
    if (val) duckdb_destroy_value(&val);

    if (!output_path) {
        output_path = decompress ? default_bgunzip_output_path(input_path) : default_bgzip_output_path(input_path);
        if (!output_path) {
            duckdb_bind_set_error(info, "Out of memory");
            duckdb_free(input_path);
            return;
        }
    }

    if (!overwrite) {
        struct stat st;
        if (stat(output_path, &st) == 0) {
            char err[512];
            snprintf(err, sizeof(err), "%s: output '%s' already exists (use overwrite := TRUE to replace)",
                     decompress ? "bgunzip" : "bgzip", output_path);
            duckdb_bind_set_error(info, err);
            duckdb_free(input_path);
            duckdb_free(output_path);
            return;
        }
    }

    if (decompress) {
        bgzf_fp = bgzf_open(input_path, "r");
        if (!bgzf_fp) {
            char err[512];
            snprintf(err, sizeof(err), "bgunzip: cannot open input %s", input_path);
            duckdb_bind_set_error(info, err);
            duckdb_free(input_path);
            duckdb_free(output_path);
            return;
        }
        if (threads > 1 && bgzf_mt(bgzf_fp, threads, 256) != 0) {
            duckdb_bind_set_error(info, "bgunzip: failed to enable multithreaded BGZF I/O");
            bgzf_close(bgzf_fp);
            duckdb_free(input_path);
            duckdb_free(output_path);
            return;
        }
        plain_fp = fopen(output_path, "wb");
        if (!plain_fp) {
            char err[512];
            snprintf(err, sizeof(err), "bgunzip: cannot open output %s: %s", output_path, strerror(errno));
            duckdb_bind_set_error(info, err);
            bgzf_close(bgzf_fp);
            duckdb_free(input_path);
            duckdb_free(output_path);
            return;
        }
        for (;;) {
            ssize_t n_read = bgzf_read(bgzf_fp, io_buf, sizeof(io_buf));
            if (n_read < 0) {
                duckdb_bind_set_error(info, "bgunzip: read error");
                fclose(plain_fp);
                bgzf_close(bgzf_fp);
                duckdb_free(input_path);
                duckdb_free(output_path);
                return;
            }
            if (n_read == 0) break;
            if (fwrite(io_buf, 1, (size_t)n_read, plain_fp) != (size_t)n_read) {
                duckdb_bind_set_error(info, "bgunzip: write error");
                fclose(plain_fp);
                bgzf_close(bgzf_fp);
                duckdb_free(input_path);
                duckdb_free(output_path);
                return;
            }
            bytes_out += n_read;
        }
        fclose(plain_fp);
        if (bgzf_close(bgzf_fp) != 0) {
            duckdb_bind_set_error(info, "bgunzip: close error");
            duckdb_free(input_path);
            duckdb_free(output_path);
            return;
        }
        if (stat_file_size(input_path, &bytes_in) != 0) bytes_in = 0;
    } else {
        char mode[8];
        plain_fp = fopen(input_path, "rb");
        if (!plain_fp) {
            char err[512];
            snprintf(err, sizeof(err), "bgzip: cannot open input %s: %s", input_path, strerror(errno));
            duckdb_bind_set_error(info, err);
            duckdb_free(input_path);
            duckdb_free(output_path);
            return;
        }
        if (level >= 0 && level <= 9) snprintf(mode, sizeof(mode), "w%d", level);
        else snprintf(mode, sizeof(mode), "w");
        bgzf_fp = bgzf_open(output_path, mode);
        if (!bgzf_fp) {
            char err[512];
            snprintf(err, sizeof(err), "bgzip: cannot open output %s", output_path);
            duckdb_bind_set_error(info, err);
            fclose(plain_fp);
            duckdb_free(input_path);
            duckdb_free(output_path);
            return;
        }
        if (threads > 1 && bgzf_mt(bgzf_fp, threads, 256) != 0) {
            duckdb_bind_set_error(info, "bgzip: failed to enable multithreaded BGZF I/O");
            fclose(plain_fp);
            bgzf_close(bgzf_fp);
            duckdb_free(input_path);
            duckdb_free(output_path);
            return;
        }
        for (;;) {
            size_t n_read = fread(io_buf, 1, sizeof(io_buf), plain_fp);
            if (n_read > 0) {
                if (bgzf_write(bgzf_fp, io_buf, n_read) < 0) {
                    duckdb_bind_set_error(info, "bgzip: write error");
                    fclose(plain_fp);
                    bgzf_close(bgzf_fp);
                    duckdb_free(input_path);
                    duckdb_free(output_path);
                    return;
                }
                bytes_in += (int64_t)n_read;
            }
            if (n_read < sizeof(io_buf)) {
                if (ferror(plain_fp)) {
                    duckdb_bind_set_error(info, "bgzip: read error");
                    fclose(plain_fp);
                    bgzf_close(bgzf_fp);
                    duckdb_free(input_path);
                    duckdb_free(output_path);
                    return;
                }
                break;
            }
        }
        fclose(plain_fp);
        if (bgzf_close(bgzf_fp) != 0) {
            duckdb_bind_set_error(info, "bgzip: close error");
            duckdb_free(input_path);
            duckdb_free(output_path);
            return;
        }
        if (stat_file_size(output_path, &bytes_out) != 0) bytes_out = 0;
    }

    if (!keep) {
        unlink(input_path);
    }

    add_bgzip_result_columns(info);
    bgzip_bind_t *bind = (bgzip_bind_t *)duckdb_malloc(sizeof(bgzip_bind_t));
    bind->path = output_path;
    bind->bytes_in = bytes_in;
    bind->bytes_out = bytes_out;
    bind->emitted = 0;
    duckdb_bind_set_bind_data(info, bind, destroy_bgzip_bind);
    duckdb_free(input_path);
}

static void bgzip_bind(duckdb_bind_info info) {
    bgzip_bind_common(info, 0);
}

static void bgunzip_bind(duckdb_bind_info info) {
    bgzip_bind_common(info, 1);
}

static void bgzip_init(duckdb_init_info info) {
    bgzip_bind_t *bind = (bgzip_bind_t *)duckdb_init_get_bind_data(info);
    bind->emitted = 0;
}

static void bgzip_scan(duckdb_function_info info, duckdb_data_chunk output) {
    bgzip_bind_t *bind = (bgzip_bind_t *)duckdb_function_get_bind_data(info);
    duckdb_vector success_vec;
    duckdb_vector path_vec;
    duckdb_vector bytes_in_vec;
    duckdb_vector bytes_out_vec;
    bool *success_data;
    int64_t *bytes_in_data;
    int64_t *bytes_out_data;

    if (bind->emitted) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    success_vec = duckdb_data_chunk_get_vector(output, 0);
    path_vec = duckdb_data_chunk_get_vector(output, 1);
    bytes_in_vec = duckdb_data_chunk_get_vector(output, 2);
    bytes_out_vec = duckdb_data_chunk_get_vector(output, 3);

    success_data = (bool *)duckdb_vector_get_data(success_vec);
    bytes_in_data = (int64_t *)duckdb_vector_get_data(bytes_in_vec);
    bytes_out_data = (int64_t *)duckdb_vector_get_data(bytes_out_vec);

    success_data[0] = true;
    duckdb_vector_assign_string_element(path_vec, 0, bind->path ? bind->path : "");
    bytes_in_data[0] = bind->bytes_in;
    bytes_out_data[0] = bind->bytes_out;

    bind->emitted = 1;
    duckdb_data_chunk_set_size(output, 1);
}

void register_bgzip_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type int_type = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);

    duckdb_table_function_set_name(tf, "bgzip");
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "output_path", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "threads", int_type);
    duckdb_table_function_add_named_parameter(tf, "level", int_type);
    duckdb_table_function_add_named_parameter(tf, "keep", bool_type);
    duckdb_table_function_add_named_parameter(tf, "overwrite", bool_type);
    duckdb_table_function_set_bind(tf, bgzip_bind);
    duckdb_table_function_set_init(tf, bgzip_init);
    duckdb_table_function_set_function(tf, bgzip_scan);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&int_type);
    duckdb_destroy_logical_type(&bool_type);
}

void register_bgunzip_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type int_type = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);

    duckdb_table_function_set_name(tf, "bgunzip");
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "output_path", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "threads", int_type);
    duckdb_table_function_add_named_parameter(tf, "keep", bool_type);
    duckdb_table_function_add_named_parameter(tf, "overwrite", bool_type);
    duckdb_table_function_set_bind(tf, bgunzip_bind);
    duckdb_table_function_set_init(tf, bgzip_init);
    duckdb_table_function_set_function(tf, bgzip_scan);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&int_type);
    duckdb_destroy_logical_type(&bool_type);
}
