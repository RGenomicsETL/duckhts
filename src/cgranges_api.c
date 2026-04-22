#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <limits.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cgranges.h"

#define DUCKHTS_CGRANGES_ERRLEN 512
#define DUCKHTS_CGRANGES_OVERLAP_CHUNK 2048

typedef enum {
    DUCKHTS_CGR_LABEL_ORDINAL = 0,
    DUCKHTS_CGR_LABEL_BIGINT,
    DUCKHTS_CGR_LABEL_DOUBLE,
    DUCKHTS_CGR_LABEL_VARCHAR,
    DUCKHTS_CGR_LABEL_BOOLEAN
} duckhts_cgr_label_kind_t;

typedef struct {
    duckhts_cgr_label_kind_t kind;
    size_t count;
    size_t cap;
    char **chroms;
    int32_t *starts;
    int32_t *ends;
    union {
        int64_t *i64;
        double *f64;
        char **str;
        uint8_t *b;
    } labels;
    uint8_t *label_valid;
} duckhts_cgranges_payload_t;

typedef struct duckhts_cgranges_entry {
    char *name;
    bool indexed;
    bool building;
    uint32_t ref_count;
    cgranges_t *cr;
    duckhts_cgranges_payload_t payload;
    struct duckhts_cgranges_entry *next;
} duckhts_cgranges_entry_t;

typedef struct {
    pthread_mutex_t mutex;
    uint32_t ref_count;
    duckdb_database database;
    duckdb_connection connection;
    duckhts_cgranges_entry_t *entries;
} duckhts_cgranges_registry_t;

typedef struct {
    duckhts_cgranges_registry_t *reg;
    char *name;
    char *chrom;
    int64_t start;
    int64_t end;
    char *mode;
    bool has_query_row_id;
    int64_t query_row_id;
    duckhts_cgranges_entry_t *entry;
    int64_t *hits;
    int64_t hits_cap;
    int64_t n_hits;
    int64_t emitted;
} duckhts_cgranges_overlaps_bind_t;

typedef struct {
    duckhts_cgranges_registry_t *reg;
    char *name;
    char *query;
    char *chrom_col;
    char *start_col;
    char *end_col;
    char *query_row_id_col;
    char *mode;
    duckhts_cgranges_entry_t *entry;
} duckhts_cgranges_overlaps_bulk_bind_t;

typedef struct {
    duckdb_result query_result;
    bool has_query_result;
    int64_t chrom_idx;
    int64_t start_idx;
    int64_t end_idx;
    int64_t query_row_id_idx;
    duckdb_type chrom_type;
    duckdb_type start_type;
    duckdb_type end_type;
    duckdb_type query_row_id_type;
    duckdb_data_chunk current_chunk;
    idx_t current_chunk_rows;
    idx_t current_chunk_row;
    int64_t next_generated_query_row_id;
    int64_t active_query_row_id;
    int64_t *hits;
    int64_t hits_cap;
    int64_t n_hits;
    int64_t emitted;
    bool done;
} duckhts_cgranges_overlaps_bulk_init_t;

static char *chunk_cell_strdup(duckdb_data_chunk input, idx_t col, idx_t row);
static int64_t chunk_cell_to_int64(duckdb_vector vec, idx_t row, duckdb_type type_id);
static duckhts_cgr_label_kind_t normalized_label_kind_from_type(duckdb_type type_id);
static int append_interval_raw(duckhts_cgranges_entry_t *entry, const char *chrom, int64_t start, int64_t end,
                               bool label_valid, duckhts_cgr_label_kind_t label_kind,
                               int64_t label_i64, double label_f64, const char *label_str, bool label_bool,
                               char *err, size_t errlen);

static void set_null(duckdb_vector vec, idx_t row) {
    duckdb_vector_ensure_validity_writable(vec);
    uint64_t *validity = duckdb_vector_get_validity(vec);
    validity[row / 64] &= ~((uint64_t)1 << (row % 64));
}

static char *dup_cstr_local(const char *s) {
    size_t n;
    char *out;
    if (!s) return NULL;
    n = strlen(s);
    out = (char *)duckdb_malloc(n + 1);
    if (!out) return NULL;
    memcpy(out, s, n + 1);
    return out;
}

static void free_payload(duckhts_cgranges_payload_t *payload) {
    size_t i;
    if (!payload) return;
    if (payload->chroms) {
        for (i = 0; i < payload->count; i++) {
            if (payload->chroms[i]) duckdb_free(payload->chroms[i]);
        }
        free(payload->chroms);
    }
    free(payload->starts);
    free(payload->ends);
    switch (payload->kind) {
        case DUCKHTS_CGR_LABEL_BIGINT: free(payload->labels.i64); break;
        case DUCKHTS_CGR_LABEL_DOUBLE: free(payload->labels.f64); break;
        case DUCKHTS_CGR_LABEL_BOOLEAN: free(payload->labels.b); break;
        case DUCKHTS_CGR_LABEL_VARCHAR:
            if (payload->labels.str) {
                for (i = 0; i < payload->count; i++) {
                    if (payload->labels.str[i]) duckdb_free(payload->labels.str[i]);
                }
                free(payload->labels.str);
            }
            break;
        default: break;
    }
    free(payload->label_valid);
    memset(payload, 0, sizeof(*payload));
}

static void destroy_entry(duckhts_cgranges_entry_t *entry) {
    if (!entry) return;
    if (entry->name) duckdb_free(entry->name);
    if (entry->cr) cr_destroy(entry->cr);
    free_payload(&entry->payload);
    free(entry);
}

static void free_registry(duckhts_cgranges_registry_t *reg) {
    duckhts_cgranges_entry_t *cur, *next;
    if (!reg) return;
    cur = reg->entries;
    while (cur) {
        next = cur->next;
        destroy_entry(cur);
        cur = next;
    }
    if (reg->connection) duckdb_disconnect(&reg->connection);
    pthread_mutex_destroy(&reg->mutex);
    free(reg);
}

static void destroy_registry(void *ptr) {
    duckhts_cgranges_registry_t *reg = (duckhts_cgranges_registry_t *)ptr;
    bool do_free = false;
    if (!reg) return;
    pthread_mutex_lock(&reg->mutex);
    if (reg->ref_count > 0) {
        reg->ref_count--;
        do_free = (reg->ref_count == 0);
    }
    pthread_mutex_unlock(&reg->mutex);
    if (do_free) free_registry(reg);
}

static void retain_registry(duckhts_cgranges_registry_t *reg) {
    pthread_mutex_lock(&reg->mutex);
    reg->ref_count++;
    pthread_mutex_unlock(&reg->mutex);
}

static duckhts_cgranges_entry_t *lookup_entry_locked(duckhts_cgranges_registry_t *reg, const char *name) {
    duckhts_cgranges_entry_t *cur = reg ? reg->entries : NULL;
    while (cur) {
        if (strcmp(cur->name, name) == 0) return cur;
        cur = cur->next;
    }
    return NULL;
}

static duckhts_cgranges_entry_t *lookup_entry(duckhts_cgranges_registry_t *reg, const char *name) {
    duckhts_cgranges_entry_t *entry;
    pthread_mutex_lock(&reg->mutex);
    entry = lookup_entry_locked(reg, name);
    pthread_mutex_unlock(&reg->mutex);
    return entry;
}

static void unlink_entry_locked(duckhts_cgranges_registry_t *reg, duckhts_cgranges_entry_t *target) {
    duckhts_cgranges_entry_t *cur = reg->entries;
    duckhts_cgranges_entry_t *prev = NULL;
    while (cur) {
        if (cur == target) {
            if (prev) prev->next = cur->next;
            else reg->entries = cur->next;
            cur->next = NULL;
            return;
        }
        prev = cur;
        cur = cur->next;
    }
}

static duckhts_cgranges_entry_t *lookup_and_pin_entry_locked(duckhts_cgranges_registry_t *reg, const char *name) {
    duckhts_cgranges_entry_t *entry = lookup_entry_locked(reg, name);
    if (entry) entry->ref_count++;
    return entry;
}

static void unpin_entry(duckhts_cgranges_registry_t *reg, duckhts_cgranges_entry_t *entry) {
    bool do_destroy = false;
    if (!reg || !entry) return;
    pthread_mutex_lock(&reg->mutex);
    if (entry->ref_count > 0) entry->ref_count--;
    do_destroy = (entry->ref_count == 0 && lookup_entry_locked(reg, entry->name) == NULL);
    pthread_mutex_unlock(&reg->mutex);
    if (do_destroy) destroy_entry(entry);
}

static int ensure_payload_capacity(duckhts_cgranges_payload_t *payload, size_t need, char *err, size_t errlen) {
    size_t old_cap = payload->cap;
    size_t new_cap;
    void *tmp;
    if (need <= payload->cap) return 0;
    new_cap = payload->cap ? payload->cap * 2 : 128;
    while (new_cap < need) new_cap *= 2;

    tmp = realloc(payload->chroms, new_cap * sizeof(char *));
    if (!tmp) goto oom;
    payload->chroms = (char **)tmp;
    tmp = realloc(payload->starts, new_cap * sizeof(int32_t));
    if (!tmp) goto oom;
    payload->starts = (int32_t *)tmp;
    tmp = realloc(payload->ends, new_cap * sizeof(int32_t));
    if (!tmp) goto oom;
    payload->ends = (int32_t *)tmp;
    tmp = realloc(payload->label_valid, new_cap * sizeof(uint8_t));
    if (!tmp) goto oom;
    payload->label_valid = (uint8_t *)tmp;

    switch (payload->kind) {
        case DUCKHTS_CGR_LABEL_BIGINT:
            tmp = realloc(payload->labels.i64, new_cap * sizeof(int64_t));
            if (!tmp) goto oom;
            payload->labels.i64 = (int64_t *)tmp;
            break;
        case DUCKHTS_CGR_LABEL_DOUBLE:
            tmp = realloc(payload->labels.f64, new_cap * sizeof(double));
            if (!tmp) goto oom;
            payload->labels.f64 = (double *)tmp;
            break;
        case DUCKHTS_CGR_LABEL_BOOLEAN:
            tmp = realloc(payload->labels.b, new_cap * sizeof(uint8_t));
            if (!tmp) goto oom;
            payload->labels.b = (uint8_t *)tmp;
            break;
        case DUCKHTS_CGR_LABEL_VARCHAR:
            tmp = realloc(payload->labels.str, new_cap * sizeof(char *));
            if (!tmp) goto oom;
            payload->labels.str = (char **)tmp;
            break;
        default:
            break;
    }
    if (new_cap > old_cap) {
        memset(payload->chroms + old_cap, 0, (new_cap - old_cap) * sizeof(char *));
        memset(payload->starts + old_cap, 0, (new_cap - old_cap) * sizeof(int32_t));
        memset(payload->ends + old_cap, 0, (new_cap - old_cap) * sizeof(int32_t));
        memset(payload->label_valid + old_cap, 0, (new_cap - old_cap) * sizeof(uint8_t));
        switch (payload->kind) {
            case DUCKHTS_CGR_LABEL_BIGINT:
                memset(payload->labels.i64 + old_cap, 0, (new_cap - old_cap) * sizeof(int64_t));
                break;
            case DUCKHTS_CGR_LABEL_DOUBLE:
                memset(payload->labels.f64 + old_cap, 0, (new_cap - old_cap) * sizeof(double));
                break;
            case DUCKHTS_CGR_LABEL_BOOLEAN:
                memset(payload->labels.b + old_cap, 0, (new_cap - old_cap) * sizeof(uint8_t));
                break;
            case DUCKHTS_CGR_LABEL_VARCHAR:
                memset(payload->labels.str + old_cap, 0, (new_cap - old_cap) * sizeof(char *));
                break;
            default:
                break;
        }
    }
    payload->cap = new_cap;
    return 0;
oom:
    snprintf(err, errlen, "duckhts_cgranges: out of memory");
    return -1;
}

static duckhts_cgr_label_kind_t label_kind_from_value(duckdb_value val, char *err, size_t errlen) {
    duckdb_logical_type type;
    duckdb_type type_id;
    if (!val || duckdb_is_null_value(val)) return DUCKHTS_CGR_LABEL_ORDINAL;
    type = duckdb_get_value_type(val);
    type_id = duckdb_get_type_id(type);
    duckdb_destroy_logical_type(&type);
    switch (type_id) {
        case DUCKDB_TYPE_BIGINT:
        case DUCKDB_TYPE_INTEGER:
        case DUCKDB_TYPE_SMALLINT:
        case DUCKDB_TYPE_TINYINT:
        case DUCKDB_TYPE_UBIGINT:
        case DUCKDB_TYPE_UINTEGER:
        case DUCKDB_TYPE_USMALLINT:
        case DUCKDB_TYPE_UTINYINT:
            return DUCKHTS_CGR_LABEL_BIGINT;
        case DUCKDB_TYPE_DOUBLE:
        case DUCKDB_TYPE_FLOAT:
            return DUCKHTS_CGR_LABEL_DOUBLE;
        case DUCKDB_TYPE_VARCHAR:
            return DUCKHTS_CGR_LABEL_VARCHAR;
        case DUCKDB_TYPE_BOOLEAN:
            return DUCKHTS_CGR_LABEL_BOOLEAN;
        default:
            snprintf(err, errlen,
                     "duckhts_cgranges: label type must be BIGINT-like, DOUBLE, VARCHAR, or BOOLEAN");
            return DUCKHTS_CGR_LABEL_ORDINAL;
    }
}

static int append_interval(duckhts_cgranges_entry_t *entry, const char *chrom, int64_t start, int64_t end,
                           duckdb_value label_val, char *err, size_t errlen) {
    size_t idx;
    int32_t ordinal;
    duckhts_cgr_label_kind_t incoming_kind;
    if (!entry || !chrom || !*chrom) {
        snprintf(err, errlen, "duckhts_cgranges_add: chrom must be non-empty");
        return -1;
    }
    if (entry->indexed) {
        snprintf(err, errlen, "duckhts_cgranges_add: index '%s' is already finalized", entry->name);
        return -1;
    }
    if (start < 0 || end < start || end > INT32_MAX) {
        snprintf(err, errlen, "duckhts_cgranges_add: interval out of int32 range");
        return -1;
    }
    if (entry->payload.count > INT32_MAX) {
        snprintf(err, errlen,
                 "duckhts_cgranges_add: too many intervals for cgranges int32 label space");
        return -1;
    }

    incoming_kind = label_kind_from_value(label_val, err, errlen);
    if (label_val && !duckdb_is_null_value(label_val)) {
        if (entry->payload.kind == DUCKHTS_CGR_LABEL_ORDINAL) {
            entry->payload.kind = incoming_kind;
        } else if (incoming_kind != entry->payload.kind) {
            snprintf(err, errlen, "duckhts_cgranges_add: label type mismatch within index '%s'", entry->name);
            return -1;
        }
    }

    if (ensure_payload_capacity(&entry->payload, entry->payload.count + 1, err, errlen) != 0) {
        return -1;
    }
    idx = entry->payload.count;
    ordinal = (int32_t)idx;
    entry->payload.chroms[idx] = dup_cstr_local(chrom);
    if (!entry->payload.chroms[idx]) {
        snprintf(err, errlen, "duckhts_cgranges_add: out of memory");
        return -1;
    }
    entry->payload.starts[idx] = (int32_t)start;
    entry->payload.ends[idx] = (int32_t)end;
    entry->payload.label_valid[idx] = (label_val && !duckdb_is_null_value(label_val)) ? 1 : 0;

    switch (entry->payload.kind) {
        case DUCKHTS_CGR_LABEL_BIGINT:
            entry->payload.labels.i64[idx] = entry->payload.label_valid[idx] ? duckdb_get_int64(label_val) : ordinal;
            break;
        case DUCKHTS_CGR_LABEL_DOUBLE:
            entry->payload.labels.f64[idx] = entry->payload.label_valid[idx] ? duckdb_get_double(label_val) : (double)ordinal;
            break;
        case DUCKHTS_CGR_LABEL_BOOLEAN:
            entry->payload.labels.b[idx] = entry->payload.label_valid[idx] ? (duckdb_get_bool(label_val) ? 1 : 0) : 0;
            break;
        case DUCKHTS_CGR_LABEL_VARCHAR:
            if (entry->payload.label_valid[idx]) {
                entry->payload.labels.str[idx] = duckdb_get_varchar(label_val);
                if (!entry->payload.labels.str[idx]) {
                    snprintf(err, errlen, "duckhts_cgranges_add: out of memory");
                    return -1;
                }
            } else {
                entry->payload.labels.str[idx] = NULL;
            }
            break;
        case DUCKHTS_CGR_LABEL_ORDINAL:
        default:
            break;
    }

    cr_add(entry->cr, chrom, (int32_t)start, (int32_t)end, ordinal);
    entry->payload.count++;
    return 0;
}

static bool row_is_null(duckdb_vector vec, idx_t row) {
    uint64_t *validity = duckdb_vector_get_validity(vec);
    if (!validity) return false;
    return ((validity[row / 64] >> (row % 64)) & 1ULL) == 0;
}

static bool is_stringlike_type(duckdb_type type_id) {
    return type_id == DUCKDB_TYPE_VARCHAR || type_id == DUCKDB_TYPE_STRING_LITERAL;
}

static bool is_numericlike_type(duckdb_type type_id) {
    switch (type_id) {
        case DUCKDB_TYPE_BIGINT:
        case DUCKDB_TYPE_INTEGER_LITERAL:
        case DUCKDB_TYPE_INTEGER:
        case DUCKDB_TYPE_SMALLINT:
        case DUCKDB_TYPE_TINYINT:
        case DUCKDB_TYPE_UBIGINT:
        case DUCKDB_TYPE_UINTEGER:
        case DUCKDB_TYPE_USMALLINT:
        case DUCKDB_TYPE_UTINYINT:
        case DUCKDB_TYPE_DOUBLE:
        case DUCKDB_TYPE_FLOAT:
            return true;
        default:
            return false;
    }
}

static duckdb_value scalar_value_from_chunk(duckdb_data_chunk input, idx_t col, idx_t row) {
    duckdb_vector vec = duckdb_data_chunk_get_vector(input, col);
    duckdb_logical_type type = duckdb_vector_get_column_type(vec);
    duckdb_type type_id = duckdb_get_type_id(type);
    duckdb_value out = NULL;
    if (row_is_null(vec, row)) {
        duckdb_destroy_logical_type(&type);
        return NULL;
    }
    switch (type_id) {
        case DUCKDB_TYPE_VARCHAR:
        case DUCKDB_TYPE_STRING_LITERAL: {
            duckdb_string_t *data = (duckdb_string_t *)duckdb_vector_get_data(vec);
            out = duckdb_create_varchar_length(duckdb_string_t_data(&data[row]), duckdb_string_t_length(data[row]));
            break;
        }
        case DUCKDB_TYPE_BOOLEAN: {
            bool *data = (bool *)duckdb_vector_get_data(vec);
            out = duckdb_create_bool(data[row]);
            break;
        }
        case DUCKDB_TYPE_DOUBLE: {
            double *data = (double *)duckdb_vector_get_data(vec);
            out = duckdb_create_double(data[row]);
            break;
        }
        case DUCKDB_TYPE_FLOAT: {
            float *data = (float *)duckdb_vector_get_data(vec);
            out = duckdb_create_double((double)data[row]);
            break;
        }
        case DUCKDB_TYPE_BIGINT:
        case DUCKDB_TYPE_INTEGER_LITERAL: {
            int64_t *data = (int64_t *)duckdb_vector_get_data(vec);
            out = duckdb_create_int64(data[row]);
            break;
        }
        case DUCKDB_TYPE_INTEGER: {
            int32_t *data = (int32_t *)duckdb_vector_get_data(vec);
            out = duckdb_create_int64((int64_t)data[row]);
            break;
        }
        case DUCKDB_TYPE_SMALLINT: {
            int16_t *data = (int16_t *)duckdb_vector_get_data(vec);
            out = duckdb_create_int64((int64_t)data[row]);
            break;
        }
        case DUCKDB_TYPE_TINYINT: {
            int8_t *data = (int8_t *)duckdb_vector_get_data(vec);
            out = duckdb_create_int64((int64_t)data[row]);
            break;
        }
        case DUCKDB_TYPE_UBIGINT: {
            uint64_t *data = (uint64_t *)duckdb_vector_get_data(vec);
            out = duckdb_create_uint64(data[row]);
            break;
        }
        case DUCKDB_TYPE_UINTEGER: {
            uint32_t *data = (uint32_t *)duckdb_vector_get_data(vec);
            out = duckdb_create_uint64((uint64_t)data[row]);
            break;
        }
        case DUCKDB_TYPE_USMALLINT: {
            uint16_t *data = (uint16_t *)duckdb_vector_get_data(vec);
            out = duckdb_create_uint64((uint64_t)data[row]);
            break;
        }
        case DUCKDB_TYPE_UTINYINT: {
            uint8_t *data = (uint8_t *)duckdb_vector_get_data(vec);
            out = duckdb_create_uint64((uint64_t)data[row]);
            break;
        }
        default:
            break;
    }
    duckdb_destroy_logical_type(&type);
    return out;
}

static void set_bool_output(duckdb_vector output, idx_t row, bool value) {
    bool *data = (bool *)duckdb_vector_get_data(output);
    data[row] = value;
}

static duckhts_cgranges_registry_t *get_registry_from_function(duckdb_function_info info) {
    return (duckhts_cgranges_registry_t *)duckdb_scalar_function_get_extra_info(info);
}

static duckhts_cgranges_registry_t *get_registry_from_bind(duckdb_bind_info info) {
    return (duckhts_cgranges_registry_t *)duckdb_bind_get_extra_info(info);
}

static void cgranges_create_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckhts_cgranges_registry_t *reg = get_registry_from_function(info);
    idx_t n = duckdb_data_chunk_get_size(input), i;
    for (i = 0; i < n; i++) {
        duckdb_value name_val = scalar_value_from_chunk(input, 0, i);
        duckhts_cgranges_entry_t *entry;
        char *name;
        if (!name_val) {
            duckdb_scalar_function_set_error(info, "duckhts_cgranges_create: name must be non-null");
            return;
        }
        name = duckdb_get_varchar(name_val);
        duckdb_destroy_value(&name_val);
        if (!name || !*name) {
            if (name) duckdb_free(name);
            duckdb_scalar_function_set_error(info, "duckhts_cgranges_create: name must be non-empty");
            return;
        }
        pthread_mutex_lock(&reg->mutex);
        if (lookup_entry_locked(reg, name)) {
            pthread_mutex_unlock(&reg->mutex);
            duckdb_free(name);
            duckdb_scalar_function_set_error(info, "duckhts_cgranges_create: index already exists in this session");
            return;
        }
        entry = (duckhts_cgranges_entry_t *)calloc(1, sizeof(*entry));
        if (!entry) {
            pthread_mutex_unlock(&reg->mutex);
            duckdb_free(name);
            duckdb_scalar_function_set_error(info, "duckhts_cgranges_create: out of memory");
            return;
        }
        entry->name = name;
        entry->cr = cr_init();
        if (!entry->cr) {
            pthread_mutex_unlock(&reg->mutex);
            destroy_entry(entry);
            duckdb_scalar_function_set_error(info, "duckhts_cgranges_create: failed to allocate cgranges index");
            return;
        }
        entry->next = reg->entries;
        reg->entries = entry;
        pthread_mutex_unlock(&reg->mutex);
        set_bool_output(output, i, true);
    }
}

static void cgranges_index_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckhts_cgranges_registry_t *reg = get_registry_from_function(info);
    idx_t n = duckdb_data_chunk_get_size(input), i;
    for (i = 0; i < n; i++) {
        duckdb_value name_val = scalar_value_from_chunk(input, 0, i);
        duckhts_cgranges_entry_t *entry;
        char *name;
        if (!name_val) {
            duckdb_scalar_function_set_error(info, "duckhts_cgranges_index: name must be non-null");
            return;
        }
        name = duckdb_get_varchar(name_val);
        duckdb_destroy_value(&name_val);
        pthread_mutex_lock(&reg->mutex);
        entry = lookup_entry_locked(reg, name);
        if (!entry) {
            pthread_mutex_unlock(&reg->mutex);
            duckdb_free(name);
            duckdb_scalar_function_set_error(info, "duckhts_cgranges_index: unknown index name");
            return;
        }
        if (entry->building) {
            pthread_mutex_unlock(&reg->mutex);
            duckdb_free(name);
            duckdb_scalar_function_set_error(info, "duckhts_cgranges_index: index is still being constructed");
            return;
        }
        if (!entry->indexed) {
            cr_index(entry->cr);
            entry->indexed = true;
        }
        pthread_mutex_unlock(&reg->mutex);
        duckdb_free(name);
        set_bool_output(output, i, true);
    }
}

static void cgranges_destroy_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckhts_cgranges_registry_t *reg = get_registry_from_function(info);
    idx_t n = duckdb_data_chunk_get_size(input), i;
    for (i = 0; i < n; i++) {
        duckdb_value name_val = scalar_value_from_chunk(input, 0, i);
        duckhts_cgranges_entry_t *cur;
        char *name;
        if (!name_val) {
            duckdb_scalar_function_set_error(info, "duckhts_cgranges_destroy: name must be non-null");
            return;
        }
        name = duckdb_get_varchar(name_val);
        duckdb_destroy_value(&name_val);
        pthread_mutex_lock(&reg->mutex);
        cur = lookup_entry_locked(reg, name);
        if (!cur) {
            pthread_mutex_unlock(&reg->mutex);
            duckdb_free(name);
            duckdb_scalar_function_set_error(info, "duckhts_cgranges_destroy: unknown index name");
            return;
        }
        if (cur->ref_count > 0) {
            pthread_mutex_unlock(&reg->mutex);
            duckdb_free(name);
            duckdb_scalar_function_set_error(info, "duckhts_cgranges_destroy: index is in use by an active scan");
            return;
        }
        unlink_entry_locked(reg, cur);
        pthread_mutex_unlock(&reg->mutex);
        destroy_entry(cur);
        duckdb_free(name);
        set_bool_output(output, i, true);
    }
}

static void cgranges_add_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckhts_cgranges_registry_t *reg = get_registry_from_function(info);
    idx_t n = duckdb_data_chunk_get_size(input), i;
    for (i = 0; i < n; i++) {
        char *name = chunk_cell_strdup(input, 0, i);
        char *chrom = chunk_cell_strdup(input, 1, i);
        duckdb_vector start_vec = duckdb_data_chunk_get_vector(input, 2);
        duckdb_vector end_vec = duckdb_data_chunk_get_vector(input, 3);
        duckdb_type start_type = duckdb_get_type_id(duckdb_vector_get_column_type(start_vec));
        duckdb_type end_type = duckdb_get_type_id(duckdb_vector_get_column_type(end_vec));
        duckhts_cgranges_entry_t *entry;
        int64_t start, end;
        bool label_valid = false;
        duckhts_cgr_label_kind_t label_kind = DUCKHTS_CGR_LABEL_ORDINAL;
        int64_t label_i64 = 0;
        double label_f64 = 0.0;
        const char *label_str = NULL;
        char *label_owned = NULL;
        bool label_bool = false;
        char err[DUCKHTS_CGRANGES_ERRLEN];
        memset(err, 0, sizeof(err));

        if (!name || !chrom || row_is_null(start_vec, i) || row_is_null(end_vec, i)) {
            if (name) duckdb_free(name);
            if (chrom) duckdb_free(chrom);
            if (label_owned) duckdb_free(label_owned);
            duckdb_scalar_function_set_error(info, "duckhts_cgranges_add: name, chrom, start, and end must be non-null");
            return;
        }
        start = chunk_cell_to_int64(start_vec, i, start_type);
        end = chunk_cell_to_int64(end_vec, i, end_type);

        if (duckdb_data_chunk_get_column_count(input) > 4) {
            duckdb_vector label_vec = duckdb_data_chunk_get_vector(input, 4);
            duckdb_type label_type = duckdb_get_type_id(duckdb_vector_get_column_type(label_vec));
            if (!row_is_null(label_vec, i)) {
                label_valid = true;
                label_kind = normalized_label_kind_from_type(label_type);
                switch (label_kind) {
                    case DUCKHTS_CGR_LABEL_BIGINT:
                        label_i64 = chunk_cell_to_int64(label_vec, i, label_type);
                        break;
                    case DUCKHTS_CGR_LABEL_DOUBLE:
                        if (label_type == DUCKDB_TYPE_FLOAT) {
                            label_f64 = (double)((float *)duckdb_vector_get_data(label_vec))[i];
                        } else {
                            label_f64 = ((double *)duckdb_vector_get_data(label_vec))[i];
                        }
                        break;
                    case DUCKHTS_CGR_LABEL_VARCHAR:
                        label_owned = chunk_cell_strdup(input, 4, i);
                        if (!label_owned) {
                            duckdb_free(name);
                            duckdb_free(chrom);
                            duckdb_scalar_function_set_error(info, "duckhts_cgranges_add: out of memory");
                            return;
                        }
                        label_str = label_owned;
                        break;
                    case DUCKHTS_CGR_LABEL_BOOLEAN:
                        label_bool = ((bool *)duckdb_vector_get_data(label_vec))[i];
                        break;
                    case DUCKHTS_CGR_LABEL_ORDINAL:
                    default:
                        duckdb_free(name);
                        duckdb_free(chrom);
                        if (label_owned) duckdb_free(label_owned);
                        duckdb_scalar_function_set_error(
                            info,
                            "duckhts_cgranges: label type must be BIGINT-like, DOUBLE, VARCHAR, or BOOLEAN"
                        );
                        return;
                }
            }
        }

        pthread_mutex_lock(&reg->mutex);
        entry = lookup_entry_locked(reg, name);
        if (!entry) {
            pthread_mutex_unlock(&reg->mutex);
            duckdb_free(name);
            duckdb_free(chrom);
            if (label_owned) duckdb_free(label_owned);
            duckdb_scalar_function_set_error(info, "duckhts_cgranges_add: unknown index name");
            return;
        }
        if (entry->building) {
            pthread_mutex_unlock(&reg->mutex);
            duckdb_free(name);
            duckdb_free(chrom);
            if (label_owned) duckdb_free(label_owned);
            duckdb_scalar_function_set_error(info, "duckhts_cgranges_add: index is still being constructed");
            return;
        }
        if (append_interval_raw(entry, chrom, start, end, label_valid, label_kind, label_i64, label_f64,
                                label_str, label_bool, err, sizeof(err)) != 0) {
            pthread_mutex_unlock(&reg->mutex);
            duckdb_free(name);
            duckdb_free(chrom);
            if (label_owned) duckdb_free(label_owned);
            duckdb_scalar_function_set_error(info, err);
            return;
        }
        pthread_mutex_unlock(&reg->mutex);
        duckdb_free(name);
        duckdb_free(chrom);
        if (label_owned) duckdb_free(label_owned);
        set_bool_output(output, i, true);
    }
}

static int make_label_value_from_text(const char *label_text, const char *label_type,
                                      duckdb_value *out, char *err, size_t errlen) {
    char *endptr = NULL;
    if (!out) return -1;
    *out = NULL;
    if (!label_type || !label_text) return 0;
    if (strcmp(label_type, "VARCHAR") == 0) {
        *out = duckdb_create_varchar(label_text);
        return *out ? 0 : -1;
    }
    if (strcmp(label_type, "BOOLEAN") == 0) {
        if (strcmp(label_text, "true") == 0 || strcmp(label_text, "TRUE") == 0) {
            *out = duckdb_create_bool(true);
            return 0;
        }
        if (strcmp(label_text, "false") == 0 || strcmp(label_text, "FALSE") == 0) {
            *out = duckdb_create_bool(false);
            return 0;
        }
        snprintf(err, errlen, "duckhts_cgranges: invalid BOOLEAN label value");
        return -1;
    }
    if (strcmp(label_type, "DOUBLE") == 0 || strcmp(label_type, "FLOAT") == 0 ||
        strncmp(label_type, "DECIMAL", 7) == 0) {
        double value = strtod(label_text, &endptr);
        if (!endptr || *endptr != '\0') {
            snprintf(err, errlen, "duckhts_cgranges: invalid DOUBLE label value");
            return -1;
        }
        *out = duckdb_create_double(value);
        return 0;
    }
    if (strcmp(label_type, "BIGINT") == 0 || strcmp(label_type, "INTEGER") == 0 ||
        strcmp(label_type, "SMALLINT") == 0 || strcmp(label_type, "TINYINT") == 0 ||
        strcmp(label_type, "UBIGINT") == 0 || strcmp(label_type, "UINTEGER") == 0 ||
        strcmp(label_type, "USMALLINT") == 0 || strcmp(label_type, "UTINYINT") == 0) {
        long long value = strtoll(label_text, &endptr, 10);
        if (!endptr || *endptr != '\0') {
            snprintf(err, errlen, "duckhts_cgranges: invalid BIGINT label value");
            return -1;
        }
        *out = duckdb_create_int64((int64_t)value);
        return 0;
    }
    snprintf(err, errlen,
             "duckhts_cgranges: label type must be BIGINT-like, DOUBLE, VARCHAR, or BOOLEAN");
    return -1;
}

static int64_t find_result_col(duckdb_result *res, const char *name) {
    idx_t cols = duckdb_column_count(res);
    idx_t i;
    for (i = 0; i < cols; i++) {
        const char *col_name = duckdb_column_name(res, i);
        if (col_name && strcmp(col_name, name) == 0) return (int64_t)i;
    }
    return -1;
}

static duckdb_value chunk_cell_to_value(duckdb_vector vec, idx_t row, duckdb_type type_id) {
    if (row_is_null(vec, row)) return NULL;
    switch (type_id) {
        case DUCKDB_TYPE_VARCHAR:
        case DUCKDB_TYPE_STRING_LITERAL: {
            duckdb_string_t *data = (duckdb_string_t *)duckdb_vector_get_data(vec);
            return duckdb_create_varchar_length(duckdb_string_t_data(&data[row]), duckdb_string_t_length(data[row]));
        }
        case DUCKDB_TYPE_BOOLEAN: {
            bool *data = (bool *)duckdb_vector_get_data(vec);
            return duckdb_create_bool(data[row]);
        }
        case DUCKDB_TYPE_DOUBLE: {
            double *data = (double *)duckdb_vector_get_data(vec);
            return duckdb_create_double(data[row]);
        }
        case DUCKDB_TYPE_FLOAT: {
            float *data = (float *)duckdb_vector_get_data(vec);
            return duckdb_create_double((double)data[row]);
        }
        case DUCKDB_TYPE_BIGINT:
        case DUCKDB_TYPE_INTEGER_LITERAL: {
            int64_t *data = (int64_t *)duckdb_vector_get_data(vec);
            return duckdb_create_int64(data[row]);
        }
        case DUCKDB_TYPE_INTEGER: {
            int32_t *data = (int32_t *)duckdb_vector_get_data(vec);
            return duckdb_create_int64((int64_t)data[row]);
        }
        case DUCKDB_TYPE_SMALLINT: {
            int16_t *data = (int16_t *)duckdb_vector_get_data(vec);
            return duckdb_create_int64((int64_t)data[row]);
        }
        case DUCKDB_TYPE_TINYINT: {
            int8_t *data = (int8_t *)duckdb_vector_get_data(vec);
            return duckdb_create_int64((int64_t)data[row]);
        }
        case DUCKDB_TYPE_UBIGINT: {
            uint64_t *data = (uint64_t *)duckdb_vector_get_data(vec);
            return duckdb_create_uint64(data[row]);
        }
        case DUCKDB_TYPE_UINTEGER: {
            uint32_t *data = (uint32_t *)duckdb_vector_get_data(vec);
            return duckdb_create_uint64((uint64_t)data[row]);
        }
        case DUCKDB_TYPE_USMALLINT: {
            uint16_t *data = (uint16_t *)duckdb_vector_get_data(vec);
            return duckdb_create_uint64((uint64_t)data[row]);
        }
        case DUCKDB_TYPE_UTINYINT: {
            uint8_t *data = (uint8_t *)duckdb_vector_get_data(vec);
            return duckdb_create_uint64((uint64_t)data[row]);
        }
        default:
            return NULL;
    }
}

static int64_t chunk_cell_to_int64(duckdb_vector vec, idx_t row, duckdb_type type_id) {
    switch (type_id) {
        case DUCKDB_TYPE_BIGINT:
        case DUCKDB_TYPE_INTEGER_LITERAL: return ((int64_t *)duckdb_vector_get_data(vec))[row];
        case DUCKDB_TYPE_INTEGER: return ((int32_t *)duckdb_vector_get_data(vec))[row];
        case DUCKDB_TYPE_SMALLINT: return ((int16_t *)duckdb_vector_get_data(vec))[row];
        case DUCKDB_TYPE_TINYINT: return ((int8_t *)duckdb_vector_get_data(vec))[row];
        case DUCKDB_TYPE_UBIGINT: return (int64_t)((uint64_t *)duckdb_vector_get_data(vec))[row];
        case DUCKDB_TYPE_UINTEGER: return (int64_t)((uint32_t *)duckdb_vector_get_data(vec))[row];
        case DUCKDB_TYPE_USMALLINT: return (int64_t)((uint16_t *)duckdb_vector_get_data(vec))[row];
        case DUCKDB_TYPE_UTINYINT: return (int64_t)((uint8_t *)duckdb_vector_get_data(vec))[row];
        case DUCKDB_TYPE_DOUBLE: return (int64_t)((double *)duckdb_vector_get_data(vec))[row];
        case DUCKDB_TYPE_FLOAT: return (int64_t)((float *)duckdb_vector_get_data(vec))[row];
        default: return 0;
    }
}

static duckhts_cgr_label_kind_t normalized_label_kind_from_type(duckdb_type type_id) {
    switch (type_id) {
        case DUCKDB_TYPE_BIGINT:
        case DUCKDB_TYPE_INTEGER_LITERAL:
        case DUCKDB_TYPE_INTEGER:
        case DUCKDB_TYPE_SMALLINT:
        case DUCKDB_TYPE_TINYINT:
        case DUCKDB_TYPE_UBIGINT:
        case DUCKDB_TYPE_UINTEGER:
        case DUCKDB_TYPE_USMALLINT:
        case DUCKDB_TYPE_UTINYINT:
            return DUCKHTS_CGR_LABEL_BIGINT;
        case DUCKDB_TYPE_DOUBLE:
        case DUCKDB_TYPE_FLOAT:
            return DUCKHTS_CGR_LABEL_DOUBLE;
        case DUCKDB_TYPE_VARCHAR:
        case DUCKDB_TYPE_STRING_LITERAL:
            return DUCKHTS_CGR_LABEL_VARCHAR;
        case DUCKDB_TYPE_BOOLEAN:
            return DUCKHTS_CGR_LABEL_BOOLEAN;
        default:
            return DUCKHTS_CGR_LABEL_ORDINAL;
    }
}

static int append_interval_raw(duckhts_cgranges_entry_t *entry, const char *chrom, int64_t start, int64_t end,
                               bool label_valid, duckhts_cgr_label_kind_t label_kind,
                               int64_t label_i64, double label_f64, const char *label_str, bool label_bool,
                               char *err, size_t errlen) {
    size_t idx;
    int32_t ordinal;
    if (!entry || !chrom || !*chrom) {
        snprintf(err, errlen, "duckhts_cgranges_add: chrom must be non-empty");
        return -1;
    }
    if (entry->indexed) {
        snprintf(err, errlen, "duckhts_cgranges_add: index '%s' is already finalized", entry->name);
        return -1;
    }
    if (start < 0 || end < start || end > INT32_MAX) {
        snprintf(err, errlen, "duckhts_cgranges_add: interval out of int32 range");
        return -1;
    }
    if (entry->payload.count > INT32_MAX) {
        snprintf(err, errlen,
                 "duckhts_cgranges_add: too many intervals for cgranges int32 label space");
        return -1;
    }
    if (label_valid) {
        if (entry->payload.kind == DUCKHTS_CGR_LABEL_ORDINAL) {
            entry->payload.kind = label_kind;
        } else if (entry->payload.kind != label_kind) {
            snprintf(err, errlen, "duckhts_cgranges_add: label type mismatch within index '%s'", entry->name);
            return -1;
        }
    }
    if (ensure_payload_capacity(&entry->payload, entry->payload.count + 1, err, errlen) != 0) {
        return -1;
    }
    idx = entry->payload.count;
    ordinal = (int32_t)idx;
    entry->payload.chroms[idx] = dup_cstr_local(chrom);
    if (!entry->payload.chroms[idx]) {
        snprintf(err, errlen, "duckhts_cgranges_add: out of memory");
        return -1;
    }
    entry->payload.starts[idx] = (int32_t)start;
    entry->payload.ends[idx] = (int32_t)end;
    entry->payload.label_valid[idx] = label_valid ? 1 : 0;

    switch (entry->payload.kind) {
        case DUCKHTS_CGR_LABEL_BIGINT:
            entry->payload.labels.i64[idx] = label_valid ? label_i64 : ordinal;
            break;
        case DUCKHTS_CGR_LABEL_DOUBLE:
            entry->payload.labels.f64[idx] = label_valid ? label_f64 : (double)ordinal;
            break;
        case DUCKHTS_CGR_LABEL_BOOLEAN:
            entry->payload.labels.b[idx] = label_valid ? (label_bool ? 1 : 0) : 0;
            break;
        case DUCKHTS_CGR_LABEL_VARCHAR:
            if (label_valid) {
                entry->payload.labels.str[idx] = dup_cstr_local(label_str);
                if (!entry->payload.labels.str[idx]) {
                    snprintf(err, errlen, "duckhts_cgranges_add: out of memory");
                    return -1;
                }
            } else {
                entry->payload.labels.str[idx] = NULL;
            }
            break;
        case DUCKHTS_CGR_LABEL_ORDINAL:
        default:
            break;
    }
    cr_add(entry->cr, chrom, (int32_t)start, (int32_t)end, ordinal);
    entry->payload.count++;
    return 0;
}

static int build_index_from_result(duckhts_cgranges_entry_t *entry, duckdb_result *res,
                                   const char *chrom_col, const char *start_col, const char *end_col,
                                   const char *label_col, char *err, size_t errlen) {
    int64_t chrom_idx = find_result_col(res, chrom_col);
    int64_t start_idx = find_result_col(res, start_col);
    int64_t end_idx = find_result_col(res, end_col);
    int64_t label_idx = (label_col && *label_col) ? find_result_col(res, label_col) : -1;
    duckdb_type label_type = label_idx >= 0 ? duckdb_column_type(res, (idx_t)label_idx) : DUCKDB_TYPE_INVALID;
    duckdb_data_chunk chunk;
    idx_t row_offset = 0;

    if (chrom_idx < 0 || start_idx < 0 || end_idx < 0) {
        snprintf(err, errlen, "duckhts_cgranges_from_query: required columns not found in query result");
        return -1;
    }
    if (label_idx >= 0 && normalized_label_kind_from_type(label_type) == DUCKHTS_CGR_LABEL_ORDINAL) {
        snprintf(err, errlen, "duckhts_cgranges_from_query: label column type is unsupported");
        return -1;
    }

    while ((chunk = duckdb_fetch_chunk(*res)) != NULL) {
        idx_t rows = duckdb_data_chunk_get_size(chunk);
        duckdb_vector chrom_vec = duckdb_data_chunk_get_vector(chunk, (idx_t)chrom_idx);
        duckdb_vector start_vec = duckdb_data_chunk_get_vector(chunk, (idx_t)start_idx);
        duckdb_vector end_vec = duckdb_data_chunk_get_vector(chunk, (idx_t)end_idx);
        duckdb_vector label_vec = (label_idx >= 0) ? duckdb_data_chunk_get_vector(chunk, (idx_t)label_idx) : NULL;
        duckdb_type start_type = duckdb_get_type_id(duckdb_vector_get_column_type(start_vec));
        duckdb_type end_type = duckdb_get_type_id(duckdb_vector_get_column_type(end_vec));
        duckdb_string_t *chrom_data = (duckdb_string_t *)duckdb_vector_get_data(chrom_vec);
        duckdb_string_t *label_str_data = NULL;
        idx_t r;

        if (label_vec && normalized_label_kind_from_type(label_type) == DUCKHTS_CGR_LABEL_VARCHAR) {
            label_str_data = (duckdb_string_t *)duckdb_vector_get_data(label_vec);
        }

        for (r = 0; r < rows; r++) {
            const char *chrom = NULL;
            int64_t start;
            int64_t end;
            bool label_valid = false;
            duckhts_cgr_label_kind_t label_kind = DUCKHTS_CGR_LABEL_ORDINAL;
            int64_t label_i64 = 0;
            double label_f64 = 0.0;
            const char *label_str = NULL;
            bool label_bool = false;

            if (row_is_null(chrom_vec, r) || row_is_null(start_vec, r) || row_is_null(end_vec, r)) {
                duckdb_destroy_data_chunk(&chunk);
                snprintf(err, errlen,
                         "duckhts_cgranges_from_query: null chrom/start/end encountered at row %llu",
                         (unsigned long long)(row_offset + r + 1ULL));
                return -1;
            }

            chrom = duckdb_string_t_data(&chrom_data[r]);
            start = chunk_cell_to_int64(start_vec, r, start_type);
            end = chunk_cell_to_int64(end_vec, r, end_type);

            if (label_vec && !row_is_null(label_vec, r)) {
                label_valid = true;
                label_kind = normalized_label_kind_from_type(label_type);
                switch (label_kind) {
                    case DUCKHTS_CGR_LABEL_BIGINT:
                        label_i64 = chunk_cell_to_int64(label_vec, r, label_type);
                        break;
                    case DUCKHTS_CGR_LABEL_DOUBLE:
                        if (label_type == DUCKDB_TYPE_FLOAT) {
                            label_f64 = (double)((float *)duckdb_vector_get_data(label_vec))[r];
                        } else {
                            label_f64 = ((double *)duckdb_vector_get_data(label_vec))[r];
                        }
                        break;
                    case DUCKHTS_CGR_LABEL_BOOLEAN:
                        label_bool = ((bool *)duckdb_vector_get_data(label_vec))[r];
                        break;
                    case DUCKHTS_CGR_LABEL_VARCHAR:
                        label_str = duckdb_string_t_data(&label_str_data[r]);
                        break;
                    case DUCKHTS_CGR_LABEL_ORDINAL:
                    default:
                        duckdb_destroy_data_chunk(&chunk);
                        snprintf(err, errlen, "duckhts_cgranges_from_query: label column type is unsupported");
                        return -1;
                }
            }

            if (append_interval_raw(entry, chrom, start, end, label_valid, label_kind, label_i64, label_f64,
                                    label_str, label_bool, err, errlen) != 0) {
                duckdb_destroy_data_chunk(&chunk);
                return -1;
            }
        }
        row_offset += rows;
        duckdb_destroy_data_chunk(&chunk);
    }
    return 0;
}

static int build_index_from_sql(duckhts_cgranges_registry_t *reg, const char *name, const char *query,
                                const char *chrom_col, const char *start_col, const char *end_col,
                                const char *label_col, char *err, size_t errlen) {
    duckhts_cgranges_entry_t *entry = NULL;
    duckdb_result res;
    duckdb_state state;

    pthread_mutex_lock(&reg->mutex);
    if (lookup_entry_locked(reg, name)) {
        pthread_mutex_unlock(&reg->mutex);
        snprintf(err, errlen, "duckhts_cgranges_from_query: index already exists in this session");
        return -1;
    }
    entry = (duckhts_cgranges_entry_t *)calloc(1, sizeof(*entry));
    if (!entry) {
        pthread_mutex_unlock(&reg->mutex);
        snprintf(err, errlen, "duckhts_cgranges_from_query: out of memory");
        return -1;
    }
    entry->name = dup_cstr_local(name);
    entry->cr = cr_init();
    entry->building = true;
    if (!entry->name || !entry->cr) {
        pthread_mutex_unlock(&reg->mutex);
        destroy_entry(entry);
        snprintf(err, errlen, "duckhts_cgranges_from_query: out of memory");
        return -1;
    }
    entry->next = reg->entries;
    reg->entries = entry;
    pthread_mutex_unlock(&reg->mutex);

    if (!reg->connection) {
        pthread_mutex_lock(&reg->mutex);
        unlink_entry_locked(reg, entry);
        pthread_mutex_unlock(&reg->mutex);
        destroy_entry(entry);
        snprintf(err, errlen, "duckhts_cgranges_from_query: query connection is unavailable");
        return -1;
    }
    pthread_mutex_lock(&reg->mutex);
    state = duckdb_query(reg->connection, query, &res);
    if (state != DuckDBSuccess) {
        unlink_entry_locked(reg, entry);
        pthread_mutex_unlock(&reg->mutex);
        snprintf(err, errlen, "duckhts_cgranges_from_query: %s", duckdb_result_error(&res));
        destroy_entry(entry);
        return -1;
    }
    if (build_index_from_result(entry, &res, chrom_col, start_col, end_col, label_col, err, errlen) != 0) {
        duckdb_destroy_result(&res);
        unlink_entry_locked(reg, entry);
        pthread_mutex_unlock(&reg->mutex);
        destroy_entry(entry);
        return -1;
    }
    duckdb_destroy_result(&res);
    cr_index(entry->cr);
    entry->indexed = true;
    entry->building = false;
    pthread_mutex_unlock(&reg->mutex);
    return 0;
}

static char *chunk_cell_strdup(duckdb_data_chunk input, idx_t col, idx_t row) {
    duckdb_vector vec = duckdb_data_chunk_get_vector(input, col);
    duckdb_logical_type type = duckdb_vector_get_column_type(vec);
    duckdb_type type_id = duckdb_get_type_id(type);
    char *out = NULL;
    if (row_is_null(vec, row)) {
        duckdb_destroy_logical_type(&type);
        return NULL;
    }
    if (type_id == DUCKDB_TYPE_VARCHAR || type_id == DUCKDB_TYPE_STRING_LITERAL) {
        duckdb_string_t *data = (duckdb_string_t *)duckdb_vector_get_data(vec);
        uint32_t len = duckdb_string_t_length(data[row]);
        out = (char *)duckdb_malloc((size_t)len + 1);
        if (out) {
            memcpy(out, duckdb_string_t_data(&data[row]), len);
            out[len] = '\0';
        }
    }
    duckdb_destroy_logical_type(&type);
    return out;
}

static void cgranges_from_query_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckhts_cgranges_registry_t *reg = get_registry_from_function(info);
    idx_t n = duckdb_data_chunk_get_size(input), i;
    for (i = 0; i < n; i++) {
        char *name = chunk_cell_strdup(input, 0, i);
        char *query = chunk_cell_strdup(input, 1, i);
        char *chrom_col = chunk_cell_strdup(input, 2, i);
        char *start_col = chunk_cell_strdup(input, 3, i);
        char *end_col = chunk_cell_strdup(input, 4, i);
        char *label_col = duckdb_data_chunk_get_column_count(input) > 5 ? chunk_cell_strdup(input, 5, i) : NULL;
        char err[DUCKHTS_CGRANGES_ERRLEN];
        memset(err, 0, sizeof(err));
        if (!name || !query || !chrom_col || !start_col || !end_col) {
            if (name) duckdb_free(name);
            if (query) duckdb_free(query);
            if (chrom_col) duckdb_free(chrom_col);
            if (start_col) duckdb_free(start_col);
            if (end_col) duckdb_free(end_col);
            if (label_col) duckdb_free(label_col);
            duckdb_scalar_function_set_error(info,
                "duckhts_cgranges_from_query: name, query, and column names must be non-null VARCHAR arguments");
            return;
        }

        if (build_index_from_sql(reg, name, query, chrom_col, start_col, end_col, label_col, err, sizeof(err)) != 0) {
            duckdb_free(name);
            duckdb_free(query);
            duckdb_free(chrom_col);
            duckdb_free(start_col);
            duckdb_free(end_col);
            if (label_col) duckdb_free(label_col);
            duckdb_scalar_function_set_error(info, err);
            return;
        }
        duckdb_free(name);
        duckdb_free(query);
        duckdb_free(chrom_col);
        duckdb_free(start_col);
        duckdb_free(end_col);
        if (label_col) duckdb_free(label_col);
        set_bool_output(output, i, true);
    }
}

static void cgranges_from_table_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    (void)input;
    (void)output;
    duckdb_scalar_function_set_error(info,
                                     "duckhts_cgranges_from_table: not implemented; use duckhts_cgranges_from_query instead");
}

static void add_overlaps_result_columns(duckdb_bind_info info, duckhts_cgranges_entry_t *entry) {
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type int_type = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_logical_type label_type = NULL;

    duckdb_bind_add_result_column(info, "query_row_id", bigint_type);
    duckdb_bind_add_result_column(info, "interval_ordinal", bigint_type);
    switch (entry->payload.kind) {
        case DUCKHTS_CGR_LABEL_DOUBLE:
            label_type = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
            break;
        case DUCKHTS_CGR_LABEL_VARCHAR:
            label_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
            break;
        case DUCKHTS_CGR_LABEL_BOOLEAN:
            label_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
            break;
        case DUCKHTS_CGR_LABEL_BIGINT:
        case DUCKHTS_CGR_LABEL_ORDINAL:
        default:
            label_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
            break;
    }
    duckdb_bind_add_result_column(info, "label", label_type);
    duckdb_bind_add_result_column(info, "interval_chrom", varchar_type);
    duckdb_bind_add_result_column(info, "interval_start", int_type);
    duckdb_bind_add_result_column(info, "interval_end", int_type);

    if (label_type) duckdb_destroy_logical_type(&label_type);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&int_type);
}

static void destroy_overlaps_bind(void *ptr) {
    duckhts_cgranges_overlaps_bind_t *bind = (duckhts_cgranges_overlaps_bind_t *)ptr;
    if (!bind) return;
    if (bind->entry) unpin_entry(bind->reg, bind->entry);
    if (bind->name) duckdb_free(bind->name);
    if (bind->chrom) duckdb_free(bind->chrom);
    if (bind->mode) duckdb_free(bind->mode);
    free(bind->hits);
    free(bind);
}

static void overlaps_bind(duckdb_bind_info info) {
    duckhts_cgranges_registry_t *reg = get_registry_from_bind(info);
    duckhts_cgranges_overlaps_bind_t *bind = NULL;
    bool ok = false;
    duckdb_value name_val = duckdb_bind_get_parameter(info, 0);
    duckdb_value chrom_val = duckdb_bind_get_parameter(info, 1);
    duckdb_value start_val = duckdb_bind_get_parameter(info, 2);
    duckdb_value end_val = duckdb_bind_get_parameter(info, 3);
    duckdb_value mode_val;
    duckdb_value qid_val;

    bind = (duckhts_cgranges_overlaps_bind_t *)calloc(1, sizeof(*bind));
    if (!bind) {
        duckdb_bind_set_error(info, "duckhts_cgranges_overlaps: out of memory");
        goto cleanup;
    }
    bind->reg = reg;
    bind->name = duckdb_get_varchar(name_val);
    bind->chrom = duckdb_get_varchar(chrom_val);
    bind->start = duckdb_get_int64(start_val);
    bind->end = duckdb_get_int64(end_val);
    bind->mode = dup_cstr_local("overlap");

    mode_val = duckdb_bind_get_named_parameter(info, "mode");
    if (mode_val && !duckdb_is_null_value(mode_val)) {
        if (bind->mode) duckdb_free(bind->mode);
        bind->mode = duckdb_get_varchar(mode_val);
    }
    if (mode_val) duckdb_destroy_value(&mode_val);

    qid_val = duckdb_bind_get_named_parameter(info, "query_row_id");
    if (qid_val && !duckdb_is_null_value(qid_val)) {
        bind->has_query_row_id = true;
        bind->query_row_id = duckdb_get_int64(qid_val);
    }
    if (qid_val) duckdb_destroy_value(&qid_val);

    if (!bind->name || !bind->chrom || bind->start < 0 || bind->end < bind->start || bind->end > INT32_MAX) {
        duckdb_bind_set_error(info, "duckhts_cgranges_overlaps: invalid name/chrom/start/end");
        goto cleanup;
    }
    if (strcmp(bind->mode, "overlap") != 0 && strcmp(bind->mode, "contain") != 0) {
        duckdb_bind_set_error(info, "duckhts_cgranges_overlaps: mode must be 'overlap' or 'contain'");
        goto cleanup;
    }
    pthread_mutex_lock(&reg->mutex);
    bind->entry = lookup_and_pin_entry_locked(reg, bind->name);
    if (!bind->entry) {
        pthread_mutex_unlock(&reg->mutex);
        duckdb_bind_set_error(info, "duckhts_cgranges_overlaps: unknown index name");
        goto cleanup;
    }
    if (bind->entry->building) {
        pthread_mutex_unlock(&reg->mutex);
        duckdb_bind_set_error(info, "duckhts_cgranges_overlaps: index is still being constructed");
        goto cleanup;
    }
    if (!bind->entry->indexed) {
        pthread_mutex_unlock(&reg->mutex);
        duckdb_bind_set_error(info, "duckhts_cgranges_overlaps: index is not finalized; call duckhts_cgranges_index first");
        goto cleanup;
    }
    pthread_mutex_unlock(&reg->mutex);
    bind->n_hits = strcmp(bind->mode, "contain") == 0
        ? cr_contain(bind->entry->cr, bind->chrom, (int32_t)bind->start, (int32_t)bind->end,
                     &bind->hits, &bind->hits_cap)
        : cr_overlap(bind->entry->cr, bind->chrom, (int32_t)bind->start, (int32_t)bind->end,
                     &bind->hits, &bind->hits_cap);
    add_overlaps_result_columns(info, bind->entry);
    duckdb_bind_set_bind_data(info, bind, destroy_overlaps_bind);
    ok = true;
cleanup:
    if (name_val) duckdb_destroy_value(&name_val);
    if (chrom_val) duckdb_destroy_value(&chrom_val);
    if (start_val) duckdb_destroy_value(&start_val);
    if (end_val) duckdb_destroy_value(&end_val);
    if (bind && !ok) destroy_overlaps_bind(bind);
}

static void overlaps_init(duckdb_init_info info) {
    duckhts_cgranges_overlaps_bind_t *bind = (duckhts_cgranges_overlaps_bind_t *)duckdb_init_get_bind_data(info);
    bind->emitted = 0;
}

static void overlaps_scan(duckdb_function_info info, duckdb_data_chunk output) {
    duckhts_cgranges_overlaps_bind_t *bind = (duckhts_cgranges_overlaps_bind_t *)duckdb_function_get_bind_data(info);
    duckdb_vector qid_vec = duckdb_data_chunk_get_vector(output, 0);
    duckdb_vector ord_vec = duckdb_data_chunk_get_vector(output, 1);
    duckdb_vector label_vec = duckdb_data_chunk_get_vector(output, 2);
    duckdb_vector chrom_vec = duckdb_data_chunk_get_vector(output, 3);
    duckdb_vector start_vec = duckdb_data_chunk_get_vector(output, 4);
    duckdb_vector end_vec = duckdb_data_chunk_get_vector(output, 5);
    int64_t *qid_data = (int64_t *)duckdb_vector_get_data(qid_vec);
    int64_t *ord_data = (int64_t *)duckdb_vector_get_data(ord_vec);
    int32_t *start_data = (int32_t *)duckdb_vector_get_data(start_vec);
    int32_t *end_data = (int32_t *)duckdb_vector_get_data(end_vec);
    idx_t out_idx = 0;
    idx_t max_rows = duckdb_vector_size();

    while (bind->emitted < bind->n_hits && out_idx < max_rows) {
        int64_t ridx = bind->hits[bind->emitted];
        int32_t ordinal = cr_label(bind->entry->cr, ridx);
        qid_data[out_idx] = bind->query_row_id;
        if (!bind->has_query_row_id) set_null(qid_vec, out_idx);
        ord_data[out_idx] = ordinal;
        switch (bind->entry->payload.kind) {
            case DUCKHTS_CGR_LABEL_BIGINT: {
                int64_t *p = (int64_t *)duckdb_vector_get_data(label_vec);
                p[out_idx] = bind->entry->payload.labels.i64[ordinal];
                if (!bind->entry->payload.label_valid[ordinal]) set_null(label_vec, out_idx);
                break;
            }
            case DUCKHTS_CGR_LABEL_DOUBLE: {
                double *p = (double *)duckdb_vector_get_data(label_vec);
                p[out_idx] = bind->entry->payload.labels.f64[ordinal];
                if (!bind->entry->payload.label_valid[ordinal]) set_null(label_vec, out_idx);
                break;
            }
            case DUCKHTS_CGR_LABEL_VARCHAR:
                if (bind->entry->payload.label_valid[ordinal]) {
                    duckdb_vector_assign_string_element(label_vec, out_idx,
                                                        bind->entry->payload.labels.str[ordinal]);
                } else {
                    set_null(label_vec, out_idx);
                }
                break;
            case DUCKHTS_CGR_LABEL_BOOLEAN: {
                bool *p = (bool *)duckdb_vector_get_data(label_vec);
                p[out_idx] = bind->entry->payload.labels.b[ordinal] ? true : false;
                if (!bind->entry->payload.label_valid[ordinal]) set_null(label_vec, out_idx);
                break;
            }
            case DUCKHTS_CGR_LABEL_ORDINAL:
            default: {
                int64_t *p = (int64_t *)duckdb_vector_get_data(label_vec);
                p[out_idx] = ordinal;
                break;
            }
        }
        duckdb_vector_assign_string_element(chrom_vec, out_idx, bind->entry->payload.chroms[ordinal]);
        start_data[out_idx] = bind->entry->payload.starts[ordinal];
        end_data[out_idx] = bind->entry->payload.ends[ordinal];
        bind->emitted++;
        out_idx++;
    }
    duckdb_data_chunk_set_size(output, out_idx);
}

static void destroy_overlaps_bulk_bind(void *ptr) {
    duckhts_cgranges_overlaps_bulk_bind_t *bind = (duckhts_cgranges_overlaps_bulk_bind_t *)ptr;
    if (!bind) return;
    if (bind->entry) unpin_entry(bind->reg, bind->entry);
    if (bind->name) duckdb_free(bind->name);
    if (bind->query) duckdb_free(bind->query);
    if (bind->chrom_col) duckdb_free(bind->chrom_col);
    if (bind->start_col) duckdb_free(bind->start_col);
    if (bind->end_col) duckdb_free(bind->end_col);
    if (bind->query_row_id_col) duckdb_free(bind->query_row_id_col);
    if (bind->mode) duckdb_free(bind->mode);
    free(bind);
}

static void destroy_overlaps_bulk_init(void *ptr) {
    duckhts_cgranges_overlaps_bulk_init_t *init = (duckhts_cgranges_overlaps_bulk_init_t *)ptr;
    if (!init) return;
    if (init->current_chunk) duckdb_destroy_data_chunk(&init->current_chunk);
    if (init->has_query_result) duckdb_destroy_result(&init->query_result);
    free(init->hits);
    free(init);
}

static void overlaps_bulk_bind(duckdb_bind_info info) {
    duckhts_cgranges_registry_t *reg = get_registry_from_bind(info);
    duckhts_cgranges_overlaps_bulk_bind_t *bind = NULL;
    bool ok = false;
    duckdb_value name_val = duckdb_bind_get_parameter(info, 0);
    duckdb_value query_val = duckdb_bind_get_parameter(info, 1);
    duckdb_value chrom_col_val = duckdb_bind_get_parameter(info, 2);
    duckdb_value start_col_val = duckdb_bind_get_parameter(info, 3);
    duckdb_value end_col_val = duckdb_bind_get_parameter(info, 4);
    duckdb_value mode_val = NULL;
    duckdb_value qid_col_val = NULL;

    bind = (duckhts_cgranges_overlaps_bulk_bind_t *)calloc(1, sizeof(*bind));
    if (!bind) {
        duckdb_bind_set_error(info, "duckhts_cgranges_overlaps_bulk: out of memory");
        goto cleanup;
    }
    bind->reg = reg;
    bind->name = duckdb_get_varchar(name_val);
    bind->query = duckdb_get_varchar(query_val);
    bind->chrom_col = duckdb_get_varchar(chrom_col_val);
    bind->start_col = duckdb_get_varchar(start_col_val);
    bind->end_col = duckdb_get_varchar(end_col_val);
    bind->mode = dup_cstr_local("overlap");

    mode_val = duckdb_bind_get_named_parameter(info, "mode");
    if (mode_val && !duckdb_is_null_value(mode_val)) {
        if (bind->mode) duckdb_free(bind->mode);
        bind->mode = duckdb_get_varchar(mode_val);
    }
    qid_col_val = duckdb_bind_get_named_parameter(info, "query_row_id_col");
    if (qid_col_val && !duckdb_is_null_value(qid_col_val)) {
        bind->query_row_id_col = duckdb_get_varchar(qid_col_val);
    }

    if (!bind->name || !bind->query || !bind->chrom_col || !bind->start_col || !bind->end_col ||
        !bind->mode || !*bind->name || !*bind->query || !*bind->chrom_col || !*bind->start_col ||
        !*bind->end_col || (bind->query_row_id_col && !*bind->query_row_id_col)) {
        duckdb_bind_set_error(
            info,
            "duckhts_cgranges_overlaps_bulk: name, query, chrom_col, start_col, and end_col must be non-empty"
        );
        goto cleanup;
    }
    if (strcmp(bind->mode, "overlap") != 0 && strcmp(bind->mode, "contain") != 0) {
        duckdb_bind_set_error(info, "duckhts_cgranges_overlaps_bulk: mode must be 'overlap' or 'contain'");
        goto cleanup;
    }

    pthread_mutex_lock(&reg->mutex);
    bind->entry = lookup_and_pin_entry_locked(reg, bind->name);
    if (!bind->entry) {
        pthread_mutex_unlock(&reg->mutex);
        duckdb_bind_set_error(info, "duckhts_cgranges_overlaps_bulk: unknown index name");
        goto cleanup;
    }
    if (bind->entry->building) {
        pthread_mutex_unlock(&reg->mutex);
        duckdb_bind_set_error(info, "duckhts_cgranges_overlaps_bulk: index is still being constructed");
        goto cleanup;
    }
    if (!bind->entry->indexed) {
        pthread_mutex_unlock(&reg->mutex);
        duckdb_bind_set_error(
            info,
            "duckhts_cgranges_overlaps_bulk: index is not finalized; call duckhts_cgranges_index first"
        );
        goto cleanup;
    }
    pthread_mutex_unlock(&reg->mutex);

    add_overlaps_result_columns(info, bind->entry);
    duckdb_bind_set_bind_data(info, bind, destroy_overlaps_bulk_bind);
    ok = true;
cleanup:
    if (name_val) duckdb_destroy_value(&name_val);
    if (query_val) duckdb_destroy_value(&query_val);
    if (chrom_col_val) duckdb_destroy_value(&chrom_col_val);
    if (start_col_val) duckdb_destroy_value(&start_col_val);
    if (end_col_val) duckdb_destroy_value(&end_col_val);
    if (mode_val) duckdb_destroy_value(&mode_val);
    if (qid_col_val) duckdb_destroy_value(&qid_col_val);
    if (bind && !ok) destroy_overlaps_bulk_bind(bind);
}

static void overlaps_bulk_init(duckdb_init_info info) {
    duckhts_cgranges_overlaps_bulk_bind_t *bind =
        (duckhts_cgranges_overlaps_bulk_bind_t *)duckdb_init_get_bind_data(info);
    duckhts_cgranges_overlaps_bulk_init_t *init =
        (duckhts_cgranges_overlaps_bulk_init_t *)calloc(1, sizeof(*init));
    duckdb_state state;

    duckdb_init_set_max_threads(info, 1);
    if (!init) {
        duckdb_init_set_error(info, "duckhts_cgranges_overlaps_bulk: out of memory");
        return;
    }
    init->chrom_idx = -1;
    init->start_idx = -1;
    init->end_idx = -1;
    init->query_row_id_idx = -1;
    init->next_generated_query_row_id = 1;

    if (!bind || !bind->reg || !bind->entry) {
        duckdb_init_set_error(info, "duckhts_cgranges_overlaps_bulk: missing bind state");
        destroy_overlaps_bulk_init(init);
        return;
    }
    if (!bind->reg->connection) {
        duckdb_init_set_error(info, "duckhts_cgranges_overlaps_bulk: query connection is unavailable");
        destroy_overlaps_bulk_init(init);
        return;
    }
    pthread_mutex_lock(&bind->reg->mutex);
    state = duckdb_query(bind->reg->connection, bind->query, &init->query_result);
    pthread_mutex_unlock(&bind->reg->mutex);
    init->has_query_result = true;
    if (state != DuckDBSuccess) {
        duckdb_init_set_error(info, duckdb_result_error(&init->query_result));
        destroy_overlaps_bulk_init(init);
        return;
    }

    init->chrom_idx = find_result_col(&init->query_result, bind->chrom_col);
    init->start_idx = find_result_col(&init->query_result, bind->start_col);
    init->end_idx = find_result_col(&init->query_result, bind->end_col);
    if (bind->query_row_id_col) {
        init->query_row_id_idx = find_result_col(&init->query_result, bind->query_row_id_col);
    }
    if (init->chrom_idx < 0 || init->start_idx < 0 || init->end_idx < 0 ||
        (bind->query_row_id_col && init->query_row_id_idx < 0)) {
        duckdb_init_set_error(info, "duckhts_cgranges_overlaps_bulk: required columns not found in probe query");
        destroy_overlaps_bulk_init(init);
        return;
    }

    init->chrom_type = duckdb_column_type(&init->query_result, (idx_t)init->chrom_idx);
    init->start_type = duckdb_column_type(&init->query_result, (idx_t)init->start_idx);
    init->end_type = duckdb_column_type(&init->query_result, (idx_t)init->end_idx);
    if (!is_stringlike_type(init->chrom_type) || !is_numericlike_type(init->start_type) ||
        !is_numericlike_type(init->end_type)) {
        duckdb_init_set_error(
            info,
            "duckhts_cgranges_overlaps_bulk: probe query columns must be VARCHAR chrom and numeric start/end"
        );
        destroy_overlaps_bulk_init(init);
        return;
    }
    if (init->query_row_id_idx >= 0) {
        init->query_row_id_type = duckdb_column_type(&init->query_result, (idx_t)init->query_row_id_idx);
        if (!is_numericlike_type(init->query_row_id_type)) {
            duckdb_init_set_error(
                info,
                "duckhts_cgranges_overlaps_bulk: query_row_id_col must be numeric when provided"
            );
            destroy_overlaps_bulk_init(init);
            return;
        }
    }
    duckdb_init_set_init_data(info, init, destroy_overlaps_bulk_init);
}

static int overlaps_bulk_prepare_next_probe(duckhts_cgranges_overlaps_bulk_bind_t *bind,
                                            duckhts_cgranges_overlaps_bulk_init_t *init,
                                            char *err, size_t errlen) {
    for (;;) {
        if (!init->current_chunk || init->current_chunk_row >= init->current_chunk_rows) {
            if (init->current_chunk) duckdb_destroy_data_chunk(&init->current_chunk);
            init->current_chunk = duckdb_fetch_chunk(init->query_result);
            if (!init->current_chunk) {
                init->done = true;
                return 0;
            }
            init->current_chunk_rows = duckdb_data_chunk_get_size(init->current_chunk);
            init->current_chunk_row = 0;
            if (init->current_chunk_rows == 0) {
                duckdb_destroy_data_chunk(&init->current_chunk);
                continue;
            }
        }

        {
            idx_t row = init->current_chunk_row;
            int64_t row_number = init->next_generated_query_row_id;
            duckdb_vector chrom_vec = duckdb_data_chunk_get_vector(init->current_chunk, (idx_t)init->chrom_idx);
            duckdb_vector start_vec = duckdb_data_chunk_get_vector(init->current_chunk, (idx_t)init->start_idx);
            duckdb_vector end_vec = duckdb_data_chunk_get_vector(init->current_chunk, (idx_t)init->end_idx);
            duckdb_vector qid_vec = init->query_row_id_idx >= 0
                ? duckdb_data_chunk_get_vector(init->current_chunk, (idx_t)init->query_row_id_idx)
                : NULL;
            duckdb_string_t *chrom_data;
            const char *chrom;
            int64_t start;
            int64_t end;
            int64_t query_row_id;

            if (row_is_null(chrom_vec, row) || row_is_null(start_vec, row) || row_is_null(end_vec, row)) {
                snprintf(err, errlen,
                         "duckhts_cgranges_overlaps_bulk: null chrom/start/end encountered at row %llu",
                         (unsigned long long)row_number);
                return -1;
            }

            chrom_data = (duckdb_string_t *)duckdb_vector_get_data(chrom_vec);
            chrom = duckdb_string_t_data(&chrom_data[row]);
            start = chunk_cell_to_int64(start_vec, row, init->start_type);
            end = chunk_cell_to_int64(end_vec, row, init->end_type);
            if (!chrom || !*chrom || start < 0 || end < start || end > INT32_MAX) {
                snprintf(err, errlen,
                         "duckhts_cgranges_overlaps_bulk: invalid chrom/start/end at row %llu",
                         (unsigned long long)row_number);
                return -1;
            }

            query_row_id = init->next_generated_query_row_id++;
            if (qid_vec && !row_is_null(qid_vec, row)) {
                query_row_id = chunk_cell_to_int64(qid_vec, row, init->query_row_id_type);
            }

            init->current_chunk_row++;
            init->active_query_row_id = query_row_id;
            init->n_hits = strcmp(bind->mode, "contain") == 0
                ? cr_contain(bind->entry->cr, chrom, (int32_t)start, (int32_t)end, &init->hits, &init->hits_cap)
                : cr_overlap(bind->entry->cr, chrom, (int32_t)start, (int32_t)end, &init->hits, &init->hits_cap);
            init->emitted = 0;
            if (init->n_hits > 0) return 1;
        }
    }
}

static void overlaps_bulk_scan(duckdb_function_info info, duckdb_data_chunk output) {
    duckhts_cgranges_overlaps_bulk_bind_t *bind =
        (duckhts_cgranges_overlaps_bulk_bind_t *)duckdb_function_get_bind_data(info);
    duckhts_cgranges_overlaps_bulk_init_t *init =
        (duckhts_cgranges_overlaps_bulk_init_t *)duckdb_function_get_init_data(info);
    duckdb_vector qid_vec = duckdb_data_chunk_get_vector(output, 0);
    duckdb_vector ord_vec = duckdb_data_chunk_get_vector(output, 1);
    duckdb_vector label_vec = duckdb_data_chunk_get_vector(output, 2);
    duckdb_vector chrom_vec = duckdb_data_chunk_get_vector(output, 3);
    duckdb_vector start_vec = duckdb_data_chunk_get_vector(output, 4);
    duckdb_vector end_vec = duckdb_data_chunk_get_vector(output, 5);
    int64_t *qid_data = (int64_t *)duckdb_vector_get_data(qid_vec);
    int64_t *ord_data = (int64_t *)duckdb_vector_get_data(ord_vec);
    int32_t *start_data = (int32_t *)duckdb_vector_get_data(start_vec);
    int32_t *end_data = (int32_t *)duckdb_vector_get_data(end_vec);
    idx_t out_idx = 0;
    idx_t max_rows = duckdb_vector_size();
    char err[DUCKHTS_CGRANGES_ERRLEN];

    if (!bind || !init) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    memset(err, 0, sizeof(err));
    while (out_idx < max_rows) {
        if (init->emitted >= init->n_hits) {
            int prep = overlaps_bulk_prepare_next_probe(bind, init, err, sizeof(err));
            if (prep < 0) {
                duckdb_function_set_error(info, err);
                duckdb_data_chunk_set_size(output, 0);
                return;
            }
            if (prep == 0) break;
        }

        while (init->emitted < init->n_hits && out_idx < max_rows) {
            int64_t ridx = init->hits[init->emitted];
            int32_t ordinal = cr_label(bind->entry->cr, ridx);
            qid_data[out_idx] = init->active_query_row_id;
            ord_data[out_idx] = ordinal;
            switch (bind->entry->payload.kind) {
                case DUCKHTS_CGR_LABEL_BIGINT: {
                    int64_t *p = (int64_t *)duckdb_vector_get_data(label_vec);
                    p[out_idx] = bind->entry->payload.labels.i64[ordinal];
                    if (!bind->entry->payload.label_valid[ordinal]) set_null(label_vec, out_idx);
                    break;
                }
                case DUCKHTS_CGR_LABEL_DOUBLE: {
                    double *p = (double *)duckdb_vector_get_data(label_vec);
                    p[out_idx] = bind->entry->payload.labels.f64[ordinal];
                    if (!bind->entry->payload.label_valid[ordinal]) set_null(label_vec, out_idx);
                    break;
                }
                case DUCKHTS_CGR_LABEL_VARCHAR:
                    if (bind->entry->payload.label_valid[ordinal]) {
                        duckdb_vector_assign_string_element(label_vec, out_idx,
                                                            bind->entry->payload.labels.str[ordinal]);
                    } else {
                        set_null(label_vec, out_idx);
                    }
                    break;
                case DUCKHTS_CGR_LABEL_BOOLEAN: {
                    bool *p = (bool *)duckdb_vector_get_data(label_vec);
                    p[out_idx] = bind->entry->payload.labels.b[ordinal] ? true : false;
                    if (!bind->entry->payload.label_valid[ordinal]) set_null(label_vec, out_idx);
                    break;
                }
                case DUCKHTS_CGR_LABEL_ORDINAL:
                default: {
                    int64_t *p = (int64_t *)duckdb_vector_get_data(label_vec);
                    p[out_idx] = ordinal;
                    break;
                }
            }
            duckdb_vector_assign_string_element(chrom_vec, out_idx, bind->entry->payload.chroms[ordinal]);
            start_data[out_idx] = bind->entry->payload.starts[ordinal];
            end_data[out_idx] = bind->entry->payload.ends[ordinal];
            init->emitted++;
            out_idx++;
        }
    }

    duckdb_data_chunk_set_size(output, out_idx);
}

static void register_scalar_bool_1(duckdb_connection connection, duckhts_cgranges_registry_t *reg,
                                   const char *name, duckdb_scalar_function_t func) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_scalar_function_set_name(fn, name);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_set_return_type(fn, bool_type);
    duckdb_scalar_function_set_volatile(fn);
    retain_registry(reg);
    duckdb_scalar_function_set_extra_info(fn, reg, destroy_registry);
    duckdb_scalar_function_set_function(fn, func);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_scalar_function(&fn);
}

void register_duckhts_cgranges_functions(duckdb_connection connection, duckdb_database database) {
    duckhts_cgranges_registry_t *reg;
    duckdb_scalar_function fn;
    duckdb_table_function tf;
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);

    reg = (duckhts_cgranges_registry_t *)calloc(1, sizeof(*reg));
    if (!reg) return;
    pthread_mutex_init(&reg->mutex, NULL);
    reg->database = database;
    reg->connection = NULL;
    reg->ref_count = 0;
    if (duckdb_connect(database, &reg->connection) == DuckDBError || !reg->connection) {
        pthread_mutex_destroy(&reg->mutex);
        free(reg);
        duckdb_destroy_logical_type(&varchar_type);
        duckdb_destroy_logical_type(&bigint_type);
        duckdb_destroy_logical_type(&bool_type);
        return;
    }

    register_scalar_bool_1(connection, reg, "duckhts_cgranges_create", cgranges_create_scalar);
    register_scalar_bool_1(connection, reg, "duckhts_cgranges_index", cgranges_index_scalar);
    register_scalar_bool_1(connection, reg, "duckhts_cgranges_destroy", cgranges_destroy_scalar);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "duckhts_cgranges_add");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_set_return_type(fn, bool_type);
    duckdb_scalar_function_set_volatile(fn);
    retain_registry(reg);
    duckdb_scalar_function_set_extra_info(fn, reg, destroy_registry);
    duckdb_scalar_function_set_function(fn, cgranges_add_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "duckhts_cgranges_add");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_set_return_type(fn, bool_type);
    duckdb_scalar_function_set_volatile(fn);
    retain_registry(reg);
    duckdb_scalar_function_set_extra_info(fn, reg, destroy_registry);
    duckdb_scalar_function_set_function(fn, cgranges_add_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "duckhts_cgranges_add");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_set_return_type(fn, bool_type);
    duckdb_scalar_function_set_volatile(fn);
    retain_registry(reg);
    duckdb_scalar_function_set_extra_info(fn, reg, destroy_registry);
    duckdb_scalar_function_set_function(fn, cgranges_add_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "duckhts_cgranges_add");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, bool_type);
    duckdb_scalar_function_set_return_type(fn, bool_type);
    duckdb_scalar_function_set_volatile(fn);
    retain_registry(reg);
    duckdb_scalar_function_set_extra_info(fn, reg, destroy_registry);
    duckdb_scalar_function_set_function(fn, cgranges_add_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    {
        duckdb_logical_type double_type = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
        fn = duckdb_create_scalar_function();
        duckdb_scalar_function_set_name(fn, "duckhts_cgranges_add");
        duckdb_scalar_function_add_parameter(fn, varchar_type);
        duckdb_scalar_function_add_parameter(fn, varchar_type);
        duckdb_scalar_function_add_parameter(fn, bigint_type);
        duckdb_scalar_function_add_parameter(fn, bigint_type);
        duckdb_scalar_function_add_parameter(fn, double_type);
        duckdb_scalar_function_set_return_type(fn, bool_type);
        duckdb_scalar_function_set_volatile(fn);
        retain_registry(reg);
        duckdb_scalar_function_set_extra_info(fn, reg, destroy_registry);
        duckdb_scalar_function_set_function(fn, cgranges_add_scalar);
        duckdb_register_scalar_function(connection, fn);
        duckdb_destroy_scalar_function(&fn);
        duckdb_destroy_logical_type(&double_type);
    }

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "duckhts_cgranges_from_query");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_set_return_type(fn, bool_type);
    duckdb_scalar_function_set_volatile(fn);
    retain_registry(reg);
    duckdb_scalar_function_set_extra_info(fn, reg, destroy_registry);
    duckdb_scalar_function_set_function(fn, cgranges_from_query_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "duckhts_cgranges_from_query");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_set_return_type(fn, bool_type);
    duckdb_scalar_function_set_volatile(fn);
    retain_registry(reg);
    duckdb_scalar_function_set_extra_info(fn, reg, destroy_registry);
    duckdb_scalar_function_set_function(fn, cgranges_from_query_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "duckhts_cgranges_from_table");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_set_return_type(fn, bool_type);
    duckdb_scalar_function_set_volatile(fn);
    retain_registry(reg);
    duckdb_scalar_function_set_extra_info(fn, reg, destroy_registry);
    duckdb_scalar_function_set_function(fn, cgranges_from_table_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "duckhts_cgranges_from_table");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_set_return_type(fn, bool_type);
    duckdb_scalar_function_set_volatile(fn);
    retain_registry(reg);
    duckdb_scalar_function_set_extra_info(fn, reg, destroy_registry);
    duckdb_scalar_function_set_function(fn, cgranges_from_table_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "duckhts_cgranges_overlaps");
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_parameter(tf, bigint_type);
    duckdb_table_function_add_parameter(tf, bigint_type);
    duckdb_table_function_add_named_parameter(tf, "mode", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "query_row_id", bigint_type);
    retain_registry(reg);
    duckdb_table_function_set_extra_info(tf, reg, destroy_registry);
    duckdb_table_function_set_bind(tf, overlaps_bind);
    duckdb_table_function_set_init(tf, overlaps_init);
    duckdb_table_function_set_function(tf, overlaps_scan);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);

    tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "duckhts_cgranges_overlaps_bulk");
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "mode", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "query_row_id_col", varchar_type);
    retain_registry(reg);
    duckdb_table_function_set_extra_info(tf, reg, destroy_registry);
    duckdb_table_function_set_bind(tf, overlaps_bulk_bind);
    duckdb_table_function_set_init(tf, overlaps_bulk_init);
    duckdb_table_function_set_function(tf, overlaps_bulk_scan);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&bool_type);
}
