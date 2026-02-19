#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "htslib/tbx.h"
#include "htslib/vcf.h"
#include "htslib/faidx.h"

// ===============================
// Helpers
// ===============================

typedef enum {
    HTS_KIND_AUTO = 0,
    HTS_KIND_VCF,
    HTS_KIND_BCF,
    HTS_KIND_SAM,
    HTS_KIND_BAM,
    HTS_KIND_CRAM,
    HTS_KIND_FASTA,
    HTS_KIND_FASTQ,
    HTS_KIND_TABIX,
    HTS_KIND_UNKNOWN
} hts_kind_t;

static const char *compression_to_string(enum htsCompression c) {
    switch (c) {
        case no_compression: return "none";
        case gzip: return "gzip";
        case bgzf: return "bgzf";
        case bzip2_compression: return "bzip2";
        case xz_compression: return "xz";
        case zstd_compression: return "zstd";
        case razf_compression: return "razf";
        case custom: return "custom";
        default: return "unknown";
    }
}

static hts_kind_t parse_format_hint(const char *s) {
    if (!s || !*s) return HTS_KIND_AUTO;
    char buf[32];
    size_t n = strlen(s);
    if (n >= sizeof(buf)) n = sizeof(buf) - 1;
    for (size_t i = 0; i < n; i++) buf[i] = (char)tolower((unsigned char)s[i]);
    buf[n] = '\0';
    if (strcmp(buf, "auto") == 0) return HTS_KIND_AUTO;
    if (strcmp(buf, "vcf") == 0) return HTS_KIND_VCF;
    if (strcmp(buf, "bcf") == 0) return HTS_KIND_BCF;
    if (strcmp(buf, "sam") == 0) return HTS_KIND_SAM;
    if (strcmp(buf, "bam") == 0) return HTS_KIND_BAM;
    if (strcmp(buf, "cram") == 0) return HTS_KIND_CRAM;
    if (strcmp(buf, "fasta") == 0) return HTS_KIND_FASTA;
    if (strcmp(buf, "fastq") == 0) return HTS_KIND_FASTQ;
    if (strcmp(buf, "tabix") == 0) return HTS_KIND_TABIX;
    return HTS_KIND_UNKNOWN;
}

static hts_kind_t kind_from_hts_format(const htsFormat *fmt) {
    if (!fmt) return HTS_KIND_UNKNOWN;
    switch (fmt->format) {
        case vcf: return HTS_KIND_VCF;
        case bcf: return HTS_KIND_BCF;
        case sam: return HTS_KIND_SAM;
        case bam: return HTS_KIND_BAM;
        case cram: return HTS_KIND_CRAM;
        case fasta_format: return HTS_KIND_FASTA;
        case fastq_format: return HTS_KIND_FASTQ;
        case tbi: return HTS_KIND_TABIX;
        default: return HTS_KIND_UNKNOWN;
    }
}

static char *dup_str(const char *s) {
    if (!s) return NULL;
    size_t n = strlen(s) + 1;
    char *out = (char *)malloc(n);
    if (out) memcpy(out, s, n);
    return out;
}

static int64_t parse_int64(const char *s) {
    if (!s || !*s) return -1;
    char *endptr = NULL;
    long long v = strtoll(s, &endptr, 10);
    if (!endptr || endptr == s) return -1;
    return (int64_t)v;
}

// ===============================
// Header table function
// ===============================

typedef struct {
    char *record_type;
    char *id;
    char *number;
    char *value_type;
    char *description;
    int64_t length;
    int64_t idx;
    char *raw;
    int n_kv;
    char **kv_keys;
    char **kv_vals;
} hts_header_entry_t;

typedef struct {
    char *file_path;
    char *format_hint;
    char *file_format;
    char *compression;
    hts_header_entry_t *entries;
    int64_t n_entries;
} hts_header_bind_t;

typedef struct {
    int64_t offset;
} hts_header_init_t;

static void free_header_entries(hts_header_entry_t *entries, int64_t n) {
    if (!entries) return;
    for (int64_t i = 0; i < n; i++) {
        free(entries[i].record_type);
        free(entries[i].id);
        free(entries[i].number);
        free(entries[i].value_type);
        free(entries[i].description);
        free(entries[i].raw);
        if (entries[i].kv_keys) {
            for (int j = 0; j < entries[i].n_kv; j++) {
                free(entries[i].kv_keys[j]);
                free(entries[i].kv_vals[j]);
            }
            free(entries[i].kv_keys);
            free(entries[i].kv_vals);
        }
    }
    free(entries);
}

static void destroy_hts_header_bind(void *data) {
    hts_header_bind_t *bind = (hts_header_bind_t *)data;
    if (!bind) return;
    free(bind->file_path);
    free(bind->format_hint);
    free(bind->file_format);
    free(bind->compression);
    free_header_entries(bind->entries, bind->n_entries);
    free(bind);
}

static void destroy_hts_header_init(void *data) {
    if (data) free(data);
}

static void header_add_kv(hts_header_entry_t *e, const char *key, const char *val) {
    int n = e->n_kv + 1;
    char **new_keys = (char **)realloc(e->kv_keys, sizeof(char *) * (size_t)n);
    char **new_vals = (char **)realloc(e->kv_vals, sizeof(char *) * (size_t)n);
    if (!new_keys || !new_vals) return;
    e->kv_keys = new_keys;
    e->kv_vals = new_vals;
    e->kv_keys[e->n_kv] = dup_str(key ? key : "");
    e->kv_vals[e->n_kv] = dup_str(val ? val : "");
    e->n_kv = n;
}

static const char *header_kv_value(const hts_header_entry_t *e, const char *key) {
    if (!e || !key) return NULL;
    for (int i = 0; i < e->n_kv; i++) {
        if (e->kv_keys[i] && strcmp(e->kv_keys[i], key) == 0) return e->kv_vals[i];
    }
    return NULL;
}

static void build_vcf_header_entries(bcf_hdr_t *hdr, hts_header_bind_t *bind) {
    if (!hdr || !bind) return;
    if (hdr->nhrec <= 0) return;

    bind->entries = (hts_header_entry_t *)calloc((size_t)hdr->nhrec, sizeof(hts_header_entry_t));
    if (!bind->entries) return;
    bind->n_entries = hdr->nhrec;

    kstring_t ks = {0, 0, NULL};
    for (int i = 0; i < hdr->nhrec; i++) {
        bcf_hrec_t *hrec = hdr->hrec[i];
        hts_header_entry_t *e = &bind->entries[i];
        e->idx = i;
        e->length = -1;

        if (hrec && hrec->key) {
            e->record_type = dup_str(hrec->key);
        }

        if (hrec && hrec->nkeys > 0 && hrec->keys && hrec->vals) {
            for (int k = 0; k < hrec->nkeys; k++) {
                header_add_kv(e, hrec->keys[k], hrec->vals[k]);
            }
        } else if (hrec && hrec->value) {
            header_add_kv(e, "value", hrec->value);
        }

        if (bcf_hrec_format(hrec, &ks) == 0 && ks.s) {
            size_t len = strlen(ks.s);
            while (len > 0 && (ks.s[len - 1] == '\n' || ks.s[len - 1] == '\r')) {
                ks.s[len - 1] = '\0';
                len--;
            }
            e->raw = dup_str(ks.s);
        }

        const char *id = header_kv_value(e, "ID");
        if (id) e->id = dup_str(id);
        const char *num = header_kv_value(e, "Number");
        if (num) e->number = dup_str(num);
        const char *typ = header_kv_value(e, "Type");
        if (typ) e->value_type = dup_str(typ);
        const char *desc = header_kv_value(e, "Description");
        if (desc) e->description = dup_str(desc);
        const char *len = header_kv_value(e, "length");
        if (!len) len = header_kv_value(e, "Length");
        if (len) {
            int64_t v = parse_int64(len);
            if (v >= 0) e->length = v;
        }
    }
    free(ks.s);
}

static void build_sam_header_entries(sam_hdr_t *hdr, hts_header_bind_t *bind) {
    if (!hdr || !bind) return;
    const char *text = sam_hdr_str(hdr);
    size_t text_len = sam_hdr_length(hdr);
    if (!text || text_len == 0) return;

    size_t cap = 16;
    hts_header_entry_t *entries = (hts_header_entry_t *)calloc(cap, sizeof(hts_header_entry_t));
    if (!entries) return;

    int64_t n_entries = 0;
    size_t i = 0;
    while (i < text_len) {
        size_t line_start = i;
        while (i < text_len && text[i] != '\n') i++;
        size_t line_len = (i > line_start) ? (i - line_start) : 0;
        i++;
        if (line_len == 0) continue;

        const char *line = text + line_start;
        if (line[0] != '@') continue;

        if (n_entries >= (int64_t)cap) {
            cap *= 2;
            hts_header_entry_t *tmp = (hts_header_entry_t *)realloc(entries, cap * sizeof(hts_header_entry_t));
            if (!tmp) break;
            entries = tmp;
            memset(entries + n_entries, 0, (cap - (size_t)n_entries) * sizeof(hts_header_entry_t));
        }

        hts_header_entry_t *e = &entries[n_entries];
        e->idx = n_entries;
        e->length = -1;

        const char *tab = memchr(line, '\t', line_len);
        size_t type_len = tab ? (size_t)(tab - line - 1) : (line_len > 1 ? line_len - 1 : 0);
        if (type_len > 0) {
            char type_buf[8];
            size_t tlen = type_len < sizeof(type_buf) - 1 ? type_len : sizeof(type_buf) - 1;
            memcpy(type_buf, line + 1, tlen);
            type_buf[tlen] = '\0';
            e->record_type = dup_str(type_buf);
        }

        // Parse tags after first tab
        size_t pos = tab ? (size_t)(tab - line + 1) : line_len;
        while (pos < line_len) {
            size_t tok_start = pos;
            while (pos < line_len && line[pos] != '\t') pos++;
            size_t tok_len = pos - tok_start;
            pos++;
            if (tok_len == 0) continue;
            const char *tok = line + tok_start;
            const char *colon = memchr(tok, ':', tok_len);
            if (!colon) continue;
            size_t key_len = (size_t)(colon - tok);
            size_t val_len = tok_len - key_len - 1;
            if (key_len == 0) continue;
            char *key = (char *)malloc(key_len + 1);
            char *val = (char *)malloc(val_len + 1);
            if (!key || !val) {
                free(key);
                free(val);
                continue;
            }
            memcpy(key, tok, key_len);
            key[key_len] = '\0';
            memcpy(val, colon + 1, val_len);
            val[val_len] = '\0';
            header_add_kv(e, key, val);
            free(key);
            free(val);
        }

        const char *id = NULL;
        if (e->record_type) {
            if (strcmp(e->record_type, "SQ") == 0) {
                id = header_kv_value(e, "SN");
                const char *ln = header_kv_value(e, "LN");
                if (ln) {
                    int64_t v = parse_int64(ln);
                    if (v >= 0) e->length = v;
                }
            } else if (strcmp(e->record_type, "RG") == 0 || strcmp(e->record_type, "PG") == 0) {
                id = header_kv_value(e, "ID");
            }
        }
        if (id) e->id = dup_str(id);

        // Raw line
        char *raw = (char *)malloc(line_len + 1);
        if (raw) {
            memcpy(raw, line, line_len);
            raw[line_len] = '\0';
            e->raw = raw;
        }

        n_entries++;
    }

    bind->entries = entries;
    bind->n_entries = n_entries;
}

static void build_tabix_header_entries(htsFile *fp, hts_header_bind_t *bind) {
    if (!fp || !bind) return;
    kstring_t line = {0, 0, NULL};
    int64_t cap = 16;
    hts_header_entry_t *entries = (hts_header_entry_t *)calloc((size_t)cap, sizeof(hts_header_entry_t));
    if (!entries) return;

    int64_t n_entries = 0;
    while (hts_getline(fp, '\n', &line) >= 0) {
        if (line.l == 0) continue;
        if (line.s[0] != '#') break;
        if (n_entries >= cap) {
            cap *= 2;
            hts_header_entry_t *tmp = (hts_header_entry_t *)realloc(entries, (size_t)cap * sizeof(hts_header_entry_t));
            if (!tmp) break;
            entries = tmp;
            memset(entries + n_entries, 0, (size_t)(cap - n_entries) * sizeof(hts_header_entry_t));
        }
        hts_header_entry_t *e = &entries[n_entries];
        e->idx = n_entries;
        e->length = -1;
        e->record_type = dup_str("META");
        e->raw = dup_str(line.s);
        n_entries++;
    }
    free(line.s);
    bind->entries = entries;
    bind->n_entries = n_entries;
}

static void read_hts_header_bind(duckdb_bind_info info) {
    duckdb_value path_val = duckdb_bind_get_parameter(info, 0);
    char *file_path = duckdb_get_varchar(path_val);
    duckdb_destroy_value(&path_val);

    if (!file_path || strlen(file_path) == 0) {
        duckdb_bind_set_error(info, "read_hts_header requires a file path");
        if (file_path) duckdb_free(file_path);
        return;
    }

    char *format_hint = NULL;
    duckdb_value fmt_val = duckdb_bind_get_named_parameter(info, "format");
    if (fmt_val && !duckdb_is_null_value(fmt_val)) {
        format_hint = duckdb_get_varchar(fmt_val);
    }
    if (fmt_val) duckdb_destroy_value(&fmt_val);

    hts_header_bind_t *bind = (hts_header_bind_t *)calloc(1, sizeof(hts_header_bind_t));
    if (!bind) {
        duckdb_bind_set_error(info, "Out of memory");
        if (file_path) duckdb_free(file_path);
        if (format_hint) duckdb_free(format_hint);
        return;
    }
    bind->file_path = dup_str(file_path);
    bind->format_hint = format_hint ? dup_str(format_hint) : NULL;

    htsFile *fp = hts_open(file_path, "r");
    if (!fp) {
        duckdb_bind_set_error(info, "Failed to open file for header reading");
        destroy_hts_header_bind(bind);
        duckdb_free(file_path);
        if (format_hint) duckdb_free(format_hint);
        return;
    }

    const htsFormat *fmt = hts_get_format(fp);
    const char *fmt_ext = fmt ? hts_format_file_extension(fmt) : "unknown";
    bind->file_format = dup_str(fmt_ext ? fmt_ext : "unknown");
    bind->compression = dup_str(fmt ? compression_to_string(fmt->compression) : "unknown");

    hts_kind_t hint_kind = parse_format_hint(format_hint);
    hts_kind_t kind = (hint_kind == HTS_KIND_AUTO) ? kind_from_hts_format(fmt) : hint_kind;

    if (kind == HTS_KIND_VCF || kind == HTS_KIND_BCF) {
        bcf_hdr_t *hdr = bcf_hdr_read(fp);
        if (!hdr) {
            hts_close(fp);
            duckdb_bind_set_error(info, "Failed to read VCF/BCF header");
            destroy_hts_header_bind(bind);
            duckdb_free(file_path);
            if (format_hint) duckdb_free(format_hint);
            return;
        }
        build_vcf_header_entries(hdr, bind);
        bcf_hdr_destroy(hdr);
    } else if (kind == HTS_KIND_SAM || kind == HTS_KIND_BAM || kind == HTS_KIND_CRAM ||
               kind == HTS_KIND_FASTA || kind == HTS_KIND_FASTQ) {
        sam_hdr_t *hdr = sam_hdr_read(fp);
        if (hdr) {
            build_sam_header_entries(hdr, bind);
            sam_hdr_destroy(hdr);
        }
    } else if (kind == HTS_KIND_TABIX) {
        build_tabix_header_entries(fp, bind);
    }

    hts_close(fp);
    duckdb_free(file_path);
    if (format_hint) duckdb_free(format_hint);

    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type key_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type val_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type map_type = duckdb_create_map_type(key_type, val_type);

    duckdb_bind_add_result_column(info, "file_format", varchar_type);
    duckdb_bind_add_result_column(info, "compression", varchar_type);
    duckdb_bind_add_result_column(info, "record_type", varchar_type);
    duckdb_bind_add_result_column(info, "id", varchar_type);
    duckdb_bind_add_result_column(info, "number", varchar_type);
    duckdb_bind_add_result_column(info, "value_type", varchar_type);
    duckdb_bind_add_result_column(info, "length", bigint_type);
    duckdb_bind_add_result_column(info, "description", varchar_type);
    duckdb_bind_add_result_column(info, "idx", bigint_type);
    duckdb_bind_add_result_column(info, "key_values", map_type);
    duckdb_bind_add_result_column(info, "raw", varchar_type);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&key_type);
    duckdb_destroy_logical_type(&val_type);
    duckdb_destroy_logical_type(&map_type);

    duckdb_bind_set_bind_data(info, bind, destroy_hts_header_bind);
}

static void read_hts_header_init(duckdb_init_info info) {
    hts_header_init_t *init = (hts_header_init_t *)calloc(1, sizeof(hts_header_init_t));
    duckdb_init_set_init_data(info, init, destroy_hts_header_init);
}

static void read_hts_header_scan(duckdb_function_info info, duckdb_data_chunk output) {
    hts_header_bind_t *bind = (hts_header_bind_t *)duckdb_function_get_bind_data(info);
    hts_header_init_t *init = (hts_header_init_t *)duckdb_function_get_init_data(info);

    if (!bind || !init || bind->n_entries <= 0) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    idx_t vector_size = duckdb_vector_size();
    idx_t row_count = 0;

    duckdb_vector vec_file_format = duckdb_data_chunk_get_vector(output, 0);
    duckdb_vector vec_compression = duckdb_data_chunk_get_vector(output, 1);
    duckdb_vector vec_record_type = duckdb_data_chunk_get_vector(output, 2);
    duckdb_vector vec_id = duckdb_data_chunk_get_vector(output, 3);
    duckdb_vector vec_number = duckdb_data_chunk_get_vector(output, 4);
    duckdb_vector vec_value_type = duckdb_data_chunk_get_vector(output, 5);
    duckdb_vector vec_length = duckdb_data_chunk_get_vector(output, 6);
    duckdb_vector vec_desc = duckdb_data_chunk_get_vector(output, 7);
    duckdb_vector vec_idx = duckdb_data_chunk_get_vector(output, 8);
    duckdb_vector vec_kv = duckdb_data_chunk_get_vector(output, 9);
    duckdb_vector vec_raw = duckdb_data_chunk_get_vector(output, 10);

    while (row_count < vector_size && init->offset < bind->n_entries) {
        hts_header_entry_t *e = &bind->entries[init->offset];

        duckdb_vector_assign_string_element(vec_file_format, row_count,
                                            bind->file_format ? bind->file_format : "unknown");
        duckdb_vector_assign_string_element(vec_compression, row_count,
                                            bind->compression ? bind->compression : "unknown");
        if (e->record_type) duckdb_vector_assign_string_element(vec_record_type, row_count, e->record_type);
        else {
            duckdb_vector_ensure_validity_writable(vec_record_type);
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vec_record_type), row_count);
        }
        if (e->id) duckdb_vector_assign_string_element(vec_id, row_count, e->id);
        else {
            duckdb_vector_ensure_validity_writable(vec_id);
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vec_id), row_count);
        }
        if (e->number) duckdb_vector_assign_string_element(vec_number, row_count, e->number);
        else {
            duckdb_vector_ensure_validity_writable(vec_number);
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vec_number), row_count);
        }
        if (e->value_type) duckdb_vector_assign_string_element(vec_value_type, row_count, e->value_type);
        else {
            duckdb_vector_ensure_validity_writable(vec_value_type);
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vec_value_type), row_count);
        }
        if (e->length >= 0) {
            int64_t *data = (int64_t *)duckdb_vector_get_data(vec_length);
            data[row_count] = e->length;
        } else {
            duckdb_vector_ensure_validity_writable(vec_length);
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vec_length), row_count);
        }
        if (e->description) duckdb_vector_assign_string_element(vec_desc, row_count, e->description);
        else {
            duckdb_vector_ensure_validity_writable(vec_desc);
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vec_desc), row_count);
        }
        {
            int64_t *data = (int64_t *)duckdb_vector_get_data(vec_idx);
            data[row_count] = e->idx;
        }

        if (e->n_kv > 0) {
            duckdb_list_entry entry;
            entry.offset = duckdb_list_vector_get_size(vec_kv);
            entry.length = e->n_kv;
            duckdb_list_vector_reserve(vec_kv, entry.offset + entry.length);
            duckdb_list_vector_set_size(vec_kv, entry.offset + entry.length);

            duckdb_vector child = duckdb_list_vector_get_child(vec_kv);
            duckdb_vector key_vec = duckdb_struct_vector_get_child(child, 0);
            duckdb_vector val_vec = duckdb_struct_vector_get_child(child, 1);

            for (int i = 0; i < e->n_kv; i++) {
                duckdb_vector_assign_string_element(key_vec, entry.offset + i, e->kv_keys[i]);
                duckdb_vector_assign_string_element(val_vec, entry.offset + i, e->kv_vals[i]);
            }

            duckdb_list_entry *list_data = (duckdb_list_entry *)duckdb_vector_get_data(vec_kv);
            list_data[row_count] = entry;
        } else {
            duckdb_vector_ensure_validity_writable(vec_kv);
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vec_kv), row_count);
            duckdb_list_entry entry = {duckdb_list_vector_get_size(vec_kv), 0};
            duckdb_list_entry *list_data = (duckdb_list_entry *)duckdb_vector_get_data(vec_kv);
            list_data[row_count] = entry;
        }

        if (e->raw) duckdb_vector_assign_string_element(vec_raw, row_count, e->raw);
        else {
            duckdb_vector_ensure_validity_writable(vec_raw);
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vec_raw), row_count);
        }

        row_count++;
        init->offset++;
    }

    duckdb_data_chunk_set_size(output, row_count);
}

// ===============================
// Index table function
// ===============================

typedef struct {
    char *seqname;
    int64_t tid;
    int64_t length;
    int has_stat;
    uint64_t mapped;
    uint64_t unmapped;
    int has_n_no_coor;
    uint64_t n_no_coor;
    char *index_type;
    char *index_path;
} hts_index_entry_t;

typedef struct {
    char *file_path;
    char *format_hint;
    char *file_format;
    char *compression;
    hts_index_entry_t *entries;
    int64_t n_entries;
    uint8_t *meta;
    uint32_t meta_len;
} hts_index_bind_t;

typedef struct {
    int64_t offset;
} hts_index_init_t;

static void free_index_entries(hts_index_entry_t *entries, int64_t n) {
    if (!entries) return;
    for (int64_t i = 0; i < n; i++) {
        free(entries[i].seqname);
        free(entries[i].index_type);
        free(entries[i].index_path);
    }
    free(entries);
}

static void destroy_hts_index_bind(void *data) {
    hts_index_bind_t *bind = (hts_index_bind_t *)data;
    if (!bind) return;
    free(bind->file_path);
    free(bind->format_hint);
    free(bind->file_format);
    free(bind->compression);
    free_index_entries(bind->entries, bind->n_entries);
    free(bind->meta);
    free(bind);
}

static void destroy_hts_index_init(void *data) {
    if (data) free(data);
}

static const char *index_fmt_to_string(int fmt) {
    switch (fmt) {
        case HTS_FMT_BAI: return "BAI";
        case HTS_FMT_CSI: return "CSI";
        case HTS_FMT_TBI: return "TBI";
        case HTS_FMT_CRAI: return "CRAI";
        default: return "UNKNOWN";
    }
}

static void add_index_entry(hts_index_bind_t *bind, const char *seqname, int64_t tid,
                            int64_t length, int has_stat, uint64_t mapped, uint64_t unmapped,
                            int has_n_no_coor, uint64_t n_no_coor, const char *index_type,
                            const char *index_path) {
    int64_t n = bind->n_entries + 1;
    hts_index_entry_t *tmp = (hts_index_entry_t *)realloc(bind->entries, (size_t)n * sizeof(hts_index_entry_t));
    if (!tmp) return;
    bind->entries = tmp;
    hts_index_entry_t *e = &bind->entries[bind->n_entries];
    memset(e, 0, sizeof(*e));
    e->seqname = dup_str(seqname);
    e->tid = tid;
    e->length = length;
    e->has_stat = has_stat;
    e->mapped = mapped;
    e->unmapped = unmapped;
    e->has_n_no_coor = has_n_no_coor;
    e->n_no_coor = n_no_coor;
    e->index_type = dup_str(index_type ? index_type : "UNKNOWN");
    e->index_path = index_path ? dup_str(index_path) : NULL;
    bind->n_entries = n;
}

static void capture_index_meta(hts_index_bind_t *bind, hts_idx_t *idx) {
    if (!bind || !idx) return;
    uint32_t meta_len = 0;
    uint8_t *meta = hts_idx_get_meta(idx, &meta_len);
    if (meta && meta_len > 0) {
        bind->meta = (uint8_t *)malloc(meta_len);
        if (bind->meta) {
            memcpy(bind->meta, meta, meta_len);
            bind->meta_len = meta_len;
        }
    }
}

static void read_hts_index_bind(duckdb_bind_info info) {
    duckdb_value path_val = duckdb_bind_get_parameter(info, 0);
    char *file_path = duckdb_get_varchar(path_val);
    duckdb_destroy_value(&path_val);

    if (!file_path || strlen(file_path) == 0) {
        duckdb_bind_set_error(info, "read_hts_index requires a file path");
        if (file_path) duckdb_free(file_path);
        return;
    }

    char *format_hint = NULL;
    duckdb_value fmt_val = duckdb_bind_get_named_parameter(info, "format");
    if (fmt_val && !duckdb_is_null_value(fmt_val)) {
        format_hint = duckdb_get_varchar(fmt_val);
    }
    if (fmt_val) duckdb_destroy_value(&fmt_val);

    char *index_path = NULL;
    duckdb_value idx_val = duckdb_bind_get_named_parameter(info, "index_path");
    if (idx_val && !duckdb_is_null_value(idx_val)) {
        index_path = duckdb_get_varchar(idx_val);
    }
    if (idx_val) duckdb_destroy_value(&idx_val);

    hts_index_bind_t *bind = (hts_index_bind_t *)calloc(1, sizeof(hts_index_bind_t));
    if (!bind) {
        duckdb_bind_set_error(info, "Out of memory");
        if (file_path) duckdb_free(file_path);
        if (format_hint) duckdb_free(format_hint);
        if (index_path) duckdb_free(index_path);
        return;
    }
    bind->file_path = dup_str(file_path);
    bind->format_hint = format_hint ? dup_str(format_hint) : NULL;

    htsFile *fp = hts_open(file_path, "r");
    if (!fp) {
        duckdb_bind_set_error(info, "Failed to open file for index reading");
        destroy_hts_index_bind(bind);
        duckdb_free(file_path);
        if (format_hint) duckdb_free(format_hint);
        if (index_path) duckdb_free(index_path);
        return;
    }

    const htsFormat *fmt = hts_get_format(fp);
    const char *fmt_ext = fmt ? hts_format_file_extension(fmt) : "unknown";
    bind->file_format = dup_str(fmt_ext ? fmt_ext : "unknown");
    bind->compression = dup_str(fmt ? compression_to_string(fmt->compression) : "unknown");

    hts_kind_t hint_kind = parse_format_hint(format_hint);
    hts_kind_t kind = (hint_kind == HTS_KIND_AUTO) ? kind_from_hts_format(fmt) : hint_kind;

    if (kind == HTS_KIND_SAM || kind == HTS_KIND_BAM || kind == HTS_KIND_CRAM) {
        sam_hdr_t *hdr = sam_hdr_read(fp);
        if (!hdr) {
            hts_close(fp);
            duckdb_bind_set_error(info, "Failed to read SAM/BAM/CRAM header");
            destroy_hts_index_bind(bind);
            duckdb_free(file_path);
            if (format_hint) duckdb_free(format_hint);
            if (index_path) duckdb_free(index_path);
            return;
        }

        int flags = HTS_IDX_SILENT_FAIL;
        hts_idx_t *idx = sam_index_load3(fp, file_path, index_path, flags);
        if (!idx) {
            sam_hdr_destroy(hdr);
            hts_close(fp);
            duckdb_bind_set_error(info, "Failed to load index for SAM/BAM/CRAM file");
            destroy_hts_index_bind(bind);
            duckdb_free(file_path);
            if (format_hint) duckdb_free(format_hint);
            if (index_path) duckdb_free(index_path);
            return;
        }

        capture_index_meta(bind, idx);
        int nseq = hts_idx_nseq(idx);
        uint64_t n_no_coor = hts_idx_get_n_no_coor(idx);
        for (int tid = 0; tid < nseq; tid++) {
            const char *name = sam_hdr_tid2name(hdr, tid);
            int64_t len = (int64_t)sam_hdr_tid2len(hdr, tid);
            uint64_t mapped = 0, unmapped = 0;
            int has_stat = (hts_idx_get_stat(idx, tid, &mapped, &unmapped) == 0);
            add_index_entry(bind, name, tid, len, has_stat, mapped, unmapped, 1, n_no_coor,
                            index_fmt_to_string(hts_idx_fmt(idx)), index_path);
        }
        hts_idx_destroy(idx);
        sam_hdr_destroy(hdr);
    } else if (kind == HTS_KIND_VCF || kind == HTS_KIND_BCF) {
        bcf_hdr_t *hdr = bcf_hdr_read(fp);
        if (!hdr) {
            hts_close(fp);
            duckdb_bind_set_error(info, "Failed to read VCF/BCF header");
            destroy_hts_index_bind(bind);
            duckdb_free(file_path);
            if (format_hint) duckdb_free(format_hint);
            if (index_path) duckdb_free(index_path);
            return;
        }

        int flags = HTS_IDX_SILENT_FAIL;
        hts_idx_t *idx = NULL;
        tbx_t *tbx = NULL;
        tbx = tbx_index_load3(file_path, index_path, flags);
        if (tbx) idx = tbx->idx;
        if (!idx) {
            idx = bcf_index_load3(file_path, index_path, flags);
        }
        if (!idx) {
            bcf_hdr_destroy(hdr);
            hts_close(fp);
            duckdb_bind_set_error(info, "Failed to load index for VCF/BCF file");
            destroy_hts_index_bind(bind);
            duckdb_free(file_path);
            if (format_hint) duckdb_free(format_hint);
            if (index_path) duckdb_free(index_path);
            return;
        }

        capture_index_meta(bind, idx);
        int nseq = hts_idx_nseq(idx);
        for (int tid = 0; tid < nseq; tid++) {
            const char *name = bcf_hdr_id2name(hdr, tid);
            uint64_t mapped = 0, unmapped = 0;
            int has_stat = (hts_idx_get_stat(idx, tid, &mapped, &unmapped) == 0);
            add_index_entry(bind, name, tid, -1, has_stat, mapped, unmapped, 0, 0,
                            index_fmt_to_string(hts_idx_fmt(idx)), index_path);
        }
        if (tbx) tbx_destroy(tbx);
        else hts_idx_destroy(idx);
        bcf_hdr_destroy(hdr);
    } else if (kind == HTS_KIND_TABIX) {
        int flags = HTS_IDX_SILENT_FAIL;
        tbx_t *tbx = tbx_index_load3(file_path, index_path, flags);
        if (!tbx) {
            hts_close(fp);
            duckdb_bind_set_error(info, "Failed to load tabix index");
            destroy_hts_index_bind(bind);
            duckdb_free(file_path);
            if (format_hint) duckdb_free(format_hint);
            if (index_path) duckdb_free(index_path);
            return;
        }
        capture_index_meta(bind, tbx->idx);
        int n = 0;
        const char **names = tbx_seqnames(tbx, &n);
        for (int tid = 0; tid < n; tid++) {
            uint64_t mapped = 0, unmapped = 0;
            int has_stat = (hts_idx_get_stat(tbx->idx, tid, &mapped, &unmapped) == 0);
            add_index_entry(bind, names[tid], tid, -1, has_stat, mapped, unmapped, 0, 0,
                            index_fmt_to_string(hts_idx_fmt(tbx->idx)), index_path);
        }
        free(names);
        tbx_destroy(tbx);
    } else if (kind == HTS_KIND_FASTA || kind == HTS_KIND_FASTQ) {
        enum fai_format_options ffmt = (kind == HTS_KIND_FASTQ) ? FAI_FASTQ : FAI_FASTA;
        faidx_t *fai = fai_load3_format(file_path, index_path, NULL, 0, ffmt);
        if (!fai) {
            hts_close(fp);
            duckdb_bind_set_error(info, "Failed to load FASTA/FASTQ index");
            destroy_hts_index_bind(bind);
            duckdb_free(file_path);
            if (format_hint) duckdb_free(format_hint);
            if (index_path) duckdb_free(index_path);
            return;
        }
        int nseq = faidx_nseq(fai);
        const char *index_type = (kind == HTS_KIND_FASTQ) ? "FQI" : "FAI";
        for (int tid = 0; tid < nseq; tid++) {
            const char *name = faidx_iseq(fai, tid);
            int64_t len = (int64_t)faidx_seq_len64(fai, name);
            add_index_entry(bind, name, tid, len, 0, 0, 0, 0, 0, index_type, index_path);
        }
        fai_destroy(fai);
    }

    hts_close(fp);
    duckdb_free(file_path);
    if (format_hint) duckdb_free(format_hint);
    if (index_path) duckdb_free(index_path);

    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type blob_type = duckdb_create_logical_type(DUCKDB_TYPE_BLOB);

    duckdb_bind_add_result_column(info, "file_format", varchar_type);
    duckdb_bind_add_result_column(info, "seqname", varchar_type);
    duckdb_bind_add_result_column(info, "tid", bigint_type);
    duckdb_bind_add_result_column(info, "length", bigint_type);
    duckdb_bind_add_result_column(info, "mapped", bigint_type);
    duckdb_bind_add_result_column(info, "unmapped", bigint_type);
    duckdb_bind_add_result_column(info, "n_no_coor", bigint_type);
    duckdb_bind_add_result_column(info, "index_type", varchar_type);
    duckdb_bind_add_result_column(info, "index_path", varchar_type);
    duckdb_bind_add_result_column(info, "meta", blob_type);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&blob_type);

    duckdb_bind_set_bind_data(info, bind, destroy_hts_index_bind);
}

static void read_hts_index_init(duckdb_init_info info) {
    hts_index_init_t *init = (hts_index_init_t *)calloc(1, sizeof(hts_index_init_t));
    duckdb_init_set_init_data(info, init, destroy_hts_index_init);
}

static void read_hts_index_scan(duckdb_function_info info, duckdb_data_chunk output) {
    hts_index_bind_t *bind = (hts_index_bind_t *)duckdb_function_get_bind_data(info);
    hts_index_init_t *init = (hts_index_init_t *)duckdb_function_get_init_data(info);

    if (!bind || !init || bind->n_entries <= 0) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    idx_t vector_size = duckdb_vector_size();
    idx_t row_count = 0;

    duckdb_vector vec_file_format = duckdb_data_chunk_get_vector(output, 0);
    duckdb_vector vec_seqname = duckdb_data_chunk_get_vector(output, 1);
    duckdb_vector vec_tid = duckdb_data_chunk_get_vector(output, 2);
    duckdb_vector vec_length = duckdb_data_chunk_get_vector(output, 3);
    duckdb_vector vec_mapped = duckdb_data_chunk_get_vector(output, 4);
    duckdb_vector vec_unmapped = duckdb_data_chunk_get_vector(output, 5);
    duckdb_vector vec_n_no_coor = duckdb_data_chunk_get_vector(output, 6);
    duckdb_vector vec_index_type = duckdb_data_chunk_get_vector(output, 7);
    duckdb_vector vec_index_path = duckdb_data_chunk_get_vector(output, 8);
    duckdb_vector vec_meta = duckdb_data_chunk_get_vector(output, 9);

    while (row_count < vector_size && init->offset < bind->n_entries) {
        hts_index_entry_t *e = &bind->entries[init->offset];

        duckdb_vector_assign_string_element(vec_file_format, row_count,
                                            bind->file_format ? bind->file_format : "unknown");
        if (e->seqname) duckdb_vector_assign_string_element(vec_seqname, row_count, e->seqname);
        else {
            duckdb_vector_ensure_validity_writable(vec_seqname);
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vec_seqname), row_count);
        }

        int64_t *tid_data = (int64_t *)duckdb_vector_get_data(vec_tid);
        tid_data[row_count] = e->tid;

        if (e->length >= 0) {
            int64_t *len_data = (int64_t *)duckdb_vector_get_data(vec_length);
            len_data[row_count] = e->length;
        } else {
            duckdb_vector_ensure_validity_writable(vec_length);
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vec_length), row_count);
        }

        if (e->has_stat) {
            int64_t *mapped_data = (int64_t *)duckdb_vector_get_data(vec_mapped);
            int64_t *unmapped_data = (int64_t *)duckdb_vector_get_data(vec_unmapped);
            mapped_data[row_count] = (int64_t)e->mapped;
            unmapped_data[row_count] = (int64_t)e->unmapped;
        } else {
            duckdb_vector_ensure_validity_writable(vec_mapped);
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vec_mapped), row_count);
            duckdb_vector_ensure_validity_writable(vec_unmapped);
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vec_unmapped), row_count);
        }

        if (e->has_n_no_coor) {
            int64_t *nn_data = (int64_t *)duckdb_vector_get_data(vec_n_no_coor);
            nn_data[row_count] = (int64_t)e->n_no_coor;
        } else {
            duckdb_vector_ensure_validity_writable(vec_n_no_coor);
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vec_n_no_coor), row_count);
        }

        if (e->index_type) duckdb_vector_assign_string_element(vec_index_type, row_count, e->index_type);
        else {
            duckdb_vector_ensure_validity_writable(vec_index_type);
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vec_index_type), row_count);
        }
        if (e->index_path) duckdb_vector_assign_string_element(vec_index_path, row_count, e->index_path);
        else {
            duckdb_vector_ensure_validity_writable(vec_index_path);
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vec_index_path), row_count);
        }

        if (bind->meta && bind->meta_len > 0) {
            duckdb_vector_assign_string_element_len(vec_meta, row_count, (const char *)bind->meta, bind->meta_len);
        } else {
            duckdb_vector_ensure_validity_writable(vec_meta);
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vec_meta), row_count);
        }

        row_count++;
        init->offset++;
    }

    duckdb_data_chunk_set_size(output, row_count);
}

// ===============================
// Registration
// ===============================

void register_read_hts_header_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "read_hts_header");

    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "format", varchar_type);

    duckdb_table_function_set_bind(tf, read_hts_header_bind);
    duckdb_table_function_set_init(tf, read_hts_header_init);
    duckdb_table_function_set_function(tf, read_hts_header_scan);

    duckdb_register_table_function(connection, tf);
    duckdb_destroy_logical_type(&varchar_type);
}

void register_read_hts_index_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "read_hts_index");

    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "format", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "index_path", varchar_type);

    duckdb_table_function_set_bind(tf, read_hts_index_bind);
    duckdb_table_function_set_init(tf, read_hts_index_init);
    duckdb_table_function_set_function(tf, read_hts_index_scan);

    duckdb_register_table_function(connection, tf);
    duckdb_destroy_logical_type(&varchar_type);
}
