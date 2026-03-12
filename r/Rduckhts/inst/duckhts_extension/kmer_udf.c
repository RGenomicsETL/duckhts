#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <ctype.h>
#include <stdint.h>
#include <string.h>

#define SAM_FLAG_PAIRED 0x1
#define SAM_FLAG_PROPER_PAIR 0x2
#define SAM_FLAG_UNMAPPED 0x4
#define SAM_FLAG_MATE_UNMAPPED 0x8
#define SAM_FLAG_REVERSE 0x10
#define SAM_FLAG_MATE_REVERSE 0x20
#define SAM_FLAG_READ1 0x40
#define SAM_FLAG_READ2 0x80
#define SAM_FLAG_SECONDARY 0x100
#define SAM_FLAG_QCFAIL 0x200
#define SAM_FLAG_DUPLICATE 0x400
#define SAM_FLAG_SUPPLEMENTARY 0x800

static inline void set_null_at(duckdb_vector vector, idx_t row) {
    duckdb_vector_ensure_validity_writable(vector);
    uint64_t *validity = duckdb_vector_get_validity(vector);
    duckdb_validity_set_row_invalid(validity, row);
}

static inline int row_is_valid(duckdb_vector vector, idx_t row) {
    uint64_t *validity = duckdb_vector_get_validity(vector);
    if (!validity) {
        return 1;
    }
    return duckdb_validity_row_is_valid(validity, row);
}

static inline char dna_complement(char c) {
    switch ((unsigned char)toupper((unsigned char)c)) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    case 'N': return 'N';
    default:  return '\0';
    }
}

static inline int dna_to_2bit(char c) {
    switch ((unsigned char)toupper((unsigned char)c)) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default:  return -1;
    }
}

static inline int iupac_to_4bit(char c) {
    switch ((unsigned char)toupper((unsigned char)c)) {
    case 'A': return 0x1;
    case 'C': return 0x2;
    case 'G': return 0x4;
    case 'T': return 0x8;
    case 'M': return 0x3;
    case 'R': return 0x5;
    case 'S': return 0x6;
    case 'V': return 0x7;
    case 'W': return 0x9;
    case 'Y': return 0xa;
    case 'H': return 0xb;
    case 'K': return 0xc;
    case 'D': return 0xd;
    case 'B': return 0xe;
    case 'N': return 0xf;
    default:  return -1;
    }
}

static inline char bit4_to_iupac(uint8_t code) {
    switch (code) {
    case 0x1: return 'A';
    case 0x2: return 'C';
    case 0x4: return 'G';
    case 0x8: return 'T';
    case 0x3: return 'M';
    case 0x5: return 'R';
    case 0x6: return 'S';
    case 0x7: return 'V';
    case 0x9: return 'W';
    case 0xa: return 'Y';
    case 0xb: return 'H';
    case 0xc: return 'K';
    case 0xd: return 'D';
    case 0xe: return 'B';
    case 0xf: return 'N';
    default:  return '\0';
    }
}

static inline const char *get_string_at(duckdb_vector vector, idx_t row, idx_t *len) {
    duckdb_string_t *data = (duckdb_string_t *)duckdb_vector_get_data(vector);
    duckdb_string_t *val = &data[row];
    *len = duckdb_string_t_length(*val);
    return duckdb_string_t_data(val);
}

static inline int64_t get_int64_at(duckdb_vector vector, idx_t row) {
    duckdb_logical_type logical_type = duckdb_vector_get_column_type(vector);
    duckdb_type type = duckdb_get_type_id(logical_type);
    void *data = duckdb_vector_get_data(vector);
    int64_t result = 0;

    switch (type) {
    case DUCKDB_TYPE_TINYINT:
        result = ((int8_t *)data)[row];
        break;
    case DUCKDB_TYPE_SMALLINT:
        result = ((int16_t *)data)[row];
        break;
    case DUCKDB_TYPE_INTEGER:
        result = ((int32_t *)data)[row];
        break;
    case DUCKDB_TYPE_BIGINT:
        result = ((int64_t *)data)[row];
        break;
    case DUCKDB_TYPE_UTINYINT:
        result = ((uint8_t *)data)[row];
        break;
    case DUCKDB_TYPE_USMALLINT:
        result = ((uint16_t *)data)[row];
        break;
    case DUCKDB_TYPE_UINTEGER:
        result = ((uint32_t *)data)[row];
        break;
    case DUCKDB_TYPE_UBIGINT:
        result = (int64_t)((uint64_t *)data)[row];
        break;
    default:
        result = 0;
        break;
    }
    duckdb_destroy_logical_type(&logical_type);
    return result;
}

static void seq_revcomp_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    (void)info;
    duckdb_vector seq_vec = duckdb_data_chunk_get_vector(input, 0);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(seq_vec, row)) {
            set_null_at(output, row);
            continue;
        }

        idx_t len = 0;
        const char *seq = get_string_at(seq_vec, row, &len);
        char *out = (char *)duckdb_malloc((size_t)len + 1);
        if (!out) {
            duckdb_scalar_function_set_error(info, "seq_revcomp: out of memory");
            return;
        }

        int valid = 1;
        for (idx_t i = 0; i < len; i++) {
            char rc = dna_complement(seq[len - 1 - i]);
            if (!rc) {
                valid = 0;
                break;
            }
            out[i] = rc;
        }

        if (!valid) {
            set_null_at(output, row);
            duckdb_free(out);
            continue;
        }

        out[len] = '\0';
        duckdb_vector_assign_string_element_len(output, row, out, len);
        duckdb_free(out);
    }
}

static void seq_canonical_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    (void)info;
    duckdb_vector seq_vec = duckdb_data_chunk_get_vector(input, 0);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(seq_vec, row)) {
            set_null_at(output, row);
            continue;
        }

        idx_t len = 0;
        const char *seq = get_string_at(seq_vec, row, &len);

        char *fwd = (char *)duckdb_malloc((size_t)len + 1);
        char *rev = (char *)duckdb_malloc((size_t)len + 1);
        if (!fwd || !rev) {
            if (fwd) duckdb_free(fwd);
            if (rev) duckdb_free(rev);
            duckdb_scalar_function_set_error(info, "seq_canonical: out of memory");
            return;
        }

        int valid = 1;
        for (idx_t i = 0; i < len; i++) {
            unsigned char c = (unsigned char)toupper((unsigned char)seq[i]);
            char rc = dna_complement(seq[len - 1 - i]);
            if (!rc || (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N')) {
                valid = 0;
                break;
            }
            fwd[i] = (char)c;
            rev[i] = rc;
        }

        if (!valid) {
            set_null_at(output, row);
            duckdb_free(fwd);
            duckdb_free(rev);
            continue;
        }

        fwd[len] = '\0';
        rev[len] = '\0';

        const char *chosen = (memcmp(fwd, rev, (size_t)len) <= 0) ? fwd : rev;
        duckdb_vector_assign_string_element_len(output, row, chosen, len);
        duckdb_free(fwd);
        duckdb_free(rev);
    }
}

static void seq_hash_2bit_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    (void)info;
    duckdb_vector seq_vec = duckdb_data_chunk_get_vector(input, 0);
    uint64_t *out_data = (uint64_t *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(seq_vec, row)) {
            set_null_at(output, row);
            continue;
        }

        idx_t len = 0;
        const char *seq = get_string_at(seq_vec, row, &len);
        if (len > 32) {
            set_null_at(output, row);
            continue;
        }

        uint64_t h = 0;
        int valid = 1;
        for (idx_t i = 0; i < len; i++) {
            int code = dna_to_2bit(seq[i]);
            if (code < 0) {
                valid = 0;
                break;
            }
            h = (h << 2) | (uint64_t)code;
        }

        if (!valid) {
            set_null_at(output, row);
            continue;
        }

        out_data[row] = h;
    }
}

static void seq_encode_4bit_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    (void)info;
    duckdb_vector seq_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_list_entry *list_data = (duckdb_list_entry *)duckdb_vector_get_data(output);
    duckdb_vector child_vec = duckdb_list_vector_get_child(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);
    idx_t child_offset = duckdb_list_vector_get_size(output);

    for (idx_t row = 0; row < row_count; row++) {
        list_data[row].offset = child_offset;
        list_data[row].length = 0;

        if (!row_is_valid(seq_vec, row)) {
            set_null_at(output, row);
            continue;
        }

        idx_t len = 0;
        const char *seq = get_string_at(seq_vec, row, &len);
        if (duckdb_list_vector_reserve(output, child_offset + len) != DuckDBSuccess) {
            duckdb_scalar_function_set_error(info, "seq_encode_4bit: failed to reserve list storage");
            return;
        }
        if (duckdb_list_vector_set_size(output, child_offset + len) != DuckDBSuccess) {
            duckdb_scalar_function_set_error(info, "seq_encode_4bit: failed to grow list storage");
            return;
        }

        uint8_t *child_data = (uint8_t *)duckdb_vector_get_data(child_vec);
        int valid = 1;
        for (idx_t i = 0; i < len; i++) {
            int code = iupac_to_4bit(seq[i]);
            if (code < 0) {
                valid = 0;
                break;
            }
            child_data[child_offset + i] = (uint8_t)code;
        }

        if (!valid) {
            set_null_at(output, row);
            if (duckdb_list_vector_set_size(output, child_offset) != DuckDBSuccess) {
                duckdb_scalar_function_set_error(info, "seq_encode_4bit: failed to roll back list storage");
                return;
            }
            continue;
        }

        list_data[row].length = len;
        child_offset += len;
    }
}

static void seq_decode_4bit_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    (void)info;
    duckdb_vector codes_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_list_entry *list_data = (duckdb_list_entry *)duckdb_vector_get_data(codes_vec);
    duckdb_vector child_vec = duckdb_list_vector_get_child(codes_vec);
    uint8_t *child_data = (uint8_t *)duckdb_vector_get_data(child_vec);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(codes_vec, row)) {
            set_null_at(output, row);
            continue;
        }

        duckdb_list_entry entry = list_data[row];
        char *decoded = (char *)duckdb_malloc((size_t)entry.length + 1);
        if (!decoded) {
            duckdb_scalar_function_set_error(info, "seq_decode_4bit: out of memory");
            return;
        }

        int valid = 1;
        for (idx_t i = 0; i < entry.length; i++) {
            idx_t child_index = entry.offset + i;
            if (!row_is_valid(child_vec, child_index)) {
                valid = 0;
                break;
            }
            char base = bit4_to_iupac(child_data[child_index]);
            if (!base) {
                valid = 0;
                break;
            }
            decoded[i] = base;
        }

        if (!valid) {
            set_null_at(output, row);
            duckdb_free(decoded);
            continue;
        }

        decoded[entry.length] = '\0';
        duckdb_vector_assign_string_element_len(output, row, decoded, entry.length);
        duckdb_free(decoded);
    }
}

static void seq_gc_content_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    (void)info;
    duckdb_vector seq_vec = duckdb_data_chunk_get_vector(input, 0);
    double *out_data = (double *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(seq_vec, row)) {
            set_null_at(output, row);
            continue;
        }

        idx_t len = 0;
        const char *seq = get_string_at(seq_vec, row, &len);
        if (len == 0) {
            set_null_at(output, row);
            continue;
        }

        idx_t gc = 0;
        idx_t called = 0;
        int valid = 1;
        for (idx_t i = 0; i < len; i++) {
            unsigned char c = (unsigned char)toupper((unsigned char)seq[i]);
            switch (c) {
            case 'G':
            case 'C':
                gc++;
                called++;
                break;
            case 'A':
            case 'T':
                called++;
                break;
            case 'N':
                break;
            default:
                valid = 0;
                break;
            }
            if (!valid) {
                break;
            }
        }

        if (!valid || called == 0) {
            set_null_at(output, row);
            continue;
        }
        out_data[row] = (double)gc / (double)called;
    }
}

static void sam_flag_scalar(duckdb_function_info info,
                            duckdb_data_chunk input,
                            duckdb_vector output,
                            uint16_t mask) {
    (void)info;
    duckdb_vector flag_vec = duckdb_data_chunk_get_vector(input, 0);
    bool *out_data = (bool *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(flag_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        int64_t flag_value = get_int64_at(flag_vec, row);
        if (flag_value < 0 || flag_value > 0xffff) {
            set_null_at(output, row);
            continue;
        }
        out_data[row] = ((((uint16_t)flag_value) & mask) != 0);
    }
}

static void sam_is_segmented_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    sam_flag_scalar(info, input, output, SAM_FLAG_PAIRED);
}

static void sam_is_properly_aligned_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    sam_flag_scalar(info, input, output, SAM_FLAG_PROPER_PAIR);
}

static void sam_is_unmapped_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    sam_flag_scalar(info, input, output, SAM_FLAG_UNMAPPED);
}

static void sam_is_mate_unmapped_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    sam_flag_scalar(info, input, output, SAM_FLAG_MATE_UNMAPPED);
}

static void sam_is_reverse_complemented_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    sam_flag_scalar(info, input, output, SAM_FLAG_REVERSE);
}

static void sam_is_mate_reverse_complemented_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    sam_flag_scalar(info, input, output, SAM_FLAG_MATE_REVERSE);
}

static void sam_is_first_segment_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    sam_flag_scalar(info, input, output, SAM_FLAG_READ1);
}

static void sam_is_last_segment_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    sam_flag_scalar(info, input, output, SAM_FLAG_READ2);
}

static void sam_is_secondary_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    sam_flag_scalar(info, input, output, SAM_FLAG_SECONDARY);
}

static void sam_is_qc_fail_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    sam_flag_scalar(info, input, output, SAM_FLAG_QCFAIL);
}

static void sam_is_duplicate_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    sam_flag_scalar(info, input, output, SAM_FLAG_DUPLICATE);
}

static void sam_is_supplementary_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    sam_flag_scalar(info, input, output, SAM_FLAG_SUPPLEMENTARY);
}

typedef struct {
    char *sequence;
    idx_t seq_len;
    idx_t k;
    int canonical;
} seq_kmers_bind_t;

typedef struct {
    idx_t next_pos;
} seq_kmers_init_t;

static void destroy_seq_kmers_bind(void *data) {
    seq_kmers_bind_t *bind = (seq_kmers_bind_t *)data;
    if (!bind) {
        return;
    }
    if (bind->sequence) {
        duckdb_free(bind->sequence);
    }
    duckdb_free(bind);
}

static void destroy_seq_kmers_init(void *data) {
    if (data) {
        duckdb_free(data);
    }
}

static void seq_kmers_bind(duckdb_bind_info info) {
    duckdb_value seq_val = duckdb_bind_get_parameter(info, 0);
    duckdb_value k_val = duckdb_bind_get_parameter(info, 1);
    duckdb_value canonical_val = duckdb_bind_get_named_parameter(info, "canonical");

    if (!seq_val || duckdb_is_null_value(seq_val)) {
        duckdb_bind_set_error(info, "seq_kmers: sequence must not be NULL");
        if (seq_val) duckdb_destroy_value(&seq_val);
        if (k_val) duckdb_destroy_value(&k_val);
        if (canonical_val) duckdb_destroy_value(&canonical_val);
        return;
    }

    if (!k_val || duckdb_is_null_value(k_val)) {
        duckdb_bind_set_error(info, "seq_kmers: k must not be NULL");
        if (seq_val) duckdb_destroy_value(&seq_val);
        if (k_val) duckdb_destroy_value(&k_val);
        if (canonical_val) duckdb_destroy_value(&canonical_val);
        return;
    }

    char *sequence = duckdb_get_varchar(seq_val);
    int64_t k = duckdb_get_int64(k_val);
    int canonical = 0;

    if (canonical_val && !duckdb_is_null_value(canonical_val)) {
        canonical = duckdb_get_bool(canonical_val) ? 1 : 0;
    }

    if (seq_val) duckdb_destroy_value(&seq_val);
    if (k_val) duckdb_destroy_value(&k_val);
    if (canonical_val) duckdb_destroy_value(&canonical_val);

    if (!sequence) {
        duckdb_bind_set_error(info, "seq_kmers: failed to read sequence");
        return;
    }
    if (k <= 0) {
        duckdb_bind_set_error(info, "seq_kmers: k must be > 0");
        duckdb_free(sequence);
        return;
    }

    seq_kmers_bind_t *bind = (seq_kmers_bind_t *)duckdb_malloc(sizeof(seq_kmers_bind_t));
    if (!bind) {
        duckdb_bind_set_error(info, "seq_kmers: out of memory");
        duckdb_free(sequence);
        return;
    }

    bind->sequence = sequence;
    bind->seq_len = (idx_t)strlen(sequence);
    bind->k = (idx_t)k;
    bind->canonical = canonical;

    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_bind_add_result_column(info, "pos", bigint_type);
    duckdb_bind_add_result_column(info, "kmer", varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&varchar_type);

    idx_t n_kmers = 0;
    if (bind->seq_len >= bind->k) {
        n_kmers = bind->seq_len - bind->k + 1;
    }
    duckdb_bind_set_cardinality(info, n_kmers, true);
    duckdb_bind_set_bind_data(info, bind, destroy_seq_kmers_bind);
}

static void seq_kmers_init(duckdb_init_info info) {
    seq_kmers_init_t *init = (seq_kmers_init_t *)duckdb_malloc(sizeof(seq_kmers_init_t));
    if (!init) {
        duckdb_init_set_error(info, "seq_kmers: out of memory");
        return;
    }
    init->next_pos = 0;
    duckdb_init_set_max_threads(info, 1);
    duckdb_init_set_init_data(info, init, destroy_seq_kmers_init);
}

static void seq_kmers_function(duckdb_function_info info, duckdb_data_chunk output) {
    seq_kmers_bind_t *bind = (seq_kmers_bind_t *)duckdb_function_get_bind_data(info);
    seq_kmers_init_t *init = (seq_kmers_init_t *)duckdb_function_get_init_data(info);

    if (!bind || !init || bind->seq_len < bind->k) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    idx_t total = bind->seq_len - bind->k + 1;
    if (init->next_pos >= total) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    idx_t vector_size = duckdb_vector_size();
    idx_t remaining = total - init->next_pos;
    idx_t emit = remaining < vector_size ? remaining : vector_size;

    duckdb_vector pos_vec = duckdb_data_chunk_get_vector(output, 0);
    duckdb_vector kmer_vec = duckdb_data_chunk_get_vector(output, 1);
    int64_t *pos_data = (int64_t *)duckdb_vector_get_data(pos_vec);

    char *fwd = NULL;
    char *rev = NULL;
    if (bind->canonical) {
        fwd = (char *)duckdb_malloc((size_t)bind->k + 1);
        rev = (char *)duckdb_malloc((size_t)bind->k + 1);
        if (!fwd || !rev) {
            if (fwd) duckdb_free(fwd);
            if (rev) duckdb_free(rev);
            duckdb_function_set_error(info, "seq_kmers: out of memory");
            return;
        }
        fwd[bind->k] = '\0';
        rev[bind->k] = '\0';
    }

    for (idx_t i = 0; i < emit; i++) {
        idx_t start = init->next_pos + i;
        pos_data[i] = (int64_t)start + 1;

        if (!bind->canonical) {
            duckdb_vector_assign_string_element_len(kmer_vec, i, bind->sequence + start, bind->k);
            continue;
        }

        int valid = 1;
        for (idx_t j = 0; j < bind->k; j++) {
            unsigned char b = (unsigned char)toupper((unsigned char)bind->sequence[start + j]);
            char rc = dna_complement(bind->sequence[start + bind->k - 1 - j]);
            if (!rc || (b != 'A' && b != 'C' && b != 'G' && b != 'T' && b != 'N')) {
                valid = 0;
                break;
            }
            fwd[j] = (char)b;
            rev[j] = rc;
        }

        if (!valid) {
            set_null_at(kmer_vec, i);
            continue;
        }

        const char *chosen = (memcmp(fwd, rev, (size_t)bind->k) <= 0) ? fwd : rev;
        duckdb_vector_assign_string_element_len(kmer_vec, i, chosen, bind->k);
    }

    if (fwd) duckdb_free(fwd);
    if (rev) duckdb_free(rev);

    init->next_pos += emit;
    duckdb_data_chunk_set_size(output, emit);
}

static void register_seq_revcomp_function(duckdb_connection connection) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "seq_revcomp");

    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_set_return_type(fn, varchar_type);
    duckdb_scalar_function_set_function(fn, seq_revcomp_scalar);

    duckdb_register_scalar_function(connection, fn);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_scalar_function(&fn);
}

static void register_seq_canonical_function(duckdb_connection connection) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "seq_canonical");

    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_set_return_type(fn, varchar_type);
    duckdb_scalar_function_set_function(fn, seq_canonical_scalar);

    duckdb_register_scalar_function(connection, fn);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_scalar_function(&fn);
}

static void register_seq_hash_2bit_function(duckdb_connection connection) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "seq_hash_2bit");

    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type ubigint_type = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_set_return_type(fn, ubigint_type);
    duckdb_scalar_function_set_function(fn, seq_hash_2bit_scalar);

    duckdb_register_scalar_function(connection, fn);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&ubigint_type);
    duckdb_destroy_scalar_function(&fn);
}

static void register_seq_encode_4bit_function(duckdb_connection connection) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "seq_encode_4bit");

    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type utinyint_type = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
    duckdb_logical_type list_type = duckdb_create_list_type(utinyint_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_set_return_type(fn, list_type);
    duckdb_scalar_function_set_function(fn, seq_encode_4bit_scalar);

    duckdb_register_scalar_function(connection, fn);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&utinyint_type);
    duckdb_destroy_logical_type(&list_type);
    duckdb_destroy_scalar_function(&fn);
}

static void register_seq_decode_4bit_function(duckdb_connection connection) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "seq_decode_4bit");

    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type utinyint_type = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
    duckdb_logical_type list_type = duckdb_create_list_type(utinyint_type);
    duckdb_scalar_function_add_parameter(fn, list_type);
    duckdb_scalar_function_set_return_type(fn, varchar_type);
    duckdb_scalar_function_set_function(fn, seq_decode_4bit_scalar);

    duckdb_register_scalar_function(connection, fn);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&utinyint_type);
    duckdb_destroy_logical_type(&list_type);
    duckdb_destroy_scalar_function(&fn);
}

static void register_seq_gc_content_function(duckdb_connection connection) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "seq_gc_content");

    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type double_type = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_set_return_type(fn, double_type);
    duckdb_scalar_function_set_function(fn, seq_gc_content_scalar);

    duckdb_register_scalar_function(connection, fn);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&double_type);
    duckdb_destroy_scalar_function(&fn);
}

static void register_sam_flag_predicate_function(duckdb_connection connection,
                                                 const char *name,
                                                 duckdb_scalar_function_t function) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, name);

    duckdb_logical_type usmallint_type = duckdb_create_logical_type(DUCKDB_TYPE_USMALLINT);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_scalar_function_add_parameter(fn, usmallint_type);
    duckdb_scalar_function_set_return_type(fn, bool_type);
    duckdb_scalar_function_set_function(fn, function);

    duckdb_register_scalar_function(connection, fn);

    duckdb_destroy_logical_type(&usmallint_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_scalar_function(&fn);
}

static void register_seq_kmers_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "seq_kmers");

    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);

    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_parameter(tf, bigint_type);
    duckdb_table_function_add_named_parameter(tf, "canonical", bool_type);

    duckdb_table_function_set_bind(tf, seq_kmers_bind);
    duckdb_table_function_set_init(tf, seq_kmers_init);
    duckdb_table_function_set_function(tf, seq_kmers_function);

    duckdb_register_table_function(connection, tf);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_table_function(&tf);
}

void register_kmer_udf_functions(duckdb_connection connection) {
    register_seq_revcomp_function(connection);
    register_seq_canonical_function(connection);
    register_seq_hash_2bit_function(connection);
    register_seq_encode_4bit_function(connection);
    register_seq_decode_4bit_function(connection);
    register_seq_gc_content_function(connection);
    register_seq_kmers_function(connection);
    register_sam_flag_predicate_function(connection, "is_segmented", sam_is_segmented_scalar);
    register_sam_flag_predicate_function(connection, "is_properly_aligned", sam_is_properly_aligned_scalar);
    register_sam_flag_predicate_function(connection, "is_properly_segmented", sam_is_properly_aligned_scalar);
    register_sam_flag_predicate_function(connection, "is_unmapped", sam_is_unmapped_scalar);
    register_sam_flag_predicate_function(connection, "is_mate_unmapped", sam_is_mate_unmapped_scalar);
    register_sam_flag_predicate_function(connection, "is_reverse_complemented", sam_is_reverse_complemented_scalar);
    register_sam_flag_predicate_function(connection, "is_mate_reverse_complemented", sam_is_mate_reverse_complemented_scalar);
    register_sam_flag_predicate_function(connection, "is_first_segment", sam_is_first_segment_scalar);
    register_sam_flag_predicate_function(connection, "is_last_segment", sam_is_last_segment_scalar);
    register_sam_flag_predicate_function(connection, "is_secondary", sam_is_secondary_scalar);
    register_sam_flag_predicate_function(connection, "is_qc_fail", sam_is_qc_fail_scalar);
    register_sam_flag_predicate_function(connection, "is_duplicate", sam_is_duplicate_scalar);
    register_sam_flag_predicate_function(connection, "is_supplementary", sam_is_supplementary_scalar);
}
