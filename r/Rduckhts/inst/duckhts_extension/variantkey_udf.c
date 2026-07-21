#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <stdint.h>
#include <string.h>

#include "include/variantkey_compat.h"

static const char *VK_RANGE_FIELD_NAMES[] = {"min", "max"};

enum {
    DECODE_VK_CHROM_CODE = 0,
    DECODE_VK_POS0 = 1,
    DECODE_VK_REFALT_CODE = 2,
    DECODE_VK_FIELD_COUNT = 3
};

static const char *DECODE_VK_FIELD_NAMES[] = {
    "chrom_code", "pos0", "refalt_code"
};

enum {
    REVERSE_VK_CHROM = 0,
    REVERSE_VK_CHROM_CODE = 1,
    REVERSE_VK_POS = 2,
    REVERSE_VK_POS0 = 3,
    REVERSE_VK_REF = 4,
    REVERSE_VK_ALT = 5,
    REVERSE_VK_REFALT_CODE = 6,
    REVERSE_VK_REVERSIBLE = 7,
    REVERSE_VK_FIELD_COUNT = 8
};

static const char *REVERSE_VK_FIELD_NAMES[] = {
    "chrom", "chrom_code", "pos", "pos0", "ref", "alt", "refalt_code", "reversible"
};

enum {
    REVERSE_RK_CHROM = 0,
    REVERSE_RK_CHROM_CODE = 1,
    REVERSE_RK_START = 2,
    REVERSE_RK_END = 3,
    REVERSE_RK_STRAND = 4,
    REVERSE_RK_STRAND_CODE = 5,
    REVERSE_RK_FIELD_COUNT = 6
};

static const char *REVERSE_RK_FIELD_NAMES[] = {
    "chrom", "chrom_code", "start", "end", "strand", "strand_code"
};

enum {
    DECODE_RK_CHROM_CODE = 0,
    DECODE_RK_START = 1,
    DECODE_RK_END = 2,
    DECODE_RK_STRAND_CODE = 3,
    DECODE_RK_FIELD_COUNT = 4
};

static const char *DECODE_RK_FIELD_NAMES[] = {
    "chrom_code", "start", "end", "strand_code"
};

static inline void set_null_at(duckdb_vector vector, idx_t row) {
    duckdb_vector_ensure_validity_writable(vector);
    duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vector), row);
}

static inline int row_is_valid(duckdb_vector vector, idx_t row) {
    uint64_t *validity = duckdb_vector_get_validity(vector);
    if (!validity) return 1;
    return duckdb_validity_row_is_valid(validity, row);
}

static inline const char *get_string_at(duckdb_vector vector, idx_t row, idx_t *len) {
    duckdb_string_t *data = (duckdb_string_t *)duckdb_vector_get_data(vector);
    duckdb_string_t *val = &data[row];
    *len = duckdb_string_t_length(*val);
    return duckdb_string_t_data(val);
}

static inline int ascii_lower(int c) {
    if (c >= 'A' && c <= 'Z') return c + ('a' - 'A');
    return c;
}

static inline int span_equals_ascii_ci(const char *value, idx_t value_len, const char *literal) {
    idx_t literal_len = (idx_t)strlen(literal);
    if (!value || value_len != literal_len) return 0;
    for (idx_t i = 0; i < value_len; i++) {
        if (ascii_lower((unsigned char)value[i]) != ascii_lower((unsigned char)literal[i])) return 0;
    }
    return 1;
}

static void contig_key_span(const char *value,
                            idx_t value_len,
                            const char **key,
                            idx_t *key_len) {
    static const char mt[] = "MT";
    static const char x[] = "X";
    static const char y[] = "Y";

    *key = value;
    *key_len = value_len;
    if (!value) return;

    if (value_len > 3 &&
        ascii_lower((unsigned char)value[0]) == 'c' &&
        ascii_lower((unsigned char)value[1]) == 'h' &&
        ascii_lower((unsigned char)value[2]) == 'r') {
        *key = value + 3;
        *key_len = value_len - 3;
    }

    if (span_equals_ascii_ci(*key, *key_len, "M") ||
        span_equals_ascii_ci(*key, *key_len, "MT")) {
        *key = mt;
        *key_len = 2;
    } else if (span_equals_ascii_ci(*key, *key_len, "X")) {
        *key = x;
        *key_len = 1;
    } else if (span_equals_ascii_ci(*key, *key_len, "Y")) {
        *key = y;
        *key_len = 1;
    }
}

static void duckhts_contig_key_scalar(duckdb_function_info info,
                                      duckdb_data_chunk input,
                                      duckdb_vector output) {
    duckdb_vector contig_vec = duckdb_data_chunk_get_vector(input, 0);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (idx_t row = 0; row < row_count; row++) {
        idx_t contig_len = 0;
        idx_t key_len = 0;
        const char *contig;
        const char *key;

        if (!row_is_valid(contig_vec, row)) {
            set_null_at(output, row);
            continue;
        }

        contig = get_string_at(contig_vec, row, &contig_len);
        contig_key_span(contig, contig_len, &key, &key_len);
        duckdb_vector_assign_string_element_len(output, row, key, key_len);
    }
}

static int get_variantkey_from_parts(const char *chrom,
                                     idx_t chrom_len,
                                     int64_t pos1,
                                     const char *ref,
                                     idx_t ref_len,
                                     const char *alt,
                                     idx_t alt_len,
                                     uint64_t *out_vk) {
    uint8_t chrom_code;
    uint32_t pos0;
    uint32_t refalt_code;

    if (!chrom || !ref || !alt || !out_vk) return 0;
    if (chrom_len == 0 || ref_len == 0 || alt_len == 0) return 0;

    chrom_code = encode_chrom(chrom, (size_t)chrom_len);
    if (chrom_code == 0) return 0;
    if (!duckhts_variantkey_pos1_to_pos0(pos1, &pos0)) return 0;

    refalt_code = encode_refalt(ref, (size_t)ref_len, alt, (size_t)alt_len);
    *out_vk = encode_variantkey(chrom_code, pos0, refalt_code);
    return 1;
}

static int get_regionkey_from_parts(const char *chrom,
                                    idx_t chrom_len,
                                    int64_t start,
                                    int64_t end,
                                    int64_t strand,
                                    uint64_t *out_rk) {
    uint8_t chrom_code;

    if (!chrom || !out_rk) return 0;
    if (chrom_len == 0) return 0;
    if (!duckhts_regionkey_pos0_valid(start) || !duckhts_regionkey_pos0_valid(end)) return 0;
    if (end < start) return 0;
    if (!(strand == -1 || strand == 0 || strand == 1)) return 0;

    chrom_code = encode_chrom(chrom, (size_t)chrom_len);
    if (chrom_code == 0) return 0;

    *out_rk = encode_regionkey(chrom_code,
                               (uint32_t)start,
                               (uint32_t)end,
                               encode_region_strand((int8_t)strand));
    return 1;
}

static void encode_variantkey_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector chrom_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector pos0_vec = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector refalt_vec = duckdb_data_chunk_get_vector(input, 2);
    uint8_t *chrom_data = (uint8_t *)duckdb_vector_get_data(chrom_vec);
    uint32_t *pos0_data = (uint32_t *)duckdb_vector_get_data(pos0_vec);
    uint32_t *refalt_data = (uint32_t *)duckdb_vector_get_data(refalt_vec);
    uint64_t *out_data = (uint64_t *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(chrom_vec, row) || !row_is_valid(pos0_vec, row) || !row_is_valid(refalt_vec, row) ||
            chrom_data[row] == 0 || pos0_data[row] > DUCKHTS_VARIANTKEY_POS0_MAX || refalt_data[row] > 0x7FFFFFFFU) {
            set_null_at(output, row);
            continue;
        }
        out_data[row] = encode_variantkey(chrom_data[row], pos0_data[row], refalt_data[row]);
    }
}

static void extract_variantkey_chrom_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector vk_vec = duckdb_data_chunk_get_vector(input, 0);
    uint64_t *vk_data = (uint64_t *)duckdb_vector_get_data(vk_vec);
    uint8_t *out_data = (uint8_t *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(vk_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        out_data[row] = extract_variantkey_chrom(vk_data[row]);
    }
}

static void extract_variantkey_pos_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector vk_vec = duckdb_data_chunk_get_vector(input, 0);
    uint64_t *vk_data = (uint64_t *)duckdb_vector_get_data(vk_vec);
    uint32_t *out_data = (uint32_t *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(vk_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        out_data[row] = extract_variantkey_pos(vk_data[row]);
    }
}

static void extract_variantkey_refalt_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector vk_vec = duckdb_data_chunk_get_vector(input, 0);
    uint64_t *vk_data = (uint64_t *)duckdb_vector_get_data(vk_vec);
    uint32_t *out_data = (uint32_t *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(vk_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        out_data[row] = extract_variantkey_refalt(vk_data[row]);
    }
}

static void decode_variantkey_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector vk_vec = duckdb_data_chunk_get_vector(input, 0);
    uint64_t *vk_data = (uint64_t *)duckdb_vector_get_data(vk_vec);
    duckdb_vector child_vecs[DECODE_VK_FIELD_COUNT];
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (int i = 0; i < DECODE_VK_FIELD_COUNT; i++) {
        child_vecs[i] = duckdb_struct_vector_get_child(output, (idx_t)i);
    }

    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(vk_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        ((uint8_t *)duckdb_vector_get_data(child_vecs[DECODE_VK_CHROM_CODE]))[row] = extract_variantkey_chrom(vk_data[row]);
        ((uint32_t *)duckdb_vector_get_data(child_vecs[DECODE_VK_POS0]))[row] = extract_variantkey_pos(vk_data[row]);
        ((uint32_t *)duckdb_vector_get_data(child_vecs[DECODE_VK_REFALT_CODE]))[row] = extract_variantkey_refalt(vk_data[row]);
    }
}

static void variantkey_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector chrom_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector pos_vec = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector ref_vec = duckdb_data_chunk_get_vector(input, 2);
    duckdb_vector alt_vec = duckdb_data_chunk_get_vector(input, 3);
    int64_t *pos_data = (int64_t *)duckdb_vector_get_data(pos_vec);
    uint64_t *out_data = (uint64_t *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (idx_t row = 0; row < row_count; row++) {
        idx_t chrom_len = 0, ref_len = 0, alt_len = 0;
        const char *chrom;
        const char *ref;
        const char *alt;

        if (!row_is_valid(chrom_vec, row) || !row_is_valid(pos_vec, row) ||
            !row_is_valid(ref_vec, row) || !row_is_valid(alt_vec, row)) {
            set_null_at(output, row);
            continue;
        }

        chrom = get_string_at(chrom_vec, row, &chrom_len);
        ref = get_string_at(ref_vec, row, &ref_len);
        alt = get_string_at(alt_vec, row, &alt_len);

        if (!get_variantkey_from_parts(chrom, chrom_len, pos_data[row], ref, ref_len, alt, alt_len, &out_data[row])) {
            set_null_at(output, row);
        }
    }
}

static void variantkey_hex_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector vk_vec = duckdb_data_chunk_get_vector(input, 0);
    uint64_t *vk_data = (uint64_t *)duckdb_vector_get_data(vk_vec);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (idx_t row = 0; row < row_count; row++) {
        char buf[17];
        if (!row_is_valid(vk_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        variantkey_hex(vk_data[row], buf);
        duckdb_vector_assign_string_element_len(output, row, buf, 16);
    }
}

static void parse_variantkey_hex_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector hex_vec = duckdb_data_chunk_get_vector(input, 0);
    uint64_t *out_data = (uint64_t *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (idx_t row = 0; row < row_count; row++) {
        idx_t len = 0;
        const char *hex;
        char buf[17];

        if (!row_is_valid(hex_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        hex = get_string_at(hex_vec, row, &len);
        if (!duckhts_variantkey_hex_is_valid(hex, (size_t)len)) {
            set_null_at(output, row);
            continue;
        }
        memcpy(buf, hex, 16);
        buf[16] = '\0';
        out_data[row] = parse_variantkey_hex(buf);
    }
}

static void reverse_variantkey_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector vk_vec = duckdb_data_chunk_get_vector(input, 0);
    uint64_t *vk_data = (uint64_t *)duckdb_vector_get_data(vk_vec);
    duckdb_vector child_vecs[REVERSE_VK_FIELD_COUNT];
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (int i = 0; i < REVERSE_VK_FIELD_COUNT; i++) {
        child_vecs[i] = duckdb_struct_vector_get_child(output, (idx_t)i);
    }

    for (idx_t row = 0; row < row_count; row++) {
        uint64_t vk;
        uint8_t chrom_code;
        uint32_t pos0;
        uint32_t refalt_code;
        char chrom_buf[8];
        char ref_buf[12];
        char alt_buf[12];
        size_t ref_len = 0;
        size_t alt_len = 0;
        size_t decoded;
        bool reversible;

        if (!row_is_valid(vk_vec, row)) {
            set_null_at(output, row);
            continue;
        }

        vk = vk_data[row];
        chrom_code = extract_variantkey_chrom(vk);
        pos0 = extract_variantkey_pos(vk);
        refalt_code = extract_variantkey_refalt(vk);
        decode_chrom(chrom_code, chrom_buf);
        decoded = decode_refalt(refalt_code, ref_buf, &ref_len, alt_buf, &alt_len);
        reversible = (decoded > 0);

        duckdb_vector_assign_string_element(child_vecs[REVERSE_VK_CHROM], row, chrom_buf);
        ((uint8_t *)duckdb_vector_get_data(child_vecs[REVERSE_VK_CHROM_CODE]))[row] = chrom_code;
        ((int64_t *)duckdb_vector_get_data(child_vecs[REVERSE_VK_POS]))[row] = (int64_t)pos0 + 1;
        ((uint32_t *)duckdb_vector_get_data(child_vecs[REVERSE_VK_POS0]))[row] = pos0;
        if (reversible) {
            duckdb_vector_assign_string_element_len(child_vecs[REVERSE_VK_REF], row, ref_buf, ref_len);
            duckdb_vector_assign_string_element_len(child_vecs[REVERSE_VK_ALT], row, alt_buf, alt_len);
        } else {
            set_null_at(child_vecs[REVERSE_VK_REF], row);
            set_null_at(child_vecs[REVERSE_VK_ALT], row);
        }
        ((uint32_t *)duckdb_vector_get_data(child_vecs[REVERSE_VK_REFALT_CODE]))[row] = refalt_code;
        ((bool *)duckdb_vector_get_data(child_vecs[REVERSE_VK_REVERSIBLE]))[row] = reversible;
    }
}

static void variantkey_range_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector chrom_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector min_vec = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector max_vec = duckdb_data_chunk_get_vector(input, 2);
    int64_t *min_data = (int64_t *)duckdb_vector_get_data(min_vec);
    int64_t *max_data = (int64_t *)duckdb_vector_get_data(max_vec);
    duckdb_vector child_min = duckdb_struct_vector_get_child(output, 0);
    duckdb_vector child_max = duckdb_struct_vector_get_child(output, 1);
    uint64_t *min_out = (uint64_t *)duckdb_vector_get_data(child_min);
    uint64_t *max_out = (uint64_t *)duckdb_vector_get_data(child_max);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (idx_t row = 0; row < row_count; row++) {
        idx_t chrom_len = 0;
        const char *chrom;
        uint8_t chrom_code;
        uint32_t pos0_min;
        uint32_t pos0_max;
        vkrange_t range;

        if (!row_is_valid(chrom_vec, row) || !row_is_valid(min_vec, row) || !row_is_valid(max_vec, row)) {
            set_null_at(output, row);
            continue;
        }

        chrom = get_string_at(chrom_vec, row, &chrom_len);
        chrom_code = encode_chrom(chrom, (size_t)chrom_len);
        if (chrom_code == 0 ||
            !duckhts_variantkey_pos1_to_pos0(min_data[row], &pos0_min) ||
            !duckhts_variantkey_pos1_to_pos0(max_data[row], &pos0_max) ||
            pos0_min > pos0_max) {
            set_null_at(output, row);
            continue;
        }

        variantkey_range(chrom_code, pos0_min, pos0_max, &range);
        min_out[row] = range.min;
        max_out[row] = range.max;
    }
}

static void encode_regionkey_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector chrom_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector start_vec = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector end_vec = duckdb_data_chunk_get_vector(input, 2);
    duckdb_vector strand_vec = duckdb_data_chunk_get_vector(input, 3);
    uint8_t *chrom_data = (uint8_t *)duckdb_vector_get_data(chrom_vec);
    uint32_t *start_data = (uint32_t *)duckdb_vector_get_data(start_vec);
    uint32_t *end_data = (uint32_t *)duckdb_vector_get_data(end_vec);
    uint8_t *strand_data = (uint8_t *)duckdb_vector_get_data(strand_vec);
    uint64_t *out_data = (uint64_t *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(chrom_vec, row) || !row_is_valid(start_vec, row) ||
            !row_is_valid(end_vec, row) || !row_is_valid(strand_vec, row) ||
            chrom_data[row] == 0 || start_data[row] > RK_MAX_POS || end_data[row] > RK_MAX_POS ||
            end_data[row] < start_data[row] || strand_data[row] > 2) {
            set_null_at(output, row);
            continue;
        }
        out_data[row] = encode_regionkey(chrom_data[row], start_data[row], end_data[row], strand_data[row]);
    }
}

static void extract_regionkey_chrom_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector rk_vec = duckdb_data_chunk_get_vector(input, 0);
    uint64_t *rk_data = (uint64_t *)duckdb_vector_get_data(rk_vec);
    uint8_t *out_data = (uint8_t *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(rk_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        out_data[row] = extract_regionkey_chrom(rk_data[row]);
    }
}

static void extract_regionkey_startpos_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector rk_vec = duckdb_data_chunk_get_vector(input, 0);
    uint64_t *rk_data = (uint64_t *)duckdb_vector_get_data(rk_vec);
    uint32_t *out_data = (uint32_t *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(rk_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        out_data[row] = extract_regionkey_startpos(rk_data[row]);
    }
}

static void extract_regionkey_endpos_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector rk_vec = duckdb_data_chunk_get_vector(input, 0);
    uint64_t *rk_data = (uint64_t *)duckdb_vector_get_data(rk_vec);
    uint32_t *out_data = (uint32_t *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(rk_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        out_data[row] = extract_regionkey_endpos(rk_data[row]);
    }
}

static void extract_regionkey_strand_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector rk_vec = duckdb_data_chunk_get_vector(input, 0);
    uint64_t *rk_data = (uint64_t *)duckdb_vector_get_data(rk_vec);
    uint8_t *out_data = (uint8_t *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(rk_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        out_data[row] = extract_regionkey_strand(rk_data[row]);
    }
}

static void decode_regionkey_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector rk_vec = duckdb_data_chunk_get_vector(input, 0);
    uint64_t *rk_data = (uint64_t *)duckdb_vector_get_data(rk_vec);
    duckdb_vector child_vecs[DECODE_RK_FIELD_COUNT];
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (int i = 0; i < DECODE_RK_FIELD_COUNT; i++) {
        child_vecs[i] = duckdb_struct_vector_get_child(output, (idx_t)i);
    }

    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(rk_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        ((uint8_t *)duckdb_vector_get_data(child_vecs[DECODE_RK_CHROM_CODE]))[row] = extract_regionkey_chrom(rk_data[row]);
        ((uint32_t *)duckdb_vector_get_data(child_vecs[DECODE_RK_START]))[row] = extract_regionkey_startpos(rk_data[row]);
        ((uint32_t *)duckdb_vector_get_data(child_vecs[DECODE_RK_END]))[row] = extract_regionkey_endpos(rk_data[row]);
        ((uint8_t *)duckdb_vector_get_data(child_vecs[DECODE_RK_STRAND_CODE]))[row] = extract_regionkey_strand(rk_data[row]);
    }
}

static void regionkey_scalar_impl(duckdb_function_info info,
                                  duckdb_data_chunk input,
                                  duckdb_vector output,
                                  int with_strand) {
    duckdb_vector chrom_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector start_vec = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector end_vec = duckdb_data_chunk_get_vector(input, 2);
    duckdb_vector strand_vec = with_strand ? duckdb_data_chunk_get_vector(input, 3) : NULL;
    int64_t *start_data = (int64_t *)duckdb_vector_get_data(start_vec);
    int64_t *end_data = (int64_t *)duckdb_vector_get_data(end_vec);
    int64_t *strand_data = with_strand ? (int64_t *)duckdb_vector_get_data(strand_vec) : NULL;
    uint64_t *out_data = (uint64_t *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (idx_t row = 0; row < row_count; row++) {
        idx_t chrom_len = 0;
        const char *chrom;
        int64_t strand = 0;

        if (!row_is_valid(chrom_vec, row) || !row_is_valid(start_vec, row) || !row_is_valid(end_vec, row) ||
            (with_strand && !row_is_valid(strand_vec, row))) {
            set_null_at(output, row);
            continue;
        }

        chrom = get_string_at(chrom_vec, row, &chrom_len);
        if (with_strand) strand = strand_data[row];
        if (!get_regionkey_from_parts(chrom, chrom_len, start_data[row], end_data[row], strand, &out_data[row])) {
            set_null_at(output, row);
        }
    }
}

static void regionkey_scalar3(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    regionkey_scalar_impl(info, input, output, 0);
}

static void regionkey_scalar4(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    regionkey_scalar_impl(info, input, output, 1);
}

static void regionkey_hex_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector rk_vec = duckdb_data_chunk_get_vector(input, 0);
    uint64_t *rk_data = (uint64_t *)duckdb_vector_get_data(rk_vec);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (idx_t row = 0; row < row_count; row++) {
        char buf[17];
        if (!row_is_valid(rk_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        regionkey_hex(rk_data[row], buf);
        duckdb_vector_assign_string_element_len(output, row, buf, 16);
    }
}

static void parse_regionkey_hex_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector hex_vec = duckdb_data_chunk_get_vector(input, 0);
    uint64_t *out_data = (uint64_t *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (idx_t row = 0; row < row_count; row++) {
        idx_t len = 0;
        const char *hex;
        char buf[17];

        if (!row_is_valid(hex_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        hex = get_string_at(hex_vec, row, &len);
        if (!duckhts_variantkey_hex_is_valid(hex, (size_t)len)) {
            set_null_at(output, row);
            continue;
        }
        memcpy(buf, hex, 16);
        buf[16] = '\0';
        out_data[row] = parse_regionkey_hex(buf);
    }
}

static void reverse_regionkey_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector rk_vec = duckdb_data_chunk_get_vector(input, 0);
    uint64_t *rk_data = (uint64_t *)duckdb_vector_get_data(rk_vec);
    duckdb_vector child_vecs[REVERSE_RK_FIELD_COUNT];
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (int i = 0; i < REVERSE_RK_FIELD_COUNT; i++) {
        child_vecs[i] = duckdb_struct_vector_get_child(output, (idx_t)i);
    }

    for (idx_t row = 0; row < row_count; row++) {
        uint64_t rk;
        regionkey_rev_t rev;
        uint8_t chrom_code;
        uint8_t strand_code;

        if (!row_is_valid(rk_vec, row)) {
            set_null_at(output, row);
            continue;
        }

        rk = rk_data[row];
        memset(&rev, 0, sizeof(rev));
        reverse_regionkey(rk, &rev);
        chrom_code = extract_regionkey_chrom(rk);
        strand_code = extract_regionkey_strand(rk);

        duckdb_vector_assign_string_element(child_vecs[REVERSE_RK_CHROM], row, rev.chrom);
        ((uint8_t *)duckdb_vector_get_data(child_vecs[REVERSE_RK_CHROM_CODE]))[row] = chrom_code;
        ((uint32_t *)duckdb_vector_get_data(child_vecs[REVERSE_RK_START]))[row] = rev.startpos;
        ((uint32_t *)duckdb_vector_get_data(child_vecs[REVERSE_RK_END]))[row] = rev.endpos;
        ((int8_t *)duckdb_vector_get_data(child_vecs[REVERSE_RK_STRAND]))[row] = rev.strand;
        ((uint8_t *)duckdb_vector_get_data(child_vecs[REVERSE_RK_STRAND_CODE]))[row] = strand_code;
    }
}

static void extend_regionkey_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector rk_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector size_vec = duckdb_data_chunk_get_vector(input, 1);
    uint64_t *rk_data = (uint64_t *)duckdb_vector_get_data(rk_vec);
    int64_t *size_data = (int64_t *)duckdb_vector_get_data(size_vec);
    uint64_t *out_data = (uint64_t *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(rk_vec, row) || !row_is_valid(size_vec, row) || size_data[row] < 0) {
            set_null_at(output, row);
            continue;
        }
        out_data[row] = extend_regionkey(rk_data[row], (uint32_t)size_data[row]);
    }
}

static void are_overlapping_regions_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector chrom_a_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector start_a_vec = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector end_a_vec = duckdb_data_chunk_get_vector(input, 2);
    duckdb_vector chrom_b_vec = duckdb_data_chunk_get_vector(input, 3);
    duckdb_vector start_b_vec = duckdb_data_chunk_get_vector(input, 4);
    duckdb_vector end_b_vec = duckdb_data_chunk_get_vector(input, 5);
    int64_t *start_a_data = (int64_t *)duckdb_vector_get_data(start_a_vec);
    int64_t *end_a_data = (int64_t *)duckdb_vector_get_data(end_a_vec);
    int64_t *start_b_data = (int64_t *)duckdb_vector_get_data(start_b_vec);
    int64_t *end_b_data = (int64_t *)duckdb_vector_get_data(end_b_vec);
    bool *out_data = (bool *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (idx_t row = 0; row < row_count; row++) {
        idx_t chrom_a_len = 0, chrom_b_len = 0;
        const char *chrom_a;
        const char *chrom_b;
        uint8_t chrom_a_code;
        uint8_t chrom_b_code;

        if (!row_is_valid(chrom_a_vec, row) || !row_is_valid(start_a_vec, row) || !row_is_valid(end_a_vec, row) ||
            !row_is_valid(chrom_b_vec, row) || !row_is_valid(start_b_vec, row) || !row_is_valid(end_b_vec, row)) {
            set_null_at(output, row);
            continue;
        }

        chrom_a = get_string_at(chrom_a_vec, row, &chrom_a_len);
        chrom_b = get_string_at(chrom_b_vec, row, &chrom_b_len);
        chrom_a_code = encode_chrom(chrom_a, (size_t)chrom_a_len);
        chrom_b_code = encode_chrom(chrom_b, (size_t)chrom_b_len);

        if (chrom_a_code == 0 || chrom_b_code == 0 ||
            !duckhts_regionkey_pos0_valid(start_a_data[row]) || !duckhts_regionkey_pos0_valid(end_a_data[row]) || end_a_data[row] < start_a_data[row] ||
            !duckhts_regionkey_pos0_valid(start_b_data[row]) || !duckhts_regionkey_pos0_valid(end_b_data[row]) || end_b_data[row] < start_b_data[row]) {
            set_null_at(output, row);
            continue;
        }

        out_data[row] = are_overlapping_regions(chrom_a_code,
                                                (uint32_t)start_a_data[row],
                                                (uint32_t)end_a_data[row],
                                                chrom_b_code,
                                                (uint32_t)start_b_data[row],
                                                (uint32_t)end_b_data[row]) != 0;
    }
}

static void are_overlapping_region_regionkey_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector chrom_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector start_vec = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector end_vec = duckdb_data_chunk_get_vector(input, 2);
    duckdb_vector rk_vec = duckdb_data_chunk_get_vector(input, 3);
    int64_t *start_data = (int64_t *)duckdb_vector_get_data(start_vec);
    int64_t *end_data = (int64_t *)duckdb_vector_get_data(end_vec);
    uint64_t *rk_data = (uint64_t *)duckdb_vector_get_data(rk_vec);
    bool *out_data = (bool *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (idx_t row = 0; row < row_count; row++) {
        idx_t chrom_len = 0;
        const char *chrom;
        uint8_t chrom_code;

        if (!row_is_valid(chrom_vec, row) || !row_is_valid(start_vec, row) ||
            !row_is_valid(end_vec, row) || !row_is_valid(rk_vec, row)) {
            set_null_at(output, row);
            continue;
        }

        chrom = get_string_at(chrom_vec, row, &chrom_len);
        chrom_code = encode_chrom(chrom, (size_t)chrom_len);
        if (chrom_code == 0 ||
            !duckhts_regionkey_pos0_valid(start_data[row]) || !duckhts_regionkey_pos0_valid(end_data[row]) ||
            end_data[row] < start_data[row]) {
            set_null_at(output, row);
            continue;
        }

        out_data[row] = are_overlapping_region_regionkey(chrom_code,
                                                         (uint32_t)start_data[row],
                                                         (uint32_t)end_data[row],
                                                         rk_data[row]) != 0;
    }
}

static void are_overlapping_regionkeys_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector a_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector b_vec = duckdb_data_chunk_get_vector(input, 1);
    uint64_t *a_data = (uint64_t *)duckdb_vector_get_data(a_vec);
    uint64_t *b_data = (uint64_t *)duckdb_vector_get_data(b_vec);
    bool *out_data = (bool *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    (void)info;
    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(a_vec, row) || !row_is_valid(b_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        out_data[row] = are_overlapping_regionkeys(a_data[row], b_data[row]) != 0;
    }
}

static void register_decode_variantkey_function(duckdb_connection connection) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_logical_type ubigint_type = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    duckdb_logical_type fields[DECODE_VK_FIELD_COUNT];
    duckdb_logical_type struct_type;

    duckdb_scalar_function_set_name(fn, "decode_variantkey");
    duckdb_scalar_function_add_parameter(fn, ubigint_type);
    fields[DECODE_VK_CHROM_CODE] = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
    fields[DECODE_VK_POS0] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
    fields[DECODE_VK_REFALT_CODE] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
    struct_type = duckdb_create_struct_type(fields, DECODE_VK_FIELD_NAMES, DECODE_VK_FIELD_COUNT);
    duckdb_scalar_function_set_return_type(fn, struct_type);
    duckdb_scalar_function_set_special_handling(fn);
    duckdb_scalar_function_set_function(fn, decode_variantkey_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    duckdb_destroy_logical_type(&ubigint_type);
    for (int i = 0; i < DECODE_VK_FIELD_COUNT; i++) duckdb_destroy_logical_type(&fields[i]);
    duckdb_destroy_logical_type(&struct_type);
}

static void register_variantkey_range_function(duckdb_connection connection) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type ubigint_type = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    duckdb_logical_type fields[2];
    duckdb_logical_type struct_type;

    duckdb_scalar_function_set_name(fn, "variantkey_range");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    fields[0] = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    fields[1] = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    struct_type = duckdb_create_struct_type(fields, VK_RANGE_FIELD_NAMES, 2);
    duckdb_scalar_function_set_return_type(fn, struct_type);
    duckdb_scalar_function_set_special_handling(fn);
    duckdb_scalar_function_set_function(fn, variantkey_range_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    duckdb_destroy_logical_type(&fields[0]);
    duckdb_destroy_logical_type(&fields[1]);
    duckdb_destroy_logical_type(&struct_type);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&ubigint_type);
}

static void register_reverse_variantkey_function(duckdb_connection connection) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_logical_type ubigint_type = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type utinyint_type = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
    duckdb_logical_type uinteger_type = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type fields[REVERSE_VK_FIELD_COUNT];
    duckdb_logical_type struct_type;

    duckdb_scalar_function_set_name(fn, "reverse_variantkey");
    duckdb_scalar_function_add_parameter(fn, ubigint_type);
    fields[REVERSE_VK_CHROM] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[REVERSE_VK_CHROM_CODE] = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
    fields[REVERSE_VK_POS] = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    fields[REVERSE_VK_POS0] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
    fields[REVERSE_VK_REF] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[REVERSE_VK_ALT] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[REVERSE_VK_REFALT_CODE] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
    fields[REVERSE_VK_REVERSIBLE] = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    struct_type = duckdb_create_struct_type(fields, REVERSE_VK_FIELD_NAMES, REVERSE_VK_FIELD_COUNT);
    duckdb_scalar_function_set_return_type(fn, struct_type);
    duckdb_scalar_function_set_special_handling(fn);
    duckdb_scalar_function_set_function(fn, reverse_variantkey_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&utinyint_type);
    duckdb_destroy_logical_type(&uinteger_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_logical_type(&ubigint_type);
    for (int i = 0; i < REVERSE_VK_FIELD_COUNT; i++) duckdb_destroy_logical_type(&fields[i]);
    duckdb_destroy_logical_type(&struct_type);
}

static void register_decode_regionkey_function(duckdb_connection connection) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_logical_type ubigint_type = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    duckdb_logical_type fields[DECODE_RK_FIELD_COUNT];
    duckdb_logical_type struct_type;

    duckdb_scalar_function_set_name(fn, "decode_regionkey");
    duckdb_scalar_function_add_parameter(fn, ubigint_type);
    fields[DECODE_RK_CHROM_CODE] = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
    fields[DECODE_RK_START] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
    fields[DECODE_RK_END] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
    fields[DECODE_RK_STRAND_CODE] = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
    struct_type = duckdb_create_struct_type(fields, DECODE_RK_FIELD_NAMES, DECODE_RK_FIELD_COUNT);
    duckdb_scalar_function_set_return_type(fn, struct_type);
    duckdb_scalar_function_set_special_handling(fn);
    duckdb_scalar_function_set_function(fn, decode_regionkey_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    duckdb_destroy_logical_type(&ubigint_type);
    for (int i = 0; i < DECODE_RK_FIELD_COUNT; i++) duckdb_destroy_logical_type(&fields[i]);
    duckdb_destroy_logical_type(&struct_type);
}

static void register_reverse_regionkey_function(duckdb_connection connection) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_logical_type ubigint_type = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    duckdb_logical_type fields[REVERSE_RK_FIELD_COUNT];
    duckdb_logical_type struct_type;

    duckdb_scalar_function_set_name(fn, "reverse_regionkey");
    duckdb_scalar_function_add_parameter(fn, ubigint_type);
    fields[REVERSE_RK_CHROM] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[REVERSE_RK_CHROM_CODE] = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
    fields[REVERSE_RK_START] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
    fields[REVERSE_RK_END] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
    fields[REVERSE_RK_STRAND] = duckdb_create_logical_type(DUCKDB_TYPE_TINYINT);
    fields[REVERSE_RK_STRAND_CODE] = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
    struct_type = duckdb_create_struct_type(fields, REVERSE_RK_FIELD_NAMES, REVERSE_RK_FIELD_COUNT);
    duckdb_scalar_function_set_return_type(fn, struct_type);
    duckdb_scalar_function_set_special_handling(fn);
    duckdb_scalar_function_set_function(fn, reverse_regionkey_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    duckdb_destroy_logical_type(&ubigint_type);
    for (int i = 0; i < REVERSE_RK_FIELD_COUNT; i++) duckdb_destroy_logical_type(&fields[i]);
    duckdb_destroy_logical_type(&struct_type);
}

void register_variantkey_functions(duckdb_connection connection) {
    duckdb_scalar_function fn;
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type ubigint_type = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type utinyint_type = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
    duckdb_logical_type uinteger_type = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "duckhts_contig_key");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_set_return_type(fn, varchar_type);
    duckdb_scalar_function_set_function(fn, duckhts_contig_key_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "encode_variantkey");
    duckdb_scalar_function_add_parameter(fn, utinyint_type);
    duckdb_scalar_function_add_parameter(fn, uinteger_type);
    duckdb_scalar_function_add_parameter(fn, uinteger_type);
    duckdb_scalar_function_set_return_type(fn, ubigint_type);
    duckdb_scalar_function_set_function(fn, encode_variantkey_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "extract_variantkey_chrom");
    duckdb_scalar_function_add_parameter(fn, ubigint_type);
    duckdb_scalar_function_set_return_type(fn, utinyint_type);
    duckdb_scalar_function_set_function(fn, extract_variantkey_chrom_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "extract_variantkey_pos");
    duckdb_scalar_function_add_parameter(fn, ubigint_type);
    duckdb_scalar_function_set_return_type(fn, uinteger_type);
    duckdb_scalar_function_set_function(fn, extract_variantkey_pos_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "extract_variantkey_refalt");
    duckdb_scalar_function_add_parameter(fn, ubigint_type);
    duckdb_scalar_function_set_return_type(fn, uinteger_type);
    duckdb_scalar_function_set_function(fn, extract_variantkey_refalt_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    register_decode_variantkey_function(connection);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "variantkey");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_set_return_type(fn, ubigint_type);
    duckdb_scalar_function_set_function(fn, variantkey_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "variantkey_hex");
    duckdb_scalar_function_add_parameter(fn, ubigint_type);
    duckdb_scalar_function_set_return_type(fn, varchar_type);
    duckdb_scalar_function_set_function(fn, variantkey_hex_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "parse_variantkey_hex");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_set_return_type(fn, ubigint_type);
    duckdb_scalar_function_set_function(fn, parse_variantkey_hex_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    register_reverse_variantkey_function(connection);
    register_variantkey_range_function(connection);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "encode_regionkey");
    duckdb_scalar_function_add_parameter(fn, utinyint_type);
    duckdb_scalar_function_add_parameter(fn, uinteger_type);
    duckdb_scalar_function_add_parameter(fn, uinteger_type);
    duckdb_scalar_function_add_parameter(fn, utinyint_type);
    duckdb_scalar_function_set_return_type(fn, ubigint_type);
    duckdb_scalar_function_set_function(fn, encode_regionkey_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "extract_regionkey_chrom");
    duckdb_scalar_function_add_parameter(fn, ubigint_type);
    duckdb_scalar_function_set_return_type(fn, utinyint_type);
    duckdb_scalar_function_set_function(fn, extract_regionkey_chrom_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "extract_regionkey_startpos");
    duckdb_scalar_function_add_parameter(fn, ubigint_type);
    duckdb_scalar_function_set_return_type(fn, uinteger_type);
    duckdb_scalar_function_set_function(fn, extract_regionkey_startpos_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "extract_regionkey_endpos");
    duckdb_scalar_function_add_parameter(fn, ubigint_type);
    duckdb_scalar_function_set_return_type(fn, uinteger_type);
    duckdb_scalar_function_set_function(fn, extract_regionkey_endpos_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "extract_regionkey_strand");
    duckdb_scalar_function_add_parameter(fn, ubigint_type);
    duckdb_scalar_function_set_return_type(fn, utinyint_type);
    duckdb_scalar_function_set_function(fn, extract_regionkey_strand_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    register_decode_regionkey_function(connection);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "regionkey");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_set_return_type(fn, ubigint_type);
    duckdb_scalar_function_set_function(fn, regionkey_scalar3);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "regionkey");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_set_return_type(fn, ubigint_type);
    duckdb_scalar_function_set_function(fn, regionkey_scalar4);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "regionkey_hex");
    duckdb_scalar_function_add_parameter(fn, ubigint_type);
    duckdb_scalar_function_set_return_type(fn, varchar_type);
    duckdb_scalar_function_set_function(fn, regionkey_hex_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "parse_regionkey_hex");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_set_return_type(fn, ubigint_type);
    duckdb_scalar_function_set_function(fn, parse_regionkey_hex_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    register_reverse_regionkey_function(connection);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "extend_regionkey");
    duckdb_scalar_function_add_parameter(fn, ubigint_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_set_return_type(fn, ubigint_type);
    duckdb_scalar_function_set_function(fn, extend_regionkey_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "are_overlapping_regions");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_set_return_type(fn, bool_type);
    duckdb_scalar_function_set_function(fn, are_overlapping_regions_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "are_overlapping_region_regionkey");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, ubigint_type);
    duckdb_scalar_function_set_return_type(fn, bool_type);
    duckdb_scalar_function_set_function(fn, are_overlapping_region_regionkey_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "are_overlapping_regionkeys");
    duckdb_scalar_function_add_parameter(fn, ubigint_type);
    duckdb_scalar_function_add_parameter(fn, ubigint_type);
    duckdb_scalar_function_set_return_type(fn, bool_type);
    duckdb_scalar_function_set_function(fn, are_overlapping_regionkeys_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&ubigint_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_logical_type(&utinyint_type);
    duckdb_destroy_logical_type(&uinteger_type);
}
