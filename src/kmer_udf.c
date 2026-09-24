#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "duckhts_cigar.h"
#include "duckhts_simd_internal.h"
#include "seq_encoding.h"

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

static const char *SAM_FLAG_FIELD_NAMES[] = {
    "is_paired",
    "is_proper_pair",
    "is_unmapped",
    "is_next_segment_unmapped",
    "is_reverse_complemented",
    "is_next_segment_reverse_complemented",
    "is_first_segment",
    "is_last_segment",
    "is_secondary",
    "is_qc_fail",
    "is_duplicate",
    "is_supplementary"
};

static const uint16_t SAM_FLAG_FIELD_MASKS[] = {
    SAM_FLAG_PAIRED,
    SAM_FLAG_PROPER_PAIR,
    SAM_FLAG_UNMAPPED,
    SAM_FLAG_MATE_UNMAPPED,
    SAM_FLAG_REVERSE,
    SAM_FLAG_MATE_REVERSE,
    SAM_FLAG_READ1,
    SAM_FLAG_READ2,
    SAM_FLAG_SECONDARY,
    SAM_FLAG_QCFAIL,
    SAM_FLAG_DUPLICATE,
    SAM_FLAG_SUPPLEMENTARY
};

enum { SAM_FLAG_FIELD_COUNT = (int)(sizeof(SAM_FLAG_FIELD_MASKS) / sizeof(SAM_FLAG_FIELD_MASKS[0])) };

typedef struct {
    uint32_t seen_ops;
    int has_soft_clip;
    int has_hard_clip;
    int64_t left_soft_clip;
    int64_t right_soft_clip;
    int64_t query_length;
    int64_t aligned_query_length;
    int64_t reference_length;
} cigar_metrics_t;

enum {
    CIGAR_METRIC_HAS_SOFT_CLIP = 1,
    CIGAR_METRIC_HAS_HARD_CLIP = 2,
    CIGAR_METRIC_LEFT_SOFT_CLIP = 3,
    CIGAR_METRIC_RIGHT_SOFT_CLIP = 4,
    CIGAR_METRIC_QUERY_LENGTH = 5,
    CIGAR_METRIC_ALIGNED_QUERY_LENGTH = 6,
    CIGAR_METRIC_REFERENCE_LENGTH = 7
};

/* Static lifetime: registration and execution borrow the same metric metadata. */
typedef struct {
    const char *name;
    int id;
    duckdb_type return_type;
} cigar_metric_definition_t;

static const cigar_metric_definition_t CIGAR_METRICS[] = {
    {"cigar_has_soft_clip", CIGAR_METRIC_HAS_SOFT_CLIP, DUCKDB_TYPE_BOOLEAN},
    {"cigar_has_hard_clip", CIGAR_METRIC_HAS_HARD_CLIP, DUCKDB_TYPE_BOOLEAN},
    {"cigar_left_soft_clip", CIGAR_METRIC_LEFT_SOFT_CLIP, DUCKDB_TYPE_BIGINT},
    {"cigar_right_soft_clip", CIGAR_METRIC_RIGHT_SOFT_CLIP, DUCKDB_TYPE_BIGINT},
    {"cigar_query_length", CIGAR_METRIC_QUERY_LENGTH, DUCKDB_TYPE_BIGINT},
    {"cigar_aligned_query_length", CIGAR_METRIC_ALIGNED_QUERY_LENGTH, DUCKDB_TYPE_BIGINT},
    {"cigar_reference_length", CIGAR_METRIC_REFERENCE_LENGTH, DUCKDB_TYPE_BIGINT}
};

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

/* iupac_to_4bit and bit4_to_iupac removed — use seq_text_to_nt16() and
   seq_nt16_to_text() from seq_encoding.h instead. */

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

/* Both representations use the same complete scan. No-CIGAR succeeds with
   seen_ops == 0 so wrappers can distinguish metric NULL from presence false.
   On failure the caller must discard all accumulated metrics. */
static inline duckhts_cigar_status_t cigar_scan_metrics(duckhts_cigar_cursor_t *cursor,
                                                      cigar_metrics_t *metrics) {
    duckhts_cigar_op_t op;
    duckhts_cigar_status_t status;

    memset(metrics, 0, sizeof(*metrics));
    while ((status = duckhts_cigar_next(cursor, &op)) > 0) {
        if (op.type == 3) {
            /* Aligned query is a subset of the already checked query span. */
            metrics->aligned_query_length += op.length;
        }
        if (metrics->seen_ops == 0 && op.code == BAM_CSOFT_CLIP) {
            metrics->left_soft_clip = op.length;
        }
        metrics->right_soft_clip = op.code == BAM_CSOFT_CLIP ? op.length : 0;
        metrics->seen_ops |= UINT32_C(1) << op.code;
    }
    if (status < 0) {
        return status;
    }
    metrics->query_length = cursor->query_span;
    metrics->reference_length = cursor->reference_span;
    metrics->has_soft_clip = (metrics->seen_ops & (UINT32_C(1) << BAM_CSOFT_CLIP)) != 0;
    metrics->has_hard_clip = (metrics->seen_ops & (UINT32_C(1) << BAM_CHARD_CLIP)) != 0;
    return DUCKHTS_CIGAR_OK;
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

/* nt16 code -> complemented nt16 code; 0xFF marks a code the text path would
   reject (dna_complement returns 0 for anything but A/C/G/T/N). A=1<->T=8,
   C=2<->G=4, N=15 self-complements; '=' (0) and IUPAC codes are invalid. */
static const uint8_t nt16_complement[16] = {
    0xFF, 8, 4, 0xFF, 2, 0xFF, 0xFF, 0xFF, 1, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 15};
/* nt16 code -> uppercase base char, for canonical's char-order comparison
   (nt16 code order is A<C<G<T<N but char order is A<C<G<N<T, so we must
   compare decoded chars, not codes). Only A/C/G/T/N are reachable. */
static const unsigned char nt16_upper_char[16] = {
    0, 'A', 'C', 0, 'G', 0, 0, 0, 'T', 0, 0, 0, 0, 0, 0, 'N'};

/* Validate that every nt16 code in [off, off+len) is A/C/G/T/N (what the
   VARCHAR path tolerates) and non-NULL. Returns 1 if valid. */
static int nt16_all_acgtn(duckdb_vector child_vec, const uint8_t *child_data,
                          idx_t off, idx_t len) {
    for (idx_t i = 0; i < len; i++) {
        idx_t ci = off + i;
        if (!row_is_valid(child_vec, ci)) return 0;
        uint8_t code = child_data[ci];
        if (code >= 16 || nt16_complement[code] == 0xFF) return 0;
    }
    return 1;
}

/* nt16 overload of seq_revcomp: reverse-complement a UTINYINT[] of nt16 codes,
   returning nt16 codes. Bit-identical to seq_revcomp over the decoded text. */
static void seq_revcomp_nt16_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector seq_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_list_entry *in_list = (duckdb_list_entry *)duckdb_vector_get_data(seq_vec);
    duckdb_vector in_child = duckdb_list_vector_get_child(seq_vec);
    uint8_t *in_data = (uint8_t *)duckdb_vector_get_data(in_child);
    duckdb_list_entry *out_list = (duckdb_list_entry *)duckdb_vector_get_data(output);
    duckdb_vector out_child = duckdb_list_vector_get_child(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);
    idx_t child_offset = duckdb_list_vector_get_size(output);

    for (idx_t row = 0; row < row_count; row++) {
        out_list[row].offset = child_offset;
        out_list[row].length = 0;
        if (!row_is_valid(seq_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        duckdb_list_entry e = in_list[row];
        if (!nt16_all_acgtn(in_child, in_data, e.offset, e.length)) {
            set_null_at(output, row);
            continue;
        }
        if (duckdb_list_vector_reserve(output, child_offset + e.length) != DuckDBSuccess ||
            duckdb_list_vector_set_size(output, child_offset + e.length) != DuckDBSuccess) {
            duckdb_scalar_function_set_error(info, "seq_revcomp: failed to grow list storage");
            return;
        }
        uint8_t *out_data = (uint8_t *)duckdb_vector_get_data(out_child);
        for (idx_t i = 0; i < e.length; i++) {
            out_data[child_offset + i] = nt16_complement[in_data[e.offset + e.length - 1 - i]];
        }
        out_list[row].length = e.length;
        child_offset += e.length;
    }
}

/* nt16 overload of seq_canonical: min(seq, revcomp(seq)) by decoded-char order,
   returning nt16 codes. Bit-identical to seq_canonical over the decoded text. */
static void seq_canonical_nt16_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector seq_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_list_entry *in_list = (duckdb_list_entry *)duckdb_vector_get_data(seq_vec);
    duckdb_vector in_child = duckdb_list_vector_get_child(seq_vec);
    uint8_t *in_data = (uint8_t *)duckdb_vector_get_data(in_child);
    duckdb_list_entry *out_list = (duckdb_list_entry *)duckdb_vector_get_data(output);
    duckdb_vector out_child = duckdb_list_vector_get_child(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);
    idx_t child_offset = duckdb_list_vector_get_size(output);

    for (idx_t row = 0; row < row_count; row++) {
        out_list[row].offset = child_offset;
        out_list[row].length = 0;
        if (!row_is_valid(seq_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        duckdb_list_entry e = in_list[row];
        if (!nt16_all_acgtn(in_child, in_data, e.offset, e.length)) {
            set_null_at(output, row);
            continue;
        }
        /* Choose forward vs reverse-complement by first differing decoded char;
           ties keep forward, matching the text path's memcmp(fwd, rev) <= 0. */
        int use_rev = 0;
        for (idx_t i = 0; i < e.length; i++) {
            unsigned char fc = nt16_upper_char[in_data[e.offset + i]];
            unsigned char rc = nt16_upper_char[nt16_complement[in_data[e.offset + e.length - 1 - i]]];
            if (fc != rc) {
                use_rev = (rc < fc);
                break;
            }
        }
        if (duckdb_list_vector_reserve(output, child_offset + e.length) != DuckDBSuccess ||
            duckdb_list_vector_set_size(output, child_offset + e.length) != DuckDBSuccess) {
            duckdb_scalar_function_set_error(info, "seq_canonical: failed to grow list storage");
            return;
        }
        uint8_t *out_data = (uint8_t *)duckdb_vector_get_data(out_child);
        for (idx_t i = 0; i < e.length; i++) {
            out_data[child_offset + i] = use_rev
                ? nt16_complement[in_data[e.offset + e.length - 1 - i]]
                : in_data[e.offset + i];
        }
        out_list[row].length = e.length;
        child_offset += e.length;
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

/* nt16 (UTINYINT[]) overload of seq_hash_2bit. Packs A/C/G/T -> 2 bits/base
   into a uint64, like the text path; sequences longer than 32 bases or with a
   non-ACGT (N/IUPAC/=) code return NULL. nt16 code -> 2-bit: A=1->0, C=2->1,
   G=4->2, T=8->3; every other code is -1 (invalid), matching dna_to_2bit. */
static void seq_hash_2bit_nt16_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    (void)info;
    static const int8_t nt16_2bit[16] = {-1, 0, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1};
    duckdb_vector seq_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_list_entry *list_data = (duckdb_list_entry *)duckdb_vector_get_data(seq_vec);
    duckdb_vector child_vec = duckdb_list_vector_get_child(seq_vec);
    uint8_t *child_data = (uint8_t *)duckdb_vector_get_data(child_vec);
    uint64_t *out_data = (uint64_t *)duckdb_vector_get_data(output);
    uint64_t *child_validity = duckdb_vector_get_validity(child_vec);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(seq_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        duckdb_list_entry entry = list_data[row];
        if (entry.length > 32) {
            set_null_at(output, row);
            continue;
        }
        uint64_t h = 0;
        int valid = 1;
        for (idx_t i = 0; i < entry.length; i++) {
            idx_t ci = entry.offset + i;
            if (child_validity && !duckdb_validity_row_is_valid(child_validity, ci)) {
                valid = 0;
                break;
            }
            uint8_t c = child_data[ci];
            int code = (c < 16) ? nt16_2bit[c] : -1;
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
        /* Permissive encoding via shared helper — unknown chars → N (15),
           matching htslib seq_nt16_table behavior exactly. */
        seq_text_to_nt16(seq, child_data + child_offset, len);

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

        /* Check for NULL child elements, then decode via shared helper.
           seq_nt16_to_text returns -1 if any code > 15. */
        int has_null_child = 0;
        for (idx_t i = 0; i < entry.length; i++) {
            if (!row_is_valid(child_vec, entry.offset + i)) {
                has_null_child = 1;
                break;
            }
        }

        if (has_null_child ||
            seq_nt16_to_text(child_data + entry.offset, decoded, entry.length) < 0) {
            set_null_at(output, row);
            duckdb_free(decoded);
            continue;
        }

        duckdb_vector_assign_string_element_len(output, row, decoded, entry.length);
        duckdb_free(decoded);
    }
}

static void seq_gc_content_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    (void)info;
    duckdb_vector seq_vec = duckdb_data_chunk_get_vector(input, 0);
    double *out_data = (double *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);
    const duckhts_simd_dispatch_table_t *simd_table = duckhts_simd_dispatch_snapshot();

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

        duckhts_simd_base_counts_t counts;
        duckhts_simd_base_counts_with_table(simd_table, seq, (size_t)len, &counts);

        if (counts.invalid || counts.called == 0) {
            set_null_at(output, row);
            continue;
        }
        out_data[row] = (double)counts.gc / (double)counts.called;
    }
}

/* nt16 overload of seq_gc_content. `read_bam(sequence_encoding := 'nt16')`
   yields SEQ as LIST(UTINYINT), one htslib nt16 code per base, so we classify
   the codes directly instead of round-tripping through text. htslib's nt16
   encoding is fixed (A=1, C=2, G=4, T=8), and `called`/`gc` follow the same
   `seq_nt16_int` authority the VARCHAR path reconciles to: a base is called
   iff its code is A/C/G/T, GC iff C/G. Result is bit-identical to
   seq_gc_content over the decoded text of the same read. */
static void seq_gc_content_nt16_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    (void)info;
    /* nt16 code (0..15) -> class, matching the VARCHAR path's alphabet exactly:
       the text kernel counts A/C/G/T, ignores N, and marks ANY other byte
       invalid (making seq_gc_content return NULL). So only A/C/G/T/N are
       tolerated here; `=` (code 0) and the IUPAC ambiguity codes are invalid.
         0 = invalid (=, M, R, S, V, W, Y, H, K, D, B)
         1 = N        (valid, not called)
         2 = A or T   (called)
         3 = C or G   (called + gc) */
    static const uint8_t nt16_class[16] = {
        /* = */ 0, /* A */ 2, /* C */ 3, /* M */ 0,
        /* G */ 3, /* R */ 0, /* S */ 0, /* V */ 0,
        /* T */ 2, /* W */ 0, /* Y */ 0, /* H */ 0,
        /* K */ 0, /* D */ 0, /* B */ 0, /* N */ 1};

    duckdb_vector seq_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_list_entry *list_data = (duckdb_list_entry *)duckdb_vector_get_data(seq_vec);
    duckdb_vector child_vec = duckdb_list_vector_get_child(seq_vec);
    uint8_t *child_data = (uint8_t *)duckdb_vector_get_data(child_vec);
    /* Fetch the child validity mask once: it is loop-invariant, and BAM SEQ has
       no NULL bases, so the common case is a NULL mask (handled by the
       vectorized kernel below) and the rare per-base branch never runs. */
    uint64_t *child_validity = duckdb_vector_get_validity(child_vec);
    double *out_data = (double *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);
    /* Snapshot the SIMD dispatch table once per chunk for a stable view. */
    const duckhts_simd_dispatch_table_t *simd = duckhts_simd_dispatch_snapshot();

    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(seq_vec, row)) {
            set_null_at(output, row);
            continue;
        }

        duckdb_list_entry entry = list_data[row];
        if (entry.length == 0) {
            set_null_at(output, row);
            continue;
        }

        uint64_t gc;
        uint64_t called;
        int invalid;

        if (!child_validity) {
            /* Common path: no NULL bases, so classify the row's contiguous nt16
               codes through the dispatched (scalar/AVX2/...) kernel. */
            duckhts_simd_base_counts_t counts;
            duckhts_simd_nt16_gc_counts_with_table(simd, child_data + entry.offset,
                                                   (size_t)entry.length, &counts);
            gc = counts.gc;
            called = counts.called;
            invalid = counts.invalid;
        } else {
            /* Rare path: a child element may be NULL. Classify scalar with a
               per-base validity check; a NULL base is invalid -> NULL result. */
            gc = 0;
            called = 0;
            invalid = 0;
            for (idx_t i = 0; i < entry.length; i++) {
                idx_t ci = entry.offset + i;
                if (!duckdb_validity_row_is_valid(child_validity, ci)) {
                    invalid = 1;
                    break;
                }
                uint8_t code = child_data[ci];
                uint8_t klass = (code < 16) ? nt16_class[code] : 0;
                if (klass == 0) {
                    invalid = 1;
                    break;
                }
                if (klass >= 2) {
                    called++;
                    if (klass == 3) gc++;
                }
            }
        }

        if (invalid || called == 0) {
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

static void sam_flag_predicate_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    uint16_t mask = (uint16_t)(uintptr_t)duckdb_scalar_function_get_extra_info(info);
    sam_flag_scalar(info, input, output, mask);
}

static void sam_is_forward_aligned_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
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
        uint16_t flag = (uint16_t)flag_value;
        if ((flag & SAM_FLAG_UNMAPPED) != 0) {
            set_null_at(output, row);
            continue;
        }
        out_data[row] = ((flag & SAM_FLAG_REVERSE) == 0);
    }
}

static void sam_flag_has_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    (void)info;
    duckdb_vector flag_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector mask_vec = duckdb_data_chunk_get_vector(input, 1);
    bool *out_data = (bool *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(flag_vec, row) || !row_is_valid(mask_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        int64_t flag_value = get_int64_at(flag_vec, row);
        int64_t mask_value = get_int64_at(mask_vec, row);
        if (flag_value < 0 || flag_value > 0xffff || mask_value < 0 || mask_value > 0xffff) {
            set_null_at(output, row);
            continue;
        }
        out_data[row] = ((((uint16_t)flag_value) & ((uint16_t)mask_value)) != 0);
    }
}

static void sam_flag_bits_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    (void)info;
    duckdb_vector flag_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector child_vecs[SAM_FLAG_FIELD_COUNT];
    bool *child_data[SAM_FLAG_FIELD_COUNT];
    idx_t row_count = duckdb_data_chunk_get_size(input);

    for (int i = 0; i < SAM_FLAG_FIELD_COUNT; i++) {
        child_vecs[i] = duckdb_struct_vector_get_child(output, (idx_t)i);
        child_data[i] = (bool *)duckdb_vector_get_data(child_vecs[i]);
    }

    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(flag_vec, row)) {
            set_null_at(output, row);
            for (int i = 0; i < SAM_FLAG_FIELD_COUNT; i++) {
                set_null_at(child_vecs[i], row);
            }
            continue;
        }

        int64_t flag_value = get_int64_at(flag_vec, row);
        if (flag_value < 0 || flag_value > 0xffff) {
            set_null_at(output, row);
            for (int i = 0; i < SAM_FLAG_FIELD_COUNT; i++) {
                set_null_at(child_vecs[i], row);
            }
            continue;
        }

        uint16_t flag_bits = (uint16_t)flag_value;
        for (int i = 0; i < SAM_FLAG_FIELD_COUNT; i++) {
            child_data[i][row] = ((flag_bits & SAM_FLAG_FIELD_MASKS[i]) != 0);
        }
    }
}

static int cigar_requested_op(const char *text, idx_t length) {
    if (length != 1) {
        return -1;
    }
    unsigned char op = (unsigned char)text[0];
    if (op >= 'a' && op <= 'z') {
        op = (unsigned char)(op - 'a' + 'A');
    }
    return bam_cigar_table[op];
}

static duckdb_vector cigar_strict_vector(duckdb_data_chunk input, idx_t index) {
    return duckdb_data_chunk_get_column_count(input) > index
        ? duckdb_data_chunk_get_vector(input, index) : NULL;
}

/* Call only after all top-level arguments, including strict, are non-NULL.
   A true return means execution must stop after setting the DuckDB error. */
static bool cigar_strict_error(duckdb_function_info info, const char *name,
                              duckdb_vector strict_vec, idx_t row,
                              duckhts_cigar_status_t status,
                              const duckhts_cigar_cursor_t *cursor) {
    if (status >= 0 || !strict_vec || !((bool *)duckdb_vector_get_data(strict_vec))[row]) {
        return false;
    }
    const char *reason;
    switch (status) {
    case DUCKHTS_CIGAR_MISSING_LENGTH:
        reason = "missing operation length";
        break;
    case DUCKHTS_CIGAR_ZERO_LENGTH:
        reason = "zero operation length";
        break;
    case DUCKHTS_CIGAR_MISSING_OP:
        reason = "trailing length without an operator";
        break;
    case DUCKHTS_CIGAR_UNSUPPORTED_OP:
        reason = "unsupported CIGAR operator";
        break;
    case DUCKHTS_CIGAR_LENGTH_OVERFLOW:
        reason = "operation length exceeds BIGINT";
        break;
    case DUCKHTS_CIGAR_QUERY_OVERFLOW:
        reason = "query span exceeds BIGINT";
        break;
    case DUCKHTS_CIGAR_REFERENCE_OVERFLOW:
        reason = "reference span exceeds BIGINT";
        break;
    case DUCKHTS_CIGAR_NULL_CHILD:
        reason = "NULL packed operation";
        break;
    case DUCKHTS_CIGAR_INPUT_TOO_LARGE:
        reason = "packed input exceeds addressable size";
        break;
    case DUCKHTS_CIGAR_INVALID_REQUESTED_OP:
        reason = "requested operator must be one of M, I, D, N, S, H, P, =, X";
        break;
    case DUCKHTS_CIGAR_POSITION_OVERFLOW:
        reason = "position plus reference span exceeds BIGINT";
        break;
    default:
        reason = "invalid CIGAR";
        break;
    }
    char message[256];
    if (cursor && status != DUCKHTS_CIGAR_INPUT_TOO_LARGE) {
        snprintf(message, sizeof(message), "%s: %s at %s %zu", name, reason,
                 cursor->text ? "operation-start byte offset" : "op index", cursor->offset + 1);
    } else {
        snprintf(message, sizeof(message), "%s: %s", name, reason);
    }
    duckdb_scalar_function_set_error(info, message);
    return true;
}

/* Validate child NULLs before borrowing packed storage. On failure offset
   identifies the NULL child within this row, independent of its list offset. */
static duckhts_cigar_status_t cigar_packed_cursor(duckdb_list_entry entry,
                                                const uint32_t *data,
                                                uint64_t *validity,
                                                duckhts_cigar_cursor_t *cursor) {
    memset(cursor, 0, sizeof(*cursor));
    if (entry.length > SIZE_MAX) {
        return DUCKHTS_CIGAR_INPUT_TOO_LARGE;
    }
    if (validity) {
        for (idx_t i = 0; i < entry.length; i++) {
            if (!duckdb_validity_row_is_valid(validity, entry.offset + i)) {
                cursor->offset = (size_t)i;
                return DUCKHTS_CIGAR_NULL_CHILD;
            }
        }
    }
    cursor->packed = entry.length ? data + entry.offset : NULL;
    cursor->count = (size_t)entry.length;
    return DUCKHTS_CIGAR_OK;
}

static void cigar_metric_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    const cigar_metric_definition_t *metric = duckdb_scalar_function_get_extra_info(info);
    int metric_id = metric->id;
    duckdb_vector cigar_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector strict_vec = cigar_strict_vector(input, 1);
    idx_t row_count = duckdb_data_chunk_get_size(input);
    cigar_metrics_t metrics;

    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(cigar_vec, row) || (strict_vec && !row_is_valid(strict_vec, row))) {
            set_null_at(output, row);
            continue;
        }

        idx_t cigar_len = 0;
        const char *cigar = get_string_at(cigar_vec, row, &cigar_len);
        duckhts_cigar_cursor_t cursor = {.text = cigar, .count = cigar_len};
        duckhts_cigar_status_t status = cigar_scan_metrics(&cursor, &metrics);
        if (status < 0 || metrics.seen_ops == 0) {
            if (cigar_strict_error(info, metric->name, strict_vec, row, status, &cursor)) {
                return;
            }
            set_null_at(output, row);
            continue;
        }

        if (metric_id == CIGAR_METRIC_HAS_SOFT_CLIP || metric_id == CIGAR_METRIC_HAS_HARD_CLIP) {
            bool *out_data = (bool *)duckdb_vector_get_data(output);
            out_data[row] = (metric_id == CIGAR_METRIC_HAS_SOFT_CLIP) ? (metrics.has_soft_clip != 0)
                                                                       : (metrics.has_hard_clip != 0);
        } else {
            int64_t *out_data = (int64_t *)duckdb_vector_get_data(output);
            switch (metric_id) {
            case CIGAR_METRIC_LEFT_SOFT_CLIP:
                out_data[row] = metrics.left_soft_clip;
                break;
            case CIGAR_METRIC_RIGHT_SOFT_CLIP:
                out_data[row] = metrics.right_soft_clip;
                break;
            case CIGAR_METRIC_QUERY_LENGTH:
                out_data[row] = metrics.query_length;
                break;
            case CIGAR_METRIC_ALIGNED_QUERY_LENGTH:
                out_data[row] = metrics.aligned_query_length;
                break;
            case CIGAR_METRIC_REFERENCE_LENGTH:
                out_data[row] = metrics.reference_length;
                break;
            default:
                set_null_at(output, row);
                break;
            }
        }
    }
}

static void cigar_has_op_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector cigar_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector op_vec = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector strict_vec = cigar_strict_vector(input, 2);
    bool *out_data = (bool *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(cigar_vec, row) || !row_is_valid(op_vec, row) ||
            (strict_vec && !row_is_valid(strict_vec, row))) {
            set_null_at(output, row);
            continue;
        }

        idx_t cigar_len = 0;
        idx_t op_len = 0;
        const char *cigar = get_string_at(cigar_vec, row, &cigar_len);
        const char *op_text = get_string_at(op_vec, row, &op_len);
        int op = cigar_requested_op(op_text, op_len);
        if (!duckhts_cigar_op_supported(op)) {
            if (cigar_strict_error(info, "cigar_has_op", strict_vec, row,
                                   DUCKHTS_CIGAR_INVALID_REQUESTED_OP, NULL)) {
                return;
            }
            set_null_at(output, row);
            continue;
        }

        duckhts_cigar_cursor_t cursor = {.text = cigar, .count = cigar_len};
        cigar_metrics_t metrics;
        duckhts_cigar_status_t status = cigar_scan_metrics(&cursor, &metrics);
        if (status < 0) {
            if (cigar_strict_error(info, "cigar_has_op", strict_vec, row, status, &cursor)) {
                return;
            }
            set_null_at(output, row);
            continue;
        }
        out_data[row] = (metrics.seen_ops & (UINT32_C(1) << op)) != 0;
    }
}

/* Binary (UINTEGER[]) overload of cigar_metric_scalar. Same metrics, computed
   from the packed BAM CIGAR list read_bam(cigar_representation := 'binary')
   emits; bit-identical to the text path for the same read. */
static void cigar_metric_binary_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    const cigar_metric_definition_t *metric = duckdb_scalar_function_get_extra_info(info);
    int metric_id = metric->id;
    duckdb_vector cigar_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector strict_vec = cigar_strict_vector(input, 1);
    duckdb_list_entry *list_data = (duckdb_list_entry *)duckdb_vector_get_data(cigar_vec);
    duckdb_vector child_vec = duckdb_list_vector_get_child(cigar_vec);
    uint32_t *child_data = (uint32_t *)duckdb_vector_get_data(child_vec);
    uint64_t *child_validity = duckdb_vector_get_validity(child_vec);
    idx_t row_count = duckdb_data_chunk_get_size(input);
    cigar_metrics_t metrics;

    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(cigar_vec, row) || (strict_vec && !row_is_valid(strict_vec, row))) {
            set_null_at(output, row);
            continue;
        }
        duckhts_cigar_cursor_t cursor;
        duckhts_cigar_status_t status = cigar_packed_cursor(list_data[row], child_data,
                                                           child_validity, &cursor);
        if (status > 0) {
            status = cigar_scan_metrics(&cursor, &metrics);
        }
        if (status < 0 || metrics.seen_ops == 0) {
            if (cigar_strict_error(info, metric->name, strict_vec, row, status, &cursor)) {
                return;
            }
            set_null_at(output, row);
            continue;
        }
        if (metric_id == CIGAR_METRIC_HAS_SOFT_CLIP || metric_id == CIGAR_METRIC_HAS_HARD_CLIP) {
            bool *out_data = (bool *)duckdb_vector_get_data(output);
            out_data[row] = (metric_id == CIGAR_METRIC_HAS_SOFT_CLIP) ? (metrics.has_soft_clip != 0)
                                                                       : (metrics.has_hard_clip != 0);
        } else {
            int64_t *out_data = (int64_t *)duckdb_vector_get_data(output);
            switch (metric_id) {
            case CIGAR_METRIC_LEFT_SOFT_CLIP:      out_data[row] = metrics.left_soft_clip; break;
            case CIGAR_METRIC_RIGHT_SOFT_CLIP:     out_data[row] = metrics.right_soft_clip; break;
            case CIGAR_METRIC_QUERY_LENGTH:        out_data[row] = metrics.query_length; break;
            case CIGAR_METRIC_ALIGNED_QUERY_LENGTH:out_data[row] = metrics.aligned_query_length; break;
            case CIGAR_METRIC_REFERENCE_LENGTH:    out_data[row] = metrics.reference_length; break;
            default: set_null_at(output, row); break;
            }
        }
    }
}

/* Binary (UINTEGER[] cigar, VARCHAR op) overload of cigar_has_op_scalar. */
static void cigar_has_op_binary_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector cigar_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector op_vec = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector strict_vec = cigar_strict_vector(input, 2);
    duckdb_list_entry *list_data = (duckdb_list_entry *)duckdb_vector_get_data(cigar_vec);
    duckdb_vector child_vec = duckdb_list_vector_get_child(cigar_vec);
    uint32_t *child_data = (uint32_t *)duckdb_vector_get_data(child_vec);
    uint64_t *child_validity = duckdb_vector_get_validity(child_vec);
    bool *out_data = (bool *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(cigar_vec, row) || !row_is_valid(op_vec, row) ||
            (strict_vec && !row_is_valid(strict_vec, row))) {
            set_null_at(output, row);
            continue;
        }
        idx_t op_len = 0;
        const char *op_text = get_string_at(op_vec, row, &op_len);
        int op = cigar_requested_op(op_text, op_len);
        if (!duckhts_cigar_op_supported(op)) {
            if (cigar_strict_error(info, "cigar_has_op", strict_vec, row,
                                   DUCKHTS_CIGAR_INVALID_REQUESTED_OP, NULL)) {
                return;
            }
            set_null_at(output, row);
            continue;
        }
        duckhts_cigar_cursor_t cursor;
        duckhts_cigar_status_t status = cigar_packed_cursor(list_data[row], child_data,
                                                           child_validity, &cursor);
        cigar_metrics_t metrics;
        if (status > 0) {
            status = cigar_scan_metrics(&cursor, &metrics);
        }
        if (status < 0) {
            if (cigar_strict_error(info, "cigar_has_op", strict_vec, row, status, &cursor)) {
                return;
            }
            set_null_at(output, row);
            continue;
        }
        out_data[row] = (metrics.seen_ops & (UINT32_C(1) << op)) != 0;
    }
}

/* ================================================================
 * cigar_aligned_blocks(cigar, pos)
 *   -> STRUCT(ref_start BIGINT[], query_start BIGINT[], width BIGINT[])
 *
 * One block per M, = or X op, in CIGAR order, never merged: pysam's
 * get_blocks() and GenomicAlignments' cigarRangesAlongReferenceSpace /
 * cigarRangesAlongQuerySpace over the aligned ops. ref_start is pos plus the
 * reference consumed before the block, so it carries whatever base pos uses;
 * query_start is the 0-based offset into the stored SEQ (S counts, H does
 * not); width is the op length. Invalid input is NULL by default and an error
 * with strict TRUE. NULL arguments and empty or '*' CIGAR remain NULL.
 * A valid CIGAR with no aligned op is three empty lists.
 *
 * Op semantics come from htslib's packed table: bam_cigar_type(op) bit 1
 * consumes query, bit 2 consumes reference, and type 3 is exactly M, =, X.
 * ================================================================ */

enum {
    CIGAR_BLOCK_REF_START = 0,
    CIGAR_BLOCK_QUERY_START = 1,
    CIGAR_BLOCK_WIDTH = 2,
    CIGAR_BLOCK_FIELD_COUNT = 3
};
static const char *CIGAR_BLOCK_FIELD_NAMES[CIGAR_BLOCK_FIELD_COUNT] = {"ref_start", "query_start", "width"};

/* Write cursor over the three block lists of one chunk. Capacity is reserved
   once from an upper bound on blocks (the chunk's op count, or its byte count
   for text) and the lists are sized to the cursor when the chunk is published,
   so nothing partial is visible if a reserve fails. */
typedef struct {
    duckdb_vector list[CIGAR_BLOCK_FIELD_COUNT]; /* struct children, LIST(BIGINT) */
    duckdb_list_entry *entry[CIGAR_BLOCK_FIELD_COUNT];
    int64_t *value[CIGAR_BLOCK_FIELD_COUNT];     /* list child data, capacity elements */
    idx_t n;                                     /* blocks written so far */
} cigar_block_sink_t;

static int cigar_block_sink_open(duckdb_function_info info, duckdb_vector output, idx_t capacity,
                                 cigar_block_sink_t *sink) {
    sink->n = 0;
    for (int f = 0; f < CIGAR_BLOCK_FIELD_COUNT; f++) {
        sink->list[f] = duckdb_struct_vector_get_child(output, (idx_t)f);
        if (duckdb_list_vector_set_size(sink->list[f], 0) != DuckDBSuccess ||
            duckdb_list_vector_reserve(sink->list[f], capacity) != DuckDBSuccess) {
            duckdb_scalar_function_set_error(info, "cigar_aligned_blocks: failed to grow list storage");
            return 0;
        }
    }
    for (int f = 0; f < CIGAR_BLOCK_FIELD_COUNT; f++) {
        sink->entry[f] = (duckdb_list_entry *)duckdb_vector_get_data(sink->list[f]);
        sink->value[f] = (int64_t *)duckdb_vector_get_data(duckdb_list_vector_get_child(sink->list[f]));
    }
    return 1;
}

static void cigar_block_sink_null_row(cigar_block_sink_t *sink, duckdb_vector output, idx_t row) {
    for (int f = 0; f < CIGAR_BLOCK_FIELD_COUNT; f++) {
        sink->entry[f][row].offset = sink->n;
        sink->entry[f][row].length = 0;
        set_null_at(sink->list[f], row);
    }
    set_null_at(output, row);
}

static void cigar_block_sink_row(cigar_block_sink_t *sink, idx_t row, idx_t start) {
    for (int f = 0; f < CIGAR_BLOCK_FIELD_COUNT; f++) {
        sink->entry[f][row].offset = start;
        sink->entry[f][row].length = sink->n - start;
    }
}

static void cigar_block_sink_publish(duckdb_function_info info, cigar_block_sink_t *sink) {
    for (int f = 0; f < CIGAR_BLOCK_FIELD_COUNT; f++) {
        if (duckdb_list_vector_set_size(sink->list[f], sink->n) != DuckDBSuccess) {
            duckdb_scalar_function_set_error(info, "cigar_aligned_blocks: failed to size list storage");
            return;
        }
    }
}

/* Validate each operation before its candidate write; the wrapper rewinds the entire row on
   failure. Reservation includes every input op, so writing a candidate for
   nonaligned ops is safe and only aligned ops advance the output cursor. */
static inline duckhts_cigar_status_t cigar_scan_blocks(cigar_block_sink_t *sink,
                                                     duckhts_cigar_cursor_t *cursor, int64_t pos) {
    duckhts_cigar_op_t op;
    duckhts_cigar_status_t status;
    size_t op_start = cursor->offset;

    while ((status = duckhts_cigar_next(cursor, &op)) > 0) {
        /* Spans are nonnegative BIGINTs. Negative pos is allowed, and adding
           a span cannot underflow. Check the reference end before addition. */
        if (pos > INT64_MAX - cursor->reference_span) {
            /* Retain the failing operation's location for adapter diagnostics. */
            cursor->offset = op_start;
            return DUCKHTS_CIGAR_POSITION_OVERFLOW;
        }
        idx_t n = sink->n;
        sink->value[CIGAR_BLOCK_REF_START][n] = pos + op.reference_start;
        sink->value[CIGAR_BLOCK_QUERY_START][n] = op.query_start;
        sink->value[CIGAR_BLOCK_WIDTH][n] = op.length;
        sink->n = n + (idx_t)(op.type == 3);
        op_start = cursor->offset;
    }
    if (status < 0) {
        return status;
    }
    return cursor->offset != 0 ? DUCKHTS_CIGAR_OK : DUCKHTS_CIGAR_END;
}

static void cigar_aligned_blocks_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector cigar_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector pos_vec = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector strict_vec = cigar_strict_vector(input, 2);
    int64_t *pos_data = (int64_t *)duckdb_vector_get_data(pos_vec);
    idx_t row_count = duckdb_data_chunk_get_size(input);
    cigar_block_sink_t sink;
    idx_t capacity = 0;

    /* A text CIGAR has no more ops than bytes. */
    for (idx_t row = 0; row < row_count; row++) {
        idx_t cigar_len = 0;
        if (!row_is_valid(cigar_vec, row) || !row_is_valid(pos_vec, row) ||
            (strict_vec && !row_is_valid(strict_vec, row))) {
            continue;
        }
        (void)get_string_at(cigar_vec, row, &cigar_len);
        if (cigar_len > (idx_t)-1 - capacity) {
            duckdb_scalar_function_set_error(info, "cigar_aligned_blocks: chunk CIGAR length overflows");
            return;
        }
        capacity += cigar_len;
    }
    if (!cigar_block_sink_open(info, output, capacity, &sink)) {
        return;
    }

    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(cigar_vec, row) || !row_is_valid(pos_vec, row) ||
            (strict_vec && !row_is_valid(strict_vec, row))) {
            cigar_block_sink_null_row(&sink, output, row);
            continue;
        }

        idx_t cigar_len = 0;
        const char *cigar = get_string_at(cigar_vec, row, &cigar_len);
        duckhts_cigar_cursor_t cursor = {.text = cigar, .count = cigar_len};
        idx_t start = sink.n;

        duckhts_cigar_status_t status = cigar_scan_blocks(&sink, &cursor, pos_data[row]);
        if (status <= 0) {
            sink.n = start;
            cigar_block_sink_null_row(&sink, output, row);
            if (cigar_strict_error(info, "cigar_aligned_blocks", strict_vec, row, status, &cursor)) {
                return;
            }
            continue;
        }
        cigar_block_sink_row(&sink, row, start);
    }
    cigar_block_sink_publish(info, &sink);
}

/* Binary (UINTEGER[]) overload: the packed BAM CIGAR read_bam(
   cigar_representation := 'binary') emits, decoded in place from the list
   child; bit-identical to the text path for the same read. */
static void cigar_aligned_blocks_binary_scalar(duckdb_function_info info, duckdb_data_chunk input,
                                               duckdb_vector output) {
    duckdb_vector cigar_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector pos_vec = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector strict_vec = cigar_strict_vector(input, 2);
    duckdb_list_entry *list_data = (duckdb_list_entry *)duckdb_vector_get_data(cigar_vec);
    duckdb_vector child_vec = duckdb_list_vector_get_child(cigar_vec);
    uint32_t *child_data = (uint32_t *)duckdb_vector_get_data(child_vec);
    uint64_t *child_validity = duckdb_vector_get_validity(child_vec);
    int64_t *pos_data = (int64_t *)duckdb_vector_get_data(pos_vec);
    idx_t row_count = duckdb_data_chunk_get_size(input);
    cigar_block_sink_t sink;
    idx_t capacity = 0;

    /* Blocks cannot outnumber ops, and every row's op count is in its entry. */
    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(cigar_vec, row) || !row_is_valid(pos_vec, row) ||
            (strict_vec && !row_is_valid(strict_vec, row))) {
            continue;
        }
        if (list_data[row].length > (idx_t)-1 - capacity) {
            duckdb_scalar_function_set_error(info, "cigar_aligned_blocks: chunk op count overflows");
            return;
        }
        capacity += list_data[row].length;
    }
    if (!cigar_block_sink_open(info, output, capacity, &sink)) {
        return;
    }

    for (idx_t row = 0; row < row_count; row++) {
        if (!row_is_valid(cigar_vec, row) || !row_is_valid(pos_vec, row) ||
            (strict_vec && !row_is_valid(strict_vec, row))) {
            cigar_block_sink_null_row(&sink, output, row);
            continue;
        }
        duckhts_cigar_cursor_t cursor;
        duckhts_cigar_status_t status = cigar_packed_cursor(list_data[row], child_data,
                                                           child_validity, &cursor);
        idx_t start = sink.n;
        if (status > 0) {
            status = cigar_scan_blocks(&sink, &cursor, pos_data[row]);
        }
        if (status <= 0) {
            sink.n = start;
            cigar_block_sink_null_row(&sink, output, row);
            if (cigar_strict_error(info, "cigar_aligned_blocks", strict_vec, row, status, &cursor)) {
                return;
            }
            continue;
        }
        cigar_block_sink_row(&sink, row, start);
    }
    cigar_block_sink_publish(info, &sink);
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
    /* Overload set: VARCHAR->VARCHAR and LIST(UTINYINT)->LIST(UTINYINT) (BAM
       nt16 codes), so read_bam(sequence_encoding := 'nt16') stays in nt16. */
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type utinyint_type = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
    duckdb_logical_type list_type = duckdb_create_list_type(utinyint_type);

    duckdb_scalar_function_set set = duckdb_create_scalar_function_set("seq_revcomp");

    duckdb_scalar_function fn_txt = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn_txt, "seq_revcomp");
    duckdb_scalar_function_add_parameter(fn_txt, varchar_type);
    duckdb_scalar_function_set_return_type(fn_txt, varchar_type);
    duckdb_scalar_function_set_function(fn_txt, seq_revcomp_scalar);
    duckdb_add_scalar_function_to_set(set, fn_txt);

    duckdb_scalar_function fn_nt16 = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn_nt16, "seq_revcomp");
    duckdb_scalar_function_add_parameter(fn_nt16, list_type);
    duckdb_scalar_function_set_return_type(fn_nt16, list_type);
    duckdb_scalar_function_set_function(fn_nt16, seq_revcomp_nt16_scalar);
    duckdb_add_scalar_function_to_set(set, fn_nt16);

    duckdb_register_scalar_function_set(connection, set);

    duckdb_destroy_scalar_function(&fn_txt);
    duckdb_destroy_scalar_function(&fn_nt16);
    duckdb_destroy_scalar_function_set(&set);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&utinyint_type);
    duckdb_destroy_logical_type(&list_type);
}

static void register_seq_canonical_function(duckdb_connection connection) {
    /* Overload set: VARCHAR->VARCHAR and LIST(UTINYINT)->LIST(UTINYINT) (BAM
       nt16 codes), so read_bam(sequence_encoding := 'nt16') stays in nt16. */
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type utinyint_type = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
    duckdb_logical_type list_type = duckdb_create_list_type(utinyint_type);

    duckdb_scalar_function_set set = duckdb_create_scalar_function_set("seq_canonical");

    duckdb_scalar_function fn_txt = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn_txt, "seq_canonical");
    duckdb_scalar_function_add_parameter(fn_txt, varchar_type);
    duckdb_scalar_function_set_return_type(fn_txt, varchar_type);
    duckdb_scalar_function_set_function(fn_txt, seq_canonical_scalar);
    duckdb_add_scalar_function_to_set(set, fn_txt);

    duckdb_scalar_function fn_nt16 = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn_nt16, "seq_canonical");
    duckdb_scalar_function_add_parameter(fn_nt16, list_type);
    duckdb_scalar_function_set_return_type(fn_nt16, list_type);
    duckdb_scalar_function_set_function(fn_nt16, seq_canonical_nt16_scalar);
    duckdb_add_scalar_function_to_set(set, fn_nt16);

    duckdb_register_scalar_function_set(connection, set);

    duckdb_destroy_scalar_function(&fn_txt);
    duckdb_destroy_scalar_function(&fn_nt16);
    duckdb_destroy_scalar_function_set(&set);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&utinyint_type);
    duckdb_destroy_logical_type(&list_type);
}

static void register_seq_hash_2bit_function(duckdb_connection connection) {
    /* Overload set: VARCHAR and LIST(UTINYINT) (BAM nt16 codes). */
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type utinyint_type = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
    duckdb_logical_type list_type = duckdb_create_list_type(utinyint_type);
    duckdb_logical_type ubigint_type = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);

    duckdb_scalar_function_set set = duckdb_create_scalar_function_set("seq_hash_2bit");

    duckdb_scalar_function fn_txt = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn_txt, "seq_hash_2bit");
    duckdb_scalar_function_add_parameter(fn_txt, varchar_type);
    duckdb_scalar_function_set_return_type(fn_txt, ubigint_type);
    duckdb_scalar_function_set_function(fn_txt, seq_hash_2bit_scalar);
    duckdb_add_scalar_function_to_set(set, fn_txt);

    duckdb_scalar_function fn_nt16 = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn_nt16, "seq_hash_2bit");
    duckdb_scalar_function_add_parameter(fn_nt16, list_type);
    duckdb_scalar_function_set_return_type(fn_nt16, ubigint_type);
    duckdb_scalar_function_set_function(fn_nt16, seq_hash_2bit_nt16_scalar);
    duckdb_add_scalar_function_to_set(set, fn_nt16);

    duckdb_register_scalar_function_set(connection, set);

    duckdb_destroy_scalar_function(&fn_txt);
    duckdb_destroy_scalar_function(&fn_nt16);
    duckdb_destroy_scalar_function_set(&set);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&utinyint_type);
    duckdb_destroy_logical_type(&list_type);
    duckdb_destroy_logical_type(&ubigint_type);
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
    /* Overload set: VARCHAR (text sequence) and LIST(UTINYINT) (BAM nt16 codes
       from read_bam(sequence_encoding := 'nt16')). DuckDB resolves by argument
       type, so decode-free BAM pipelines skip the text round-trip. */
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type utinyint_type = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
    duckdb_logical_type list_type = duckdb_create_list_type(utinyint_type);
    duckdb_logical_type double_type = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);

    duckdb_scalar_function_set set = duckdb_create_scalar_function_set("seq_gc_content");

    duckdb_scalar_function fn_txt = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn_txt, "seq_gc_content");
    duckdb_scalar_function_add_parameter(fn_txt, varchar_type);
    duckdb_scalar_function_set_return_type(fn_txt, double_type);
    duckdb_scalar_function_set_function(fn_txt, seq_gc_content_scalar);
    duckdb_add_scalar_function_to_set(set, fn_txt);

    duckdb_scalar_function fn_nt16 = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn_nt16, "seq_gc_content");
    duckdb_scalar_function_add_parameter(fn_nt16, list_type);
    duckdb_scalar_function_set_return_type(fn_nt16, double_type);
    duckdb_scalar_function_set_function(fn_nt16, seq_gc_content_nt16_scalar);
    duckdb_add_scalar_function_to_set(set, fn_nt16);

    duckdb_register_scalar_function_set(connection, set);

    duckdb_destroy_scalar_function(&fn_txt);
    duckdb_destroy_scalar_function(&fn_nt16);
    duckdb_destroy_scalar_function_set(&set);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&utinyint_type);
    duckdb_destroy_logical_type(&list_type);
    duckdb_destroy_logical_type(&double_type);
}

static void register_sam_flag_predicate_function(duckdb_connection connection,
                                                 const char *name,
                                                 uint16_t mask) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, name);

    duckdb_logical_type usmallint_type = duckdb_create_logical_type(DUCKDB_TYPE_USMALLINT);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_scalar_function_add_parameter(fn, usmallint_type);
    duckdb_scalar_function_set_return_type(fn, bool_type);
    duckdb_scalar_function_set_function(fn, sam_flag_predicate_scalar);
    duckdb_scalar_function_set_extra_info(fn, (void *)(uintptr_t)mask, NULL);

    duckdb_register_scalar_function(connection, fn);

    duckdb_destroy_logical_type(&usmallint_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_scalar_function(&fn);
}

static void register_sam_flag_has_function(duckdb_connection connection) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "sam_flag_has");

    duckdb_logical_type usmallint_type = duckdb_create_logical_type(DUCKDB_TYPE_USMALLINT);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_scalar_function_add_parameter(fn, usmallint_type);
    duckdb_scalar_function_add_parameter(fn, usmallint_type);
    duckdb_scalar_function_set_return_type(fn, bool_type);
    duckdb_scalar_function_set_function(fn, sam_flag_has_scalar);

    duckdb_register_scalar_function(connection, fn);

    duckdb_destroy_logical_type(&usmallint_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_scalar_function(&fn);
}

static void register_sam_is_forward_aligned_function(duckdb_connection connection) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "is_forward_aligned");

    duckdb_logical_type int_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_scalar_function_add_parameter(fn, int_type);
    duckdb_scalar_function_set_return_type(fn, bool_type);
    duckdb_scalar_function_set_function(fn, sam_is_forward_aligned_scalar);

    duckdb_register_scalar_function(connection, fn);

    duckdb_destroy_logical_type(&int_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_scalar_function(&fn);
}

static void register_sam_flag_bits_function(duckdb_connection connection) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "sam_flag_bits");

    duckdb_logical_type usmallint_type = duckdb_create_logical_type(DUCKDB_TYPE_USMALLINT);
    duckdb_logical_type bool_types[SAM_FLAG_FIELD_COUNT];
    duckdb_logical_type struct_type;

    for (int i = 0; i < SAM_FLAG_FIELD_COUNT; i++) {
        bool_types[i] = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    }
    struct_type = duckdb_create_struct_type(bool_types, SAM_FLAG_FIELD_NAMES, SAM_FLAG_FIELD_COUNT);

    duckdb_scalar_function_add_parameter(fn, usmallint_type);
    duckdb_scalar_function_set_return_type(fn, struct_type);
    duckdb_scalar_function_set_function(fn, sam_flag_bits_scalar);

    duckdb_register_scalar_function(connection, fn);

    duckdb_destroy_logical_type(&usmallint_type);
    for (int i = 0; i < SAM_FLAG_FIELD_COUNT; i++) {
        duckdb_destroy_logical_type(&bool_types[i]);
    }
    duckdb_destroy_logical_type(&struct_type);
    duckdb_destroy_scalar_function(&fn);
}

static void register_cigar_metric_function(duckdb_connection connection,
                                           const cigar_metric_definition_t *metric) {
    /* Overload set: VARCHAR (text CIGAR) and UINTEGER[] (binary CIGAR as
       read_bam(cigar_representation := 'binary') emits). Both carry the same
       metadata in extra_info; results are bit-identical for the same read. */
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type uinteger_type = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
    duckdb_logical_type list_type = duckdb_create_list_type(uinteger_type);
    duckdb_logical_type return_type = duckdb_create_logical_type(metric->return_type);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);

    duckdb_scalar_function_set set = duckdb_create_scalar_function_set(metric->name);

    for (int binary = 0; binary < 2; binary++) {
        for (int strict = 0; strict < 2; strict++) {
            duckdb_scalar_function fn = duckdb_create_scalar_function();
            duckdb_scalar_function_set_name(fn, metric->name);
            duckdb_scalar_function_add_parameter(fn, binary ? list_type : varchar_type);
            if (strict) {
                duckdb_scalar_function_add_parameter(fn, bool_type);
            }
            duckdb_scalar_function_set_return_type(fn, return_type);
            duckdb_scalar_function_set_function(fn, binary ? cigar_metric_binary_scalar
                                                           : cigar_metric_scalar);
            duckdb_scalar_function_set_extra_info(fn, (void *)metric, NULL);
            duckdb_add_scalar_function_to_set(set, fn);
            duckdb_destroy_scalar_function(&fn);
        }
    }

    duckdb_register_scalar_function_set(connection, set);

    duckdb_destroy_scalar_function_set(&set);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&uinteger_type);
    duckdb_destroy_logical_type(&list_type);
    duckdb_destroy_logical_type(&return_type);
    duckdb_destroy_logical_type(&bool_type);
}

static void register_cigar_has_op_function(duckdb_connection connection) {
    /* Text/packed inputs, each with an optional final BOOLEAN strict argument. */
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type uinteger_type = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
    duckdb_logical_type list_type = duckdb_create_list_type(uinteger_type);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);

    duckdb_scalar_function_set set = duckdb_create_scalar_function_set("cigar_has_op");

    for (int binary = 0; binary < 2; binary++) {
        for (int strict = 0; strict < 2; strict++) {
            duckdb_scalar_function fn = duckdb_create_scalar_function();
            duckdb_scalar_function_set_name(fn, "cigar_has_op");
            duckdb_scalar_function_add_parameter(fn, binary ? list_type : varchar_type);
            duckdb_scalar_function_add_parameter(fn, varchar_type);
            if (strict) {
                duckdb_scalar_function_add_parameter(fn, bool_type);
            }
            duckdb_scalar_function_set_return_type(fn, bool_type);
            duckdb_scalar_function_set_function(fn, binary ? cigar_has_op_binary_scalar
                                                           : cigar_has_op_scalar);
            duckdb_add_scalar_function_to_set(set, fn);
            duckdb_destroy_scalar_function(&fn);
        }
    }

    duckdb_register_scalar_function_set(connection, set);

    duckdb_destroy_scalar_function_set(&set);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&uinteger_type);
    duckdb_destroy_logical_type(&list_type);
    duckdb_destroy_logical_type(&bool_type);
}

static void register_cigar_aligned_blocks_function(duckdb_connection connection) {
    /* Overload set as register_cigar_metric_function: VARCHAR or UINTEGER[]
       cigar with a BIGINT pos, returning one STRUCT of three BIGINT lists. */
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type uinteger_type = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
    duckdb_logical_type list_type = duckdb_create_list_type(uinteger_type);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type field_types[CIGAR_BLOCK_FIELD_COUNT];
    duckdb_logical_type struct_type;

    for (int f = 0; f < CIGAR_BLOCK_FIELD_COUNT; f++) {
        field_types[f] = duckdb_create_list_type(bigint_type);
    }
    struct_type = duckdb_create_struct_type(field_types, CIGAR_BLOCK_FIELD_NAMES, CIGAR_BLOCK_FIELD_COUNT);

    duckdb_scalar_function_set set = duckdb_create_scalar_function_set("cigar_aligned_blocks");

    for (int binary = 0; binary < 2; binary++) {
        for (int strict = 0; strict < 2; strict++) {
            duckdb_scalar_function fn = duckdb_create_scalar_function();
            duckdb_scalar_function_set_name(fn, "cigar_aligned_blocks");
            duckdb_scalar_function_add_parameter(fn, binary ? list_type : varchar_type);
            duckdb_scalar_function_add_parameter(fn, bigint_type);
            if (strict) {
                duckdb_scalar_function_add_parameter(fn, bool_type);
            }
            duckdb_scalar_function_set_return_type(fn, struct_type);
            duckdb_scalar_function_set_function(fn, binary ? cigar_aligned_blocks_binary_scalar
                                                           : cigar_aligned_blocks_scalar);
            duckdb_add_scalar_function_to_set(set, fn);
            duckdb_destroy_scalar_function(&fn);
        }
    }

    duckdb_register_scalar_function_set(connection, set);

    duckdb_destroy_scalar_function_set(&set);
    for (int f = 0; f < CIGAR_BLOCK_FIELD_COUNT; f++) {
        duckdb_destroy_logical_type(&field_types[f]);
    }
    duckdb_destroy_logical_type(&struct_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_logical_type(&list_type);
    duckdb_destroy_logical_type(&uinteger_type);
    duckdb_destroy_logical_type(&varchar_type);
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
    for (size_t i = 0; i < sizeof(CIGAR_METRICS) / sizeof(CIGAR_METRICS[0]); i++) {
        register_cigar_metric_function(connection, &CIGAR_METRICS[i]);
    }
    register_cigar_has_op_function(connection);
    register_cigar_aligned_blocks_function(connection);
    register_sam_flag_bits_function(connection);
    register_sam_flag_has_function(connection);
    register_sam_is_forward_aligned_function(connection);
    register_sam_flag_predicate_function(connection, "is_paired", SAM_FLAG_PAIRED);
    register_sam_flag_predicate_function(connection, "is_proper_pair", SAM_FLAG_PROPER_PAIR);
    register_sam_flag_predicate_function(connection, "is_unmapped", SAM_FLAG_UNMAPPED);
    register_sam_flag_predicate_function(connection, "is_next_segment_unmapped", SAM_FLAG_MATE_UNMAPPED);
    register_sam_flag_predicate_function(connection, "is_reverse_complemented", SAM_FLAG_REVERSE);
    register_sam_flag_predicate_function(connection, "is_next_segment_reverse_complemented", SAM_FLAG_MATE_REVERSE);
    register_sam_flag_predicate_function(connection, "is_first_segment", SAM_FLAG_READ1);
    register_sam_flag_predicate_function(connection, "is_last_segment", SAM_FLAG_READ2);
    register_sam_flag_predicate_function(connection, "is_secondary", SAM_FLAG_SECONDARY);
    register_sam_flag_predicate_function(connection, "is_qc_fail", SAM_FLAG_QCFAIL);
    register_sam_flag_predicate_function(connection, "is_duplicate", SAM_FLAG_DUPLICATE);
    register_sam_flag_predicate_function(connection, "is_supplementary", SAM_FLAG_SUPPLEMENTARY);
}
