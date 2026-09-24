#ifndef DUCKHTS_CIGAR_H
#define DUCKHTS_CIGAR_H

#include <stddef.h>
#include <stdint.h>

#include <htslib/sam.h>

_Static_assert(BAM_CMATCH == 0 && BAM_CINS == 1 && BAM_CDEL == 2 && BAM_CREF_SKIP == 3 &&
               BAM_CSOFT_CLIP == 4 && BAM_CHARD_CLIP == 5 && BAM_CPAD == 6 && BAM_CEQUAL == 7 &&
               BAM_CDIFF == 8 && BAM_CBACK == 9,
               "htslib BAM CIGAR op codes");
_Static_assert(bam_cigar_type(BAM_CMATCH) == 3 && bam_cigar_type(BAM_CEQUAL) == 3 &&
               bam_cigar_type(BAM_CDIFF) == 3 && bam_cigar_type(BAM_CINS) == 1 &&
               bam_cigar_type(BAM_CSOFT_CLIP) == 1 && bam_cigar_type(BAM_CDEL) == 2 &&
               bam_cigar_type(BAM_CREF_SKIP) == 2 && bam_cigar_type(BAM_CHARD_CLIP) == 0 &&
               bam_cigar_type(BAM_CPAD) == 0,
               "aligned ops consume both query and reference");

/* Positive means success, zero means end/no-CIGAR, negative means failure.
   Adapters also use these statuses for input validity and coordinate checks. */
typedef enum {
    DUCKHTS_CIGAR_OK = 1,
    DUCKHTS_CIGAR_END = 0,
    DUCKHTS_CIGAR_MISSING_LENGTH = -1,
    DUCKHTS_CIGAR_ZERO_LENGTH = -2,
    DUCKHTS_CIGAR_MISSING_OP = -3,
    DUCKHTS_CIGAR_UNSUPPORTED_OP = -4,
    DUCKHTS_CIGAR_LENGTH_OVERFLOW = -5,
    DUCKHTS_CIGAR_QUERY_OVERFLOW = -6,
    DUCKHTS_CIGAR_REFERENCE_OVERFLOW = -7,
    DUCKHTS_CIGAR_NULL_CHILD = -8,
    DUCKHTS_CIGAR_INPUT_TOO_LARGE = -9,
    DUCKHTS_CIGAR_INVALID_REQUESTED_OP = -10,
    DUCKHTS_CIGAR_POSITION_OVERFLOW = -11
} duckhts_cigar_status_t;

/* Borrowed text bytes or packed BAM ops, with zero-initialized cursor/spans.
   text != NULL selects text; count is bytes for text, ops otherwise. Empty
   input and the single text byte '*' are no-CIGAR representations. Packed
   elements must all be present; host adapters own child NULL validation. */
typedef struct {
    const char *text;
    const uint32_t *packed;
    size_t count;
    size_t offset;
    int64_t query_span;
    int64_t reference_span;
} duckhts_cigar_cursor_t;

typedef struct {
    int code;
    int type; /* htslib query/reference consumption bits */
    int64_t length;
    int64_t query_start;     /* stored-query offset before this op */
    int64_t reference_start; /* reference offset before this op */
} duckhts_cigar_op_t;

static inline int duckhts_cigar_op_supported(int code) {
    return code >= BAM_CMATCH && code <= BAM_CDIFF;
}

/* Return OK for a checked op, END at end/no-CIGAR, a negative status for
   malformed input or overflow. A successful prefix is not whole-input validation: callers must
   pull through end and discard row results on error. On error, neither the
   cursor nor the output op advances. All lengths and consumed spans fit in
   nonnegative int64_t; H and P lengths are checked even though not consumed. */
static inline duckhts_cigar_status_t duckhts_cigar_next(duckhts_cigar_cursor_t *cursor,
                                                     duckhts_cigar_op_t *op) {
    size_t next = cursor->offset;
    int64_t length = 0;
    int code;

    if (next == cursor->count ||
        (cursor->text && cursor->count == 1 && cursor->text[0] == '*')) {
        return DUCKHTS_CIGAR_END;
    }
    if (cursor->text) {
        while (next < cursor->count) {
            unsigned char c = (unsigned char)cursor->text[next];
            if (c < '0' || c > '9') {
                break;
            }
            int digit = c - '0';
            if (length > (INT64_MAX - digit) / 10) {
                return DUCKHTS_CIGAR_LENGTH_OVERFLOW;
            }
            length = length * 10 + digit;
            next++;
        }
        if (next == cursor->count) {
            return DUCKHTS_CIGAR_MISSING_OP; /* Includes zero-valued trailing digits. */
        }
        if (next == cursor->offset) {
            return DUCKHTS_CIGAR_MISSING_LENGTH;
        }
        code = bam_cigar_table[(unsigned char)cursor->text[next++]];
    } else {
        uint32_t packed = cursor->packed[next++];
        code = bam_cigar_op(packed);
        length = bam_cigar_oplen(packed);
    }
    if (length == 0) {
        return DUCKHTS_CIGAR_ZERO_LENGTH;
    }
    if (!duckhts_cigar_op_supported(code)) {
        return DUCKHTS_CIGAR_UNSUPPORTED_OP;
    }

    int type = bam_cigar_type(code);
    int64_t query_step = (type & 1) ? length : 0;
    int64_t reference_step = (type & 2) ? length : 0;
    if (cursor->query_span > INT64_MAX - query_step) {
        return DUCKHTS_CIGAR_QUERY_OVERFLOW;
    }
    if (cursor->reference_span > INT64_MAX - reference_step) {
        return DUCKHTS_CIGAR_REFERENCE_OVERFLOW;
    }

    op->code = code;
    op->type = type;
    op->length = length;
    op->query_start = cursor->query_span;
    op->reference_start = cursor->reference_span;
    cursor->query_span += query_step;
    cursor->reference_span += reference_step;
    cursor->offset = next;
    return DUCKHTS_CIGAR_OK;
}

#endif
