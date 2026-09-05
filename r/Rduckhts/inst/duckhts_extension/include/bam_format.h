#ifndef DUCKHTS_BAM_FORMAT_H
#define DUCKHTS_BAM_FORMAT_H

#include <htslib/sam.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>

/* Host-neutral formatters over borrowed, HTSlib-validated record fields.
 * No allocation. Output is NUL-terminated on success; on failure it must not
 * be published. Capacity includes the terminator; length excludes it.
 * Preserve DuckHTS's value-only AUX representation and %g precision. */
static inline int duckhts_bam_cigar_capacity(uint32_t count, size_t *capacity) {
    /* Nine decimal digits for the 28-bit length, one operation character. */
    if (count && (SIZE_MAX - 1) / count < 10) return 0;
    *capacity = (size_t)count * 10 + 1;
    return 1;
}

static inline int duckhts_bam_aux_capacity(const uint8_t *aux, const uint8_t *end,
                                          size_t *capacity) {
    if (aux >= end) return 0;
    switch (bam_aux_type(aux)) {
        case 'A': *capacity = 2; return 1;
        case 'c': case 'C': case 's': case 'S': case 'i': case 'I':
        case 'f': case 'd': *capacity = 32; return 1;
        case 'Z': case 'H': {
            const uint8_t *nul = memchr(aux + 1, 0, (size_t)(end - aux - 1));
            if (!nul) return 0;
            *capacity = (size_t)(nul - aux);
            return 1;
        }
        case 'B': {
            if (end - aux < 6) return 0;
            size_t width;
            switch (aux[1]) {
                case 'c': case 'C': width = 1; break;
                case 's': case 'S': width = 2; break;
                case 'i': case 'I': case 'f': width = 4; break;
                default: return 0;
            }
            uint32_t count = bam_auxB_len(aux);
            if (count > (size_t)(end - aux - 6) / width) return 0;
            /* Comma plus signed integer or six-significant-digit %g text;
             * reserve 32 bytes/value, plus subtype and NUL. */
            if (count && (SIZE_MAX - 2) / count < 32) return 0;
            *capacity = (size_t)count * 32 + 2;
            return 1;
        }
        default: return 0;
    }
}

static inline int duckhts_bam_cigar_format(const uint32_t *cigar, uint32_t count,
                                          char *data, size_t capacity, size_t *length) {
    *length = 0;
    if (!capacity) return 0;
    data[0] = '\0';
    size_t used = 0;
    for (uint32_t i = 0; i < count; i++) {
        if (bam_cigar_op(cigar[i]) >= 10) return 0;
        /* At most nine digits; avoid general printf parsing per operation. */
        char digits[9];
        size_t count_digits = 0;
        uint32_t value = bam_cigar_oplen(cigar[i]);
        do {
            digits[count_digits++] = (char)('0' + value % 10);
            value /= 10;
        } while (value);
        if (count_digits + 1 >= capacity - used) return 0;
        while (count_digits) data[used++] = digits[--count_digits];
        data[used++] = bam_cigar_opchr(cigar[i]);
        data[used] = '\0';
    }
    *length = used;
    return 1;
}

static inline int duckhts_bam_aux_format(const uint8_t *aux, char *data,
                                        size_t capacity, size_t *length) {
    *length = 0;
    if (!capacity) return 0;
    data[0] = '\0';
    int written;
    switch (bam_aux_type(aux)) {
        case 'A': written = snprintf(data, capacity, "%c", bam_aux2A(aux)); break;
        case 'c': case 'C': case 's': case 'S': case 'i': case 'I':
            written = snprintf(data, capacity, "%" PRId64, bam_aux2i(aux)); break;
        case 'f': case 'd': written = snprintf(data, capacity, "%g", bam_aux2f(aux)); break;
        case 'Z': case 'H': {
            size_t size = strlen(bam_aux2Z(aux));
            if (size >= capacity) return 0;
            memcpy(data, bam_aux2Z(aux), size + 1);
            *length = size;
            return 1;
        }
        case 'B': {
            if (capacity < 2) return 0;
            data[0] = aux[1];
            data[1] = '\0';
            size_t used = 1;
            uint32_t count = bam_auxB_len(aux);
            for (uint32_t i = 0; i < count; i++) {
                written = aux[1] == 'f'
                    ? snprintf(data + used, capacity - used, ",%g", bam_auxB2f(aux, i))
                    : snprintf(data + used, capacity - used, ",%" PRId64, bam_auxB2i(aux, i));
                if (written < 0 || (size_t)written >= capacity - used) return 0;
                used += (size_t)written;
            }
            *length = used;
            return 1;
        }
        default: return 0;
    }
    if (written < 0 || (size_t)written >= capacity) return 0;
    *length = (size_t)written;
    return 1;
}

#endif
