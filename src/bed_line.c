#include "include/bed_line.h"

#include <string.h>

int duckhts_bed_line_is_meta(const char *line) {
    if (!line || !*line) return 1;
    if (line[0] == '#') return 1;
    if (strncmp(line, "track", 5) == 0) return 1;
    if (strncmp(line, "browser", 7) == 0) return 1;
    return 0;
}

const char *duckhts_bed_line_field(const char *line, int index, int *len) {
    int cur = 0;
    const char *start = line;
    const char *p = line;
    while (*p) {
        if (*p == '\t') {
            if (cur == index) {
                *len = (int)(p - start);
                return start;
            }
            cur++;
            start = p + 1;
        }
        p++;
    }
    if (cur == index) {
        *len = (int)(p - start);
        return start;
    }
    *len = 0;
    return NULL;
}

int duckhts_bed_line_split(const char *line, int max_fields,
                           const char **starts, int *lens) {
    int count = 0;
    const char *start = line;
    const char *p = line;

    if (!line || max_fields <= 0) return 0;
    for (;; p++) {
        if (*p != '\t' && *p != '\0') continue;
        starts[count] = start;
        lens[count] = (int)(p - start);
        count++;
        if (*p == '\0' || count == max_fields) return count;
        start = p + 1;
    }
}

/* The whitespace strtoll() skips in the C locale. */
static int bed_is_space(char c) {
    return c == ' ' || (c >= '\t' && c <= '\r');
}

int duckhts_bed_parse_int64_span(const char *s, int len, int64_t *out) {
    const uint64_t min_magnitude = (uint64_t)INT64_MAX + 1;
    const char *p = s;
    const char *end;
    uint64_t magnitude = 0;
    uint64_t limit;
    int negative = 0;

    if (!s || len <= 0 || !out) return 0;
    end = s + len;
    while (p < end && bed_is_space(*p)) p++;
    if (p < end && (*p == '+' || *p == '-')) {
        negative = *p == '-';
        p++;
    }
    if (p == end) return 0;
    limit = negative ? min_magnitude : (uint64_t)INT64_MAX;
    for (; p < end; p++) {
        unsigned digit = (unsigned)(unsigned char)*p - (unsigned)'0';
        if (digit > 9) return 0;
        /* Saturation is sticky: once at the limit, later digits keep it there. */
        if (magnitude > (limit - digit) / 10) {
            magnitude = limit;
        } else {
            magnitude = magnitude * 10 + digit;
        }
    }
    if (!negative) {
        *out = (int64_t)magnitude;
    } else if (magnitude == min_magnitude) {
        *out = INT64_MIN;
    } else {
        *out = -(int64_t)magnitude;
    }
    return 1;
}
