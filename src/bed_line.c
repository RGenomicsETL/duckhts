#include "include/bed_line.h"

#include <stdlib.h>
#include <string.h>

/* Spans shorter than this parse from the stack; longer spans use the heap. */
#define BED_PARSE_STACK_BYTES 64

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

int duckhts_bed_parse_int64_span(const char *s, int len, int64_t *out) {
    char stack_buf[BED_PARSE_STACK_BYTES];
    char *buf = stack_buf;
    char *end = NULL;
    long long value;
    int ok;

    if (!s || len <= 0 || !out) return 0;
    if ((size_t)len >= sizeof(stack_buf)) {
        buf = (char *)malloc((size_t)len + 1);
        if (!buf) return 0;
    }
    memcpy(buf, s, (size_t)len);
    buf[len] = '\0';
    value = strtoll(buf, &end, 10);
    ok = end && *end == '\0';
    if (ok) *out = (int64_t)value;
    if (buf != stack_buf) free(buf);
    return ok;
}
