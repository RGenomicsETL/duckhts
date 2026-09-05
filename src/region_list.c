#include "include/region_list.h"

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int region_space(char c) {
    return c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f' || c == '\v';
}

/* HTSlib permits {quoted,reference}:start-end. A comma inside the initial
 * quoted reference is not a list separator. HTSlib validates the quote syntax. */
static const char *region_separator(const char *start) {
    const char *p = start;
    while (region_space(*p)) p++;
    if (*p == '{') {
        const char *close = strchr(p, '}');
        if (!close) return NULL;
        p = close + 1;
    }
    return strchr(p, ',');
}

int duckhts_region_list_parse(const char *text, char ***items, unsigned int *count,
                             char *error, size_t error_size) {
    *items = NULL;
    *count = 0;
    if (!text || !text[0]) return 1;

    size_t length = strlen(text);
    unsigned int n = 0;
    const char *start = text;
    for (;;) {
        if (n == INT_MAX) {
            snprintf(error, error_size, "region list: too many items for HTSlib");
            return 0;
        }
        const char *comma = region_separator(start);
        const char *end = comma ? comma : text + length;
        while (start < end && region_space(*start)) start++;
        while (end > start && region_space(end[-1])) end--;
        if (start == end) {
            snprintf(error, error_size, "region list: empty item at position %u", n + 1);
            return 0;
        }
        n++;
        if (!comma) break;
        start = comma + 1;
    }
    if (length == SIZE_MAX || (size_t)n > (SIZE_MAX - length - 1) / sizeof(char *)) {
        snprintf(error, error_size, "region list: allocation size overflow");
        return 0;
    }
    char **list = malloc((size_t)n * sizeof(*list) + length + 1);
    if (!list) {
        snprintf(error, error_size, "region list: out of memory");
        return 0;
    }
    char *storage = (char *)(list + n);
    start = text;
    for (unsigned int i = 0; i < n; i++) {
        const char *comma = region_separator(start);
        const char *end = comma ? comma : text + length;
        while (start < end && region_space(*start)) start++;
        while (end > start && region_space(end[-1])) end--;
        size_t bytes = (size_t)(end - start);
        list[i] = storage;
        memcpy(storage, start, bytes);
        storage[bytes] = '\0';
        storage += bytes + 1;
        if (comma) start = comma + 1;
    }
    *items = list;
    *count = n;
    return 1;
}

typedef struct {
    hts_name2id_f name2id;
    void *header; /* borrowed */
    int looked_up;
    int matched;
} region_name_lookup_t;

static int region_lookup(void *data, const char *name) {
    region_name_lookup_t *lookup = data;
    int tid = lookup->name2id(lookup->header, name);
    lookup->looked_up = 1;
    if (tid >= 0) lookup->matched = 1;
    return tid;
}

int duckhts_region_list_validate(char *const *items, unsigned int count,
                                hts_name2id_f name2id, void *header,
                                char *error, size_t error_size) {
    for (unsigned int i = 0; i < count; i++) {
        if (strcmp(items[i], ".") == 0 || strcmp(items[i], "*") == 0) continue;
        region_name_lookup_t lookup = {name2id, header, 0, 0};
        int tid = -1;
        hts_pos_t beg, end;
        const char *parsed = hts_parse_region(items[i], &tid, &beg, &end,
                                              region_lookup, &lookup,
                                              HTS_PARSE_THOUSANDS_SEP);
        if (!parsed && (tid != -1 || lookup.matched || !lookup.looked_up)) {
            snprintf(error, error_size, "region list: invalid item %u: %.160s", i + 1, items[i]);
            return 0;
        }
    }
    return 1;
}
