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

/* ---- Typed interval plans ---------------------------------------------- */

void duckhts_interval_plan_destroy(duckhts_interval_plan_t *plan) {
    if (!plan) return;
    free(plan->items);
    free(plan->names);
    memset(plan, 0, sizeof(*plan));
}

static int interval_grow(void **buffer, size_t *capacity, size_t needed, size_t element) {
    if (needed <= *capacity) return 1;
    size_t next = *capacity ? *capacity : 16;
    while (next < needed) {
        if (next > SIZE_MAX / 2) return 0;
        next *= 2;
    }
    if (next > SIZE_MAX / element) return 0;
    void *grown = realloc(*buffer, next * element);
    if (!grown) return 0;
    *buffer = grown;
    *capacity = next;
    return 1;
}

int duckhts_interval_plan_add(duckhts_interval_plan_t *plan, const char *chrom,
                              int64_t start, int64_t end, size_t item,
                              char *error, size_t error_size) {
    if (plan->finished) {
        snprintf(error, error_size, "regions: interval plan is already finished");
        return 0;
    }
    if (!chrom) {
        snprintf(error, error_size, "regions: interval %zu has a NULL chrom", item);
        return 0;
    }
    if (start < 0) {
        snprintf(error, error_size, "regions: interval %zu has a negative start (%lld)",
                 item, (long long)start);
        return 0;
    }
    if (end <= start) {
        snprintf(error, error_size,
                 "regions: interval %zu has end <= start ([%lld, %lld) is empty or inverted)",
                 item, (long long)start, (long long)end);
        return 0;
    }
    if (plan->added >= DUCKHTS_INTERVAL_PLAN_MAX_INTERVALS) {
        snprintf(error, error_size,
                 "regions: more than %u intervals; the cap is %u intervals per query",
                 DUCKHTS_INTERVAL_PLAN_MAX_INTERVALS, DUCKHTS_INTERVAL_PLAN_MAX_INTERVALS);
        return 0;
    }
    size_t length = strlen(chrom);
    uint64_t cost = (uint64_t)length + 1 + (uint64_t)sizeof(duckhts_interval_t);
    if (cost < length || plan->payload > UINT64_MAX - cost ||
        plan->payload + cost > DUCKHTS_INTERVAL_PLAN_MAX_BYTES) {
        snprintf(error, error_size,
                 "regions: native payload exceeds the cap of %llu MiB",
                 (unsigned long long)(DUCKHTS_INTERVAL_PLAN_MAX_BYTES >> 20));
        return 0;
    }
    if (!interval_grow((void **)&plan->items, &plan->capacity, plan->count + 1,
                       sizeof(*plan->items))) goto oom;
    /* Consecutive intervals on one contig share a single stored name. */
    size_t offset;
    if (plan->count &&
        strcmp(plan->names + plan->items[plan->count - 1].name_offset, chrom) == 0) {
        offset = plan->items[plan->count - 1].name_offset;
    } else {
        if (length + 1 > SIZE_MAX - plan->names_length) goto oom;
        if (!interval_grow((void **)&plan->names, &plan->names_capacity,
                           plan->names_length + length + 1, 1)) goto oom;
        offset = plan->names_length;
        memcpy(plan->names + offset, chrom, length + 1);
        plan->names_length += length + 1;
    }
    plan->items[plan->count++] = (duckhts_interval_t){NULL, start, end, offset};
    plan->added++;
    plan->payload += cost;
    return 1;
oom:
    snprintf(error, error_size, "regions: out of memory storing intervals");
    return 0;
}

static int interval_compare(const void *a, const void *b) {
    const duckhts_interval_t *x = a, *y = b;
    int c = strcmp(x->chrom, y->chrom);
    if (c) return c;
    if (x->start != y->start) return x->start < y->start ? -1 : 1;
    if (x->end != y->end) return x->end < y->end ? -1 : 1;
    return 0;
}

int duckhts_interval_plan_finish(duckhts_interval_plan_t *plan, char *error, size_t error_size) {
    (void)error;
    (void)error_size;
    if (plan->finished) return 1;
    /* The name buffer no longer grows, so pointers into it are stable. */
    for (size_t i = 0; i < plan->count; i++)
        plan->items[i].chrom = plan->names + plan->items[i].name_offset;
    if (plan->count > 1) {
        qsort(plan->items, plan->count, sizeof(*plan->items), interval_compare);
        size_t out = 0;
        for (size_t i = 1; i < plan->count; i++) {
            duckhts_interval_t *cur = &plan->items[out];
            const duckhts_interval_t *next = &plan->items[i];
            if (strcmp(cur->chrom, next->chrom) == 0 &&
                next->start <= cur->end) {
                if (next->end > cur->end) cur->end = next->end;
            } else {
                plan->items[++out] = *next;
            }
        }
        plan->count = out + 1;
    }
    plan->finished = 1;
    return 1;
}
