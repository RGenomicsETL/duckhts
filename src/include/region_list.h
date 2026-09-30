#ifndef DUCKHTS_REGION_LIST_H
#define DUCKHTS_REGION_LIST_H

#include <stddef.h>
#include <stdint.h>
#include <htslib/hts.h>

/* Parse comma-separated requests, trimming surrounding ASCII whitespace.
 * NULL/empty strings preserve the no-filter API; an empty list item is an error.
 * On success, items and their strings occupy ONE malloc-owned block: free only
 * *items, never individual strings. On failure, *items=NULL and *count=0.
 * No sorting or deduplication: FASTA retains one row per requested interval. */
int duckhts_region_list_parse(const char *text, char ***items, unsigned int *count,
                             char *error, size_t error_size);

/* HTSlib is the coordinate/name authority. Truly unknown contigs remain allowed
 * for the iterator's existing skip policy; malformed known-coordinate requests,
 * ambiguous names, header failures and parser allocation errors are rejected. */
int duckhts_region_list_validate(char *const *items, unsigned int count,
                                hts_name2id_f name2id, void *header,
                                char *error, size_t error_size);

/* ---- Typed interval plans -------------------------------------------------
 * DuckDB-independent, immutable-after-finish description of literal contig
 * names with 0-based half-open [start, end) intervals. Chromosome names are
 * never parsed as HTSlib region expressions. Zero-initialize, add every
 * interval, then finish (validated count and payload caps are enforced while
 * adding). Finish sorts by (name bytes, start) and coalesces overlapping or
 * adjacent intervals per name without broadening any of them. A finished plan
 * with zero intervals is a valid, empty selection. Destroy is always safe. */
#define DUCKHTS_INTERVAL_PLAN_MAX_INTERVALS 1000000u
#define DUCKHTS_INTERVAL_PLAN_MAX_BYTES ((uint64_t)128 * 1024 * 1024)

typedef struct {
    const char *chrom; /* borrowed from the plan; valid once finished */
    int64_t start, end;
    size_t name_offset; /* private until finish */
} duckhts_interval_t;

typedef struct {
    duckhts_interval_t *items;
    size_t count, capacity;
    char *names;
    size_t names_length, names_capacity;
    uint64_t payload;
    size_t added; /* intervals accepted before coalescing */
    int finished;
} duckhts_interval_plan_t;

/* chrom NULL, start < 0 or end <= start are errors; item is the 1-based
 * caller position used in messages. Returns 0 with an explicit error. */
int duckhts_interval_plan_add(duckhts_interval_plan_t *plan, const char *chrom,
                              int64_t start, int64_t end, size_t item,
                              char *error, size_t error_size);
int duckhts_interval_plan_finish(duckhts_interval_plan_t *plan, char *error, size_t error_size);
void duckhts_interval_plan_destroy(duckhts_interval_plan_t *plan);

#endif
