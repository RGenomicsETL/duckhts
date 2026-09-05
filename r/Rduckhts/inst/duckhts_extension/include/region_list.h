#ifndef DUCKHTS_REGION_LIST_H
#define DUCKHTS_REGION_LIST_H

#include <stddef.h>
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

#endif
