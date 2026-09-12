#ifndef DUCKHTS_BCF_FORMAT_H
#define DUCKHTS_BCF_FORMAT_H

#include "bcf_scan.h"

/* Zero-initialize one buffer per selected tag and worker. Keep its tag and
 * selected sample set fixed. Reset loaded when advancing the input record;
 * output chunks may borrow the same decoded record. All storage is malloc-owned. */
typedef struct {
    void *data;       /* int32_t or float values, including missing/vector-end sentinels. */
    char **strings;   /* HTSlib's sample pointers and contiguous string storage. */
    int capacity;
    int count;       /* Numeric values or string bytes; zero for absent/invalid fields. */
    int stride;      /* Numeric values per selected sample. */
    int loaded;
} duckhts_bcf_format_t;

void duckhts_bcf_format_destroy(duckhts_bcf_format_t *values);

/* HTSlib owns decoding. This adapter checks encoded types, capacity arithmetic
 * and numeric scalar cardinality before publishing a view.
 * GT uses its integer payload despite its String header. Mismatches obey policy;
 * physical unpacking, capacity and allocation failures always return zero. */
int duckhts_bcf_format_decode(duckhts_bcf_format_t *values, bcf_hdr_t *header,
                              bcf1_t *record, const char *tag, int header_type,
                              duckhts_bcf_decode_policy_t policy, const char *reader_name,
                              char *error, size_t error_size);

#endif
