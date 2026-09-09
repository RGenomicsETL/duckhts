#ifndef DUCKHTS_BCF_FIELD_VECTOR_H
#define DUCKHTS_BCF_FIELD_VECTOR_H

#include "duckdb_extension.h"
#include "bcf_format.h"
#include <stdint.h>

duckdb_logical_type duckhts_bcf_field_type(int type, int is_list);
void duckhts_bcf_field_null(duckdb_vector vector, idx_t row, int is_list);
/* Values borrow one HTSlib decode. Missing items retain their ordinal; vector-end
 * padding terminates a sample slice. Zero returns before accessing a failed reserve. */
int duckhts_bcf_field_int32(duckdb_vector vector, idx_t row,
                           const int32_t *values, int count, int is_list);
int duckhts_bcf_field_float(duckdb_vector vector, idx_t row,
                           const float *values, int count, int is_list);
int duckhts_bcf_field_string(duckdb_vector vector, idx_t row,
                            const char *value, int is_list);

/* The sample index addresses the fixed selected-sample set used for decoding. */
int duckhts_bcf_format_write(duckdb_vector vector, idx_t row,
                             const duckhts_bcf_format_t *values, int header_type,
                             int is_list, int sample);

#endif
