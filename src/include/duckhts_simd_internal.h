#ifndef DUCKHTS_SIMD_INTERNAL_H
#define DUCKHTS_SIMD_INTERNAL_H

#include <stddef.h>

#include "duckhts_simd.h"

typedef void (*duckhts_base_counts_fn)(const char *seq, size_t len,
                                       duckhts_simd_base_counts_t *out);

typedef struct duckhts_simd_ops {
    const char *name;
    duckhts_base_counts_fn base_counts;
} duckhts_simd_ops_t;

const duckhts_simd_ops_t *duckhts_simd_scalar_ops(void);
const duckhts_simd_ops_t *duckhts_simd_avx2_ops_if_available(void);
int duckhts_simd_avx2_available(void);

#endif /* DUCKHTS_SIMD_INTERNAL_H */
