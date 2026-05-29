#ifndef DUCKHTS_SIMD_INTERNAL_H
#define DUCKHTS_SIMD_INTERNAL_H

#include <stddef.h>

#include "duckhts_simd.h"

typedef void (*duckhts_base_counts_fn)(const char *seq, size_t len,
                                       duckhts_simd_base_counts_t *out);

typedef struct duckhts_simd_ops {
    const char *name;
    /* Non-scalar backends may leave individual kernel slots NULL; public
     * dispatch wrappers must fall back to the scalar implementation per op.
     * The scalar backend is the only backend required to fill every slot. */
    duckhts_base_counts_fn base_counts;
} duckhts_simd_ops_t;

const duckhts_simd_ops_t *duckhts_simd_scalar_ops(void);
const duckhts_simd_ops_t *duckhts_simd_avx2_ops_if_available(void);
int duckhts_simd_avx2_compiled(void);
int duckhts_simd_avx2_cpu_supported(void);
int duckhts_simd_avx2_available(void);
const duckhts_simd_ops_t *duckhts_simd_avx512_ops_if_available(void);
int duckhts_simd_avx512_compiled(void);
int duckhts_simd_avx512_cpu_supported(void);
int duckhts_simd_avx512_available(void);
const duckhts_simd_ops_t *duckhts_simd_neon_ops_if_available(void);
int duckhts_simd_neon_compiled(void);
int duckhts_simd_neon_cpu_supported(void);
int duckhts_simd_neon_available(void);
const duckhts_simd_ops_t *duckhts_simd_wasm_simd128_ops_if_available(void);
int duckhts_simd_wasm_simd128_compiled(void);
int duckhts_simd_wasm_simd128_cpu_supported(void);
int duckhts_simd_wasm_simd128_available(void);

#endif /* DUCKHTS_SIMD_INTERNAL_H */
