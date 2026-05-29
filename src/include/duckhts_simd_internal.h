#ifndef DUCKHTS_SIMD_INTERNAL_H
#define DUCKHTS_SIMD_INTERNAL_H

#include <stddef.h>
#include <stdint.h>

#include "duckhts_simd.h"

typedef uint64_t duckhts_simd_cap_t;

#define DUCKHTS_SIMD_CAP_NONE         ((duckhts_simd_cap_t)0)
#define DUCKHTS_SIMD_CAP_SCALAR       ((duckhts_simd_cap_t)1 << 0)
#define DUCKHTS_SIMD_CAP_SSE2         ((duckhts_simd_cap_t)1 << 1)
#define DUCKHTS_SIMD_CAP_SSE41        ((duckhts_simd_cap_t)1 << 2)
#define DUCKHTS_SIMD_CAP_AVX2         ((duckhts_simd_cap_t)1 << 3)
#define DUCKHTS_SIMD_CAP_AVX512       ((duckhts_simd_cap_t)1 << 4)
#define DUCKHTS_SIMD_CAP_NEON         ((duckhts_simd_cap_t)1 << 5)
#define DUCKHTS_SIMD_CAP_WASM_SIMD128 ((duckhts_simd_cap_t)1 << 6)

typedef void (*duckhts_base_counts_fn)(const char *seq, size_t len,
                                       duckhts_simd_base_counts_t *out);

#define DUCKHTS_SIMD_FN_BASE_COUNTS duckhts_base_counts_fn

typedef enum {
#define DUCKHTS_SIMD_KERNEL(id, field, sig, name) DUCKHTS_KERNEL_##id,
#include "duckhts_simd_kernels.def"
#undef DUCKHTS_SIMD_KERNEL
    DUCKHTS_KERNEL_COUNT
} duckhts_kernel_kind_t;

typedef struct duckhts_simd_kernel_slot {
    const char *kernel;
    const char *backend;
    duckhts_simd_cap_t cap;
    int priority;
    int scalar_fallback;
} duckhts_simd_kernel_slot_t;

typedef struct duckhts_simd_dispatch_table {
    const char *requested_backend;
    const char *selected_backend;
    duckhts_simd_cap_t viable_caps;
    duckhts_simd_kernel_slot_t slots[DUCKHTS_KERNEL_COUNT];
#define DUCKHTS_SIMD_KERNEL(id, field, sig, name) DUCKHTS_SIMD_FN_##sig field;
#include "duckhts_simd_kernels.def"
#undef DUCKHTS_SIMD_KERNEL
} duckhts_simd_dispatch_table_t;

typedef struct duckhts_simd_builder {
    duckhts_simd_dispatch_table_t *table;
    int best_priority[DUCKHTS_KERNEL_COUNT];
} duckhts_simd_builder_t;

const char *duckhts_simd_kernel_name(duckhts_kernel_kind_t kernel);
const char *duckhts_simd_cap_name(duckhts_simd_cap_t cap);
const duckhts_simd_dispatch_table_t *duckhts_simd_dispatch_snapshot(void);
void duckhts_simd_base_counts_with_table(const duckhts_simd_dispatch_table_t *table,
                                         const char *seq,
                                         size_t len,
                                         duckhts_simd_base_counts_t *out);

void duckhts_simd_builder_consider_base_counts(duckhts_simd_builder_t *builder,
                                               duckhts_simd_cap_t cap,
                                               const char *backend,
                                               int priority,
                                               duckhts_base_counts_fn fn);

void duckhts_simd_scalar_register(duckhts_simd_builder_t *builder);
void duckhts_simd_avx2_register(duckhts_simd_builder_t *builder);
int duckhts_simd_avx2_compiled(void);
int duckhts_simd_avx2_cpu_supported(void);
int duckhts_simd_avx2_available(void);
void duckhts_simd_avx512_register(duckhts_simd_builder_t *builder);
int duckhts_simd_avx512_compiled(void);
int duckhts_simd_avx512_cpu_supported(void);
int duckhts_simd_avx512_available(void);
void duckhts_simd_neon_register(duckhts_simd_builder_t *builder);
int duckhts_simd_neon_compiled(void);
int duckhts_simd_neon_cpu_supported(void);
int duckhts_simd_neon_available(void);
void duckhts_simd_wasm_simd128_register(duckhts_simd_builder_t *builder);
int duckhts_simd_wasm_simd128_compiled(void);
int duckhts_simd_wasm_simd128_cpu_supported(void);
int duckhts_simd_wasm_simd128_available(void);

#endif /* DUCKHTS_SIMD_INTERNAL_H */
