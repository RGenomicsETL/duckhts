#ifndef DUCKHTS_SIMD_H
#define DUCKHTS_SIMD_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    uint64_t gc;
    uint64_t called;
    int invalid;
} duckhts_simd_base_counts_t;

void duckhts_simd_init(void);
const char *duckhts_simd_selected_backend(void);
const char *duckhts_simd_requested_backend(void);
int duckhts_simd_backend_compiled(const char *backend);
int duckhts_simd_backend_cpu_supported(const char *backend);
int duckhts_simd_backend_available(const char *backend);
int duckhts_simd_set_backend(const char *backend, char *err, size_t err_len);
void duckhts_simd_base_counts(const char *seq, size_t len,
                              duckhts_simd_base_counts_t *out);

#endif /* DUCKHTS_SIMD_H */
