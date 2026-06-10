#ifndef DUCKHTS_BCFTOOLS_SHIM_H
#define DUCKHTS_BCFTOOLS_SHIM_H

#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <htslib/hts_defs.h>
#include <htslib/kfunc.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>

static const uint64_t bcf_double_missing = 0x7ff0000000000001;
static const uint64_t bcf_double_vector_end = 0x7ff0000000000002;
static inline void bcf_double_set(double *ptr, uint64_t value) {
    union {
        uint64_t i;
        double d;
    } u;
    u.i = value;
    *ptr = u.d;
}
static inline int bcf_double_test(double d, uint64_t value) {
    union {
        uint64_t i;
        double d;
    } u;
    u.d = d;
    return u.i == value ? 1 : 0;
}
#define bcf_double_set_vector_end(x) bcf_double_set(&(x), bcf_double_vector_end)
#define bcf_double_set_missing(x) bcf_double_set(&(x), bcf_double_missing)
#define bcf_double_is_vector_end(x) bcf_double_test((x), bcf_double_vector_end)
#define bcf_double_is_missing(x) bcf_double_test((x), bcf_double_missing)
#define bcf_double_is_missing_or_vector_end(x) (bcf_double_test((x), bcf_double_missing) || bcf_double_test((x), bcf_double_vector_end))

static inline double calc_binom_two_sided(int na, int nb, double aprob) {
    double prob;
    if (!na && !nb) return -1;
    if (na == nb) return 1;
    prob = na > nb ? 2 * kf_betai(na, nb + 1, aprob) : 2 * kf_betai(nb, na + 1, aprob);
    if (prob > 1) prob = 1;
    return prob;
}

char *bcftools_version(void);
void duckhts_bcftools_error(const char *format, ...) HTS_FORMAT(HTS_PRINTF_FMT, 1, 2);
void duckhts_bcftools_error_errno(const char *format, ...) HTS_FORMAT(HTS_PRINTF_FMT, 1, 2);

int duckhts_filter_try_begin_impl(void);
jmp_buf *duckhts_filter_try_current_env(void);
void duckhts_filter_try_end(void);
const char *duckhts_filter_last_error(void);

#ifndef DUCKHTS_BCFTOOLS_SHIM_IMPLEMENTATION
#define duckhts_filter_try_begin() \
    (duckhts_filter_try_begin_impl() ? 1 : setjmp(*duckhts_filter_try_current_env()))
#endif

#endif
