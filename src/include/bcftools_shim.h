#ifndef DUCKHTS_BCFTOOLS_SHIM_H
#define DUCKHTS_BCFTOOLS_SHIM_H

#include <stdint.h>
#include <htslib/kfunc.h>
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
#define bcf_double_is_missing_or_vector_end(x) \
    (bcf_double_is_missing(x) || bcf_double_is_vector_end(x))

static inline double calc_binom_two_sided(int na, int nb, double aprob) {
    double prob;
    if (!na && !nb) return -1;
    if (na == nb) return 1;
    prob = na > nb ? 2 * kf_betai(na, nb + 1, aprob) : 2 * kf_betai(nb, na + 1, aprob);
    return prob > 1 ? 1 : prob;
}

#endif
