#include "count_mixture_model.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static void sum_check(int condition, int *checks) {
    assert(condition);
    (*checks)++;
}

static void test_exact_sums(void) {
    int checks = 0;
    duckhts_mixture_sum_t sum = {{0}};
    sum_check(sizeof(sum) == 152, &checks);
    sum_check(duckhts_mixture_sum_value(&sum) == 0, &checks);
    sum_check(duckhts_mixture_sum_add(&sum, 1), &checks);
    sum_check(duckhts_mixture_sum_add(&sum, ldexp(1, -53)), &checks);
    sum_check(duckhts_mixture_sum_value(&sum) == 1, &checks);
    sum_check(duckhts_mixture_sum_add(&sum, ldexp(1, -1074)), &checks);
    sum_check(duckhts_mixture_sum_value(&sum) == nextafter(1, 2), &checks);
    memset(&sum, 0, sizeof(sum));
    sum_check(duckhts_mixture_sum_add(&sum, nextafter(1, 2)), &checks);
    sum_check(duckhts_mixture_sum_add(&sum, ldexp(1, -53)), &checks);
    sum_check(duckhts_mixture_sum_value(&sum) == 1 + ldexp(1, -51), &checks);
    sum_check(!duckhts_mixture_sum_add(&sum, -1), &checks);
    sum_check(!duckhts_mixture_sum_add(&sum, INFINITY), &checks);
    sum_check(!duckhts_mixture_sum_add(&sum, NAN), &checks);
    sum_check(!duckhts_mixture_sum_add(&sum, 2147483649.0), &checks);
    sum_check(duckhts_mixture_sum_add(&sum, 2147483648.0), &checks);
    memset(&sum, 0, sizeof(sum));
    sum_check(duckhts_mixture_sum_add(&sum, 32768), &checks);
    for (int i = 0; i < 64; i++) {
        duckhts_mixture_sum_t copy = sum;
        sum_check(duckhts_mixture_sum_merge(&sum, &copy), &checks);
    }
    sum_check(duckhts_mixture_sum_value(&sum) == ldexp(1, 79), &checks);
    for (int exponent = -1074; exponent <= 0; exponent += 358) {
        uint64_t units = 0, random = 332;
        duckhts_mixture_sum_t parts[8] = {{{0}}};
        memset(&sum, 0, sizeof(sum));
        for (int i = 0; i < 10000; i++) {
            random = random * UINT64_C(6364136223846793005) + 1;
            uint64_t term = (random >> 32) % 1024;
            units += term;
            sum_check(duckhts_mixture_sum_add(&parts[i % 8], ldexp((double)term, exponent)), &checks);
        }
        for (int i = 7; i >= 0; i--)
            sum_check(duckhts_mixture_sum_merge(&sum, &parts[i]), &checks);
        sum_check(duckhts_mixture_sum_value(&sum) == ldexp((double)units, exponent), &checks);
    }
    sum.limbs[DUCKHTS_MIXTURE_SUM_LIMBS - 1] = UINT64_MAX;
    duckhts_mixture_sum_t high = {{0}};
    high.limbs[DUCKHTS_MIXTURE_SUM_LIMBS - 1] = 1;
    sum_check(!duckhts_mixture_sum_merge(&sum, &high), &checks);
    printf("count mixture exact sums: %d checks passed\n", checks);
}

int main(void) {
    test_exact_sums();
    duckhts_mixture_block_t b;
    memset(&b, 0, sizeof(b));
    b.matched = 1;
    assert(duckhts_mixture_add(&b, 2, 1, 4, 0.25, 0.75, 1, 0));
    assert(b.counts[MIX_READS] == 3 && b.counts[MIX_RR] == 1);
    assert(b.counts[MIX_RA] == 2 && b.counts[MIX_AA] == 0);
    assert(b.counts[MIX_OTHER] == 4 && b.denominator == 0.5625);
    assert(b.s[0][0] == 3 && b.t[0][0] == 3 && b.z[0] == 2);
    assert(b.y[1] == 0.25 && b.t[4][4] == 3 * 0.1875 * 0.1875);
    assert(!b.matched);
    assert(duckhts_mixture_add(&b, 0, 0, 0, 0, 1, 1, 0));
    assert(duckhts_mixture_add(&b, 9, 0, 0, 0.5, 0.5, 1, 1));
    assert(b.counts[MIX_ZERO] == 1 && b.counts[MIX_EXCLUDED] == 1);
    assert(b.counts[MIX_READS] == 3);
    duckhts_mixture_calibration_t c = {0, 0, 0, 30};
    double beta[3], gamma[6];
    duckhts_mixture_coefficients(0, &c, 0, beta, gamma);
    assert(beta[1] == 1 && gamma[1] == 1 && gamma[3] == -1);
    c.artefact_weight = 0.9;
    duckhts_mixture_coefficients(0, &c, 0, beta, gamma);
    assert(fabs(gamma[0] - 0.3) < 1e-15);
    c.artefact_weight = 0;
    duckhts_mixture_coefficients(0.2, &c, 1, beta, gamma);
    double ibd = gamma[1] + gamma[2] + gamma[3] + gamma[4] + gamma[5];
    assert(fabs(ibd) < 1e-15);
    double f = 0.3, d = gamma[1] * f + gamma[2] * f + (gamma[3] + gamma[4] + gamma[5]) * f * f;
    duckhts_mixture_coefficients((1 - sqrt(0.68)) / 2, &c, 0, beta, gamma);
    double unrelated = gamma[1] * f + gamma[2] * f + (gamma[3] + gamma[4] + gamma[5]) * f * f;
    assert(fabs(d - unrelated) < 1e-15);
    duckhts_mixture_fit_t fit;
    duckhts_mixture_fit(NULL, 0, &c, 0, &fit);
    assert(strcmp(fit.status, "no_sites") == 0 && isnan(fit.share));
    duckhts_mixture_block_t read_only;
    memset(&read_only, 0, sizeof(read_only));
    read_only.matched = 1;
    assert(duckhts_mixture_add(&read_only, 1, 0, 0, .5, .5, 1, 0));
    duckhts_mixture_fit(&read_only, 1, &c, 0, &fit);
    assert(strcmp(fit.status, "no_population_contrast") == 0);
    assert(isnan(fit.share));
    /* Read-only integer moments identify shares just inside both bounds. */
    for (int upper = 0; upper < 2; upper++) {
        duckhts_mixture_block_t near_bound;
        memset(&near_bound, 0, sizeof(near_bound));
        double u[3] = {1, .2, .8};
        double target = upper ? .4999 : .0001;
        double alt = upper ? 24997 : 10003;
        near_bound.counts[MIX_SITES] = 50000;
        near_bound.counts[MIX_USABLE] = 50000;
        near_bound.counts[MIX_SINGLETON] = 50000;
        near_bound.counts[MIX_READS] = 50000;
        for (int i = 0; i < 3; i++) {
            near_bound.y[i] = alt * u[i];
            for (int j = 0; j < 3; j++)
                near_bound.s[i][j] = 50000 * u[i] * u[j];
        }
        duckhts_mixture_fit(&near_bound, 1, &c, 0, &fit);
        assert(fabs(fit.share - target) < 1e-8);
        assert(strcmp(fit.status, "no_pairs") == 0);
    }
    b.counts[MIX_SITES] = INT64_MAX;
    assert(!duckhts_mixture_add(&b, 0, 0, 0, 0.5, 0.5, 1, 0));
    puts("count mixture model: 24 assertions passed");
    return 0;
}
