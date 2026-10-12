#include "count_mixture_model.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

int main(void) {
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
    b.counts[MIX_SITES] = INT64_MAX;
    assert(!duckhts_mixture_add(&b, 0, 0, 0, 0.5, 0.5, 1, 0));
    puts("count mixture model: 17 assertions passed");
    return 0;
}
