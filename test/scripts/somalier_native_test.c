#include "duckhts_somalier.h"

#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static unsigned failures;
static uint64_t next_random(void);

static const uint8_t test_panel_sha256[32] = {
    0x00u, 0x01u, 0x02u, 0x03u, 0x04u, 0x05u, 0x06u, 0x07u,
    0x08u, 0x09u, 0x0au, 0x0bu, 0x0cu, 0x0du, 0x0eu, 0x0fu,
    0x10u, 0x11u, 0x12u, 0x13u, 0x14u, 0x15u, 0x16u, 0x17u,
    0x18u, 0x19u, 0x1au, 0x1bu, 0x1cu, 0x1du, 0x1eu, 0x1fu
};

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "%s:%d: check failed: %s\n", __FILE__, __LINE__, #condition); \
        failures++; \
    } \
} while (0)

static int close_enough(double a, double b, double tolerance) {
    return fabs(a - b) <= tolerance;
}

static duckhts_somalier_counts_t counts_for_genotype(int genotype) {
    duckhts_somalier_counts_t counts = {0u, 0u, 0u, 1u};
    if (genotype == DUCKHTS_SOMALIER_HOM_A) {
        counts.allele_a = 20u;
    } else if (genotype == DUCKHTS_SOMALIER_HET) {
        counts.allele_a = 10u;
        counts.allele_b = 10u;
    } else if (genotype == DUCKHTS_SOMALIER_HOM_B) {
        counts.allele_b = 20u;
    }
    return counts;
}

static duckhts_somalier_masks_t mask_view(
    uint64_t *hom_a, uint64_t *het, uint64_t *hom_b,
    size_t site_count, size_t word_count, uint64_t middling) {
    duckhts_somalier_relatedness_settings_t settings;
    duckhts_somalier_masks_t masks;
    memset(&masks, 0, sizeof(masks));
    duckhts_somalier_relatedness_settings_default(&settings);
    masks.hom_a = hom_a;
    masks.het = het;
    masks.hom_b = hom_b;
    masks.site_count = site_count;
    masks.word_count = word_count;
    masks.middling_balance_count = middling;
    memcpy(masks.identity.panel_sha256, test_panel_sha256,
           sizeof(masks.identity.panel_sha256));
    masks.identity.min_depth = settings.min_depth;
    masks.identity.min_het_balance = settings.min_het_balance;
    masks.identity.hom_balance_cutoff = settings.hom_balance_cutoff;
    return masks;
}

static int profile_genotype(unsigned profile, size_t site) {
    unsigned state = (profile >> (2u * (unsigned)site)) & 3u;
    return state == 0u ? DUCKHTS_SOMALIER_UNKNOWN : (int)state - 1;
}

static double clamp_unit_reference(double value) {
    if (value < 0.0) return 0.0;
    if (value > 1.0) return 1.0;
    return value;
}

static void reference_pair(const int *a, const int *b, size_t count,
                           uint64_t middling_a, uint64_t middling_b,
                           duckhts_somalier_pair_stats_t *result) {
    size_t i;
    memset(result, 0, sizeof(*result));
    for (i = 0u; i < count; i++) {
        int known_a = a[i] >= 0;
        int known_b = b[i] >= 0;
        if (known_a && known_b) {
            result->jointly_called++;
            if ((a[i] == DUCKHTS_SOMALIER_HOM_A && b[i] == DUCKHTS_SOMALIER_HOM_B) ||
                (a[i] == DUCKHTS_SOMALIER_HOM_B && b[i] == DUCKHTS_SOMALIER_HOM_A)) {
                result->ibs0++;
            }
            if (a[i] == b[i]) result->ibs2++;
            if (a[i] == DUCKHTS_SOMALIER_HET && b[i] == DUCKHTS_SOMALIER_HET) {
                result->shared_hets++;
            }
            if (a[i] == DUCKHTS_SOMALIER_HOM_B && b[i] == DUCKHTS_SOMALIER_HOM_B) {
                result->shared_hom_b++;
            }
            if (a[i] == DUCKHTS_SOMALIER_HET) result->het_ab++;
            if (b[i] == DUCKHTS_SOMALIER_HET) result->het_ab++;
        }
        if (a[i] == DUCKHTS_SOMALIER_HET) result->het_count_a++;
        if (b[i] == DUCKHTS_SOMALIER_HET) result->het_count_b++;
        if (a[i] == DUCKHTS_SOMALIER_HOM_B) result->hom_b_count_a++;
        if (b[i] == DUCKHTS_SOMALIER_HOM_B) result->hom_b_count_b++;
        if ((a[i] == DUCKHTS_SOMALIER_HOM_A || a[i] == DUCKHTS_SOMALIER_HOM_B) &&
            known_b) {
            result->callable_hom_count_a++;
        }
        if ((b[i] == DUCKHTS_SOMALIER_HOM_A || b[i] == DUCKHTS_SOMALIER_HOM_B) &&
            known_a) {
            result->callable_hom_count_b++;
        }
        if (known_a && known_b && a[i] == b[i] &&
            (a[i] == DUCKHTS_SOMALIER_HOM_A || a[i] == DUCKHTS_SOMALIER_HOM_B)) {
            result->matching_hom_count++;
        }
    }
    result->relatedness = 2.0 * ((double)result->shared_hets -
        2.0 * (double)result->ibs0) /
        (result->het_ab == 0u ? 1.0 : (double)result->het_ab);
    {
        uint64_t hom_denominator = result->callable_hom_count_a <
            result->callable_hom_count_b ? result->callable_hom_count_a :
            result->callable_hom_count_b;
        uint64_t hom_b_denominator = result->hom_b_count_a < result->hom_b_count_b
            ? result->hom_b_count_a : result->hom_b_count_b;
        double pm_a = (double)middling_a / (double)(count + middling_a);
        double pm_b = (double)middling_b / (double)(count + middling_b);
        double penalty;
        double base;
        result->inferred_hom_concordance = (double)result->matching_hom_count /
            (hom_denominator == 0u ? 1.0 : (double)hom_denominator);
        result->raw_hom_b_concordance =
            ((double)result->shared_hom_b - 2.0 * (double)result->ibs0) /
            (hom_b_denominator == 0u ? 1.0 : (double)hom_b_denominator);
        result->p_middling_a = pm_a;
        result->p_middling_b = pm_b;
        penalty = fmax(0.0, 0.70 - clamp_unit_reference(result->raw_hom_b_concordance) +
            pm_a + pm_b);
        base = clamp_unit_reference(result->inferred_hom_concordance - penalty);
        result->adjusted_concordance = base <= 0.4 ? base :
            0.4 + 0.6 * cbrt((base - 0.4) / 0.6);
    }
}

static void compare_pair_stats(const duckhts_somalier_pair_stats_t *actual,
                               const duckhts_somalier_pair_stats_t *expected) {
    CHECK(actual->ibs0 == expected->ibs0);
    CHECK(actual->ibs2 == expected->ibs2);
    CHECK(actual->jointly_called == expected->jointly_called);
    CHECK(actual->shared_hets == expected->shared_hets);
    CHECK(actual->shared_hom_b == expected->shared_hom_b);
    CHECK(actual->het_ab == expected->het_ab);
    CHECK(actual->het_count_a == expected->het_count_a);
    CHECK(actual->het_count_b == expected->het_count_b);
    CHECK(actual->hom_b_count_a == expected->hom_b_count_a);
    CHECK(actual->hom_b_count_b == expected->hom_b_count_b);
    CHECK(actual->callable_hom_count_a == expected->callable_hom_count_a);
    CHECK(actual->callable_hom_count_b == expected->callable_hom_count_b);
    CHECK(actual->matching_hom_count == expected->matching_hom_count);
    CHECK(close_enough(actual->relatedness, expected->relatedness, 1e-14));
    CHECK(close_enough(actual->inferred_hom_concordance,
                       expected->inferred_hom_concordance, 1e-14));
    CHECK(close_enough(actual->raw_hom_b_concordance,
                       expected->raw_hom_b_concordance, 1e-14));
    CHECK(close_enough(actual->p_middling_a, expected->p_middling_a, 1e-14));
    CHECK(close_enough(actual->p_middling_b, expected->p_middling_b, 1e-14));
    CHECK(close_enough(actual->adjusted_concordance, expected->adjusted_concordance, 1e-14));
}

static void test_exhaustive_three_site_pairs(void) {
    duckhts_somalier_relatedness_settings_t settings;
    unsigned pa;
    unsigned pb;
    duckhts_somalier_relatedness_settings_default(&settings);
    for (pa = 0u; pa < 64u; pa++) {
        duckhts_somalier_counts_t counts_a[3];
        int genotypes_a[3];
        uint64_t hom_a_a = 0u;
        uint64_t het_a = 0u;
        uint64_t hom_b_a = 0u;
        duckhts_somalier_masks_t masks_a =
            mask_view(&hom_a_a, &het_a, &hom_b_a, 3u, 1u, 0u);
        size_t site;
        for (site = 0u; site < 3u; site++) {
            genotypes_a[site] = profile_genotype(pa, site);
            counts_a[site] = counts_for_genotype(genotypes_a[site]);
        }
        CHECK(duckhts_somalier_prepare_masks(counts_a, 3u, test_panel_sha256,
              &settings, &masks_a) ==
              DUCKHTS_SOMALIER_OK);
        for (pb = 0u; pb < 64u; pb++) {
            duckhts_somalier_counts_t counts_b[3];
            int genotypes_b[3];
            uint64_t hom_a_b = 0u;
            uint64_t het_b = 0u;
            uint64_t hom_b_b = 0u;
            duckhts_somalier_masks_t masks_b =
                mask_view(&hom_a_b, &het_b, &hom_b_b, 3u, 1u, 0u);
            duckhts_somalier_pair_stats_t actual;
            duckhts_somalier_pair_stats_t expected;
            for (site = 0u; site < 3u; site++) {
                genotypes_b[site] = profile_genotype(pb, site);
                counts_b[site] = counts_for_genotype(genotypes_b[site]);
            }
            CHECK(duckhts_somalier_prepare_masks(counts_b, 3u, test_panel_sha256,
                  &settings, &masks_b) ==
                  DUCKHTS_SOMALIER_OK);
            CHECK(duckhts_somalier_pair_stats(&masks_a, &masks_b, &actual) ==
                  DUCKHTS_SOMALIER_OK);
            reference_pair(genotypes_a, genotypes_b, 3u, 0u, 0u, &expected);
            compare_pair_stats(&actual, &expected);
        }
    }
}

static void test_classification_contracts(void) {
    duckhts_somalier_relatedness_settings_t related;
    duckhts_somalier_contamination_settings_t contamination;
    duckhts_somalier_genotype_t genotype;
    duckhts_somalier_counts_t counts;
    duckhts_somalier_relatedness_settings_default(&related);
    duckhts_somalier_charr_settings_default(&contamination);
    related.min_depth = 1u;
    contamination.min_depth = 1u;

    counts = (duckhts_somalier_counts_t){9u, 0u, 1u, 1u};
    CHECK(duckhts_somalier_classify_relatedness(&counts, &related, &genotype) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(genotype == DUCKHTS_SOMALIER_HOM_A);
    counts.other = 2u;
    CHECK(duckhts_somalier_classify_relatedness(&counts, &related, &genotype) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(genotype == DUCKHTS_SOMALIER_UNKNOWN);

    counts = (duckhts_somalier_counts_t){20u, 0u, 0u, 0u};
    CHECK(duckhts_somalier_classify_relatedness(&counts, &related, &genotype) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(genotype == DUCKHTS_SOMALIER_UNKNOWN);
    counts.available = 2u;
    CHECK(duckhts_somalier_classify_relatedness(&counts, &related, &genotype) ==
          DUCKHTS_SOMALIER_INVALID_ARGUMENT);

    counts = (duckhts_somalier_counts_t){24u, 0u, 1u, 1u};
    CHECK(duckhts_somalier_classify_contamination(&counts, &contamination, &genotype) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(genotype == DUCKHTS_SOMALIER_HOM_A);
    counts.other = 2u;
    CHECK(duckhts_somalier_classify_contamination(&counts, &contamination, &genotype) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(genotype == DUCKHTS_SOMALIER_UNKNOWN);

    counts = (duckhts_somalier_counts_t){94u, 6u, 0u, 1u};
    related.min_het_balance = 0.30;
    CHECK(duckhts_somalier_classify_relatedness(&counts, &related, &genotype) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(genotype == DUCKHTS_SOMALIER_UNKNOWN);

    counts = (duckhts_somalier_counts_t){UINT32_MAX, UINT32_MAX, 0u, 1u};
    CHECK(duckhts_somalier_classify_relatedness(&counts, &related, &genotype) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(genotype == DUCKHTS_SOMALIER_HET);
    contamination.max_depth = UINT32_MAX;
    counts = (duckhts_somalier_counts_t){UINT32_MAX, 0u, UINT32_MAX, 1u};
    CHECK(duckhts_somalier_classify_contamination(&counts, &contamination, &genotype) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(genotype == DUCKHTS_SOMALIER_UNKNOWN);

    {
        duckhts_somalier_counts_t middling[3] = {
            {94u, 6u, 0u, 1u}, {6u, 94u, 0u, 1u}, {99u, 1u, 0u, 1u}
        };
        uint64_t hom_a = 0u;
        uint64_t het = 0u;
        uint64_t hom_b = 0u;
        duckhts_somalier_masks_t masks =
            mask_view(&hom_a, &het, &hom_b, 3u, 1u, UINT64_MAX);
        CHECK(duckhts_somalier_prepare_masks(middling, 3u, test_panel_sha256,
              &related, &masks) ==
              DUCKHTS_SOMALIER_OK);
        CHECK(masks.middling_balance_count == 2u);
    }

    {
        duckhts_somalier_counts_t observations[2] = {
            {20u, 0u, 0u, 1u}, {0u, 0u, 0u, 0u}
        };
        uint64_t hom_a = UINT64_MAX;
        uint64_t het = UINT64_MAX;
        uint64_t hom_b = UINT64_MAX;
        duckhts_somalier_masks_t masks =
            mask_view(&hom_a, &het, &hom_b, 2u, 1u, 7u);
        CHECK(duckhts_somalier_prepare_masks(observations, 2u, test_panel_sha256,
              &related, &masks) == DUCKHTS_SOMALIER_OK);
        CHECK(masks.unavailable_count == 1u);
        CHECK(hom_a == 1u && het == 0u && hom_b == 0u);
        observations[1].available = 2u;
        hom_a = UINT64_C(0x55);
        masks.middling_balance_count = 9u;
        CHECK(duckhts_somalier_prepare_masks(observations, 2u, test_panel_sha256,
              &related, &masks) == DUCKHTS_SOMALIER_INVALID_ARGUMENT);
        CHECK(hom_a == UINT64_C(0x55));
        CHECK(masks.middling_balance_count == 9u);
    }
}

static void test_fixed_concordance_witness(void) {
    static const int genotypes_a[5] = {
        DUCKHTS_SOMALIER_HOM_A, DUCKHTS_SOMALIER_HOM_B, DUCKHTS_SOMALIER_HET,
        DUCKHTS_SOMALIER_HOM_A, DUCKHTS_SOMALIER_UNKNOWN
    };
    static const int genotypes_b[5] = {
        DUCKHTS_SOMALIER_HOM_A, DUCKHTS_SOMALIER_HOM_B, DUCKHTS_SOMALIER_HET,
        DUCKHTS_SOMALIER_HOM_A, DUCKHTS_SOMALIER_HOM_B
    };
    uint64_t hom_a_a = UINT64_C(0x09);
    uint64_t het_a = UINT64_C(0x04);
    uint64_t hom_b_a = UINT64_C(0x02);
    uint64_t hom_a_b = UINT64_C(0x09);
    uint64_t het_b = UINT64_C(0x04);
    uint64_t hom_b_b = UINT64_C(0x12);
    duckhts_somalier_masks_t a =
        mask_view(&hom_a_a, &het_a, &hom_b_a, 5u, 1u, 1u);
    duckhts_somalier_masks_t b =
        mask_view(&hom_a_b, &het_b, &hom_b_b, 5u, 1u, 2u);
    duckhts_somalier_pair_stats_t actual;
    duckhts_somalier_pair_stats_t expected;
    CHECK(duckhts_somalier_seal_masks(&a) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_seal_masks(&b) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_pair_stats(&a, &b, &actual) == DUCKHTS_SOMALIER_OK);
    reference_pair(genotypes_a, genotypes_b, 5u, 1u, 2u, &expected);
    compare_pair_stats(&actual, &expected);
    CHECK(actual.ibs0 == 0u && actual.ibs2 == 4u && actual.jointly_called == 4u);
    CHECK(close_enough(actual.inferred_hom_concordance, 1.0, 1e-14));
    CHECK(close_enough(actual.raw_hom_b_concordance, 1.0, 1e-14));
    CHECK(close_enough(actual.p_middling_a, 1.0 / 6.0, 1e-14));
    CHECK(close_enough(actual.p_middling_b, 2.0 / 7.0, 1e-14));
    CHECK(close_enough(actual.adjusted_concordance, 0.94417296, 1e-7));

    b.identity.panel_sha256[0] ^= 1u;
    CHECK(duckhts_somalier_validate_masks(&b) == DUCKHTS_SOMALIER_CORRUPT_MASK);
    CHECK(duckhts_somalier_seal_masks(&b) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_pair_stats(&a, &b, &actual) ==
          DUCKHTS_SOMALIER_IDENTITY_MISMATCH);
    b.identity.panel_sha256[0] ^= 1u;
    b.identity.min_depth++;
    CHECK(duckhts_somalier_seal_masks(&b) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_pair_stats(&a, &b, &actual) ==
          DUCKHTS_SOMALIER_IDENTITY_MISMATCH);
}

static void test_content_identity_and_pair_verifier(void) {
    duckhts_somalier_relatedness_settings_t settings;
    duckhts_somalier_counts_t evidence_a[3] = {
        {20u, 0u, 0u, 1u}, {10u, 10u, 0u, 1u}, {0u, 0u, 0u, 0u}
    };
    duckhts_somalier_counts_t evidence_b[3] = {
        {20u, 0u, 0u, 1u}, {0u, 20u, 0u, 1u}, {0u, 0u, 0u, 0u}
    };
    duckhts_somalier_counts_t unavailable[3] = {
        {0u, 0u, 0u, 0u}, {0u, 0u, 0u, 0u}, {0u, 0u, 0u, 0u}
    };
    uint64_t hom_a_a = 0u, het_a = 0u, hom_b_a = 0u;
    uint64_t hom_a_b = 0u, het_b = 0u, hom_b_b = 0u;
    duckhts_somalier_masks_t a =
        mask_view(&hom_a_a, &het_a, &hom_b_a, 3u, 1u, 0u);
    duckhts_somalier_masks_t b =
        mask_view(&hom_a_b, &het_b, &hom_b_b, 3u, 1u, 0u);
    duckhts_somalier_pair_stats_t actual;
    duckhts_somalier_pair_stats_t corrupted;
    duckhts_somalier_relatedness_settings_default(&settings);
    CHECK(duckhts_somalier_prepare_masks(evidence_a, 3u, test_panel_sha256,
          &settings, &a) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_prepare_masks(evidence_b, 3u, test_panel_sha256,
          &settings, &b) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_validate_masks(&a) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_validate_masks(&b) == DUCKHTS_SOMALIER_OK);

    a.middling_balance_count = 1u;
    CHECK(duckhts_somalier_validate_masks(&a) == DUCKHTS_SOMALIER_CORRUPT_MASK);
    a.middling_balance_count = 0u;
    a.unavailable_count = 0u;
    CHECK(duckhts_somalier_validate_masks(&a) == DUCKHTS_SOMALIER_CORRUPT_MASK);
    a.unavailable_count = 1u;
    a.identity.panel_sha256[0] ^= 1u;
    CHECK(duckhts_somalier_validate_masks(&a) == DUCKHTS_SOMALIER_CORRUPT_MASK);
    a.identity.panel_sha256[0] ^= 1u;
    a.identity.min_depth++;
    CHECK(duckhts_somalier_validate_masks(&a) == DUCKHTS_SOMALIER_CORRUPT_MASK);
    a.identity.min_depth--;
    a.identity.min_het_balance = 0.25;
    CHECK(duckhts_somalier_validate_masks(&a) == DUCKHTS_SOMALIER_CORRUPT_MASK);
    a.identity.min_het_balance = settings.min_het_balance;
    a.identity.hom_balance_cutoff = 0.005;
    CHECK(duckhts_somalier_validate_masks(&a) == DUCKHTS_SOMALIER_CORRUPT_MASK);
    a.identity.hom_balance_cutoff = settings.hom_balance_cutoff;
    hom_a_a |= UINT64_C(1) << 2u;
    CHECK(duckhts_somalier_validate_masks(&a) == DUCKHTS_SOMALIER_CORRUPT_MASK);
    hom_a_a &= ~(UINT64_C(1) << 2u);
    a.content_digest[0] ^= 1u;
    CHECK(duckhts_somalier_validate_masks(&a) == DUCKHTS_SOMALIER_CORRUPT_MASK);
    a.content_digest[0] ^= 1u;
    CHECK(duckhts_somalier_validate_masks(&a) == DUCKHTS_SOMALIER_OK);
    a.count_digest[1] ^= 1u;
    CHECK(duckhts_somalier_validate_masks(&a) == DUCKHTS_SOMALIER_CORRUPT_MASK);
    a.count_digest[1] ^= 1u;
    CHECK(duckhts_somalier_validate_masks(&a) == DUCKHTS_SOMALIER_OK);

    CHECK(duckhts_somalier_pair_stats(&a, &b, &actual) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_verify_pair_result(&a, &b, &actual,
          DUCKHTS_SOMALIER_OK) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_verify_pair_result(&a, &b, &actual,
          DUCKHTS_SOMALIER_NO_EVIDENCE) == DUCKHTS_SOMALIER_CORRUPT_RESULT);
#define CHECK_CORRUPT_INTEGER(field) do { \
    corrupted = actual; \
    corrupted.field++; \
    CHECK(duckhts_somalier_verify_pair_result(&a, &b, &corrupted, \
          DUCKHTS_SOMALIER_OK) == DUCKHTS_SOMALIER_CORRUPT_RESULT); \
} while (0)
#define CHECK_CORRUPT_FLOAT(field) do { \
    corrupted = actual; \
    corrupted.field += 0.125; \
    CHECK(duckhts_somalier_verify_pair_result(&a, &b, &corrupted, \
          DUCKHTS_SOMALIER_OK) == DUCKHTS_SOMALIER_CORRUPT_RESULT); \
} while (0)
    CHECK_CORRUPT_INTEGER(ibs0);
    CHECK_CORRUPT_INTEGER(ibs2);
    CHECK_CORRUPT_INTEGER(jointly_called);
    CHECK_CORRUPT_INTEGER(shared_hets);
    CHECK_CORRUPT_INTEGER(shared_hom_b);
    CHECK_CORRUPT_INTEGER(het_ab);
    CHECK_CORRUPT_INTEGER(het_count_a);
    CHECK_CORRUPT_INTEGER(het_count_b);
    CHECK_CORRUPT_INTEGER(hom_b_count_a);
    CHECK_CORRUPT_INTEGER(hom_b_count_b);
    CHECK_CORRUPT_INTEGER(callable_hom_count_a);
    CHECK_CORRUPT_INTEGER(callable_hom_count_b);
    CHECK_CORRUPT_INTEGER(matching_hom_count);
    CHECK_CORRUPT_FLOAT(relatedness);
    CHECK_CORRUPT_FLOAT(inferred_hom_concordance);
    CHECK_CORRUPT_FLOAT(raw_hom_b_concordance);
    CHECK_CORRUPT_FLOAT(p_middling_a);
    CHECK_CORRUPT_FLOAT(p_middling_b);
    CHECK_CORRUPT_FLOAT(adjusted_concordance);
#undef CHECK_CORRUPT_INTEGER
#undef CHECK_CORRUPT_FLOAT

    CHECK(duckhts_somalier_prepare_masks(unavailable, 3u, test_panel_sha256,
          &settings, &a) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_prepare_masks(unavailable, 3u, test_panel_sha256,
          &settings, &b) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_pair_stats(&a, &b, &actual) == DUCKHTS_SOMALIER_OK);
    CHECK(actual.jointly_called == 0u);
    CHECK(duckhts_somalier_verify_pair_result(&a, &b, &actual,
          DUCKHTS_SOMALIER_NO_EVIDENCE) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_verify_pair_result(&a, &b, &actual,
          DUCKHTS_SOMALIER_OK) == DUCKHTS_SOMALIER_CORRUPT_RESULT);
}

static void test_count_digest_receipt(void) {
    duckhts_somalier_counts_t observations[2] = {
        {20u, 0u, 0u, 1u}, {10u, 10u, 1u, 1u}
    };
    duckhts_somalier_counts_t same_genotype = {19u, 0u, 0u, 1u};
    duckhts_somalier_counts_t unavailable = {0u, 0u, 0u, 0u};
    duckhts_somalier_counts_t measured_zero = {0u, 0u, 0u, 1u};
    duckhts_somalier_counts_t invalid = {0u, 0u, 0u, 2u};
    duckhts_somalier_relatedness_settings_t settings;
    uint64_t direct[2] = {0u, 0u};
    uint64_t reverse[2] = {0u, 0u};
    uint64_t combined[2] = {0u, 0u};
    uint64_t right[2] = {0u, 0u};
    uint64_t unavailable_digest[2] = {0u, 0u};
    uint64_t measured_zero_digest[2] = {0u, 0u};
    uint64_t hom_a_a = 0u, het_a = 0u, hom_b_a = 0u;
    uint64_t hom_a_b = 0u, het_b = 0u, hom_b_b = 0u;
    duckhts_somalier_masks_t a =
        mask_view(&hom_a_a, &het_a, &hom_b_a, 1u, 1u, 0u);
    duckhts_somalier_masks_t b =
        mask_view(&hom_a_b, &het_b, &hom_b_b, 1u, 1u, 0u);

    CHECK(duckhts_somalier_count_digest_observe(direct, 0u,
          &observations[0]) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_count_digest_observe(direct, 1u,
          &observations[1]) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_count_digest_observe(reverse, 1u,
          &observations[1]) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_count_digest_observe(reverse, 0u,
          &observations[0]) == DUCKHTS_SOMALIER_OK);
    CHECK(memcmp(direct, reverse, sizeof(direct)) == 0);

    CHECK(duckhts_somalier_count_digest_observe(combined, 0u,
          &observations[0]) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_count_digest_observe(right, 1u,
          &observations[1]) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_count_digest_combine(combined, right) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(memcmp(direct, combined, sizeof(direct)) == 0);
    CHECK(duckhts_somalier_count_digest_observe(NULL, 0u,
          &observations[0]) == DUCKHTS_SOMALIER_INVALID_ARGUMENT);
    CHECK(duckhts_somalier_count_digest_observe(direct, 0u, &invalid) ==
          DUCKHTS_SOMALIER_INVALID_ARGUMENT);

    CHECK(duckhts_somalier_count_digest_observe(unavailable_digest, 0u,
          &unavailable) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_count_digest_observe(measured_zero_digest, 0u,
          &measured_zero) == DUCKHTS_SOMALIER_OK);
    CHECK(memcmp(unavailable_digest, measured_zero_digest,
                 sizeof(unavailable_digest)) != 0);

    duckhts_somalier_relatedness_settings_default(&settings);
    CHECK(duckhts_somalier_prepare_masks(observations, 1u,
          test_panel_sha256, &settings, &a) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_prepare_masks(&same_genotype, 1u,
          test_panel_sha256, &settings, &b) == DUCKHTS_SOMALIER_OK);
    CHECK(hom_a_a == 1u && hom_a_b == 1u);
    CHECK(memcmp(a.count_digest, b.count_digest,
                 sizeof(a.count_digest)) != 0);
    CHECK(memcmp(a.content_digest, b.content_digest,
                 sizeof(a.content_digest)) != 0);
}

static void test_mask_shapes(void) {
    static const size_t counts[] = {63u, 64u, 65u};
    duckhts_somalier_relatedness_settings_t settings;
    size_t case_index;
    duckhts_somalier_relatedness_settings_default(&settings);
    for (case_index = 0u; case_index < sizeof(counts) / sizeof(counts[0]); case_index++) {
        duckhts_somalier_counts_t input[65];
        uint64_t hom_a[2] = {UINT64_MAX, UINT64_MAX};
        uint64_t het[2] = {UINT64_MAX, UINT64_MAX};
        uint64_t hom_b[2] = {UINT64_MAX, UINT64_MAX};
        size_t words;
        size_t i;
        duckhts_somalier_masks_t masks;
        CHECK(duckhts_somalier_mask_word_count(counts[case_index], &words) ==
              DUCKHTS_SOMALIER_OK);
        CHECK(words == (counts[case_index] == 65u ? 2u : 1u));
        masks = mask_view(hom_a, het, hom_b, counts[case_index], words, UINT64_MAX);
        for (i = 0u; i < counts[case_index]; i++) {
            input[i] = counts_for_genotype(DUCKHTS_SOMALIER_HOM_A);
        }
        CHECK(duckhts_somalier_prepare_masks(input, counts[case_index],
              test_panel_sha256, &settings, &masks) ==
              DUCKHTS_SOMALIER_OK);
        CHECK(masks.middling_balance_count == 0u);
        CHECK(duckhts_somalier_validate_masks(&masks) == DUCKHTS_SOMALIER_OK);
        if (counts[case_index] % 64u != 0u) {
            size_t last = words - 1u;
            hom_a[last] |= UINT64_C(1) << (counts[case_index] % 64u);
            CHECK(duckhts_somalier_validate_masks(&masks) ==
                  DUCKHTS_SOMALIER_CORRUPT_MASK);
            hom_a[last] &= ~(UINT64_C(1) << (counts[case_index] % 64u));
        }
        het[0] = hom_a[0] & 1u;
        CHECK(duckhts_somalier_validate_masks(&masks) == DUCKHTS_SOMALIER_CORRUPT_MASK);
    }
    {
        duckhts_somalier_counts_t input[65];
        uint64_t storage[6] = {0u, 0u, 0u, 0u, 0u, 0u};
        duckhts_somalier_masks_t adjacent =
            mask_view(storage, storage + 2, storage + 4, 65u, 2u, 0u);
        duckhts_somalier_masks_t overlapping =
            mask_view(storage, storage + 1, storage + 3, 65u, 2u, 0u);
        size_t i;
        for (i = 0u; i < 65u; i++) input[i] = counts_for_genotype(DUCKHTS_SOMALIER_HOM_A);
        CHECK(duckhts_somalier_seal_masks(&adjacent) == DUCKHTS_SOMALIER_OK);
        CHECK(duckhts_somalier_validate_masks(&adjacent) == DUCKHTS_SOMALIER_OK);
        CHECK(duckhts_somalier_validate_masks(&overlapping) ==
              DUCKHTS_SOMALIER_CORRUPT_MASK);
        CHECK(duckhts_somalier_prepare_masks(input, 65u, test_panel_sha256,
              &settings, &overlapping) ==
              DUCKHTS_SOMALIER_CORRUPT_MASK);
    }
}

static void test_binomial_and_charr(void) {
    duckhts_somalier_contamination_settings_t settings;
    duckhts_somalier_genotype_t genotype;
    uint64_t threshold;
    duckhts_somalier_counts_t one;
    double af;
    duckhts_somalier_charr_result_t result;
    duckhts_somalier_charr_settings_default(&settings);

    CHECK(duckhts_somalier_binomial_max_minor(1u, 0.05, 0.05, 1000000u, &threshold) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(threshold == 1u);
    CHECK(duckhts_somalier_binomial_max_minor(
          1u, 0.05, nextafter(0.05, 0.0), 1000000u, &threshold) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(threshold == 1u);
    CHECK(duckhts_somalier_binomial_max_minor(
          1u, 0.05, nextafter(0.05, 1.0), 1000000u, &threshold) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(threshold == 0u);
    CHECK(duckhts_somalier_binomial_max_minor(
          3u, 0.25, nextafter(0.578125, 0.0), 1000000u, &threshold) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(threshold == 1u);
    CHECK(duckhts_somalier_binomial_max_minor(
          3u, 0.25, 0.578125, 1000000u, &threshold) == DUCKHTS_SOMALIER_OK);
    CHECK(threshold == 1u);
    CHECK(duckhts_somalier_binomial_max_minor(
          3u, 0.25, nextafter(0.578125, 1.0), 1000000u, &threshold) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(threshold == 0u);
    settings.min_depth = 1u;
    settings.hom_minor_rate = 0.25;
    settings.hom_tail_alpha = nextafter(0.578125, 1.0);
    one = (duckhts_somalier_counts_t){2u, 1u, 0u, 1u};
    CHECK(duckhts_somalier_classify_contamination(&one, &settings, &genotype) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(genotype == DUCKHTS_SOMALIER_UNKNOWN);
    settings.hom_tail_alpha = 0.578125;
    CHECK(duckhts_somalier_classify_contamination(&one, &settings, &genotype) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(genotype == DUCKHTS_SOMALIER_HOM_A);
    CHECK(duckhts_somalier_binomial_max_minor(
          1000000u, 0.0, 0.002, 1000000u, &threshold) == DUCKHTS_SOMALIER_OK);
    CHECK(threshold == 0u);
    CHECK(duckhts_somalier_binomial_max_minor(
          1000000u, 1.0, 0.002, 1000000u, &threshold) == DUCKHTS_SOMALIER_OK);
    CHECK(threshold == 1000000u);

    CHECK(duckhts_somalier_binomial_max_minor(100u, 0.12, 0.002, 6000u, &threshold) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(threshold == 22u);
    CHECK(duckhts_somalier_binomial_max_minor(200u, 0.12, 0.002, 6000u, &threshold) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(threshold == 38u);
    CHECK(duckhts_somalier_binomial_max_minor(6000u, 0.12, 0.002, 6000u, &threshold) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(threshold == 793u);
    CHECK(duckhts_somalier_binomial_max_minor(
          750000u, 0.49, 0.002, 1000000u, &threshold) == DUCKHTS_SOMALIER_OK);
    CHECK(threshold == 368746u);
    CHECK(duckhts_somalier_binomial_max_minor(
          1000000u, 0.12, 0.5, 1000000u, &threshold) == DUCKHTS_SOMALIER_OK);
    CHECK(threshold == 120000u);
    CHECK(duckhts_somalier_binomial_max_minor(
          1000000u, 0.49, 0.002, 1000000u, &threshold) == DUCKHTS_SOMALIER_OK);
    CHECK(threshold == 491439u);
    CHECK(duckhts_somalier_binomial_max_minor(
          1000000u, 0.01, 0.000001, 1000000u, &threshold) == DUCKHTS_SOMALIER_OK);
    CHECK(threshold == 10476u);
    CHECK(duckhts_somalier_binomial_max_minor(
          DUCKHTS_SOMALIER_MAX_BINOMIAL_DEPTH + 1u, 0.12, 0.002,
          DUCKHTS_SOMALIER_MAX_BINOMIAL_DEPTH + 1u, &threshold) ==
          DUCKHTS_SOMALIER_LIMIT_EXCEEDED);
    CHECK(duckhts_somalier_binomial_max_minor(
          100u, 0.12, 1.0, 100u, &threshold) ==
          DUCKHTS_SOMALIER_INVALID_ARGUMENT);

    settings.min_depth = 7u;
    settings.hom_minor_rate = 0.10;
    settings.hom_tail_alpha = 0.001;
    one = (duckhts_somalier_counts_t){199u, 1u, 0u, 1u};
    af = 0.25;
    CHECK(duckhts_somalier_charr(&one, &af, 1u, &settings, &result) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(close_enough(result.estimate, 0.02, 1e-12));
    CHECK(result.usable_sites == 1u && result.usable_hom_a == 1u);

    duckhts_somalier_charr_settings_default(&settings);
    one = (duckhts_somalier_counts_t){4000u, 2000u, 0u, 1u};
    af = 0.5;
    CHECK(duckhts_somalier_charr(&one, &af, 1u, &settings, &result) ==
          DUCKHTS_SOMALIER_NO_EVIDENCE);
    CHECK(result.status == DUCKHTS_SOMALIER_NO_EVIDENCE && isnan(result.estimate));

    one = (duckhts_somalier_counts_t){50u, 50u, 0u, 1u};
    CHECK(duckhts_somalier_charr(&one, &af, 1u, &settings, &result) ==
          DUCKHTS_SOMALIER_NO_EVIDENCE);

    {
        duckhts_somalier_counts_t two[2] = {
            {199u, 1u, 0u, 1u}, {199u, 1u, 0u, 1u}
        };
        double frequencies[2] = {0.25, NAN};
        CHECK(duckhts_somalier_charr(two, frequencies, 2u, &settings, &result) ==
              DUCKHTS_SOMALIER_INVALID_ARGUMENT);
        CHECK(result.status == DUCKHTS_SOMALIER_INVALID_ARGUMENT);
        CHECK(isnan(result.estimate) && result.usable_sites == 0u);
    }
}

static uint64_t reference_binomial_max_minor(unsigned depth, double probability,
                                             double alpha) {
    long double p = (long double)probability;
    long double q = 1.0L - p;
    long double pmf = powl(q, (long double)depth);
    long double tail = 1.0L;
    uint64_t result = 0u;
    unsigned k;
    for (k = 0u; k <= depth; k++) {
        if (tail >= (long double)alpha) {
            result = k;
        } else {
            break;
        }
        if (k == depth) break;
        tail -= pmf;
        if (tail < 0.0L) tail = 0.0L;
        pmf *= ((long double)(depth - k) / (long double)(k + 1u)) * p / q;
    }
    return result;
}

static void test_random_binomial_differential(void) {
    static const double alphas[] = {0.001, 0.002, 0.01, 0.05};
    unsigned trial;
    for (trial = 0u; trial < 1000u; trial++) {
        unsigned depth = 1u + (unsigned)(next_random() % 300u);
        double probability = (double)(1u + next_random() % 49u) / 100.0;
        double alpha = alphas[next_random() % 4u];
        uint64_t expected = reference_binomial_max_minor(depth, probability, alpha);
        uint64_t actual = UINT64_MAX;
        CHECK(duckhts_somalier_binomial_max_minor(depth, probability, alpha,
              300u, &actual) == DUCKHTS_SOMALIER_OK);
        CHECK(actual == expected);
    }
}

static void test_charr_orientation(void) {
    duckhts_somalier_contamination_settings_t settings;
    duckhts_somalier_counts_t original[2] = {
        {199u, 1u, 0u, 1u}, {2u, 198u, 0u, 1u}
    };
    duckhts_somalier_counts_t flipped[2] = {
        {1u, 199u, 0u, 1u}, {198u, 2u, 0u, 1u}
    };
    double af[2] = {0.25, 0.80};
    double flipped_af[2] = {0.75, 0.20};
    duckhts_somalier_charr_result_t a;
    duckhts_somalier_charr_result_t b;
    duckhts_somalier_charr_settings_default(&settings);
    settings.min_depth = 7u;
    settings.hom_minor_rate = 0.10;
    settings.hom_tail_alpha = 0.001;
    CHECK(duckhts_somalier_charr(original, af, 2u, &settings, &a) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_charr(flipped, flipped_af, 2u, &settings, &b) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(close_enough(a.estimate, b.estimate, 1e-14));
    CHECK(a.usable_sites == b.usable_sites);
}

static void test_charr_reduction_order(void) {
    enum { SITE_COUNT = 31 };
    duckhts_somalier_contamination_settings_t settings;
    duckhts_somalier_charr_accumulator_t forward = {0};
    duckhts_somalier_charr_accumulator_t reverse = {0};
    duckhts_somalier_charr_accumulator_t even = {0};
    duckhts_somalier_charr_accumulator_t odd = {0};
    duckhts_somalier_charr_accumulator_t combined_left;
    duckhts_somalier_charr_accumulator_t combined_right;
    duckhts_somalier_charr_result_t result[4];
    unsigned i;

    duckhts_somalier_charr_settings_default(&settings);
    settings.min_depth = 7u;
    for (i = 0u; i < SITE_COUNT; i++) {
        duckhts_somalier_counts_t counts = {
            199u - i % 3u, 1u + i % 3u, 0u, 1u
        };
        double frequency = (double)(10u + i) / 100.0;
        CHECK(duckhts_somalier_charr_observe(&forward, &counts, frequency,
              &settings) == DUCKHTS_SOMALIER_OK);
        CHECK(duckhts_somalier_charr_observe(i % 2u == 0u ? &even : &odd,
              &counts, frequency, &settings) == DUCKHTS_SOMALIER_OK);
    }
    for (i = SITE_COUNT; i > 0u; i--) {
        unsigned site = i - 1u;
        duckhts_somalier_counts_t counts = {
            199u - site % 3u, 1u + site % 3u, 0u, 1u
        };
        double frequency = (double)(10u + site) / 100.0;
        CHECK(duckhts_somalier_charr_observe(&reverse, &counts, frequency,
              &settings) == DUCKHTS_SOMALIER_OK);
    }
    combined_left = even;
    combined_right = odd;
    CHECK(duckhts_somalier_charr_combine(&combined_left, &odd) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_charr_combine(&combined_right, &even) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_charr_finish(&forward, &result[0]) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_charr_finish(&reverse, &result[1]) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_charr_finish(&combined_left, &result[2]) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_charr_finish(&combined_right, &result[3]) ==
          DUCKHTS_SOMALIER_OK);
    for (i = 1u; i < 4u; i++) {
        CHECK(result[i].estimate == result[0].estimate);
        CHECK(result[i].usable_sites == result[0].usable_sites);
        CHECK(result[i].usable_hom_a == result[0].usable_hom_a);
        CHECK(result[i].usable_hom_b == result[0].usable_hom_b);
    }
    {
        duckhts_somalier_charr_accumulator_t maximum = {0};
        duckhts_somalier_charr_accumulator_t one = {0};
        maximum.contribution_scaled_high = UINT64_MAX;
        one.contribution_scaled_high = 1u;
        CHECK(duckhts_somalier_charr_combine(&maximum, &one) ==
              DUCKHTS_SOMALIER_LIMIT_EXCEEDED);
        CHECK(maximum.contribution_scaled_high == UINT64_MAX &&
              maximum.contribution_scaled_low == 0u);
    }
}

static duckhts_somalier_matched_view_t matched_view(
    const uint32_t *a, const uint32_t *b, const int8_t *anchor,
    const double *frequency, size_t count) {
    static const uint8_t all_usable[3] = {1u, 1u, 1u};
    CHECK(count <= sizeof(all_usable));
    duckhts_somalier_matched_view_t view = {
        a, b, all_usable, anchor, frequency, count
    };
    return view;
}

static void test_matched_anchor(void) {
    static const uint32_t receiver_a[2] = {160u, 40u};
    static const uint32_t receiver_b[2] = {40u, 160u};
    static const int8_t anchor[2] = {
        DUCKHTS_SOMALIER_HOM_A, DUCKHTS_SOMALIER_HOM_B
    };
    static const double af[2] = {0.5, 0.5};
    duckhts_somalier_matched_settings_t settings;
    duckhts_somalier_matched_view_t view =
        matched_view(receiver_a, receiver_b, anchor, af, 2u);
    duckhts_somalier_matched_view_t empty = {
        NULL, NULL, NULL, NULL, NULL, 0u
    };
    duckhts_somalier_matched_result_t result;
    double zero_likelihood;
    double one_likelihood;
    double upstream_likelihood;
    double grid_likelihood;
    duckhts_somalier_matched_settings_default(&settings);
    CHECK(duckhts_somalier_matched_anchor(&view, &settings, &result) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(close_enough(result.alpha, 0.3975904, 2e-6));
    CHECK(result.evaluations >= 101u && result.usable_sites == 2u);
    CHECK(duckhts_somalier_matched_log_likelihood(&view, 0.4409830056250525,
          &settings, &upstream_likelihood) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_matched_log_likelihood(&view, 0.4, &settings,
          &grid_likelihood) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_matched_log_likelihood(&view, 0.0, &settings,
          &zero_likelihood) == DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_matched_log_likelihood(&view, 1.0, &settings,
          &one_likelihood) == DUCKHTS_SOMALIER_OK);
    /* Independently evaluated likelihoods for the retained two-site witness. */
    CHECK(close_enough(grid_likelihood, -201.54905839206708, 1e-10));
    CHECK(close_enough(upstream_likelihood, -202.10199747075765, 1e-10));
    CHECK(close_enough(zero_likelihood, -497.80928872839064, 1e-10));
    CHECK(close_enough(one_likelihood, -278.64516658509802, 1e-10));
    CHECK(result.log_likelihood >= grid_likelihood);
    CHECK(result.log_likelihood > upstream_likelihood + 0.5);
    {
        unsigned i;
        for (i = 0u; i <= 100u; i++) {
            double likelihood;
            CHECK(duckhts_somalier_matched_log_likelihood(&view, (double)i / 100.0,
                  &settings, &likelihood) == DUCKHTS_SOMALIER_OK);
            CHECK(result.log_likelihood + 1e-12 >= likelihood);
        }
    }
    CHECK(duckhts_somalier_matched_anchor(&empty, &settings, &result) ==
          DUCKHTS_SOMALIER_NO_EVIDENCE);
    CHECK(result.status == DUCKHTS_SOMALIER_NO_EVIDENCE && result.evaluations == 0u);

    settings.max_evaluations = 10u;
    CHECK(duckhts_somalier_matched_anchor(&view, &settings, &result) ==
          DUCKHTS_SOMALIER_LIMIT_EXCEEDED);
    CHECK(isnan(result.alpha) && isinf(result.log_likelihood) &&
          result.log_likelihood < 0.0 && result.usable_sites == 2u &&
          result.evaluations == 0u);

    duckhts_somalier_matched_settings_default(&settings);
    settings.max_evaluations = 105u;
    CHECK(duckhts_somalier_matched_anchor(&view, &settings, &result) ==
          DUCKHTS_SOMALIER_LIMIT_EXCEEDED);
    CHECK(isnan(result.alpha) && isinf(result.log_likelihood) &&
          result.log_likelihood < 0.0 && result.usable_sites == 2u &&
          result.evaluations == 105u);

    duckhts_somalier_matched_settings_default(&settings);
    settings.grid_step = ldexp(1.0, -64);
    CHECK(duckhts_somalier_matched_anchor(&view, &settings, &result) ==
          DUCKHTS_SOMALIER_LIMIT_EXCEEDED);
    CHECK(isnan(result.alpha) && result.usable_sites == 2u);
}

static void test_matched_orientation(void) {
    static const uint32_t receiver_a[3] = {90u, 30u, 140u};
    static const uint32_t receiver_b[3] = {10u, 70u, 60u};
    static const int8_t anchor[3] = {
        DUCKHTS_SOMALIER_HOM_A, DUCKHTS_SOMALIER_HOM_B, DUCKHTS_SOMALIER_HOM_A
    };
    static const double af[3] = {0.2, 0.8, 0.35};
    uint32_t flipped_a[3];
    uint32_t flipped_b[3];
    int8_t flipped_anchor[3];
    double flipped_af[3];
    duckhts_somalier_matched_settings_t settings;
    duckhts_somalier_matched_result_t original;
    duckhts_somalier_matched_result_t flipped;
    duckhts_somalier_matched_view_t a;
    duckhts_somalier_matched_view_t b;
    size_t i;
    for (i = 0u; i < 3u; i++) {
        flipped_a[i] = receiver_b[i];
        flipped_b[i] = receiver_a[i];
        flipped_anchor[i] = anchor[i] == DUCKHTS_SOMALIER_HOM_A
            ? DUCKHTS_SOMALIER_HOM_B : DUCKHTS_SOMALIER_HOM_A;
        flipped_af[i] = 1.0 - af[i];
    }
    a = matched_view(receiver_a, receiver_b, anchor, af, 3u);
    b = matched_view(flipped_a, flipped_b, flipped_anchor, flipped_af, 3u);
    duckhts_somalier_matched_settings_default(&settings);
    CHECK(duckhts_somalier_matched_anchor(&a, &settings, &original) ==
          DUCKHTS_SOMALIER_OK);
    CHECK(duckhts_somalier_matched_anchor(&b, &settings, &flipped) ==
          DUCKHTS_SOMALIER_OK);
    /* An independent R optimize() calculation gives 0.2869014364168983.
     * The mirrored likelihoods agree to 6e-14; their maximizers differ by
     * 8.3e-9 because the objective is flat at double precision. */
    CHECK(close_enough(original.alpha, 0.2869014364168983, 1e-8));
    CHECK(close_enough(flipped.alpha, 0.2869014364168983, 1e-8));
    CHECK(close_enough(original.alpha, flipped.alpha, 1e-8));
    CHECK(close_enough(original.log_likelihood, flipped.log_likelihood, 1e-10));
}

static uint64_t random_state = UINT64_C(0x91e10da5c79e7b1d);

static uint64_t next_random(void) {
    random_state ^= random_state >> 12u;
    random_state ^= random_state << 25u;
    random_state ^= random_state >> 27u;
    return random_state * UINT64_C(2685821657736338717);
}

static void test_random_pair_differential(void) {
    enum { SITE_COUNT = 65, TRIALS = 2000 };
    unsigned trial;
    for (trial = 0u; trial < TRIALS; trial++) {
        int genotypes_a[SITE_COUNT];
        int genotypes_b[SITE_COUNT];
        uint64_t hom_a_a[2] = {0u, 0u};
        uint64_t het_a[2] = {0u, 0u};
        uint64_t hom_b_a[2] = {0u, 0u};
        uint64_t hom_a_b[2] = {0u, 0u};
        uint64_t het_b[2] = {0u, 0u};
        uint64_t hom_b_b[2] = {0u, 0u};
        duckhts_somalier_masks_t masks_a = mask_view(
            hom_a_a, het_a, hom_b_a, SITE_COUNT, 2u, next_random() % 8u);
        duckhts_somalier_masks_t masks_b = mask_view(
            hom_a_b, het_b, hom_b_b, SITE_COUNT, 2u, next_random() % 8u);
        duckhts_somalier_pair_stats_t actual;
        duckhts_somalier_pair_stats_t expected;
        size_t site;
        for (site = 0u; site < SITE_COUNT; site++) {
            uint64_t bit = UINT64_C(1) << (site % 64u);
            genotypes_a[site] = (int)(next_random() & 3u) - 1;
            genotypes_b[site] = (int)(next_random() & 3u) - 1;
            if (genotypes_a[site] == DUCKHTS_SOMALIER_HOM_A) hom_a_a[site / 64u] |= bit;
            if (genotypes_a[site] == DUCKHTS_SOMALIER_HET) het_a[site / 64u] |= bit;
            if (genotypes_a[site] == DUCKHTS_SOMALIER_HOM_B) hom_b_a[site / 64u] |= bit;
            if (genotypes_b[site] == DUCKHTS_SOMALIER_HOM_A) hom_a_b[site / 64u] |= bit;
            if (genotypes_b[site] == DUCKHTS_SOMALIER_HET) het_b[site / 64u] |= bit;
            if (genotypes_b[site] == DUCKHTS_SOMALIER_HOM_B) hom_b_b[site / 64u] |= bit;
        }
        CHECK(duckhts_somalier_seal_masks(&masks_a) == DUCKHTS_SOMALIER_OK);
        CHECK(duckhts_somalier_seal_masks(&masks_b) == DUCKHTS_SOMALIER_OK);
        CHECK(duckhts_somalier_pair_stats(&masks_a, &masks_b, &actual) ==
              DUCKHTS_SOMALIER_OK);
        reference_pair(genotypes_a, genotypes_b, SITE_COUNT,
                       masks_a.middling_balance_count,
                       masks_b.middling_balance_count, &expected);
        compare_pair_stats(&actual, &expected);
    }
}

/* Statistical campaigns pass named, fixed-width TSV rows through the same
 * production kernel as the ordinary native tests. No input is skipped. */
static int parse_unsigned(const char *text, uint64_t *value) {
    char *end;
    unsigned long long parsed;
    errno = 0;
    parsed = strtoull(text, &end, 10);
    if (text[0] == '-' || text == end || *end != '\0' || errno != 0) return 0;
    *value = (uint64_t)parsed;
    return 1;
}

static int parse_signed(const char *text, int *value) {
    char *end;
    long parsed;
    errno = 0;
    parsed = strtol(text, &end, 10);
    if (text == end || *end != '\0' || errno != 0 || parsed < -1 || parsed > 2) {
        return 0;
    }
    *value = (int)parsed;
    return 1;
}

static int parse_real(const char *text, double *value) {
    char *end;
    errno = 0;
    *value = strtod(text, &end);
    return text != end && *end == '\0' && errno != ERANGE;
}

static size_t split_tsv(char *line, char **fields, size_t capacity) {
    size_t count = 0u;
    char *cursor = line;
    while (count < capacity) {
        char *tab = strchr(cursor, '\t');
        fields[count++] = cursor;
        if (tab == NULL) break;
        if (count == capacity) return capacity + 1u;
        *tab = '\0';
        cursor = tab + 1;
    }
    return count;
}

static void print_real(double value) {
    if (isnan(value)) {
        fputs("NaN", stdout);
    } else if (isinf(value)) {
        fputs(value < 0.0 ? "-Inf" : "Inf", stdout);
    } else {
        printf("%.17g", value);
    }
}

static int campaign_binomial(char **fields) {
    uint64_t depth, max_depth, result = 0u;
    double rate, tail;
    duckhts_somalier_status_t status;
    if (!parse_unsigned(fields[1], &depth) || !parse_real(fields[2], &rate) ||
        !parse_real(fields[3], &tail) || !parse_unsigned(fields[4], &max_depth)) {
        return 0;
    }
    status = duckhts_somalier_binomial_max_minor(depth, rate, tail, max_depth, &result);
    printf("%s\t%s\t", fields[0], duckhts_somalier_status_string(status));
    if (status == DUCKHTS_SOMALIER_OK) printf("%" PRIu64, result);
    else fputs("NA", stdout);
    putchar('\n');
    return 1;
}

static int campaign_eligibility(char **fields) {
    uint64_t a, b, other, available, min_depth, max_depth;
    double frequency, rate, tail;
    duckhts_somalier_counts_t counts;
    duckhts_somalier_contamination_settings_t settings;
    duckhts_somalier_charr_result_t charr;
    duckhts_somalier_genotype_t genotype = DUCKHTS_SOMALIER_UNKNOWN;
    duckhts_somalier_status_t classification_status, receiver_status, charr_status;
    uint8_t receiver_usable = 0u;
    if (!parse_unsigned(fields[1], &a) || !parse_unsigned(fields[2], &b) ||
        !parse_unsigned(fields[3], &other) || !parse_unsigned(fields[4], &available) ||
        !parse_real(fields[5], &frequency) || !parse_unsigned(fields[6], &min_depth) ||
        !parse_unsigned(fields[7], &max_depth) || !parse_real(fields[8], &rate) ||
        !parse_real(fields[9], &tail) || a > UINT32_MAX || b > UINT32_MAX ||
        other > UINT32_MAX || available > UINT8_MAX) {
        return 0;
    }
    counts = (duckhts_somalier_counts_t){(uint32_t)a, (uint32_t)b,
        (uint32_t)other, (uint8_t)available};
    duckhts_somalier_charr_settings_default(&settings);
    settings.min_depth = min_depth;
    settings.max_depth = max_depth;
    settings.hom_minor_rate = rate;
    settings.hom_tail_alpha = tail;
    classification_status = duckhts_somalier_classify_contamination(
        &counts, &settings, &genotype);
    receiver_status = duckhts_somalier_contamination_usable(
        &counts, &settings, &receiver_usable);
    charr_status = duckhts_somalier_charr(&counts, &frequency, 1u, &settings, &charr);
    printf("%s\t%s\t%d\t%s\t%u\t%s\t", fields[0],
        duckhts_somalier_status_string(classification_status), (int)genotype,
        duckhts_somalier_status_string(receiver_status), (unsigned)receiver_usable,
        duckhts_somalier_status_string(charr_status));
    print_real(charr.estimate);
    printf("\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\n",
        charr.usable_sites, charr.usable_hom_a, charr.usable_hom_b);
    return 1;
}

static int campaign_matched(char **fields) {
    uint64_t site_count, raw_a[2], raw_b[2], raw_usable[2];
    uint64_t max_depth, max_evaluations;
    int anchor[2];
    double frequency[2], alpha, error_rate, min_probability, min_prior_frequency;
    double alpha_min, alpha_max, grid_step, refine_tolerance;
    uint32_t receiver_a[2], receiver_b[2];
    uint8_t receiver_usable[2];
    int8_t anchor_genotype[2];
    duckhts_somalier_matched_view_t view;
    duckhts_somalier_matched_settings_t settings;
    duckhts_somalier_matched_result_t fitted;
    duckhts_somalier_status_t status;
    double likelihood = NAN;
    uint64_t usable_sites = 0u;
    size_t evaluations = 0u;
    size_t i;
    if ((strcmp(fields[1], "likelihood") != 0 && strcmp(fields[1], "fit") != 0) ||
        !parse_unsigned(fields[2], &site_count) ||
        !parse_unsigned(fields[3], &raw_a[0]) ||
        !parse_unsigned(fields[4], &raw_b[0]) ||
        !parse_unsigned(fields[5], &raw_usable[0]) ||
        !parse_signed(fields[6], &anchor[0]) ||
        !parse_real(fields[7], &frequency[0]) ||
        !parse_unsigned(fields[8], &raw_a[1]) ||
        !parse_unsigned(fields[9], &raw_b[1]) ||
        !parse_unsigned(fields[10], &raw_usable[1]) ||
        !parse_signed(fields[11], &anchor[1]) ||
        !parse_real(fields[12], &frequency[1]) ||
        !parse_real(fields[13], &alpha) ||
        !parse_unsigned(fields[14], &max_depth) ||
        !parse_real(fields[15], &error_rate) ||
        !parse_real(fields[16], &min_probability) ||
        !parse_real(fields[17], &min_prior_frequency) ||
        !parse_real(fields[18], &alpha_min) ||
        !parse_real(fields[19], &alpha_max) ||
        !parse_real(fields[20], &grid_step) ||
        !parse_real(fields[21], &refine_tolerance) ||
        !parse_unsigned(fields[22], &max_evaluations) ||
        site_count == 0u || site_count > 2u || max_evaluations > SIZE_MAX) {
        return 0;
    }
    for (i = 0u; i < 2u; i++) {
        if (raw_a[i] > UINT32_MAX || raw_b[i] > UINT32_MAX ||
            raw_usable[i] > UINT8_MAX) return 0;
        receiver_a[i] = (uint32_t)raw_a[i];
        receiver_b[i] = (uint32_t)raw_b[i];
        receiver_usable[i] = (uint8_t)raw_usable[i];
        anchor_genotype[i] = (int8_t)anchor[i];
    }
    view = (duckhts_somalier_matched_view_t){receiver_a, receiver_b,
        receiver_usable, anchor_genotype, frequency, (size_t)site_count};
    duckhts_somalier_matched_settings_default(&settings);
    settings.max_depth = max_depth;
    settings.max_evaluations = (size_t)max_evaluations;
    settings.error_rate = error_rate;
    settings.min_probability = min_probability;
    settings.min_prior_frequency = min_prior_frequency;
    settings.alpha_min = alpha_min;
    settings.alpha_max = alpha_max;
    settings.grid_step = grid_step;
    settings.refine_tolerance = refine_tolerance;
    for (i = 0u; i < site_count; i++) {
        if (receiver_usable[i] &&
            (anchor_genotype[i] == DUCKHTS_SOMALIER_HOM_A ||
             anchor_genotype[i] == DUCKHTS_SOMALIER_HOM_B)) {
            usable_sites++;
        }
    }
    if (strcmp(fields[1], "fit") == 0) {
        status = duckhts_somalier_matched_anchor(&view, &settings, &fitted);
        alpha = fitted.alpha;
        likelihood = fitted.log_likelihood;
        usable_sites = fitted.usable_sites;
        evaluations = fitted.evaluations;
    } else {
        status = duckhts_somalier_matched_log_likelihood(
            &view, alpha, &settings, &likelihood);
    }
    printf("%s\t%s\t%s\t", fields[0], fields[1],
        duckhts_somalier_status_string(status));
    print_real(alpha);
    putchar('\t');
    print_real(likelihood);
    printf("\t%" PRIu64 "\t%zu\n", usable_sites, evaluations);
    return 1;
}

static void pinned_row(const char *case_id, const char *metric, double value,
                       duckhts_somalier_status_t status, const char *usable_sites,
                       const char *hom_a, const char *hom_b, size_t evaluations) {
    const char *status_name = status == DUCKHTS_SOMALIER_NO_EVIDENCE
        ? "no_evidence" : duckhts_somalier_status_string(status);
    printf("%s\t%s\t", case_id, metric);
    print_real(value);
    printf("\t%s\t%s\t%s\t%s\t%zu\n", status_name, usable_sites,
        hom_a, hom_b, evaluations);
}

static void pinned_charr_row(const char *case_id, const char *metric,
                             double value, duckhts_somalier_status_t status,
                             const duckhts_somalier_charr_result_t *result) {
    char usable[32], hom_a[32], hom_b[32];
    snprintf(usable, sizeof(usable), "%" PRIu64, result->usable_sites);
    snprintf(hom_a, sizeof(hom_a), "%" PRIu64, result->usable_hom_a);
    snprintf(hom_b, sizeof(hom_b), "%" PRIu64, result->usable_hom_b);
    pinned_row(case_id, metric, value, status, usable, hom_a, hom_b, 0u);
}

static int campaign_pinned(void) {
    static const uint32_t receiver_a[2] = {160u, 40u};
    static const uint32_t receiver_b[2] = {40u, 160u};
    static const uint8_t receiver_usable[2] = {1u, 1u};
    static const int8_t anchor[2] = {
        DUCKHTS_SOMALIER_HOM_A, DUCKHTS_SOMALIER_HOM_B
    };
    static const double pair_frequency[2] = {0.5, 0.5};
    static const char *likelihood_cases[5] = {
        "pair_likelihood_0", "pair_likelihood_02", "pair_likelihood_04",
        "pair_likelihood_05", "pair_likelihood_1"
    };
    static const double likelihood_alpha[5] = {0.0, 0.2, 0.4, 0.5, 1.0};
    duckhts_somalier_counts_t deep = {4000u, 2000u, 0u, 1u};
    duckhts_somalier_counts_t control = {199u, 1u, 0u, 1u};
    duckhts_somalier_contamination_settings_t charr_settings;
    duckhts_somalier_matched_settings_t matched_settings;
    duckhts_somalier_matched_view_t view = {
        receiver_a, receiver_b, receiver_usable, anchor, pair_frequency, 2u
    };
    duckhts_somalier_charr_result_t charr;
    duckhts_somalier_matched_result_t fitted;
    duckhts_somalier_status_t status;
    uint64_t threshold;
    double deep_frequency = 0.5;
    double control_frequency = 0.25;
    double likelihood;
    size_t i;
    puts("case_id\tmetric\tvalue\tstatus\tusable_sites\tusable_hom_a\t"
         "usable_hom_b\tevaluations");
    for (i = 0u; i < 3u; i++) {
        static const uint64_t depths[3] = {100u, 200u, 6000u};
        static const char *cases[3] = {
            "charr_depth_100", "charr_depth_200", "charr_depth_6000"
        };
        status = duckhts_somalier_binomial_max_minor(depths[i], 0.12, 0.002,
            DUCKHTS_SOMALIER_MAX_BINOMIAL_DEPTH, &threshold);
        if (status != DUCKHTS_SOMALIER_OK) return 1;
        pinned_row(cases[i], "charr_threshold", (double)threshold,
            status, "NA", "NA", "NA", 0u);
    }
    duckhts_somalier_charr_settings_default(&charr_settings);
    status = duckhts_somalier_charr(&deep, &deep_frequency, 1u,
        &charr_settings, &charr);
    pinned_charr_row("charr_deep_usable", "charr_usable",
        (double)charr.usable_sites, status, &charr);
    pinned_charr_row("charr_deep_estimate", "charr_estimate",
        charr.estimate, status, &charr);
    status = duckhts_somalier_charr(&control, &control_frequency, 1u,
        &charr_settings, &charr);
    if (status != DUCKHTS_SOMALIER_OK) return 1;
    pinned_charr_row("charr_control_estimate", "charr_estimate",
        charr.estimate, status, &charr);
    duckhts_somalier_matched_settings_default(&matched_settings);
    status = duckhts_somalier_matched_anchor(&view, &matched_settings, &fitted);
    if (status != DUCKHTS_SOMALIER_OK) return 1;
    pinned_row("pair_two_sites_usable", "pair_usable", (double)fitted.usable_sites,
        status, "2", "NA", "NA", fitted.evaluations);
    pinned_row("pair_search_alpha", "pair_alpha", fitted.alpha,
        status, "2", "NA", "NA", fitted.evaluations);
    status = duckhts_somalier_matched_log_likelihood(&view,
        0.4409830056250525, &matched_settings, &likelihood);
    if (status != DUCKHTS_SOMALIER_OK) return 1;
    pinned_row("pair_search_likelihood", "pair_log_likelihood", likelihood,
        status, "2", "NA", "NA", 0u);
    for (i = 0u; i < 5u; i++) {
        status = duckhts_somalier_matched_log_likelihood(&view,
            likelihood_alpha[i], &matched_settings, &likelihood);
        if (status != DUCKHTS_SOMALIER_OK) return 1;
        pinned_row(likelihood_cases[i], "pair_log_likelihood", likelihood,
            status, "2", "NA", "NA", 0u);
    }
    return 0;
}

static int run_campaign(const char *mode) {
    static const char binomial_header[] =
        "case_id\tdepth\tminor_rate\ttail_alpha\tmax_depth";
    static const char eligibility_header[] =
        "case_id\tallele_a\tallele_b\tother\tavailable\tpopulation_b_af\t"
        "min_depth\tmax_depth\thom_minor_rate\thom_tail_alpha";
    static const char matched_header[] =
        "case_id\toperation\tsite_count\treceiver_a_1\treceiver_b_1\t"
        "receiver_usable_1\tanchor_gt_1\tpopulation_b_af_1\treceiver_a_2\t"
        "receiver_b_2\treceiver_usable_2\tanchor_gt_2\tpopulation_b_af_2\t"
        "alpha\tmax_depth\terror_rate\tmin_probability\tmin_prior_frequency\t"
        "alpha_min\talpha_max\tgrid_step\trefine_tolerance\tmax_evaluations";
    const char *header;
    size_t expected_fields;
    char line[2048];
    size_t row = 0u;
    if (strcmp(mode, "pinned") == 0) return campaign_pinned();
    if (strcmp(mode, "binomial") == 0) {
        header = binomial_header;
        expected_fields = 5u;
        puts("case_id\tstatus\tmax_minor");
    } else if (strcmp(mode, "eligibility") == 0) {
        header = eligibility_header;
        expected_fields = 10u;
        puts("case_id\tclassification_status\tgenotype\treceiver_status\t"
             "receiver_usable\tcharr_status\tcharr_estimate\tcharr_usable_sites\t"
             "charr_hom_a\tcharr_hom_b");
    } else if (strcmp(mode, "matched") == 0) {
        header = matched_header;
        expected_fields = 23u;
        puts("case_id\toperation\tstatus\talpha\tlog_likelihood\t"
             "usable_sites\tevaluations");
    } else {
        fprintf(stderr, "unknown Somalier campaign mode: %s\n", mode);
        return 1;
    }
    if (fgets(line, sizeof(line), stdin) == NULL) {
        fputs("Somalier campaign input has no header\n", stderr);
        return 1;
    }
    line[strcspn(line, "\r\n")] = '\0';
    if (strcmp(line, header) != 0) {
        fprintf(stderr, "Somalier campaign header mismatch for %s\n", mode);
        return 1;
    }
    while (fgets(line, sizeof(line), stdin) != NULL) {
        char *fields[23];
        size_t count;
        int valid;
        row++;
        if (strchr(line, '\n') == NULL && strlen(line) == sizeof(line) - 1u) {
            fprintf(stderr, "Somalier campaign row %zu exceeds input capacity\n", row);
            return 1;
        }
        line[strcspn(line, "\r\n")] = '\0';
        count = split_tsv(line, fields, 23u);
        if (count != expected_fields || fields[0][0] == '\0') {
            fprintf(stderr, "Somalier campaign row %zu has wrong field count\n", row);
            return 1;
        }
        valid = expected_fields == 5u ? campaign_binomial(fields) :
            expected_fields == 10u ? campaign_eligibility(fields) :
            campaign_matched(fields);
        if (!valid) {
            fprintf(stderr, "Somalier campaign row %zu has invalid TSV values\n", row);
            return 1;
        }
    }
    if (row == 0u || ferror(stdin)) {
        fputs("Somalier campaign input ended before a valid row\n", stderr);
        return 1;
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc == 3 && strcmp(argv[1], "--campaign") == 0) {
        return run_campaign(argv[2]);
    }
    if (argc != 1) {
        fputs("usage: somalier_native_test [--campaign pinned|binomial|eligibility|matched]\n",
              stderr);
        return 1;
    }
    test_exhaustive_three_site_pairs();
    test_classification_contracts();
    test_fixed_concordance_witness();
    test_content_identity_and_pair_verifier();
    test_count_digest_receipt();
    test_mask_shapes();
    test_binomial_and_charr();
    test_random_binomial_differential();
    test_charr_orientation();
    test_charr_reduction_order();
    test_matched_anchor();
    test_matched_orientation();
    test_random_pair_differential();
    if (failures != 0u) {
        fprintf(stderr, "somalier native test: %u failures\n", failures);
        return 1;
    }
    printf("somalier native test: OK (4096 exhaustive pairs, 2000 random pairs, "
           "1000 binomial differentials)\n");
    return 0;
}
