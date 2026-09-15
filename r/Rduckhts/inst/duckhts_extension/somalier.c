#include "duckhts_somalier.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#define SOMALIER_MIN_CONTAMINANT_AF 1e-5
#define SOMALIER_BETA_MAX_ITERATIONS 10000u
#define SOMALIER_BETA_EPSILON (8.0 * DBL_EPSILON)

#if defined(__FAST_MATH__)
#error "Somalier numerical kernels require strict floating-point semantics"
#endif

_Static_assert(CHAR_BIT == 8, "Somalier masks require eight-bit bytes");
_Static_assert(sizeof(uint64_t) * CHAR_BIT == 64, "Somalier masks require uint64_t words");
_Static_assert(sizeof(size_t) <= sizeof(uint64_t), "Somalier counters must represent size_t");
_Static_assert(sizeof(size_t) <= sizeof(uintptr_t), "Mask spans must fit uintptr_t");
_Static_assert(sizeof(double) == sizeof(uint64_t), "Somalier content digest requires binary64");
_Static_assert(FLT_RADIX == 2, "Somalier numerical kernels require binary floating point");
_Static_assert(DBL_MANT_DIG == 53 && DBL_MAX_EXP == 1024,
               "Somalier content digest requires IEEE binary64 semantics");

static int relatedness_settings_valid(const duckhts_somalier_relatedness_settings_t *settings) {
    return settings != NULL && settings->min_depth > 0u && settings->max_sites > 0u &&
        isfinite(settings->min_het_balance) && isfinite(settings->hom_balance_cutoff) &&
        settings->hom_balance_cutoff >= 0.0 &&
        settings->hom_balance_cutoff <= settings->min_het_balance &&
        settings->min_het_balance <= 0.5;
}

static int contamination_settings_valid(
    const duckhts_somalier_contamination_settings_t *settings) {
    return settings != NULL && settings->min_depth > 0u &&
        settings->max_depth >= settings->min_depth && settings->max_sites > 0u &&
        isfinite(settings->hom_minor_rate) && settings->hom_minor_rate > 0.0 &&
        settings->hom_minor_rate < 0.5 && isfinite(settings->hom_tail_alpha) &&
        settings->hom_tail_alpha > 0.0 && settings->hom_tail_alpha < 1.0;
}

static int matched_settings_valid(const duckhts_somalier_matched_settings_t *settings) {
    return settings != NULL && settings->max_depth > 0u && settings->max_sites > 0u &&
        settings->max_evaluations > 0u && isfinite(settings->error_rate) &&
        settings->error_rate >= 0.0 && settings->error_rate < 0.5 &&
        isfinite(settings->min_probability) && settings->min_probability > 0.0 &&
        settings->min_probability < 0.5 && isfinite(settings->min_prior_frequency) &&
        settings->min_prior_frequency > 0.0 && settings->min_prior_frequency < 0.5 &&
        isfinite(settings->alpha_min) && isfinite(settings->alpha_max) &&
        settings->alpha_min >= 0.0 && settings->alpha_max <= 1.0 &&
        settings->alpha_min < settings->alpha_max && isfinite(settings->grid_step) &&
        settings->grid_step > 0.0 && settings->grid_step <= 1.0 &&
        isfinite(settings->refine_tolerance) && settings->refine_tolerance > 0.0;
}

static uint64_t counts_depth(const duckhts_somalier_counts_t *counts) {
    return (uint64_t)counts->allele_a + (uint64_t)counts->allele_b;
}

static int relatedness_other_is_high(const duckhts_somalier_counts_t *counts) {
    uint64_t total = counts_depth(counts) + (uint64_t)counts->other;
    return (uint64_t)counts->other * 10u > total;
}

static int contamination_other_is_high(const duckhts_somalier_counts_t *counts) {
    uint64_t total = counts_depth(counts) + (uint64_t)counts->other;
    return (uint64_t)counts->other * 25u > total;
}

static uint64_t popcount64(uint64_t value) {
#if defined(__GNUC__) || defined(__clang__)
    return (uint64_t)__builtin_popcountll((unsigned long long)value);
#else
    uint64_t count = 0u;
    while (value != 0u) {
        value &= value - 1u;
        count++;
    }
    return count;
#endif
}

static double clamp_unit(double value) {
    if (value < 0.0) return 0.0;
    if (value > 1.0) return 1.0;
    return value;
}

static uint64_t count_digest_mix(uint64_t value) {
    value ^= value >> 30u;
    value *= UINT64_C(0xbf58476d1ce4e5b9);
    value ^= value >> 27u;
    value *= UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31u);
}

static uint64_t count_digest_tuple(uint64_t seed, uint64_t site_index,
                                   const duckhts_somalier_counts_t *counts) {
    uint64_t value = count_digest_mix(seed ^ site_index);
    value = count_digest_mix(value ^ (uint64_t)counts->available);
    value = count_digest_mix(value ^ (uint64_t)counts->allele_a);
    value = count_digest_mix(value ^ (uint64_t)counts->allele_b);
    return count_digest_mix(value ^ (uint64_t)counts->other);
}

duckhts_somalier_status_t duckhts_somalier_count_digest_observe(
    uint64_t digest[2], uint64_t site_index,
    const duckhts_somalier_counts_t *counts) {
    if (digest == NULL || counts == NULL || counts->available > 1u) {
        return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    }
    digest[0] += count_digest_tuple(UINT64_C(0x243f6a8885a308d3),
                                    site_index, counts);
    digest[1] += count_digest_tuple(UINT64_C(0x13198a2e03707344),
                                    site_index, counts);
    return DUCKHTS_SOMALIER_OK;
}

duckhts_somalier_status_t duckhts_somalier_count_digest_combine(
    uint64_t target[2], const uint64_t source[2]) {
    if (target == NULL || source == NULL) {
        return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    }
    target[0] += source[0];
    target[1] += source[1];
    return DUCKHTS_SOMALIER_OK;
}

void duckhts_somalier_relatedness_settings_default(
    duckhts_somalier_relatedness_settings_t *settings) {
    if (settings == NULL) return;
    settings->min_depth = 7u;
    settings->max_sites = 100000000u;
    settings->min_het_balance = 0.3;
    settings->hom_balance_cutoff = 0.01;
}

duckhts_somalier_status_t duckhts_somalier_classify_relatedness(
    const duckhts_somalier_counts_t *counts,
    const duckhts_somalier_relatedness_settings_t *settings,
    duckhts_somalier_genotype_t *genotype) {
    uint64_t depth;
    double balance;

    if (counts == NULL || genotype == NULL || counts->available > 1u ||
        !relatedness_settings_valid(settings)) {
        return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    }
    *genotype = DUCKHTS_SOMALIER_UNKNOWN;
    if (!counts->available) return DUCKHTS_SOMALIER_OK;
    depth = counts_depth(counts);
    if (depth < settings->min_depth || relatedness_other_is_high(counts)) {
        return DUCKHTS_SOMALIER_OK;
    }
    balance = (double)counts->allele_b / (double)depth;
    if (balance < settings->hom_balance_cutoff) {
        *genotype = DUCKHTS_SOMALIER_HOM_A;
    } else if (balance > 1.0 - settings->hom_balance_cutoff) {
        *genotype = DUCKHTS_SOMALIER_HOM_B;
    } else if (balance >= settings->min_het_balance &&
               balance <= 1.0 - settings->min_het_balance) {
        *genotype = DUCKHTS_SOMALIER_HET;
    }
    return DUCKHTS_SOMALIER_OK;
}

duckhts_somalier_status_t duckhts_somalier_classify_relatedness_evidence(
    const duckhts_somalier_counts_t *counts,
    const duckhts_somalier_relatedness_settings_t *settings,
    duckhts_somalier_classification_t *classification) {
    duckhts_somalier_status_t status;
    uint64_t depth;
    uint64_t total;
    double balance;

    if (classification == NULL) return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    memset(classification, 0, sizeof(*classification));
    classification->genotype = DUCKHTS_SOMALIER_UNKNOWN;
    status = duckhts_somalier_classify_relatedness(
        counts, settings, &classification->genotype);
    if (status != DUCKHTS_SOMALIER_OK) return status;
    if (!counts->available) {
        classification->unavailable = 1u;
        return DUCKHTS_SOMALIER_OK;
    }
    depth = counts_depth(counts);
    total = depth + counts->other;
    if (depth < settings->min_depth || (uint64_t)counts->other * 10u > total) {
        return DUCKHTS_SOMALIER_OK;
    }
    balance = (double)counts->allele_b / (double)depth;
    classification->middling_balance =
        balance > 0.02 && balance < 0.98 &&
        (balance < 0.10 || balance > 0.90);
    return DUCKHTS_SOMALIER_OK;
}

duckhts_somalier_status_t duckhts_somalier_mask_word_count(
    size_t site_count, size_t *word_count) {
    if (word_count == NULL) return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    *word_count = site_count / 64u + (site_count % 64u != 0u);
    return DUCKHTS_SOMALIER_OK;
}

static duckhts_somalier_status_t mask_shape(const duckhts_somalier_masks_t *masks) {
    size_t required;
    size_t byte_count;
    uintptr_t start[3];
    uintptr_t end[3];
    unsigned i;
    unsigned j;
    if (masks == NULL) return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    (void)duckhts_somalier_mask_word_count(masks->site_count, &required);
    if (masks->word_count != required) return DUCKHTS_SOMALIER_CORRUPT_MASK;
    if (required == 0u) {
        return masks->hom_a == NULL && masks->het == NULL && masks->hom_b == NULL
            ? DUCKHTS_SOMALIER_OK
            : DUCKHTS_SOMALIER_CORRUPT_MASK;
    }
    if (masks->hom_a == NULL || masks->het == NULL || masks->hom_b == NULL ||
        masks->word_count > SIZE_MAX / sizeof(uint64_t)) {
        return DUCKHTS_SOMALIER_CORRUPT_MASK;
    }
    byte_count = masks->word_count * sizeof(uint64_t);
    start[0] = (uintptr_t)(void *)masks->hom_a;
    start[1] = (uintptr_t)(void *)masks->het;
    start[2] = (uintptr_t)(void *)masks->hom_b;
    for (i = 0u; i < 3u; i++) {
        if (start[i] % _Alignof(uint64_t) != 0u ||
            start[i] > UINTPTR_MAX - (uintptr_t)byte_count) {
            return DUCKHTS_SOMALIER_CORRUPT_MASK;
        }
        end[i] = start[i] + (uintptr_t)byte_count;
    }
    for (i = 0u; i < 3u; i++) {
        for (j = i + 1u; j < 3u; j++) {
            if (start[i] < end[j] && start[j] < end[i]) {
                return DUCKHTS_SOMALIER_CORRUPT_MASK;
            }
        }
    }
    return DUCKHTS_SOMALIER_OK;
}

static int sketch_identity_valid(const duckhts_somalier_sketch_identity_t *identity) {
    return identity != NULL && identity->min_depth > 0u &&
        isfinite(identity->min_het_balance) &&
        isfinite(identity->hom_balance_cutoff) &&
        identity->hom_balance_cutoff >= 0.0 &&
        identity->hom_balance_cutoff <= identity->min_het_balance &&
        identity->min_het_balance <= 0.5;
}

static int sketch_identity_equal(const duckhts_somalier_sketch_identity_t *a,
                                 const duckhts_somalier_sketch_identity_t *b) {
    return memcmp(a->panel_sha256, b->panel_sha256,
                  sizeof(a->panel_sha256)) == 0 &&
        a->min_depth == b->min_depth &&
        a->min_het_balance == b->min_het_balance &&
        a->hom_balance_cutoff == b->hom_balance_cutoff;
}

static duckhts_somalier_status_t validate_mask_content(
    const duckhts_somalier_masks_t *masks) {
    duckhts_somalier_status_t status = mask_shape(masks);
    size_t i;
    if (status != DUCKHTS_SOMALIER_OK) return status;
    if (masks->middling_balance_count > (uint64_t)masks->site_count ||
        masks->unavailable_count > (uint64_t)masks->site_count) {
        return DUCKHTS_SOMALIER_CORRUPT_MASK;
    }
    if (!sketch_identity_valid(&masks->identity)) {
        return DUCKHTS_SOMALIER_CORRUPT_MASK;
    }
    for (i = 0u; i < masks->word_count; i++) {
        uint64_t overlap = (masks->hom_a[i] & masks->het[i]) |
            (masks->hom_a[i] & masks->hom_b[i]) | (masks->het[i] & masks->hom_b[i]);
        if (overlap != 0u) return DUCKHTS_SOMALIER_CORRUPT_MASK;
    }
    if (masks->word_count > 0u && masks->site_count % 64u != 0u) {
        unsigned used = (unsigned)(masks->site_count % 64u);
        uint64_t tail_mask = UINT64_MAX << used;
        size_t last = masks->word_count - 1u;
        if (((masks->hom_a[last] | masks->het[last] | masks->hom_b[last]) & tail_mask) != 0u) {
            return DUCKHTS_SOMALIER_CORRUPT_MASK;
        }
    }
    return DUCKHTS_SOMALIER_OK;
}

typedef struct somalier_digest_state {
    uint64_t word[2];
} somalier_digest_state_t;

static void digest_byte(somalier_digest_state_t *state, uint8_t value) {
    state->word[0] = (state->word[0] ^ value) * UINT64_C(0x100000001b3);
    state->word[1] = (state->word[1] ^ value) * UINT64_C(0x100000001b7);
}

static void digest_u64(somalier_digest_state_t *state, uint64_t value) {
    for (unsigned shift = 0u; shift < 64u; shift += 8u) {
        digest_byte(state, (uint8_t)(value >> shift));
    }
}

static void digest_double(somalier_digest_state_t *state, double value) {
    uint64_t bits;
    memcpy(&bits, &value, sizeof(bits));
    digest_u64(state, bits);
}

static void mask_content_digest(const duckhts_somalier_masks_t *masks,
                                uint64_t digest[2]) {
    static const uint8_t domain[] = "duckhts-somalier-sketch-content-v2";
    somalier_digest_state_t state = {
        {UINT64_C(0xcbf29ce484222325), UINT64_C(0x84222325cbf29ce4)}
    };
    for (size_t i = 0u; i < sizeof(domain); i++) digest_byte(&state, domain[i]);
    digest_u64(&state, (uint64_t)masks->site_count);
    digest_u64(&state, (uint64_t)masks->word_count);
    digest_u64(&state, masks->middling_balance_count);
    digest_u64(&state, masks->unavailable_count);
    digest_u64(&state, masks->count_digest[0]);
    digest_u64(&state, masks->count_digest[1]);
    for (size_t i = 0u; i < sizeof(masks->identity.panel_sha256); i++) {
        digest_byte(&state, masks->identity.panel_sha256[i]);
    }
    digest_u64(&state, masks->identity.min_depth);
    digest_double(&state, masks->identity.min_het_balance);
    digest_double(&state, masks->identity.hom_balance_cutoff);
    for (size_t i = 0u; i < masks->word_count; i++) {
        digest_u64(&state, masks->hom_a[i]);
        digest_u64(&state, masks->het[i]);
        digest_u64(&state, masks->hom_b[i]);
    }
    digest[0] = state.word[0];
    digest[1] = state.word[1];
}

duckhts_somalier_status_t duckhts_somalier_seal_masks(
    duckhts_somalier_masks_t *masks) {
    duckhts_somalier_status_t status = validate_mask_content(masks);
    if (status != DUCKHTS_SOMALIER_OK) return status;
    mask_content_digest(masks, masks->content_digest);
    return DUCKHTS_SOMALIER_OK;
}

duckhts_somalier_status_t duckhts_somalier_validate_masks(
    const duckhts_somalier_masks_t *masks) {
    uint64_t expected[2];
    duckhts_somalier_status_t status = validate_mask_content(masks);
    if (status != DUCKHTS_SOMALIER_OK) return status;
    mask_content_digest(masks, expected);
    return expected[0] == masks->content_digest[0] &&
        expected[1] == masks->content_digest[1]
        ? DUCKHTS_SOMALIER_OK : DUCKHTS_SOMALIER_CORRUPT_MASK;
}

duckhts_somalier_status_t duckhts_somalier_prepare_masks(
    const duckhts_somalier_counts_t *counts,
    size_t site_count,
    const uint8_t panel_sha256[32],
    const duckhts_somalier_relatedness_settings_t *settings,
    duckhts_somalier_masks_t *masks) {
    duckhts_somalier_status_t status;
    size_t i;

    if (!relatedness_settings_valid(settings) || masks == NULL || panel_sha256 == NULL ||
        (site_count > 0u && counts == NULL)) {
        return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    }
    if (site_count > settings->max_sites) return DUCKHTS_SOMALIER_LIMIT_EXCEEDED;
    if (masks->site_count != site_count) return DUCKHTS_SOMALIER_CORRUPT_MASK;
    status = mask_shape(masks);
    if (status != DUCKHTS_SOMALIER_OK) return status;
    for (i = 0u; i < site_count; i++) {
        if (counts[i].available > 1u) return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    }
    masks->middling_balance_count = 0u;
    masks->unavailable_count = 0u;
    masks->count_digest[0] = 0u;
    masks->count_digest[1] = 0u;
    memcpy(masks->identity.panel_sha256, panel_sha256,
           sizeof(masks->identity.panel_sha256));
    masks->identity.min_depth = settings->min_depth;
    masks->identity.min_het_balance = settings->min_het_balance == 0.0
        ? 0.0 : settings->min_het_balance;
    masks->identity.hom_balance_cutoff = settings->hom_balance_cutoff == 0.0
        ? 0.0 : settings->hom_balance_cutoff;
    if (masks->word_count > 0u) {
        memset(masks->hom_a, 0, masks->word_count * sizeof(*masks->hom_a));
        memset(masks->het, 0, masks->word_count * sizeof(*masks->het));
        memset(masks->hom_b, 0, masks->word_count * sizeof(*masks->hom_b));
    }
    for (i = 0u; i < site_count; i++) {
        duckhts_somalier_classification_t classification;
        uint64_t bit = UINT64_C(1) << (i % 64u);
        status = duckhts_somalier_classify_relatedness_evidence(
            &counts[i], settings, &classification);
        if (status != DUCKHTS_SOMALIER_OK) return status;
        status = duckhts_somalier_count_digest_observe(
            masks->count_digest, (uint64_t)i, &counts[i]);
        if (status != DUCKHTS_SOMALIER_OK) return status;
        masks->unavailable_count += classification.unavailable;
        masks->middling_balance_count += classification.middling_balance;
        if (classification.genotype == DUCKHTS_SOMALIER_HOM_A) {
            masks->hom_a[i / 64u] |= bit;
        } else if (classification.genotype == DUCKHTS_SOMALIER_HET) {
            masks->het[i / 64u] |= bit;
        } else if (classification.genotype == DUCKHTS_SOMALIER_HOM_B) {
            masks->hom_b[i / 64u] |= bit;
        }
    }
    return duckhts_somalier_seal_masks(masks);
}

duckhts_somalier_status_t duckhts_somalier_pair_stats(
    const duckhts_somalier_masks_t *a,
    const duckhts_somalier_masks_t *b,
    duckhts_somalier_pair_stats_t *result) {
    duckhts_somalier_status_t status;
    size_t i;

    if (result == NULL || a == NULL || b == NULL) {
        return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    }
    memset(result, 0, sizeof(*result));
    if (a->site_count != b->site_count || a->word_count != b->word_count) {
        return DUCKHTS_SOMALIER_CORRUPT_MASK;
    }
    if ((uint64_t)a->site_count > UINT64_MAX / 2u) {
        return DUCKHTS_SOMALIER_LIMIT_EXCEEDED;
    }
    status = duckhts_somalier_validate_masks(a);
    if (status != DUCKHTS_SOMALIER_OK) return status;
    status = duckhts_somalier_validate_masks(b);
    if (status != DUCKHTS_SOMALIER_OK) return status;
    if (!sketch_identity_equal(&a->identity, &b->identity)) {
        return DUCKHTS_SOMALIER_IDENTITY_MISMATCH;
    }

    for (i = 0u; i < a->word_count; i++) {
        uint64_t known_a = a->hom_a[i] | a->het[i] | a->hom_b[i];
        uint64_t known_b = b->hom_a[i] | b->het[i] | b->hom_b[i];
        uint64_t jointly_called = known_a & known_b;
        uint64_t shared_hets = a->het[i] & b->het[i];
        uint64_t shared_hom_b = a->hom_b[i] & b->hom_b[i];
        uint64_t hom_a_callable = popcount64((a->hom_a[i] | a->hom_b[i]) & known_b);
        uint64_t hom_b_callable = popcount64((b->hom_a[i] | b->hom_b[i]) & known_a);
        uint64_t hom_matches = popcount64((a->hom_a[i] & b->hom_a[i]) | shared_hom_b);
        result->ibs0 += popcount64((a->hom_a[i] & b->hom_b[i]) |
                                  (a->hom_b[i] & b->hom_a[i]));
        result->ibs2 += popcount64((a->hom_a[i] & b->hom_a[i]) |
                                  shared_hets | shared_hom_b);
        result->jointly_called += popcount64(jointly_called);
        result->shared_hets += popcount64(shared_hets);
        result->shared_hom_b += popcount64(shared_hom_b);
        result->het_count_a += popcount64(a->het[i]);
        result->het_count_b += popcount64(b->het[i]);
        result->hom_b_count_a += popcount64(a->hom_b[i]);
        result->hom_b_count_b += popcount64(b->hom_b[i]);
        result->callable_hom_count_a += hom_a_callable;
        result->callable_hom_count_b += hom_b_callable;
        result->matching_hom_count += hom_matches;
        result->het_ab += popcount64(a->het[i] & jointly_called) +
            popcount64(b->het[i] & jointly_called);
    }
    result->relatedness = 2.0 * ((double)result->shared_hets -
        2.0 * (double)result->ibs0) /
        (result->het_ab == 0u ? 1.0 : (double)result->het_ab);
    {
        uint64_t hom_denominator_count = result->callable_hom_count_a <
            result->callable_hom_count_b ? result->callable_hom_count_a :
            result->callable_hom_count_b;
        double hom_b_denominator = (double)(result->hom_b_count_a < result->hom_b_count_b
            ? result->hom_b_count_a : result->hom_b_count_b);
        double middling_denominator_a = (double)a->site_count +
            (double)a->middling_balance_count;
        double middling_denominator_b = (double)b->site_count +
            (double)b->middling_balance_count;
        double p_middling_a = (double)a->middling_balance_count /
            (middling_denominator_a < 1.0 ? 1.0 : middling_denominator_a);
        double p_middling_b = (double)b->middling_balance_count /
            (middling_denominator_b < 1.0 ? 1.0 : middling_denominator_b);
        double low_hom_b_penalty;
        double base;

        result->inferred_hom_concordance = (double)result->matching_hom_count /
            (hom_denominator_count == 0u ? 1.0 : (double)hom_denominator_count);
        result->raw_hom_b_concordance =
            ((double)result->shared_hom_b - 2.0 * (double)result->ibs0) /
            (hom_b_denominator < 1.0 ? 1.0 : hom_b_denominator);
        result->p_middling_a = p_middling_a;
        result->p_middling_b = p_middling_b;
        low_hom_b_penalty = fmax(0.0, 0.70 - clamp_unit(result->raw_hom_b_concordance) +
            p_middling_a + p_middling_b);
        base = clamp_unit(result->inferred_hom_concordance - low_hom_b_penalty);
        if (base > 0.4) {
            base = 0.4 + 0.6 * cbrt((base - 0.4) / 0.6);
        }
        result->adjusted_concordance = base;
    }
    return DUCKHTS_SOMALIER_OK;
}

static int same_double_bits(double left, double right) {
    return memcmp(&left, &right, sizeof(left)) == 0;
}

static int pair_stats_equal(const duckhts_somalier_pair_stats_t *left,
                            const duckhts_somalier_pair_stats_t *right) {
    return left->ibs0 == right->ibs0 && left->ibs2 == right->ibs2 &&
        left->jointly_called == right->jointly_called &&
        left->shared_hets == right->shared_hets &&
        left->shared_hom_b == right->shared_hom_b && left->het_ab == right->het_ab &&
        left->het_count_a == right->het_count_a &&
        left->het_count_b == right->het_count_b &&
        left->hom_b_count_a == right->hom_b_count_a &&
        left->hom_b_count_b == right->hom_b_count_b &&
        left->callable_hom_count_a == right->callable_hom_count_a &&
        left->callable_hom_count_b == right->callable_hom_count_b &&
        left->matching_hom_count == right->matching_hom_count &&
        same_double_bits(left->relatedness, right->relatedness) &&
        same_double_bits(left->inferred_hom_concordance,
                         right->inferred_hom_concordance) &&
        same_double_bits(left->raw_hom_b_concordance,
                         right->raw_hom_b_concordance) &&
        same_double_bits(left->p_middling_a, right->p_middling_a) &&
        same_double_bits(left->p_middling_b, right->p_middling_b) &&
        same_double_bits(left->adjusted_concordance,
                         right->adjusted_concordance);
}

duckhts_somalier_status_t duckhts_somalier_verify_pair_result(
    const duckhts_somalier_masks_t *a,
    const duckhts_somalier_masks_t *b,
    const duckhts_somalier_pair_stats_t *reported,
    duckhts_somalier_status_t reported_status) {
    duckhts_somalier_pair_stats_t expected;
    duckhts_somalier_status_t status;
    duckhts_somalier_status_t expected_status;
    if (reported == NULL) return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    status = duckhts_somalier_pair_stats(a, b, &expected);
    if (status != DUCKHTS_SOMALIER_OK) return status;
    expected_status = expected.jointly_called == 0u
        ? DUCKHTS_SOMALIER_NO_EVIDENCE : DUCKHTS_SOMALIER_OK;
    return reported_status == expected_status && pair_stats_equal(reported, &expected)
        ? DUCKHTS_SOMALIER_OK : DUCKHTS_SOMALIER_CORRUPT_RESULT;
}

void duckhts_somalier_charr_settings_default(
    duckhts_somalier_contamination_settings_t *settings) {
    if (settings == NULL) return;
    settings->min_depth = 15u;
    settings->max_depth = DUCKHTS_SOMALIER_MAX_BINOMIAL_DEPTH;
    settings->max_sites = 100000000u;
    settings->hom_minor_rate = 0.12;
    settings->hom_tail_alpha = 0.002;
}

void duckhts_somalier_matched_anchor_settings_default(
    duckhts_somalier_contamination_settings_t *settings) {
    if (settings == NULL) return;
    settings->min_depth = 15u;
    settings->max_depth = DUCKHTS_SOMALIER_MAX_BINOMIAL_DEPTH;
    settings->max_sites = 100000000u;
    settings->hom_minor_rate = 0.05;
    settings->hom_tail_alpha = 0.001;
}

/* Modified Lentz continued fraction for the regularized incomplete beta. */
static int beta_fraction(double a, double b, double x, double *value) {
    double qab = a + b;
    double qap = a + 1.0;
    double qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    double h;
    unsigned m;

    if (fabs(d) < DBL_MIN / SOMALIER_BETA_EPSILON) {
        d = DBL_MIN / SOMALIER_BETA_EPSILON;
    }
    d = 1.0 / d;
    h = d;
    for (m = 1u; m <= SOMALIER_BETA_MAX_ITERATIONS; m++) {
        double m2 = 2.0 * (double)m;
        double aa = (double)m * (b - (double)m) * x /
            ((qam + m2) * (a + m2));
        double delta;
        d = 1.0 + aa * d;
        if (fabs(d) < DBL_MIN / SOMALIER_BETA_EPSILON) {
            d = DBL_MIN / SOMALIER_BETA_EPSILON;
        }
        c = 1.0 + aa / c;
        if (fabs(c) < DBL_MIN / SOMALIER_BETA_EPSILON) {
            c = DBL_MIN / SOMALIER_BETA_EPSILON;
        }
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + (double)m) * (qab + (double)m) * x /
            ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < DBL_MIN / SOMALIER_BETA_EPSILON) {
            d = DBL_MIN / SOMALIER_BETA_EPSILON;
        }
        c = 1.0 + aa / c;
        if (fabs(c) < DBL_MIN / SOMALIER_BETA_EPSILON) {
            c = DBL_MIN / SOMALIER_BETA_EPSILON;
        }
        d = 1.0 / d;
        delta = d * c;
        h *= delta;
        if (!isfinite(h)) return 0;
        if (fabs(delta - 1.0) <= SOMALIER_BETA_EPSILON) {
            *value = h;
            return 1;
        }
    }
    return 0;
}

static int regularized_beta(double x, double a, double b, double *value) {
    double fraction;
    double front;
    if (x <= 0.0) {
        *value = 0.0;
        return 1;
    }
    if (x >= 1.0) {
        *value = 1.0;
        return 1;
    }
    front = exp(lgamma(a + b) - lgamma(a) - lgamma(b) +
                a * log(x) + b * log1p(-x));
    if (!isfinite(front)) return 0;
    if (x < (a + 1.0) / (a + b + 2.0)) {
        if (!beta_fraction(a, b, x, &fraction)) return 0;
        *value = front * fraction / a;
    } else {
        if (!beta_fraction(b, a, 1.0 - x, &fraction)) return 0;
        *value = 1.0 - front * fraction / b;
    }
    if (*value < 0.0 && *value > -32.0 * DBL_EPSILON) *value = 0.0;
    if (*value > 1.0 && *value < 1.0 + 32.0 * DBL_EPSILON) *value = 1.0;
    return isfinite(*value) && *value >= 0.0 && *value <= 1.0;
}

static int binomial_survival(uint64_t depth, uint64_t k, double probability,
                             double *survival) {
    if (k == 0u) {
        *survival = 1.0;
        return 1;
    }
    if (k > depth) {
        *survival = 0.0;
        return 1;
    }
    return regularized_beta(probability, (double)k, (double)(depth - k + 1u), survival);
}

typedef struct binomial_interval {
    long double lower;
    long double upper;
} binomial_interval_t;

static long double binomial_round_down(long double value) {
    volatile long double rounded = value;
    if (rounded <= 0.0L) return 0.0L;
    return nextafterl(rounded, 0.0L);
}

static long double binomial_round_up(long double value) {
    volatile long double rounded = value;
    return nextafterl(rounded, INFINITY);
}

static binomial_interval_t binomial_interval_add(binomial_interval_t left,
                                                  binomial_interval_t right) {
    binomial_interval_t result;
    result.lower = binomial_round_down(left.lower + right.lower);
    result.upper = binomial_round_up(left.upper + right.upper);
    return result;
}

static binomial_interval_t binomial_interval_multiply(
    binomial_interval_t left,
    binomial_interval_t right) {
    binomial_interval_t result;
    result.lower = binomial_round_down(left.lower * right.lower);
    result.upper = binomial_round_up(left.upper * right.upper);
    return result;
}

static int binomial_interval_divide(binomial_interval_t numerator,
                                    binomial_interval_t denominator,
                                    binomial_interval_t *result) {
    if (result == NULL || denominator.lower <= 0.0) return 0;
    result->lower = binomial_round_down(numerator.lower / denominator.upper);
    result->upper = binomial_round_up(numerator.upper / denominator.lower);
    return isfinite(result->lower) && isfinite(result->upper) &&
        result->lower >= 0.0 && result->upper >= result->lower;
}

static int binomial_interval_remainder_upper(
    binomial_interval_t weight,
    binomial_interval_t next_ratio,
    long double *upper) {
    long double numerator;
    long double denominator;
    if (upper == NULL || next_ratio.upper <= 0.0L ||
        next_ratio.upper >= 1.0L) {
        return 0;
    }
    numerator = binomial_round_up(weight.upper * next_ratio.upper);
    denominator = binomial_round_down(1.0L - next_ratio.upper);
    if (denominator <= 0.0L) return 0;
    *upper = binomial_round_up(numerator / denominator);
    return isfinite(*upper) && *upper >= 0.0;
}

static unsigned binomial_trailing_zeroes(uint64_t value) {
    unsigned count = 0u;
    while ((value & UINT64_C(1)) == 0u) {
        value >>= 1u;
        count++;
    }
    return count;
}

/* Decode an exact binary64 probability as numerator / 2^denominator_power,
 * reduced to an odd numerator. */
static int binomial_binary_fraction(double value, uint64_t *numerator,
                                    unsigned *denominator_power) {
    double fraction;
    int exponent;
    unsigned trailing;

    if (numerator == NULL || denominator_power == NULL ||
        !isfinite(value) || value <= 0.0 || value >= 1.0) {
        return 0;
    }
    fraction = frexp(value, &exponent);
    *numerator = (uint64_t)ldexp(fraction, DBL_MANT_DIG);
    *denominator_power = (unsigned)(DBL_MANT_DIG - exponent);
    trailing = binomial_trailing_zeroes(*numerator);
    *numerator >>= trailing;
    *denominator_power -= trailing;
    return 1;
}

/* When the complete binomial denominator fits in binary64's significand,
 * polynomial convolution gives the exact tail using only uint64 arithmetic. */
static int binomial_tail_exact_small(uint64_t depth, uint64_t k,
                                     double probability, double tail_alpha,
                                     int *at_least) {
    uint64_t coefficients[54] = {0u};
    uint64_t numerator;
    uint64_t denominator;
    uint64_t complement;
    uint64_t tail = 0u;
    uint64_t power;
    unsigned denominator_power;

    if (at_least == NULL ||
        !binomial_binary_fraction(probability, &numerator,
                                  &denominator_power) ||
        denominator_power == 0u || depth > 53u ||
        denominator_power > 53u / (unsigned)depth) {
        return 0;
    }
    power = denominator_power * (unsigned)depth;
    denominator = UINT64_C(1) << denominator_power;
    complement = denominator - numerator;
    coefficients[0] = 1u;
    for (uint64_t trial = 0u; trial < depth; trial++) {
        for (uint64_t successes = trial + 1u; successes > 0u; successes--) {
            uint64_t from_failure = successes <= trial
                ? coefficients[successes] * complement : 0u;
            uint64_t from_success = coefficients[successes - 1u] * numerator;
            coefficients[successes] = from_failure + from_success;
        }
        coefficients[0] *= complement;
    }
    for (uint64_t successes = k; successes <= depth; successes++) {
        tail += coefficients[successes];
    }
    *at_least = ldexp((double)tail, -(int)power) >= tail_alpha;
    return 1;
}

/* Enclose the exact tail using outward-rounded floating-point operations. The
 * incomplete beta locates a candidate only; this comparison certifies the
 * final k/k+1 decision or reports that binary64 intervals cannot decide it. */
static int binomial_tail_interval(uint64_t depth, uint64_t k,
                                  double probability,
                                  binomial_interval_t *survival) {
    binomial_interval_t one = {1.0L, 1.0L};
    binomial_interval_t probability_interval = {
        (long double)probability, (long double)probability
    };
    binomial_interval_t complement;
    binomial_interval_t total = {1.0L, 1.0L};
    binomial_interval_t side = {0.0L, 0.0L};
    binomial_interval_t weight = {1.0L, 1.0L};
    binomial_interval_t probability_ratio;
    uint64_t mode;
    volatile long double rounded_complement =
        1.0L - (long double)probability;

    if (survival == NULL) return 0;
    complement.lower = binomial_round_down(rounded_complement);
    complement.upper = binomial_round_up(rounded_complement);
    if (!binomial_interval_divide(complement, probability_interval,
                                  &probability_ratio)) {
        return 0;
    }
    mode = (uint64_t)floor(((double)depth + 1.0) * probability);
    if (mode > depth) mode = depth;

    for (uint64_t i = mode; i > 0u; i--) {
        binomial_interval_t integer_ratio = {
            (long double)i / (long double)(depth - i + 1u),
            (long double)i / (long double)(depth - i + 1u)
        };
        integer_ratio.lower = binomial_round_down(integer_ratio.lower);
        integer_ratio.upper = binomial_round_up(integer_ratio.upper);
        weight = binomial_interval_multiply(weight,
            binomial_interval_multiply(integer_ratio, probability_ratio));
        total = binomial_interval_add(total, weight);
        if (i - 1u < k && k <= mode) {
            side = binomial_interval_add(side, weight);
        }
        if (i > 1u) {
            binomial_interval_t next_integer_ratio = {
                (long double)(i - 1u) / (long double)(depth - i + 2u),
                (long double)(i - 1u) / (long double)(depth - i + 2u)
            };
            binomial_interval_t next_ratio;
            long double remainder;
            next_integer_ratio.lower =
                binomial_round_down(next_integer_ratio.lower);
            next_integer_ratio.upper =
                binomial_round_up(next_integer_ratio.upper);
            next_ratio = binomial_interval_multiply(
                next_integer_ratio, probability_ratio);
            if (binomial_interval_remainder_upper(weight, next_ratio,
                    &remainder) && remainder <= LDBL_EPSILON / 64.0L) {
                total.upper = binomial_round_up(total.upper + remainder);
                if (k <= mode) {
                    side.upper = binomial_round_up(side.upper + remainder);
                }
                break;
            }
        }
    }

    weight = one;
    if (!binomial_interval_divide(probability_interval, complement,
                                  &probability_ratio)) {
        return 0;
    }
    for (uint64_t i = mode; i < depth; i++) {
        binomial_interval_t integer_ratio = {
            (long double)(depth - i) / (long double)(i + 1u),
            (long double)(depth - i) / (long double)(i + 1u)
        };
        integer_ratio.lower = binomial_round_down(integer_ratio.lower);
        integer_ratio.upper = binomial_round_up(integer_ratio.upper);
        weight = binomial_interval_multiply(weight,
            binomial_interval_multiply(integer_ratio, probability_ratio));
        total = binomial_interval_add(total, weight);
        if (i + 1u >= k && k > mode) {
            side = binomial_interval_add(side, weight);
        }
        if (i + 1u < depth) {
            binomial_interval_t next_integer_ratio = {
                (long double)(depth - i - 1u) / (long double)(i + 2u),
                (long double)(depth - i - 1u) / (long double)(i + 2u)
            };
            binomial_interval_t next_ratio;
            long double remainder;
            next_integer_ratio.lower =
                binomial_round_down(next_integer_ratio.lower);
            next_integer_ratio.upper =
                binomial_round_up(next_integer_ratio.upper);
            next_ratio = binomial_interval_multiply(
                next_integer_ratio, probability_ratio);
            if (binomial_interval_remainder_upper(weight, next_ratio,
                    &remainder) && remainder <= LDBL_EPSILON / 64.0L) {
                total.upper = binomial_round_up(total.upper + remainder);
                if (k > mode) {
                    side.upper = binomial_round_up(side.upper + remainder);
                }
                break;
            }
        }
    }
    if (k <= mode) {
        binomial_interval_t excluded;
        if (!binomial_interval_divide(side, total, &excluded)) return 0;
        survival->lower = binomial_round_down(1.0L - excluded.upper);
        survival->upper = binomial_round_up(1.0L - excluded.lower);
    } else if (!binomial_interval_divide(side, total, survival)) {
        return 0;
    }
    if (survival->lower < 0.0L) survival->lower = 0.0L;
    if (survival->upper > 1.0L) survival->upper = 1.0L;
    return survival->lower <= survival->upper;
}

static int binomial_tail_compare(uint64_t depth, uint64_t k,
                                 double probability, double tail_alpha,
                                 int *at_least) {
    binomial_interval_t survival;
    if (at_least == NULL) return 0;
    if (k == 0u) {
        *at_least = 1;
        return 1;
    }
    if (k > depth || probability == 0.0) {
        *at_least = 0;
        return 1;
    }
    if (probability == 1.0) {
        *at_least = 1;
        return 1;
    }
    if (depth == 1u) {
        *at_least = probability >= tail_alpha;
        return 1;
    }
    if (binomial_tail_exact_small(depth, k, probability, tail_alpha,
                                  at_least)) {
        return 1;
    }
    if (!binomial_tail_interval(depth, k, probability, &survival)) return 0;
    if (survival.lower >= (long double)tail_alpha) {
        *at_least = 1;
        return 1;
    }
    if (survival.upper < (long double)tail_alpha) {
        *at_least = 0;
        return 1;
    }
    return 0;
}

duckhts_somalier_status_t duckhts_somalier_binomial_max_minor(
    uint64_t depth,
    double minor_rate,
    double tail_alpha,
    uint64_t max_depth,
    uint64_t *max_minor) {
    uint64_t lo = 0u;
    uint64_t hi;

    if (max_minor == NULL || max_depth == 0u || depth > max_depth ||
        depth > DUCKHTS_SOMALIER_MAX_BINOMIAL_DEPTH ||
        !isfinite(minor_rate) || minor_rate < 0.0 || minor_rate > 1.0 ||
        !isfinite(tail_alpha) || tail_alpha <= 0.0 || tail_alpha >= 1.0) {
        return (depth > max_depth || depth > DUCKHTS_SOMALIER_MAX_BINOMIAL_DEPTH) &&
            max_minor != NULL
            ? DUCKHTS_SOMALIER_LIMIT_EXCEEDED
            : DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    }
    if (depth == UINT64_MAX) return DUCKHTS_SOMALIER_LIMIT_EXCEEDED;
    hi = depth + 1u;
    while (hi - lo > 1u) {
        uint64_t mid = lo + (hi - lo) / 2u;
        double survival;
        if (!binomial_survival(depth, mid, minor_rate, &survival)) {
            return DUCKHTS_SOMALIER_NUMERIC_FAILURE;
        }
        if (survival >= tail_alpha) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    while (lo > 0u) {
        int at_least;
        if (!binomial_tail_compare(depth, lo, minor_rate, tail_alpha,
                                   &at_least)) {
            return DUCKHTS_SOMALIER_NUMERIC_FAILURE;
        }
        if (at_least) break;
        lo--;
    }
    while (lo < depth) {
        int at_least;
        if (!binomial_tail_compare(depth, lo + 1u, minor_rate, tail_alpha,
                                   &at_least)) {
            return DUCKHTS_SOMALIER_NUMERIC_FAILURE;
        }
        if (!at_least) break;
        lo++;
    }
    *max_minor = lo;
    return DUCKHTS_SOMALIER_OK;
}

duckhts_somalier_status_t duckhts_somalier_classify_contamination(
    const duckhts_somalier_counts_t *counts,
    const duckhts_somalier_contamination_settings_t *settings,
    duckhts_somalier_genotype_t *genotype) {
    duckhts_somalier_status_t status;
    uint64_t depth;
    uint64_t threshold;
    uint64_t minor;

    if (counts == NULL || genotype == NULL || counts->available > 1u ||
        !contamination_settings_valid(settings)) {
        return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    }
    *genotype = DUCKHTS_SOMALIER_UNKNOWN;
    if (!counts->available) return DUCKHTS_SOMALIER_OK;
    depth = counts_depth(counts);
    if (depth > settings->max_depth) return DUCKHTS_SOMALIER_LIMIT_EXCEEDED;
    if (depth < settings->min_depth || contamination_other_is_high(counts)) {
        return DUCKHTS_SOMALIER_OK;
    }
    status = duckhts_somalier_binomial_max_minor(depth, settings->hom_minor_rate,
        settings->hom_tail_alpha, settings->max_depth, &threshold);
    if (status != DUCKHTS_SOMALIER_OK) return status;
    minor = counts->allele_a < counts->allele_b ? counts->allele_a : counts->allele_b;
    if (minor > threshold || counts->allele_a == counts->allele_b) {
        return DUCKHTS_SOMALIER_OK;
    }
    *genotype = counts->allele_a > counts->allele_b
        ? DUCKHTS_SOMALIER_HOM_A
        : DUCKHTS_SOMALIER_HOM_B;
    return DUCKHTS_SOMALIER_OK;
}

duckhts_somalier_status_t duckhts_somalier_contamination_usable(
    const duckhts_somalier_counts_t *counts,
    const duckhts_somalier_contamination_settings_t *settings,
    uint8_t *usable) {
    uint64_t depth;
    if (counts == NULL || usable == NULL || counts->available > 1u ||
        !contamination_settings_valid(settings)) {
        return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    }
    *usable = 0u;
    if (!counts->available) return DUCKHTS_SOMALIER_OK;
    depth = counts_depth(counts);
    if (depth > settings->max_depth) return DUCKHTS_SOMALIER_LIMIT_EXCEEDED;
    *usable = depth >= settings->min_depth && !contamination_other_is_high(counts);
    return DUCKHTS_SOMALIER_OK;
}

duckhts_somalier_status_t duckhts_somalier_charr_observe(
    duckhts_somalier_charr_accumulator_t *accumulator,
    const duckhts_somalier_counts_t *counts,
    double population_b_frequency,
    const duckhts_somalier_contamination_settings_t *settings) {
    duckhts_somalier_genotype_t genotype;
    duckhts_somalier_status_t status;
    uint64_t depth;
    double contaminant_frequency;
    double infiltrating;
    double contribution;

    if (accumulator == NULL || !isfinite(population_b_frequency) ||
        population_b_frequency < 0.0 ||
        population_b_frequency > 1.0) {
        return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    }
    status = duckhts_somalier_classify_contamination(counts, settings, &genotype);
    if (status != DUCKHTS_SOMALIER_OK) return status;
    if (genotype == DUCKHTS_SOMALIER_UNKNOWN) return DUCKHTS_SOMALIER_OK;
    if (accumulator->usable_sites == UINT64_MAX ||
        (genotype == DUCKHTS_SOMALIER_HOM_A &&
         accumulator->usable_hom_a == UINT64_MAX) ||
        (genotype == DUCKHTS_SOMALIER_HOM_B &&
         accumulator->usable_hom_b == UINT64_MAX)) {
        return DUCKHTS_SOMALIER_LIMIT_EXCEEDED;
    }
    depth = counts_depth(counts);
    contaminant_frequency = genotype == DUCKHTS_SOMALIER_HOM_A
        ? population_b_frequency : 1.0 - population_b_frequency;
    if (contaminant_frequency <= SOMALIER_MIN_CONTAMINANT_AF) {
        return DUCKHTS_SOMALIER_OK;
    }
    infiltrating = genotype == DUCKHTS_SOMALIER_HOM_A
        ? (double)counts->allele_b : (double)counts->allele_a;
    contribution = infiltrating / (contaminant_frequency * (double)depth);
    if (contribution > 0.0) {
        double fraction;
        uint64_t significand;
        uint64_t low;
        uint64_t high;
        uint64_t combined_low;
        uint64_t carry;
        int exponent;
        unsigned shift;

        /* Accepted counts make each contribution a positive finite double in
         * [1e-6, 1e6). Scaling by 2^72 therefore produces an exact integer;
         * 1e8 accepted sites fit in 119 bits. Two limbs make aggregate
         * combination independent of DuckDB's parallel reduction order. */
        fraction = frexp(contribution, &exponent);
        if (!isfinite(contribution) || exponent < -19 || exponent > 20) {
            return DUCKHTS_SOMALIER_NUMERIC_FAILURE;
        }
        significand = (uint64_t)ldexp(fraction, 53);
        shift = (unsigned)(exponent + 19);
        low = significand << shift;
        high = shift == 0u ? 0u : significand >> (64u - shift);
        combined_low = accumulator->contribution_scaled_low + low;
        carry = combined_low < accumulator->contribution_scaled_low;
        if (high > UINT64_MAX - carry ||
            accumulator->contribution_scaled_high > UINT64_MAX - high - carry) {
            return DUCKHTS_SOMALIER_LIMIT_EXCEEDED;
        }
        accumulator->contribution_scaled_low = combined_low;
        accumulator->contribution_scaled_high += high + carry;
    }
    accumulator->usable_sites++;
    if (genotype == DUCKHTS_SOMALIER_HOM_A) accumulator->usable_hom_a++;
    else accumulator->usable_hom_b++;
    return DUCKHTS_SOMALIER_OK;
}

duckhts_somalier_status_t duckhts_somalier_charr_combine(
    duckhts_somalier_charr_accumulator_t *target,
    const duckhts_somalier_charr_accumulator_t *source) {
    uint64_t low;
    uint64_t carry;
    uint64_t high;
    if (target == NULL || source == NULL) {
        return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    }
    if (source->usable_sites > UINT64_MAX - target->usable_sites ||
        source->usable_hom_a > UINT64_MAX - target->usable_hom_a ||
        source->usable_hom_b > UINT64_MAX - target->usable_hom_b) {
        return DUCKHTS_SOMALIER_LIMIT_EXCEEDED;
    }
    low = target->contribution_scaled_low + source->contribution_scaled_low;
    carry = low < target->contribution_scaled_low;
    if (source->contribution_scaled_high > UINT64_MAX - carry ||
        target->contribution_scaled_high >
            UINT64_MAX - source->contribution_scaled_high - carry) {
        return DUCKHTS_SOMALIER_LIMIT_EXCEEDED;
    }
    high = target->contribution_scaled_high +
        source->contribution_scaled_high + carry;
    target->contribution_scaled_low = low;
    target->contribution_scaled_high = high;
    target->usable_sites += source->usable_sites;
    target->usable_hom_a += source->usable_hom_a;
    target->usable_hom_b += source->usable_hom_b;
    return DUCKHTS_SOMALIER_OK;
}

duckhts_somalier_status_t duckhts_somalier_charr_finish(
    const duckhts_somalier_charr_accumulator_t *accumulator,
    duckhts_somalier_charr_result_t *result) {
    double contribution_sum;
    if (result == NULL || accumulator == NULL ||
        (accumulator->usable_sites == 0u &&
         (accumulator->contribution_scaled_low != 0u ||
          accumulator->contribution_scaled_high != 0u)) ||
        accumulator->usable_hom_a > accumulator->usable_sites ||
        accumulator->usable_hom_b >
            accumulator->usable_sites - accumulator->usable_hom_a ||
        accumulator->usable_hom_a + accumulator->usable_hom_b !=
            accumulator->usable_sites) {
        return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    }
    contribution_sum = ldexp((double)accumulator->contribution_scaled_high, -8) +
        ldexp((double)accumulator->contribution_scaled_low, -72);
    if (!isfinite(contribution_sum)) return DUCKHTS_SOMALIER_NUMERIC_FAILURE;
    memset(result, 0, sizeof(*result));
    result->estimate = NAN;
    if (accumulator->usable_sites == 0u) {
        result->status = DUCKHTS_SOMALIER_NO_EVIDENCE;
        return result->status;
    }
    result->status = DUCKHTS_SOMALIER_OK;
    result->estimate = contribution_sum / (double)accumulator->usable_sites;
    result->usable_sites = accumulator->usable_sites;
    result->usable_hom_a = accumulator->usable_hom_a;
    result->usable_hom_b = accumulator->usable_hom_b;
    return DUCKHTS_SOMALIER_OK;
}

duckhts_somalier_status_t duckhts_somalier_charr(
    const duckhts_somalier_counts_t *counts,
    const double *population_b_frequency,
    size_t site_count,
    const duckhts_somalier_contamination_settings_t *settings,
    duckhts_somalier_charr_result_t *result) {
    duckhts_somalier_charr_accumulator_t accumulator = {0};
    size_t i;

    if (result == NULL) return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    memset(result, 0, sizeof(*result));
    result->status = DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    result->estimate = NAN;
    if (!contamination_settings_valid(settings) ||
        (site_count > 0u && (counts == NULL || population_b_frequency == NULL))) {
        return result->status;
    }
    if (site_count > settings->max_sites) {
        result->status = DUCKHTS_SOMALIER_LIMIT_EXCEEDED;
        return result->status;
    }
    for (i = 0u; i < site_count; i++) {
        duckhts_somalier_status_t status = duckhts_somalier_charr_observe(
            &accumulator, &counts[i], population_b_frequency[i], settings);
        if (status != DUCKHTS_SOMALIER_OK) {
            result->status = status;
            return status;
        }
    }
    return duckhts_somalier_charr_finish(&accumulator, result);
}

void duckhts_somalier_matched_settings_default(
    duckhts_somalier_matched_settings_t *settings) {
    if (settings == NULL) return;
    settings->max_depth = DUCKHTS_SOMALIER_MAX_BINOMIAL_DEPTH;
    settings->max_sites = 100000000u;
    settings->max_evaluations = 4096u;
    settings->error_rate = 0.002;
    settings->min_probability = 1e-10;
    settings->min_prior_frequency = 1e-6;
    settings->alpha_min = 0.0;
    settings->alpha_max = 1.0;
    settings->grid_step = 0.01;
    settings->refine_tolerance = 1e-10;
}

static duckhts_somalier_status_t matched_view_valid(
    const duckhts_somalier_matched_view_t *view,
    const duckhts_somalier_matched_settings_t *settings,
    uint64_t *usable_sites) {
    size_t i;
    uint64_t usable = 0u;
    if (view == NULL || usable_sites == NULL || !matched_settings_valid(settings)) {
        return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    }
    if (view->site_count > settings->max_sites) return DUCKHTS_SOMALIER_LIMIT_EXCEEDED;
    if (view->site_count == 0u) return DUCKHTS_SOMALIER_NO_EVIDENCE;
    if (view->receiver_a == NULL || view->receiver_b == NULL ||
        view->receiver_usable == NULL ||
        view->anchor_genotype == NULL || view->population_b_frequency == NULL) {
        return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    }
    for (i = 0u; i < view->site_count; i++) {
        uint64_t depth = (uint64_t)view->receiver_a[i] + (uint64_t)view->receiver_b[i];
        double frequency = view->population_b_frequency[i];
        int anchor = view->anchor_genotype[i];
        if (view->receiver_usable[i] > 1u ||
            (anchor != DUCKHTS_SOMALIER_UNKNOWN &&
             anchor != DUCKHTS_SOMALIER_HOM_A &&
             anchor != DUCKHTS_SOMALIER_HET &&
             anchor != DUCKHTS_SOMALIER_HOM_B)) {
            return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
        }
        if (view->receiver_usable[i] && (depth == 0u || depth > settings->max_depth)) {
            return depth > settings->max_depth
                ? DUCKHTS_SOMALIER_LIMIT_EXCEEDED
                : DUCKHTS_SOMALIER_INVALID_ARGUMENT;
        }
        if (!isfinite(frequency) || frequency < 0.0 || frequency > 1.0) {
            return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
        }
        if (view->receiver_usable[i] &&
            (anchor == DUCKHTS_SOMALIER_HOM_A ||
             anchor == DUCKHTS_SOMALIER_HOM_B)) {
            usable++;
        }
    }
    *usable_sites = usable;
    return usable == 0u ? DUCKHTS_SOMALIER_NO_EVIDENCE : DUCKHTS_SOMALIER_OK;
}

static double clamp_double(double value, double lo, double hi) {
    if (value < lo) return lo;
    if (value > hi) return hi;
    return value;
}

static double logsumexp3(const double values[3]) {
    double maximum = values[0];
    double sum;
    if (values[1] > maximum) maximum = values[1];
    if (values[2] > maximum) maximum = values[2];
    sum = exp(values[0] - maximum) + exp(values[1] - maximum) +
        exp(values[2] - maximum);
    return maximum + log(sum);
}

static duckhts_somalier_status_t matched_log_likelihood_unchecked(
    const duckhts_somalier_matched_view_t *view,
    double alpha,
    const duckhts_somalier_matched_settings_t *settings,
    double *log_likelihood) {
    double total = 0.0;
    size_t i;

    for (i = 0u; i < view->site_count; i++) {
        if (!view->receiver_usable[i] ||
            (view->anchor_genotype[i] != DUCKHTS_SOMALIER_HOM_A &&
             view->anchor_genotype[i] != DUCKHTS_SOMALIER_HOM_B)) {
            continue;
        }
        double frequency = clamp_double(view->population_b_frequency[i],
            settings->min_prior_frequency, 1.0 - settings->min_prior_frequency);
        double clean_b = view->anchor_genotype[i] == DUCKHTS_SOMALIER_HOM_B ? 1.0 : 0.0;
        double priors[3] = {
            (1.0 - frequency) * (1.0 - frequency),
            2.0 * frequency * (1.0 - frequency),
            frequency * frequency
        };
        double terms[3];
        unsigned genotype;
        for (genotype = 0u; genotype < 3u; genotype++) {
            double contaminant_b = 0.5 * (double)genotype;
            double latent_b = (1.0 - alpha) * clean_b + alpha * contaminant_b;
            double observed_b = latent_b * (1.0 - settings->error_rate) +
                (1.0 - latent_b) * settings->error_rate;
            observed_b = clamp_double(observed_b, settings->min_probability,
                1.0 - settings->min_probability);
            terms[genotype] = log(priors[genotype]) +
                (double)view->receiver_b[i] * log(observed_b) +
                (double)view->receiver_a[i] * log1p(-observed_b);
        }
        total += logsumexp3(terms);
        if (!isfinite(total)) return DUCKHTS_SOMALIER_NUMERIC_FAILURE;
    }
    *log_likelihood = total;
    return DUCKHTS_SOMALIER_OK;
}

duckhts_somalier_status_t duckhts_somalier_matched_log_likelihood(
    const duckhts_somalier_matched_view_t *view,
    double alpha,
    const duckhts_somalier_matched_settings_t *settings,
    double *log_likelihood) {
    duckhts_somalier_status_t status;
    uint64_t usable_sites;
    if (log_likelihood == NULL || !isfinite(alpha)) {
        return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    }
    status = matched_view_valid(view, settings, &usable_sites);
    if (status != DUCKHTS_SOMALIER_OK) return status;
    if (alpha < settings->alpha_min || alpha > settings->alpha_max) {
        return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    }
    return matched_log_likelihood_unchecked(view, alpha, settings, log_likelihood);
}

static duckhts_somalier_status_t matched_evaluate(
    const duckhts_somalier_matched_view_t *view,
    const duckhts_somalier_matched_settings_t *settings,
    double alpha,
    duckhts_somalier_matched_result_t *result,
    double *likelihood_out) {
    double likelihood;
    duckhts_somalier_status_t status;
    if (result->evaluations >= settings->max_evaluations) {
        return DUCKHTS_SOMALIER_LIMIT_EXCEEDED;
    }
    status = matched_log_likelihood_unchecked(view, alpha, settings, &likelihood);
    if (status != DUCKHTS_SOMALIER_OK) return status;
    result->evaluations++;
    if (likelihood_out != NULL) *likelihood_out = likelihood;
    if (likelihood > result->log_likelihood ||
        (likelihood == result->log_likelihood && alpha < result->alpha)) {
        result->alpha = alpha;
        result->log_likelihood = likelihood;
    }
    return DUCKHTS_SOMALIER_OK;
}

duckhts_somalier_status_t duckhts_somalier_matched_anchor(
    const duckhts_somalier_matched_view_t *view,
    const duckhts_somalier_matched_settings_t *settings,
    duckhts_somalier_matched_result_t *result) {
    const double golden = 0.6180339887498949;
    duckhts_somalier_status_t status;
    duckhts_somalier_matched_result_t working;
    uint64_t usable_sites;
    double range;
    double ratio;
    double best_grid_alpha;
    double lo;
    double hi;
    double c;
    double d;
    double fc;
    double fd;
    size_t whole_steps;
    size_t i;
    double last_alpha = -1.0;

    if (result == NULL) return DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    result->status = DUCKHTS_SOMALIER_INVALID_ARGUMENT;
    result->alpha = NAN;
    result->log_likelihood = -INFINITY;
    result->usable_sites = 0u;
    result->evaluations = 0u;
    status = matched_view_valid(view, settings, &usable_sites);
    if (status != DUCKHTS_SOMALIER_OK) {
        result->status = status;
        return status;
    }
    working = *result;
    working.usable_sites = usable_sites;
    result->usable_sites = usable_sites;
    range = settings->alpha_max - settings->alpha_min;
    ratio = floor(range / settings->grid_step);
    if (!isfinite(ratio) || ratio >= (double)SIZE_MAX) {
        result->status = DUCKHTS_SOMALIER_LIMIT_EXCEEDED;
        return result->status;
    }
    whole_steps = (size_t)ratio;
    if (whole_steps > SIZE_MAX - 2u ||
        whole_steps + 2u > settings->max_evaluations) {
        result->status = DUCKHTS_SOMALIER_LIMIT_EXCEEDED;
        return result->status;
    }
    for (i = 0u; i <= whole_steps; i++) {
        double alpha = settings->alpha_min + (double)i * settings->grid_step;
        if (alpha > settings->alpha_max) alpha = settings->alpha_max;
        status = matched_evaluate(view, settings, alpha, &working, NULL);
        if (status != DUCKHTS_SOMALIER_OK) goto fail;
        last_alpha = alpha;
    }
    if (last_alpha < settings->alpha_max) {
        status = matched_evaluate(view, settings, settings->alpha_max, &working, NULL);
        if (status != DUCKHTS_SOMALIER_OK) goto fail;
    }
    best_grid_alpha = working.alpha;
    lo = fmax(settings->alpha_min, best_grid_alpha - settings->grid_step);
    hi = fmin(settings->alpha_max, best_grid_alpha + settings->grid_step);
    if (hi - lo <= settings->refine_tolerance) {
        working.status = DUCKHTS_SOMALIER_OK;
        *result = working;
        return DUCKHTS_SOMALIER_OK;
    }
    c = hi - golden * (hi - lo);
    d = lo + golden * (hi - lo);
    status = matched_evaluate(view, settings, c, &working, &fc);
    if (status != DUCKHTS_SOMALIER_OK) goto fail;
    status = matched_evaluate(view, settings, d, &working, &fd);
    if (status != DUCKHTS_SOMALIER_OK) goto fail;
    while (hi - lo > settings->refine_tolerance) {
        if (fc < fd) {
            lo = c;
            c = d;
            fc = fd;
            d = lo + golden * (hi - lo);
            status = matched_evaluate(view, settings, d, &working, &fd);
        } else {
            hi = d;
            d = c;
            fd = fc;
            c = hi - golden * (hi - lo);
            status = matched_evaluate(view, settings, c, &working, &fc);
        }
        if (status != DUCKHTS_SOMALIER_OK) goto fail;
    }
    working.status = DUCKHTS_SOMALIER_OK;
    *result = working;
    return DUCKHTS_SOMALIER_OK;

fail:
    result->status = status;
    result->alpha = NAN;
    result->log_likelihood = -INFINITY;
    result->usable_sites = working.usable_sites;
    result->evaluations = working.evaluations;
    return status;
}

const char *duckhts_somalier_status_string(duckhts_somalier_status_t status) {
    switch (status) {
    case DUCKHTS_SOMALIER_OK: return "ok";
    case DUCKHTS_SOMALIER_NO_EVIDENCE: return "no evidence";
    case DUCKHTS_SOMALIER_INVALID_ARGUMENT: return "invalid argument";
    case DUCKHTS_SOMALIER_LIMIT_EXCEEDED: return "limit exceeded";
    case DUCKHTS_SOMALIER_IDENTITY_MISMATCH: return "identity mismatch";
    case DUCKHTS_SOMALIER_CORRUPT_MASK: return "corrupt mask";
    case DUCKHTS_SOMALIER_CORRUPT_RESULT: return "corrupt result";
    case DUCKHTS_SOMALIER_NUMERIC_FAILURE: return "numeric failure";
    default: return "unknown status";
    }
}
