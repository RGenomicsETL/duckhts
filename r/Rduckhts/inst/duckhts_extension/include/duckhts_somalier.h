#ifndef DUCKHTS_SOMALIER_H
#define DUCKHTS_SOMALIER_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Behavioral reference: Somalier v0.3.4 at
 * ff58fdade8f4f8293d904f10e0a4a13f1fac808d (MIT, Brent Pedersen).
 * This API uses caller-owned storage and does not read Somalier files. */

#define DUCKHTS_SOMALIER_RELATEDNESS_MAX_OTHER_FRACTION 0.10
#define DUCKHTS_SOMALIER_CONTAMINATION_MAX_OTHER_FRACTION 0.04
#define DUCKHTS_SOMALIER_MAX_BINOMIAL_DEPTH UINT64_C(1000000)
#define DUCKHTS_SOMALIER_MAX_IDENTITY_BYTES UINT32_C(1024)

typedef enum duckhts_somalier_status {
    DUCKHTS_SOMALIER_OK = 0,
    DUCKHTS_SOMALIER_NO_EVIDENCE,
    DUCKHTS_SOMALIER_INVALID_ARGUMENT,
    DUCKHTS_SOMALIER_LIMIT_EXCEEDED,
    DUCKHTS_SOMALIER_IDENTITY_MISMATCH,
    DUCKHTS_SOMALIER_CORRUPT_MASK,
    DUCKHTS_SOMALIER_NUMERIC_FAILURE,
    DUCKHTS_SOMALIER_CORRUPT_RESULT
} duckhts_somalier_status_t;

typedef enum duckhts_somalier_genotype {
    DUCKHTS_SOMALIER_UNKNOWN = -1,
    DUCKHTS_SOMALIER_HOM_A = 0,
    DUCKHTS_SOMALIER_HET = 1,
    DUCKHTS_SOMALIER_HOM_B = 2
} duckhts_somalier_genotype_t;

/* Counts are aligned to the declared panel's ordered A/B alleles. Other reads
 * are measured reads matching neither allele, not inferred from total depth. */
typedef struct duckhts_somalier_counts {
    uint32_t allele_a;
    uint32_t allele_b;
    uint32_t other;
    uint8_t available;
} duckhts_somalier_counts_t;

typedef struct duckhts_somalier_relatedness_settings {
    uint64_t min_depth;
    size_t max_sites;
    double min_het_balance;
    double hom_balance_cutoff;
} duckhts_somalier_relatedness_settings_t;

typedef struct duckhts_somalier_classification {
    duckhts_somalier_genotype_t genotype;
    uint8_t middling_balance;
    uint8_t unavailable;
} duckhts_somalier_classification_t;

void duckhts_somalier_relatedness_settings_default(
    duckhts_somalier_relatedness_settings_t *settings);

/* Autosomal relatedness classification follows Somalier's permissive 0.10
 * other-read ceiling. UNKNOWN is a valid classification, not an API error. */
duckhts_somalier_status_t duckhts_somalier_classify_relatedness(
    const duckhts_somalier_counts_t *counts,
    const duckhts_somalier_relatedness_settings_t *settings,
    duckhts_somalier_genotype_t *genotype);

duckhts_somalier_status_t duckhts_somalier_classify_relatedness_evidence(
    const duckhts_somalier_counts_t *counts,
    const duckhts_somalier_relatedness_settings_t *settings,
    duckhts_somalier_classification_t *classification);

/* panel_sha256 identifies the assembly and every ordered (region, position,
 * allele A, allele B) tuple in the autosomal, diploid, biallelic-SNP domain.
 * Public panel validation rejects the exact X/Y aliases excluded by pinned
 * Somalier v0.3.4 and requires lexical A < B. Other contig aliases cannot be
 * classified biologically from the region string alone. Diploidy is the
 * three-genotype calculation model, not inferred from count evidence. A
 * sketch also binds the settings that determined its genotype masks. */
typedef struct duckhts_somalier_sketch_identity {
    uint8_t panel_sha256[32];
    uint64_t min_depth;
    double min_het_balance;
    double hom_balance_cutoff;
} duckhts_somalier_sketch_identity_t;

/* Each mask has exactly ceil(site_count / 64) words. Its three aligned byte
 * spans do not overlap. The unused high bits of the final word are zero, and
 * no site may occur in more than one mask. count_digest is an order-independent
 * receipt over each ordinal, availability bit, and raw A/B/other count tuple.
 * content_digest is a deterministic, unkeyed integrity checksum of that
 * receipt, the panel digest, settings, counters, shape, and mask words. These
 * digests detect accidental changes but do not authenticate a sketch against
 * an actor who can reseal it. */
typedef struct duckhts_somalier_masks {
    uint64_t *hom_a;
    uint64_t *het;
    uint64_t *hom_b;
    size_t site_count;
    size_t word_count;
    uint64_t middling_balance_count;
    uint64_t unavailable_count;
    duckhts_somalier_sketch_identity_t identity;
    uint64_t count_digest[2];
    uint64_t content_digest[2];
} duckhts_somalier_masks_t;

/* Add one raw observation to a count receipt. Addition modulo 2^64 makes
 * receipts independent of input order and associatively composable across
 * parallel aggregate states. A fresh receipt is two zero words. */
duckhts_somalier_status_t duckhts_somalier_count_digest_observe(
    uint64_t digest[2], uint64_t site_index,
    const duckhts_somalier_counts_t *counts);

duckhts_somalier_status_t duckhts_somalier_count_digest_combine(
    uint64_t target[2], const uint64_t source[2]);

duckhts_somalier_status_t duckhts_somalier_mask_word_count(
    size_t site_count, size_t *word_count);

duckhts_somalier_status_t duckhts_somalier_prepare_masks(
    const duckhts_somalier_counts_t *counts,
    size_t site_count,
    const uint8_t panel_sha256[32],
    const duckhts_somalier_relatedness_settings_t *settings,
    duckhts_somalier_masks_t *masks);

duckhts_somalier_status_t duckhts_somalier_validate_masks(
    const duckhts_somalier_masks_t *masks);

/* Call after filling caller-owned masks. prepare_masks seals automatically.
 * validate_masks checks the stored digest before a pair calculation. */
duckhts_somalier_status_t duckhts_somalier_seal_masks(
    duckhts_somalier_masks_t *masks);

typedef struct duckhts_somalier_pair_stats {
    uint64_t ibs0;
    uint64_t ibs2;
    uint64_t jointly_called;
    uint64_t shared_hets;
    uint64_t shared_hom_b;
    uint64_t het_ab;
    uint64_t het_count_a;
    uint64_t het_count_b;
    uint64_t hom_b_count_a;
    uint64_t hom_b_count_b;
    uint64_t callable_hom_count_a;
    uint64_t callable_hom_count_b;
    uint64_t matching_hom_count;
    double relatedness;
    double inferred_hom_concordance;
    double raw_hom_b_concordance;
    double p_middling_a;
    double p_middling_b;
    double adjusted_concordance;
} duckhts_somalier_pair_stats_t;

/* Both sketches must have the same panel digest and classification settings.
 * relatedness = 2 * (shared_hets - 2 * ibs0) / max(1, het_ab).
 * Concordance follows v0.3.4's formula but is evaluated and reported in
 * double; this API does not promise its float32 output rounding. */
duckhts_somalier_status_t duckhts_somalier_pair_stats(
    const duckhts_somalier_masks_t *a,
    const duckhts_somalier_masks_t *b,
    duckhts_somalier_pair_stats_t *result);

/* Recompute all integer and floating metrics without allocation. The reported
 * status is OK when jointly_called > 0, NO_EVIDENCE otherwise. */
duckhts_somalier_status_t duckhts_somalier_verify_pair_result(
    const duckhts_somalier_masks_t *a,
    const duckhts_somalier_masks_t *b,
    const duckhts_somalier_pair_stats_t *reported,
    duckhts_somalier_status_t reported_status);

typedef struct duckhts_somalier_contamination_settings {
    uint64_t min_depth;
    uint64_t max_depth;
    size_t max_sites;
    double hom_minor_rate;
    double hom_tail_alpha;
} duckhts_somalier_contamination_settings_t;

void duckhts_somalier_charr_settings_default(
    duckhts_somalier_contamination_settings_t *settings);

void duckhts_somalier_matched_anchor_settings_default(
    duckhts_somalier_contamination_settings_t *settings);

/* Returns the largest k for which P[X >= k] >= tail_alpha for
 * X ~ Binomial(depth, minor_rate). max_depth bounds the incomplete-beta
 * evaluations and binary search domain. DuckHTS certifies this implementation
 * through depth 1,000,000; larger depths return LIMIT_EXCEEDED. A final
 * outward-rounded comparison returns NUMERIC_FAILURE rather than guessing
 * when the exact tail and cutoff cannot be separated. */
duckhts_somalier_status_t duckhts_somalier_binomial_max_minor(
    uint64_t depth,
    double minor_rate,
    double tail_alpha,
    uint64_t max_depth,
    uint64_t *max_minor);

/* Contamination eligibility uses the stricter 0.04 other-read ceiling and a
 * depth-aware homozygous-like binomial test. */
duckhts_somalier_status_t duckhts_somalier_classify_contamination(
    const duckhts_somalier_counts_t *counts,
    const duckhts_somalier_contamination_settings_t *settings,
    duckhts_somalier_genotype_t *genotype);

/* Receiver eligibility does not require a homozygous-like call. An
 * unavailable site is not a measured zero-depth site. */
duckhts_somalier_status_t duckhts_somalier_contamination_usable(
    const duckhts_somalier_counts_t *counts,
    const duckhts_somalier_contamination_settings_t *settings,
    uint8_t *usable);

typedef struct duckhts_somalier_charr_result {
    duckhts_somalier_status_t status;
    double estimate;
    uint64_t usable_sites;
    uint64_t usable_hom_a;
    uint64_t usable_hom_b;
} duckhts_somalier_charr_result_t;

typedef struct duckhts_somalier_charr_accumulator {
    uint64_t contribution_scaled_low;
    uint64_t contribution_scaled_high;
    uint64_t usable_sites;
    uint64_t usable_hom_a;
    uint64_t usable_hom_b;
} duckhts_somalier_charr_accumulator_t;

duckhts_somalier_status_t duckhts_somalier_charr_observe(
    duckhts_somalier_charr_accumulator_t *accumulator,
    const duckhts_somalier_counts_t *counts,
    double population_b_frequency,
    const duckhts_somalier_contamination_settings_t *settings);

duckhts_somalier_status_t duckhts_somalier_charr_combine(
    duckhts_somalier_charr_accumulator_t *target,
    const duckhts_somalier_charr_accumulator_t *source);

duckhts_somalier_status_t duckhts_somalier_charr_finish(
    const duckhts_somalier_charr_accumulator_t *accumulator,
    duckhts_somalier_charr_result_t *result);

/* population_b_frequency is aligned to the same ordered B allele as counts. */
duckhts_somalier_status_t duckhts_somalier_charr(
    const duckhts_somalier_counts_t *counts,
    const double *population_b_frequency,
    size_t site_count,
    const duckhts_somalier_contamination_settings_t *settings,
    duckhts_somalier_charr_result_t *result);

/* A prepared view is aligned to the full ordered panel. receiver_usable marks
 * sites accepted by the receiver filter; the anchor genotype is HOM_A/HOM_B
 * only where its independent homozygous-like filter accepts the site. Counts
 * and population B frequencies use the same orientation. All arrays are
 * borrowed and the pair reduction allocates no workspace. */
typedef struct duckhts_somalier_matched_view {
    const uint32_t *receiver_a;
    const uint32_t *receiver_b;
    const uint8_t *receiver_usable;
    const int8_t *anchor_genotype;
    const double *population_b_frequency;
    size_t site_count;
} duckhts_somalier_matched_view_t;

typedef struct duckhts_somalier_matched_settings {
    uint64_t max_depth;
    size_t max_sites;
    size_t max_evaluations;
    double error_rate;
    double min_probability;
    double min_prior_frequency;
    double alpha_min;
    double alpha_max;
    double grid_step;
    double refine_tolerance;
} duckhts_somalier_matched_settings_t;

void duckhts_somalier_matched_settings_default(
    duckhts_somalier_matched_settings_t *settings);

/* The score omits the alpha-independent binomial coefficient for each site.
 * It is suitable for comparing alpha values on the same prepared view, but is
 * not an absolute likelihood comparable across different observations. */
duckhts_somalier_status_t duckhts_somalier_matched_log_likelihood(
    const duckhts_somalier_matched_view_t *view,
    double alpha,
    const duckhts_somalier_matched_settings_t *settings,
    double *log_likelihood);

typedef struct duckhts_somalier_matched_result {
    duckhts_somalier_status_t status;
    double alpha;
    double log_likelihood;
    uint64_t usable_sites;
    size_t evaluations;
} duckhts_somalier_matched_result_t;

/* The reported point is the best of every declared grid and local-refinement
 * candidate evaluated by the call. Equal likelihoods select the lower alpha. */
duckhts_somalier_status_t duckhts_somalier_matched_anchor(
    const duckhts_somalier_matched_view_t *view,
    const duckhts_somalier_matched_settings_t *settings,
    duckhts_somalier_matched_result_t *result);

const char *duckhts_somalier_status_string(duckhts_somalier_status_t status);

#ifdef __cplusplus
}
#endif

#endif
