/* GT/PS phase interpretation (INTERNAL). No allocation or host types.
 * Observe every decoded allele of one call into a zero-initialized summary,
 * then resolve those same slots against the completed summary. Allele -1 is
 * missing, 0 is REF, positive values are source ALT ordinals. phase_before is
 * the decoded per-allele flag, including HTSlib's leading-slot normalization.
 * The caller retains sample, chromosome, raw GT and nullable PS provenance.
 */
#ifndef DUCKVEP_PHASE_H
#define DUCKVEP_PHASE_H

#include <stdint.h>
#include <stddef.h>

typedef enum {
    DUCKVEP_RAW_GT_OK = 0,
    DUCKVEP_RAW_GT_INVALID_ARG,
    DUCKVEP_RAW_GT_INVALID_SYNTAX,
    DUCKVEP_RAW_GT_ALLELE_OUT_OF_RANGE,
    DUCKVEP_RAW_GT_PLOIDY_LIMIT
} duckvep_raw_gt_status_t;

typedef enum {
    DUCKVEP_RAW_GT_INVALID = 0,
    DUCKVEP_RAW_GT_OMITTED_REFERENCE,
    DUCKVEP_RAW_GT_OMITTED_EMPTY,
    DUCKVEP_RAW_GT_RETAINED
} duckvep_raw_gt_disposition_t;

typedef struct {
    uint32_t allele_index[2]; /* UINT32_MAX is an undefined container slot, not REF. */
    uint32_t parsed_slots;    /* Includes slots beyond the file profile's two lanes. */
    uint16_t source_ploidy;
    uint8_t source_has_missing;
    duckvep_raw_gt_disposition_t disposition;
} duckvep_raw_gt_t;

/* VEP-116 BaseVCF4::get_samples_genotypes(non_ref_only=1), followed by
 * Haplosaurus's file-input diploid slots. Input is original VCF spelling, not
 * reconstructed typed GT. Source REF/ALTs must be nonempty allele strings;
 * source_alt_count excludes REF. Valid VCF numeric/dot GT grammar is required,
 * including one optional leading phase separator. Malformed syntax and allele
 * indices outside the source record fail instead of emulating Perl warnings.
 * OMITTED_EMPTY differs from a retained call with an undefined second slot:
 * the latter can cause an upstream conditional deletion. This routine only
 * resolves source ordinals; it does not project or apply that interpretation.
 * Constant native space, no allocation. Errors leave the output zeroed. */
duckvep_raw_gt_status_t duckvep_phase_parse_vep116_raw(
    const uint8_t *gt, size_t length, uint32_t source_alt_count, duckvep_raw_gt_t *out);

typedef enum {
    DUCKVEP_PHASE_STRICT,
    DUCKVEP_PHASE_VEP116_COMPAT,
    DUCKVEP_PHASE_VEP116_RAW
} duckvep_phase_policy_t;

typedef enum {
    DUCKVEP_PHASE_OK,
    DUCKVEP_PHASE_INVALID_ARG,
    DUCKVEP_PHASE_PLOIDY_LIMIT
} duckvep_phase_status_t;

typedef enum {
    DUCKVEP_PHASE_UNRESOLVED,
    DUCKVEP_PHASE_SET,
    DUCKVEP_PHASE_ALL_SETS,
    DUCKVEP_PHASE_ALLELE_SLOT
} duckvep_phase_scope_t;

typedef enum {
    DUCKVEP_PHASE_CALLED,
    DUCKVEP_PHASE_MISSING,
    DUCKVEP_PHASE_UNPHASED
} duckvep_phase_call_status_t;

typedef struct {
    int32_t first_allele, unphased_allele;
    uint16_t ploidy, unphased_count;
    uint8_t homozygous, unphased_equal;
} duckvep_phase_summary_t;

typedef struct {
    uint16_t lane; /* One-based; zero means no justified lane, never REF. */
    duckvep_phase_scope_t scope;
    duckvep_phase_call_status_t status;
} duckvep_phase_assignment_t;

/* Invalid observations and exceeding 65,535 slots leave the summary unchanged. */
duckvep_phase_status_t duckvep_phase_observe(
    duckvep_phase_summary_t *summary, int32_t allele, uint8_t phase_before);

/* slot1 must identify the same allele/phase observation in the completed call.
 * SET uses PS within sample/chromosome; absent PS is a distinct default set.
 * ALL_SETS is phase-invariant: apply to every phase set, not only NULL PS.
 * ALLELE_SLOT ignores separators and PS and ranks only called alleles: the
 * decoded-call subset, not VEP-116 raw GT parsing or file ploidy inference. Missing
 * entries retain provenance but have no assigned compatibility lane.
 * UNRESOLVED must remain in provenance; never insert lane=0 into a carrier index
 * or omit it when deciding whether a completed sequence is fully known.
 * Missing alleles remain MISSING even when strict mode determines their lane.
 * called_before is the number of non-missing observations before slot1.
 */
duckvep_phase_status_t duckvep_phase_assign(
    const duckvep_phase_summary_t *summary, uint16_t slot1, uint16_t called_before,
    int32_t allele, uint8_t phase_before, duckvep_phase_policy_t policy,
    duckvep_phase_assignment_t *out);

#endif
