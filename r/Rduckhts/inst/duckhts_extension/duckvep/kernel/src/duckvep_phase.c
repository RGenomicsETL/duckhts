#include "duckvep_phase.h"

#include <stddef.h>
#include <string.h>

duckvep_phase_status_t duckvep_phase_observe(
    duckvep_phase_summary_t *s, int32_t allele, uint8_t phase_before) {
    if (!s || allele < -1 || phase_before > 1u) return DUCKVEP_PHASE_INVALID_ARG;
    if (s->ploidy == UINT16_MAX) return DUCKVEP_PHASE_PLOIDY_LIMIT;
    if (!s->ploidy) {
        s->first_allele = allele;
        s->homozygous = allele >= 0;
    } else if (allele < 0 || allele != s->first_allele) {
        s->homozygous = 0u;
    }
    if (!phase_before) {
        if (!s->unphased_count) {
            s->unphased_allele = allele;
            s->unphased_equal = allele >= 0;
        } else if (allele < 0 || allele != s->unphased_allele) {
            s->unphased_equal = 0u;
        }
        s->unphased_count++;
    }
    s->ploidy++;
    return DUCKVEP_PHASE_OK;
}

duckvep_phase_status_t duckvep_phase_assign(
    const duckvep_phase_summary_t *s, uint16_t slot1, uint16_t called_before,
    int32_t allele, uint8_t phase_before, duckvep_phase_policy_t policy,
    duckvep_phase_assignment_t *out) {
    if (!out) return DUCKVEP_PHASE_INVALID_ARG;
    memset(out, 0, sizeof(*out));
    if (!s || !s->ploidy || !slot1 || slot1 > s->ploidy || called_before >= slot1 ||
        allele < -1 || phase_before > 1u ||
        (policy != DUCKVEP_PHASE_STRICT && policy != DUCKVEP_PHASE_VEP116_COMPAT))
        return DUCKVEP_PHASE_INVALID_ARG;
    out->lane = slot1;
    out->status = allele < 0 ? DUCKVEP_PHASE_MISSING : DUCKVEP_PHASE_CALLED;
    if (policy == DUCKVEP_PHASE_VEP116_COMPAT) {
        /* BaseVCF4::get_samples_genotypes omits missing entries before
         * Transcript::get_genotypes assigns allele slots and phased=1. */
        out->lane = allele < 0 ? 0u : (uint16_t)(called_before + 1u);
        out->scope = DUCKVEP_PHASE_ALLELE_SLOT;
    } else if (s->ploidy == 1u || s->homozygous) {
        out->scope = DUCKVEP_PHASE_ALL_SETS;
    } else if (phase_before || s->unphased_count == 1u || s->unphased_equal) {
        /* Permuting the remaining unphased slots cannot change their allele
         * assignment if there is only one, or all contain the same known allele. */
        out->scope = DUCKVEP_PHASE_SET;
    } else {
        out->lane = 0u;
        out->scope = DUCKVEP_PHASE_UNRESOLVED;
        if (allele >= 0) out->status = DUCKVEP_PHASE_UNPHASED;
    }
    return DUCKVEP_PHASE_OK;
}
