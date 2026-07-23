/*
 * duckvep_compat.h — versioned executable-compatibility policy (INTERNAL).
 *
 * DuckVEP's release target is Ensembl VEP 116.  Most of the kernel implements
 * that declared semantics directly.  A smaller set of outputs instead depends
 * on implementation details of the VEP 116 Perl/BioPerl execution path.  Keep
 * those details behind named flags so they cannot silently become universal
 * biological rules or leak into a later VEP compatibility profile.
 */
#ifndef DUCKVEP_COMPAT_H
#define DUCKVEP_COMPAT_H

#include "duckvep_codon.h"

#include <stdint.h>

typedef enum duckvep_compat_flag {
    /* hgvs_protein::_trim_incomplete_codon assigns instead of compares its
     * keep-length guard, retaining every alternate CDS with length >= 3. */
    DUCKVEP_COMPAT_HGVS_INCOMPLETE_CODON_ASSIGNMENT = UINT32_C(1) << 0,
    /* _stop_loss_extra_AA translates without the transcript codon-table
     * argument, so BioPerl silently selects NCBI table 1. */
    DUCKVEP_COMPAT_HGVS_LATE_STOP_STANDARD_TABLE = UINT32_C(1) << 1,
    /* A pure insertion in an incomplete terminal codon uses distinct
     * consequence and protein-HGVS peptide views. */
    DUCKVEP_COMPAT_HGVS_TERMINAL_PARTIAL_INSERTION = UINT32_C(1) << 2,
    /* Perl substr with a negative start can produce VEP's position-zero
     * start-of-peptide insertion representation. */
    DUCKVEP_COMPAT_HGVS_NEGATIVE_SUBSTR = UINT32_C(1) << 3,
    /* VEP's protein formatter conflates BioPerl Xaa with Ter. */
    DUCKVEP_COMPAT_HGVS_XAA_AS_TER = UINT32_C(1) << 4
} duckvep_compat_flag_t;

typedef struct duckvep_compat_policy {
    uint32_t flags;
} duckvep_compat_policy_t;

static inline int duckvep_compat_profile_valid(
    duckvep_compat_profile_t profile) {

    return profile == DUCKVEP_COMPAT_VEP_116 ||
           profile == DUCKVEP_COMPAT_STRICT;
}

/* This table is the sole executable authority for release-specific language
 * leaks.  Algorithms query the policy; they do not embed VEP-version tests or
 * substitute codon tables locally. */
static inline duckvep_compat_policy_t duckvep_compat_policy(
    duckvep_compat_profile_t profile) {

    duckvep_compat_policy_t policy;

    policy.flags = UINT32_C(0);
    if (profile == DUCKVEP_COMPAT_VEP_116) {
        policy.flags =
            (uint32_t)(DUCKVEP_COMPAT_HGVS_INCOMPLETE_CODON_ASSIGNMENT |
                       DUCKVEP_COMPAT_HGVS_LATE_STOP_STANDARD_TABLE |
                       DUCKVEP_COMPAT_HGVS_TERMINAL_PARTIAL_INSERTION |
                       DUCKVEP_COMPAT_HGVS_NEGATIVE_SUBSTR |
                       DUCKVEP_COMPAT_HGVS_XAA_AS_TER);
    }
    return policy;
}

static inline int duckvep_compat_enabled(
    duckvep_compat_profile_t profile,
    duckvep_compat_flag_t    flag) {

    return (duckvep_compat_policy(profile).flags & (uint32_t)flag) != 0u;
}

static inline duckvep_codon_table_t duckvep_compat_late_stop_codon_table(
    duckvep_compat_profile_t profile,
    duckvep_codon_table_t    transcript_table) {

    return duckvep_compat_enabled(
            profile, DUCKVEP_COMPAT_HGVS_LATE_STOP_STANDARD_TABLE)
        ? DUCKVEP_CODON_TABLE_STANDARD : transcript_table;
}

#endif /* DUCKVEP_COMPAT_H */
