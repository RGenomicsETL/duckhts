/*
 * duckvep_sv.h — structural-variant/CNV predicate facts (INTERNAL).
 *
 * The adapter supplies an explicit structural operation and, for generic CNVs, a
 * sample-aware copy-change direction. This module maps those semantics plus
 * exact span/topology facts to the predicates used by VEP's consequence program.
 */
#ifndef DUCKVEP_SV_H
#define DUCKVEP_SV_H

#include "duckvep_classify.h"
#include "duckvep_event.h"

#include <stdint.h>

typedef enum duckvep_breakend_status {
    DUCKVEP_BREAKEND_OK = 0,
    DUCKVEP_BREAKEND_NOT_BREAKEND,
    DUCKVEP_BREAKEND_INVALID,
    DUCKVEP_BREAKEND_POSITION_OVERFLOW
} duckvep_breakend_status_t;

/* Borrowed spans into one ALT (not a comma-separated allele list). VCF 4.5
 * section 5.4 defines the four bracket forms and the two single-breakend forms.
 * Replacement sequence includes the local retained bases: it is not necessarily
 * an inserted sequence. Mate identity/reciprocity needs INFO/MATEID and a join.
 * Coordinates are one-based, with zero retained for a virtual telomeric mate.
 * A single breakend has no mate. Parsing never validates against a reference,
 * guesses chromosome aliases, changes case, or allocates storage. */
typedef struct duckvep_breakend {
    const uint8_t *mate_chrom;
    size_t mate_chrom_length;
    uint64_t mate_position;
    const uint8_t *replacement;
    size_t replacement_length;
    uint8_t local_join_after;
    uint8_t mate_extends_right;
    uint8_t has_mate;
} duckvep_breakend_t;

/* Non-BND alleles return NOT_BREAKEND; malformed bracket or terminal-dot forms fail.
 * `out` is zeroed on failure, and borrows `alternate` only on success. */
duckvep_breakend_status_t duckvep_breakend_parse(
    const uint8_t *alternate, size_t length, duckvep_breakend_t *out);

typedef struct duckvep_sv_effect {
    uint8_t deletion;
    uint8_t insertion;
    uint8_t copy_number_gain;
    uint8_t copy_number_loss;
    uint8_t chromosome_breakpoint;
    uint8_t feature_ablation;
    uint8_t feature_amplification;
    uint8_t feature_elongation;
    uint8_t feature_truncation;
    uint8_t start_lost;
    uint8_t start_retained;
    uint8_t stop_lost;
    uint8_t frameshift;
    uint8_t inframe_deletion;
} duckvep_sv_effect_t;

duckvep_sv_effect_t duckvep_sv_effect_fill(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          translateable_cds_length,
    const duckvep_event_t            *event,
    const duckvep_region_state_t     *region);

#endif /* DUCKVEP_SV_H */
